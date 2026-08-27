// Floating streaming Bernstein separator for the d=4 local-quadratic dual.
//
// On each projected threshold cell, g is a quadratic polynomial in the three
// other centered random coordinates and centered b.  After multiplication by
// b, the dual gap has degree four, hence 126 Bernstein coefficients on each
// exact five-simplex.  This scanner streams all 40,990 simplex orbits without
// materializing the 5.16-million-row LP.  It is a candidate finder only.

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Real = long double;
constexpr int NV = 6;
constexpr int RANDOM = 4;
constexpr int BASIS = 15;
constexpr Real M = 0.2L;
constexpr Real Z = 0.09500753725790037273535596965L;

struct Simplex {
  int b_box{};
  std::array<int, 4> box{};
  std::array<int, 4> key{};
  std::array<std::array<int, 3>, 4> axes{};
  std::array<int, NV> vertex{};
};

struct Row {
  Real value{};
  int simplex{};
  std::array<int, NV> alpha{};
  bool operator<(const Row& other) const { return value < other.value; }
};

struct MinimalIndex {
  std::array<int, NV> alpha{};
  std::array<int, 4> slot{};
};

static Real parse_rat(const std::string& token) {
  auto slash = token.find('/');
  if (slash == std::string::npos) return std::stold(token);
  return std::stold(token.substr(0, slash)) / std::stold(token.substr(slash + 1));
}

static int popcount(int q) { return __builtin_popcount((unsigned)q); }

static long long choose_int(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n-k);
  long long answer = 1;
  for (int j = 1; j <= k; ++j) answer = answer * (n-k+j) / j;
  return answer;
}

static void compositions_rec(int remaining, int axis, std::array<int, NV>& alpha,
                             std::vector<std::array<int, NV>>& output) {
  if (axis == NV-1) { alpha[axis] = remaining; output.push_back(alpha); return; }
  for (int q = 0; q <= remaining; ++q) {
    alpha[axis] = q;
    compositions_rec(remaining-q, axis+1, alpha, output);
  }
}

static std::vector<std::array<int, NV>> compositions(int degree) {
  std::vector<std::array<int, NV>> output;
  std::array<int, NV> alpha{};
  compositions_rec(degree, 0, alpha, output);
  return output;
}

static std::vector<MinimalIndex> minimal_indices() {
  std::vector<MinimalIndex> output;
  for (int a = 0; a < NV; ++a) for (int b = a; b < NV; ++b)
   for (int c = b; c < NV; ++c) for (int d = c; d < NV; ++d) {
    MinimalIndex item; item.slot = {a,b,c,d};
    for (int q : item.slot) ++item.alpha[q];
    output.push_back(item);
  }
  assert(output.size() == 126);
  return output;
}

static Real basis_g(int basis, const std::array<std::array<Real, NV>, 4>& h,
                    int r, int s) {
  if (basis == 0) return 1;
  if (basis <= 4) {
    int j = basis - 1;
    return (h[j][r] + h[j][s]) / 2;
  }
  int at = 5;
  for (int j = 0; j < 4; ++j) for (int k = j; k < 4; ++k, ++at) {
    if (at != basis) continue;
    if (j == k) return h[j][r] * h[j][s];
    return (h[j][r] * h[k][s] + h[j][s] * h[k][r]) / 2;
  }
  throw std::runtime_error("bad basis index");
}

template<class G>
static Real product_blossom(const std::array<Real, NV>& b,
                            const std::array<Real, NV>& u,
                            const std::array<int, 4>& slot, G g) {
  Real answer = 0;
  for (int a = 0; a < 4; ++a) for (int c = 0; c < 4; ++c) if (a != c) {
    int rem[2], q = 0;
    for (int j = 0; j < 4; ++j) if (j != a && j != c) rem[q++] = j;
    answer += b[slot[a]] * u[slot[c]] * g(slot[rem[0]], slot[rem[1]]);
  }
  return answer / 12;
}

static void point_data(const Simplex& simplex,
                       const std::vector<std::array<Real, 5>>& local,
                       std::array<std::array<Real, 5>, NV>& point) {
  for (int v = 0; v < NV; ++v) {
    auto y = local[simplex.vertex[v]];
    for (int j = 0; j < 4; ++j) point[v][j] = (simplex.box[j] + y[j]) / 5;
    point[v][4] = (simplex.b_box + y[4]) / 5;
  }
}

static std::pair<Real, Real> kernel_terms(
    const std::array<std::array<Real, 5>, NV>& point) {
  std::array<Real, 5> center{};
  for (int j = 0; j < 5; ++j)
    for (int v = 0; v < NV; ++v) center[j] += point[v][j] / NV;
  Real sum0 = 0, sum1 = 0;
  for (int state = 0; state < 16; ++state) {
    Real subtotal = 0;
    for (int j = 0; j < 4; ++j) if ((state >> j) & 1) subtotal += center[j];
    Real weight = std::pow(Z, 4 - popcount(state));
    if (subtotal > 1 + 1e-13L) sum0 += weight;
    if (subtotal + center[4] > 1 + 1e-13L) sum1 += weight;
  }
  return {(1 + Z) * sum0, M * (sum1 - sum0)};
}

static void local_arrays(const Simplex& simplex,
                         const std::array<std::array<Real, 5>, NV>& point,
                         int omitted, std::array<Real, NV>& b,
                         std::array<Real, NV>& u,
                         std::array<std::array<Real, NV>, 4>& h) {
  for (int v = 0; v < NV; ++v) {
    b[v] = point[v][4];
    u[v] = point[v][omitted] - M;
    for (int j = 0; j < 3; ++j)
      h[j][v] = point[v][simplex.axes[omitted][j]] - M;
    h[3][v] = point[v][4] - M;
  }
}

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cerr << "usage: scanner geometry_flat coefficients top elevation_degree\n";
    return 2;
  }
  std::ifstream in(argv[1]);
  if (!in) throw std::runtime_error("cannot open geometry");
  std::string tag;
  int vertex_count;
  in >> tag >> vertex_count;
  if (tag != "VERTICES") throw std::runtime_error("bad vertex header");
  std::vector<std::array<Real, 5>> local(vertex_count);
  for (auto& point : local) for (Real& q : point) {
    std::string token; in >> token; q = parse_rat(token);
  }
  int key_count;
  in >> tag >> key_count;
  if (tag != "PROJECTED_KEYS") throw std::runtime_error("bad key header");
  std::string line;
  std::getline(in, line);
  for (int i = 0; i < key_count; ++i) std::getline(in, line);
  int simplex_count;
  in >> tag >> simplex_count;
  if (tag != "SIMPLICES") throw std::runtime_error("bad simplex header");
  std::vector<Simplex> simplices(simplex_count);
  for (auto& simplex : simplices) {
    in >> simplex.b_box;
    for (int& q : simplex.box) in >> q;
    for (int omitted = 0; omitted < 4; ++omitted) {
      in >> simplex.key[omitted];
      for (int& q : simplex.axes[omitted]) in >> q;
    }
    for (int& q : simplex.vertex) in >> q;
  }
  in >> tag;
  if (tag != "END") throw std::runtime_error("missing end marker");

  std::vector<Real> coefficient(key_count * BASIS);
  if (std::string(argv[2]) != "ZERO") {
    std::ifstream cinfile(argv[2]);
    if (!cinfile) throw std::runtime_error("cannot open coefficients");
    for (Real& q : coefficient) if (!(cinfile >> q))
      throw std::runtime_error("too few coefficients");
  }
  int top = std::stoi(argv[3]);
  int elevation_degree = std::stoi(argv[4]);
  if (elevation_degree < 4 || elevation_degree > 12)
    throw std::runtime_error("elevation degree must lie in [4,12]");
  auto minimal = minimal_indices();
  auto elevated = compositions(elevation_degree);
  std::vector<std::vector<std::pair<int, Real>>> elevation(elevated.size());
  Real elevation_denominator = choose_int(elevation_degree, 4);
  for (int eid = 0; eid < (int)elevated.size(); ++eid) {
    for (int mid = 0; mid < (int)minimal.size(); ++mid) {
      Real numerator = 1;
      bool contained = true;
      for (int axis = 0; axis < NV; ++axis) {
        if (minimal[mid].alpha[axis] > elevated[eid][axis]) { contained = false; break; }
        numerator *= choose_int(elevated[eid][axis], minimal[mid].alpha[axis]);
      }
      if (contained && numerator) elevation[eid].push_back({mid, numerator/elevation_denominator});
    }
  }
  std::priority_queue<Row> worst;
  long long row_count = 0;
  Real global_minimum = 1e100L;

  for (int sid = 0; sid < simplex_count; ++sid) {
    const auto& simplex = simplices[sid];
    std::array<std::array<Real, 5>, NV> point{};
    point_data(simplex, local, point);
    auto [slope, constant] = kernel_terms(point);
    std::array<Real, NV> b{};
    for (int v = 0; v < NV; ++v) b[v] = point[v][4];
    std::array<std::array<Real, NV>, RANDOM> u{};
    std::array<std::array<std::array<Real, NV>, 4>, RANDOM> h{};
    std::array<std::array<std::array<Real, NV>, NV>, RANDOM> candidate_g{};
    for (int omitted = 0; omitted < RANDOM; ++omitted) {
      local_arrays(simplex, point, omitted, b, u[omitted], h[omitted]);
      int offset = BASIS * simplex.key[omitted];
      for (int r = 0; r < NV; ++r) for (int s = 0; s < NV; ++s) {
        Real value = 0;
        for (int basis = 0; basis < BASIS; ++basis)
          value += coefficient[offset + basis] * basis_g(basis, h[omitted], r, s);
        candidate_g[omitted][r][s] = value;
      }
    }
    std::array<Real, 126> minimal_value{};
    for (int mid = 0; mid < (int)minimal.size(); ++mid) {
      const auto& slot = minimal[mid].slot;
      int a=slot[0],c=slot[1],d=slot[2],e=slot[3];
      Real b_average = (b[a] + b[c] + b[d] + b[e]) / 4;
      Real value = (1 - slope) * b_average - constant;
      for (int omitted = 0; omitted < RANDOM; ++omitted)
        value += product_blossom(b, u[omitted], slot,
          [&](int r, int s) { return candidate_g[omitted][r][s]; });
      // Closed singleton corrections on top faces.
      for (int axis = 0; axis < RANDOM; ++axis) {
        bool top_face = true;
        for (int q : slot) top_face = top_face && std::abs(point[q][axis] - 1) < 1e-14L;
        if (top_face) value -= std::pow(Z, 3) * (b_average * (1 + Z) - M);
      }
      bool b_top = true;
      for (int q : slot) b_top = b_top && std::abs(point[q][4] - 1) < 1e-14L;
      if (b_top) value -= M * std::pow(Z, 4);
      minimal_value[mid] = value;
    }
    for (int eid = 0; eid < (int)elevated.size(); ++eid) {
      Real value = 0;
      for (auto [mid, weight] : elevation[eid]) value += weight * minimal_value[mid];
      ++row_count;
      global_minimum = std::min(global_minimum, value);
      Row row{value, sid, elevated[eid]};
      if ((int)worst.size() < top) worst.push(row);
      else if (value < worst.top().value) { worst.pop(); worst.push(row); }
    }
    if ((sid + 1) % 5000 == 0)
      std::cerr << "processed " << sid + 1 << " simplices\n";
  }

  std::vector<Row> output;
  while (!worst.empty()) { output.push_back(worst.top()); worst.pop(); }
  std::sort(output.begin(), output.end(), [](const Row& x, const Row& y) {
    if (x.value != y.value) return x.value < y.value;
    if (x.simplex != y.simplex) return x.simplex < y.simplex;
    return x.alpha < y.alpha;
  });
  std::cout << std::setprecision(21);
  std::cout << "SUMMARY rows " << row_count << " minimum " << global_minimum
            << " returned " << output.size() << "\n";
  for (const Row& row : output) {
    const auto& simplex = simplices[row.simplex];
    std::array<std::array<Real, 5>, NV> point{};
    point_data(simplex, local, point);
    auto [slope, constant] = kernel_terms(point);
    std::array<Real, NV> b{};
    for (int v = 0; v < NV; ++v) b[v] = point[v][4];
    int eid = int(std::lower_bound(elevated.begin(), elevated.end(), row.alpha) - elevated.begin());
    assert(eid < (int)elevated.size() && elevated[eid] == row.alpha);
    Real base = 0;
    for (auto [mid, weight] : elevation[eid]) {
      const auto& slot = minimal[mid].slot;
      Real b_average = 0;
      for (int q : slot) b_average += b[q] / 4;
      Real local_base = (1 - slope) * b_average - constant;
      for (int axis = 0; axis < RANDOM; ++axis) {
        bool top_face = true;
        for (int q : slot) top_face = top_face && std::abs(point[q][axis] - 1) < 1e-14L;
        if (top_face) local_base -= std::pow(Z, 3) * (b_average * (1 + Z) - M);
      }
      bool b_top = true;
      for (int q : slot) b_top = b_top && std::abs(point[q][4] - 1) < 1e-14L;
      if (b_top) local_base -= M * std::pow(Z, 4);
      base += weight * local_base;
    }
    std::cout << "ROW " << row.value << ' ' << base << ' ' << row.simplex;
    std::cout << " A";
    for (int q : row.alpha) std::cout << ' ' << q;
    for (int omitted = 0; omitted < RANDOM; ++omitted) {
      std::array<Real, NV> u{};
      std::array<std::array<Real, NV>, 4> h{};
      local_arrays(simplex, point, omitted, b, u, h);
      std::cout << " P " << simplex.key[omitted];
      for (int basis = 0; basis < BASIS; ++basis) {
        Real feature = 0;
        for (auto [mid, weight] : elevation[eid])
          feature += weight * product_blossom(b, u, minimal[mid].slot,
            [&](int r, int s) { return basis_g(basis, h, r, s); });
        std::cout << ' ' << feature;
      }
    }
    std::cout << '\n';
  }
  std::cout << "FLOATING STREAMING DIAGNOSTIC ONLY\n";
}
