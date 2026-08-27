// Exact interval replay for the d=4 projected-cell quadratic ZS5 certificate.
//
// The repaired dual coefficients belong to Q(z), where z is the unique root
// in (19/200,12/125) of
//   3125 z^4 + 11250 z^3 + 15250 z^2 + 9225 z - 1024.
// The companion certificate file gives outward dyadic enclosures for every
// coefficient and z^0,...,z^3.  This verifier rebuilds the exact rational
// arrangement, blossoms every degree-four gap polynomial, and checks all
// 40,990 * C(9,5) = 5,164,740 Bernstein rows.  Rows on the exact active face
// are separately proved to vanish by r5_zs5_d4_quadratic_face_repair_exact.py;
// every other row must have a strictly positive interval lower endpoint here.

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

using i64 = long long;
using i128 = __int128_t;
using u64 = unsigned long long;

constexpr int NV = 6;
constexpr int RANDOM = 4;
constexpr int BASIS = 15;
constexpr int CUSP_DEGREE = 4;

static i128 abs128(i128 q) { return q < 0 ? -q : q; }

static i128 gcd128(i128 a, i128 b) {
  a = abs128(a); b = abs128(b);
  while (b) { i128 r = a % b; a = b; b = r; }
  return a;
}

static i128 checked_add(i128 a, i128 b) {
  i128 answer;
  if (__builtin_add_overflow(a, b, &answer))
    throw std::overflow_error("signed 128-bit addition overflow");
  return answer;
}

static i128 checked_subtract(i128 a, i128 b) {
  i128 answer;
  if (__builtin_sub_overflow(a, b, &answer))
    throw std::overflow_error("signed 128-bit subtraction overflow");
  return answer;
}

static i128 checked_multiply(i128 a, i128 b) {
  i128 answer;
  if (__builtin_mul_overflow(a, b, &answer))
    throw std::overflow_error("signed 128-bit multiplication overflow");
  return answer;
}

static std::string str128(i128 q) {
  if (q == 0) return "0";
  bool negative = q < 0;
  i128 a = negative ? -q : q;
  std::string output;
  while (a) { output.push_back(char('0' + a % 10)); a /= 10; }
  if (negative) output.push_back('-');
  std::reverse(output.begin(), output.end());
  return output;
}

static i128 parse128(const std::string& token) {
  if (token.empty()) throw std::runtime_error("empty integer token");
  bool negative = token[0] == '-';
  size_t at = negative ? 1 : 0;
  i128 answer = 0;
  for (; at < token.size(); ++at) {
    if (token[at] < '0' || token[at] > '9')
      throw std::runtime_error("bad integer token: " + token);
    answer = checked_add(checked_multiply(answer, 10), token[at] - '0');
  }
  return negative ? -answer : answer;
}

struct Rat {
  i64 n = 0, d = 1;
  Rat() = default;
  Rat(i64 a): n(a), d(1) {}
  Rat(i64 a, i64 b) { set128(a, b); }
  void set128(i128 a, i128 b) {
    if (!b) throw std::runtime_error("zero rational denominator");
    if (b < 0) { a = -a; b = -b; }
    i128 g = gcd128(a, b); a /= g; b /= g;
    if (a > std::numeric_limits<i64>::max() ||
        a < std::numeric_limits<i64>::min() ||
        b > std::numeric_limits<i64>::max())
      throw std::overflow_error("rational overflow");
    n = i64(a); d = i64(b);
  }
};

static Rat operator+(Rat a, Rat b) {
  i128 left = checked_multiply(a.n, b.d);
  i128 right = checked_multiply(b.n, a.d);
  Rat q; q.set128(checked_add(left, right), checked_multiply(a.d, b.d)); return q;
}
static Rat operator-(Rat a, Rat b) {
  i128 left = checked_multiply(a.n, b.d);
  i128 right = checked_multiply(b.n, a.d);
  Rat q; q.set128(checked_subtract(left, right), checked_multiply(a.d, b.d)); return q;
}
static Rat operator*(Rat a, Rat b) {
  // Cross cancellation keeps all exact geometry intermediates small.
  i64 g1 = std::gcd(a.n < 0 ? -a.n : a.n, b.d);
  i64 g2 = std::gcd(b.n < 0 ? -b.n : b.n, a.d);
  Rat q; q.set128(checked_multiply(a.n/g1, b.n/g2),
                   checked_multiply(a.d/g2, b.d/g1)); return q;
}
static Rat operator/(Rat a, Rat b) {
  if (!b.n) throw std::runtime_error("division by zero rational");
  return a * Rat(b.d, b.n);
}
static Rat operator-(Rat a) { return Rat(-a.n, a.d); }
static bool operator==(Rat a, Rat b) { return a.n == b.n && a.d == b.d; }
static bool operator<(Rat a, Rat b) {
  return i128(a.n)*b.d < i128(b.n)*a.d;
}
static bool operator>(Rat a, Rat b) { return b < a; }

static Rat parse_rat(const std::string& token) {
  size_t slash = token.find('/');
  if (slash == std::string::npos) return Rat(std::stoll(token));
  return Rat(std::stoll(token.substr(0, slash)),
             std::stoll(token.substr(slash + 1)));
}

// An interval [lo/2^scale_bits, hi/2^scale_bits].  All operations below keep
// the same scale and round outwards exactly.
struct Interval { i128 lo = 0, hi = 0; };

static Interval operator+(Interval a, Interval b) {
  return {checked_add(a.lo, b.lo), checked_add(a.hi, b.hi)};
}

static i128 floor_div(i128 a, i128 b) {
  assert(b > 0);
  i128 q = a / b, r = a % b;
  if (r && a < 0) --q;
  return q;
}

static i128 ceil_div(i128 a, i128 b) {
  assert(b > 0);
  i128 q = a / b, r = a % b;
  if (r && a > 0) ++q;
  return q;
}

static Interval multiply(Rat q, Interval x) {
  if (!q.n) return {};
  i128 left = checked_multiply(i128(q.n), x.lo);
  i128 right = checked_multiply(i128(q.n), x.hi);
  if (left > right) std::swap(left, right);
  return {floor_div(left, q.d), ceil_div(right, q.d)};
}

struct Simplex {
  int b_box{};
  std::array<int, 4> box{};
  std::array<int, 4> key{};
  std::array<std::array<int, 3>, 4> axes{};
  std::array<int, NV> vertex{};
};

struct MinimalIndex {
  std::array<int, NV> alpha{};
  std::array<int, 4> slot{};
};

using Point = std::array<Rat, 5>;
using Poly = std::array<Rat, 6>;

static int popcount(int q) { return __builtin_popcount(unsigned(q)); }

static std::vector<MinimalIndex> minimal_indices() {
  std::vector<MinimalIndex> output;
  for (int a = 0; a < NV; ++a) for (int b = a; b < NV; ++b)
   for (int c = b; c < NV; ++c) for (int d = c; d < NV; ++d) {
    MinimalIndex item; item.slot = {a,b,c,d};
    for (int q : item.slot) ++item.alpha[q];
    output.push_back(item);
  }
  if (output.size() != 126) throw std::runtime_error("bad composition count");
  return output;
}

static Rat basis_g(int basis, const std::array<std::array<Rat, NV>, 4>& h,
                   int r, int s) {
  if (basis == 0) return Rat(1);
  if (basis <= 4) {
    int j = basis - 1;
    return (h[j][r] + h[j][s]) / Rat(2);
  }
  int at = 5;
  for (int j = 0; j < 4; ++j) for (int k = j; k < 4; ++k, ++at) {
    if (at != basis) continue;
    if (j == k) return h[j][r] * h[j][s];
    return (h[j][r] * h[k][s] + h[j][s] * h[k][r]) / Rat(2);
  }
  throw std::runtime_error("bad basis index");
}

static void point_data(const Simplex& simplex,
                       const std::vector<Point>& local,
                       std::array<Point, NV>& point) {
  for (int v = 0; v < NV; ++v) {
    Point y = local[simplex.vertex[v]];
    for (int j = 0; j < 4; ++j)
      point[v][j] = (Rat(simplex.box[j]) + y[j]) / Rat(5);
    point[v][4] = (Rat(simplex.b_box) + y[4]) / Rat(5);
  }
}

static void local_arrays(const Simplex& simplex,
                         const std::array<Point, NV>& point,
                         int omitted, std::array<Rat, NV>& b,
                         std::array<Rat, NV>& u,
                         std::array<std::array<Rat, NV>, 4>& h) {
  const Rat m(1,5);
  for (int v = 0; v < NV; ++v) {
    b[v] = point[v][4];
    u[v] = point[v][omitted] - m;
    for (int j = 0; j < 3; ++j)
      h[j][v] = point[v][simplex.axes[omitted][j]] - m;
    h[3][v] = point[v][4] - m;
  }
}

static Poly base_polynomial(const std::array<Point, NV>& point,
                            const std::array<int, 4>& slot) {
  const Rat m(1,5);
  std::array<Rat, 5> center{};
  for (int j = 0; j < 5; ++j)
    for (int v = 0; v < NV; ++v) center[j] = center[j] + point[v][j] / Rat(NV);

  std::array<int, 5> sum_zero{}, sum_one{};
  for (int state = 0; state < 16; ++state) {
    Rat subtotal;
    for (int j = 0; j < 4; ++j)
      if ((state >> j) & 1) subtotal = subtotal + center[j];
    int degree = 4 - popcount(state);
    if (subtotal > Rat(1)) ++sum_zero[degree];
    if (subtotal + center[4] > Rat(1)) ++sum_one[degree];
  }

  Rat b_average;
  for (int q : slot) b_average = b_average + point[q][4] / Rat(4);
  Poly p{};
  p[0] = p[0] + b_average;
  for (int degree = 0; degree <= 4; ++degree) {
    p[degree] = p[degree] + (m - b_average) * Rat(sum_zero[degree]);
    p[degree] = p[degree] - m * Rat(sum_one[degree]);
    p[degree + 1] = p[degree + 1] - b_average * Rat(sum_zero[degree]);
  }

  // Closed singleton corrections on all random-coordinate top faces.
  for (int axis = 0; axis < RANDOM; ++axis) {
    bool top_face = true;
    for (int q : slot) top_face = top_face && point[q][axis] == Rat(1);
    if (top_face) {
      p[3] = p[3] + (m - b_average);
      p[4] = p[4] - b_average;
    }
  }
  bool b_top = true;
  for (int q : slot) b_top = b_top && point[q][4] == Rat(1);
  if (b_top) p[4] = p[4] - m;

  // Reduce exactly modulo the defining quartic of z.
  const std::array<Rat,5> cusp = {
      Rat(-1024), Rat(9225), Rat(15250), Rat(11250), Rat(3125)};
  for (int degree = 5; degree >= CUSP_DEGREE; --degree) {
    Rat lead = p[degree];
    if (!lead.n) continue;
    int shift = degree - CUSP_DEGREE;
    for (int j = 0; j <= CUSP_DEGREE; ++j)
      p[shift+j] = p[shift+j] - lead * cusp[j] / cusp[CUSP_DEGREE];
    if (p[degree].n)
      throw std::runtime_error("quartic reduction failed");
  }
  return p;
}

static Interval evaluate_reduced(const Poly& p,
                                 const std::array<Interval,4>& powers) {
  Interval answer{};
  for (int degree = 0; degree < 4; ++degree)
    answer = answer + multiply(p[degree], powers[degree]);
  if (p[4].n || p[5].n) throw std::runtime_error("unreduced polynomial");
  return answer;
}

static u64 row_code(int simplex, const std::array<int,NV>& alpha) {
  u64 code = 0, place = 1;
  for (int q : alpha) { code += u64(q) * place; place *= 5; }
  return u64(simplex) * place + code;
}

static void hash_word(u64& state, u64 word) {
  // Deterministic FNV-1a audit digest.  This is not part of the inequality
  // proof; it lets independent replays confirm the identical complete stream.
  for (int j = 0; j < 8; ++j) {
    state ^= (word >> (8*j)) & 255;
    state *= 1099511628211ULL;
  }
}

static void hash_i128(u64& state, i128 value) {
  __uint128_t u = static_cast<__uint128_t>(value);
  hash_word(state, static_cast<u64>(u));
  hash_word(state, static_cast<u64>(u >> 64));
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: exact_verify geometry_exact_flat exact_interval_flat\n";
    return 2;
  }

  // Parse the exact geometry.
  std::ifstream geometry(argv[1]);
  if (!geometry) throw std::runtime_error("cannot open geometry");
  std::string tag, token;
  int vertex_count;
  geometry >> tag >> vertex_count;
  if (tag != "VERTICES") throw std::runtime_error("bad vertex header");
  std::vector<Point> local(vertex_count);
  for (Point& point : local) for (Rat& q : point) {
    geometry >> token; q = parse_rat(token);
  }
  int key_count;
  geometry >> tag >> key_count;
  if (tag != "PROJECTED_KEYS") throw std::runtime_error("bad key header");
  std::string line; std::getline(geometry, line);
  for (int key = 0; key < key_count; ++key) std::getline(geometry, line);
  int simplex_count;
  geometry >> tag >> simplex_count;
  if (tag != "SIMPLICES") throw std::runtime_error("bad simplex header");
  std::vector<Simplex> simplices(simplex_count);
  for (Simplex& simplex : simplices) {
    geometry >> simplex.b_box;
    for (int& q : simplex.box) geometry >> q;
    for (int omitted = 0; omitted < RANDOM; ++omitted) {
      geometry >> simplex.key[omitted];
      for (int& q : simplex.axes[omitted]) geometry >> q;
    }
    for (int& q : simplex.vertex) geometry >> q;
  }
  geometry >> tag;
  if (tag != "END") throw std::runtime_error("missing geometry end marker");
  if (vertex_count != 587 || key_count != 545 || simplex_count != 40990)
    throw std::runtime_error("unexpected exact geometry cardinalities");

  // Parse the dyadic coefficient enclosure and the exact-active row manifest.
  std::ifstream certificate(argv[2]);
  if (!certificate) throw std::runtime_error("cannot open certificate");
  certificate >> tag;
  if (tag != "R5_ZS5_D4_QUADRATIC_EXACT_INTERVAL_V1")
    throw std::runtime_error("bad certificate header");
  int scale_bits;
  certificate >> tag >> scale_bits;
  if (tag != "SCALE_BITS" || scale_bits < 80 || scale_bits > 100)
    throw std::runtime_error("bad dyadic scale");
  // Hashes and the rational cusp interval are retained in the certificate for
  // the independent manifest.  The verifier consumes but does not trust them.
  for (const std::string expected : {"REPAIR_SHA256", "ACTIVE_SHA256", "GEOMETRY_SHA256"}) {
    certificate >> tag >> token;
    if (tag != expected) throw std::runtime_error("bad hash header");
  }
  std::string cusp_lower, cusp_upper;
  certificate >> tag >> cusp_lower >> cusp_upper;
  if (tag != "CUSP_INTERVAL") throw std::runtime_error("bad cusp interval header");
  int power_count;
  certificate >> tag >> power_count;
  if (tag != "POWERS" || power_count != 4) throw std::runtime_error("bad power header");
  std::array<Interval,4> powers{};
  for (Interval& q : powers) {
    certificate >> token; q.lo = parse128(token);
    certificate >> token; q.hi = parse128(token);
    if (q.lo > q.hi) throw std::runtime_error("inverted power interval");
  }
  int coefficient_count;
  certificate >> tag >> coefficient_count;
  if (tag != "COEFFICIENT_INTERVALS" || coefficient_count != 6759)
    throw std::runtime_error("bad coefficient interval count");
  std::vector<Interval> collapsed(coefficient_count);
  for (Interval& q : collapsed) {
    certificate >> token; q.lo = parse128(token);
    certificate >> token; q.hi = parse128(token);
    if (q.lo > q.hi) throw std::runtime_error("inverted coefficient interval");
  }
  int mapping_count;
  certificate >> tag >> mapping_count;
  if (tag != "ORIGINAL_TO_COLLAPSED" || mapping_count != key_count*BASIS)
    throw std::runtime_error("bad coefficient map count");
  std::vector<int> mapping(mapping_count);
  for (int& q : mapping) {
    certificate >> q;
    if (q < 0 || q >= coefficient_count) throw std::runtime_error("bad map entry");
  }
  int active_count;
  certificate >> tag >> active_count;
  if (tag != "ACTIVE_ROWS" || active_count != 12794)
    throw std::runtime_error("bad active row count");
  std::unordered_set<u64> active;
  active.reserve(active_count * 2);
  for (int row = 0; row < active_count; ++row) {
    int sid; std::array<int,NV> alpha{};
    certificate >> sid;
    for (int& q : alpha) certificate >> q;
    if (sid < 0 || sid >= simplex_count ||
        std::accumulate(alpha.begin(), alpha.end(), 0) != 4)
      throw std::runtime_error("bad active row metadata");
    active.insert(row_code(sid, alpha));
  }
  certificate >> tag;
  if (tag != "END" || int(active.size()) != active_count)
    throw std::runtime_error("bad certificate end or duplicate active row");

  const auto indices = minimal_indices();
  long long rows = 0, strict_rows = 0, active_rows = 0;
  i128 minimum_strict_lower = 0;
  i128 largest_active_width = 0;
  u64 stream_hash = 14695981039346656037ULL;

  for (int sid = 0; sid < simplex_count; ++sid) {
    const Simplex& simplex = simplices[sid];
    std::array<Point,NV> point{};
    point_data(simplex, local, point);
    std::array<Rat,NV> b{};
    for (int v = 0; v < NV; ++v) b[v] = point[v][4];

    // G_i(r,s) is the quadratic blossom of g_i at the remaining two slots.
    std::array<std::array<std::array<Interval,NV>,NV>,RANDOM> G{};
    std::array<std::array<Rat,NV>,RANDOM> u{};
    for (int omitted = 0; omitted < RANDOM; ++omitted) {
      std::array<std::array<Rat,NV>,4> h{};
      local_arrays(simplex, point, omitted, b, u[omitted], h);
      int offset = BASIS * simplex.key[omitted];
      for (int r = 0; r < NV; ++r) for (int s = 0; s < NV; ++s) {
        Interval value{};
        for (int basis = 0; basis < BASIS; ++basis)
          value = value + multiply(basis_g(basis, h, r, s),
                                   collapsed[mapping[offset+basis]]);
        G[omitted][r][s] = value;
      }
    }

    for (const MinimalIndex& index : indices) {
      const auto& slot = index.slot;
      Interval value = evaluate_reduced(base_polynomial(point, slot), powers);
      for (int omitted = 0; omitted < RANDOM; ++omitted) {
        for (int a = 0; a < 4; ++a) for (int c = 0; c < 4; ++c) if (a != c) {
          int rem[2], at = 0;
          for (int j = 0; j < 4; ++j) if (j != a && j != c) rem[at++] = j;
          Rat factor = b[slot[a]] * u[omitted][slot[c]] / Rat(12);
          value = value + multiply(factor, G[omitted][slot[rem[0]]][slot[rem[1]]]);
        }
      }

      ++rows;
      u64 code = row_code(sid, index.alpha);
      bool is_active = active.find(code) != active.end();
      hash_word(stream_hash, code);
      hash_i128(stream_hash, value.lo);
      hash_i128(stream_hash, value.hi);
      hash_word(stream_hash, is_active ? 1 : 0);
      if (is_active) {
        ++active_rows;
        if (value.lo > 0 || value.hi < 0) {
          std::cerr << "active interval misses zero at simplex " << sid << " alpha";
          for (int q : index.alpha) std::cerr << ' ' << q;
          std::cerr << " interval " << str128(value.lo) << ' ' << str128(value.hi) << '\n';
          return 1;
        }
        largest_active_width = std::max(largest_active_width, value.hi-value.lo);
      } else {
        if (value.lo <= 0) {
          std::cerr << "nonactive row is not certified positive at simplex " << sid << " alpha";
          for (int q : index.alpha) std::cerr << ' ' << q;
          std::cerr << " interval " << str128(value.lo) << ' ' << str128(value.hi) << '\n';
          return 1;
        }
        ++strict_rows;
        if (!minimum_strict_lower || value.lo < minimum_strict_lower)
          minimum_strict_lower = value.lo;
      }
    }
    if ((sid+1) % 5000 == 0)
      std::cerr << "exact interval replay processed " << sid+1 << " simplices\n";
  }

  if (rows != 5164740LL || active_rows != active_count ||
      strict_rows != rows-active_rows)
    throw std::runtime_error("final row counts do not match manifest");

  std::cout << "R5_ZS5_D4_QUADRATIC_EXACT_INTERVAL_PASS\n";
  std::cout << "vertices " << vertex_count << " projected_keys " << key_count
            << " simplices " << simplex_count << " rows " << rows << '\n';
  std::cout << "strict_rows " << strict_rows << " exact_active_rows " << active_rows << '\n';
  std::cout << "minimum_strict_lower_scaled " << str128(minimum_strict_lower)
            << " scale_bits " << scale_bits << '\n';
  std::cout << "largest_active_interval_width_scaled "
            << str128(largest_active_width) << '\n';
  std::cout << "row_stream_fnv1a64 " << std::hex << std::setfill('0')
            << std::setw(16) << stream_hash << std::dec << '\n';
  std::cout << "closed_top_face_corrections INCLUDED\n";
  std::cout << "EXACT ACTIVE IDENTITIES REQUIRED FROM FACE-REPAIR REPLAY\n";
  return 0;
}
