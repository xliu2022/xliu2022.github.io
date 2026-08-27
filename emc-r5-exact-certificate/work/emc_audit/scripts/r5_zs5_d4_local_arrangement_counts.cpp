// Exact-count diagnostic for the five-dimensional local subset-sum arrangement.
// Uses reduced int64 rationals; all intermediate values are guarded through
// __int128.  This first-stage program reports vertices/chambers only and is not
// yet the final triangulation certificate.

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

using i64 = long long;
using i128 = __int128_t;
constexpr int N = 5;

struct Rat {
  i64 n = 0, d = 1;
  Rat() = default;
  Rat(i64 a): n(a), d(1) {}
  Rat(i64 a, i64 b) { set(a,b); }
  void set128(i128 a, i128 b) {
    if (b < 0) { a=-a; b=-b; }
    if (a > INT64_MAX || a < INT64_MIN || b > INT64_MAX)
      throw std::overflow_error("rational overflow");
    set((i64)a,(i64)b);
  }
  void set(i64 a, i64 b) {
    assert(b != 0); if (b < 0) { a=-a; b=-b; }
    i64 g=std::gcd(a<0?-a:a,b); n=a/g; d=b/g;
  }
};
Rat operator+(Rat a,Rat b){Rat q;q.set128((i128)a.n*b.d+(i128)b.n*a.d,(i128)a.d*b.d);return q;}
Rat operator-(Rat a,Rat b){Rat q;q.set128((i128)a.n*b.d-(i128)b.n*a.d,(i128)a.d*b.d);return q;}
Rat operator*(Rat a,Rat b){Rat q;q.set128((i128)a.n*b.n,(i128)a.d*b.d);return q;}
Rat operator/(Rat a,Rat b){assert(b.n);Rat q;q.set128((i128)a.n*b.d,(i128)a.d*b.n);return q;}
Rat operator-(Rat a){return Rat(-a.n,a.d);}
bool operator==(Rat a,Rat b){return a.n==b.n&&a.d==b.d;}
bool operator<(Rat a,Rat b){return (i128)a.n*b.d<(i128)b.n*a.d;}
bool operator>(Rat a,Rat b){return b<a;}
bool operator<=(Rat a,Rat b){return !(b<a);}

struct Plane { std::array<int,N> a{}; int rhs=0; int mask=0; };
using Point = std::array<Rat,N>;

static int popcount(int mask){return __builtin_popcount((unsigned)mask);}

static std::vector<Plane> all_planes(){
  std::vector<Plane> out;
  for(int mask=1;mask<(1<<N);++mask){
    int s=popcount(mask),lo=s==1?0:1,hi=s==1?1:s-1;
    for(int k=lo;k<=hi;++k){Plane p;p.mask=mask;p.rhs=k;for(int i=0;i<N;++i)p.a[i]=(mask>>i)&1;out.push_back(p);}
  }
  assert(out.size()==59);return out;
}
static std::vector<Plane> walls(){
  std::vector<Plane> out;
  for(int s=2;s<=N;++s)for(int mask=1;mask<(1<<N);++mask)if(popcount(mask)==s)
    for(int k=1;k<s;++k){Plane p;p.mask=mask;p.rhs=k;for(int i=0;i<N;++i)p.a[i]=(mask>>i)&1;out.push_back(p);}
  assert(out.size()==49);return out;
}

static bool solve(const std::array<Plane,N>& eq,Point& x){
  std::array<std::array<Rat,N+1>,N> a{};
  for(int i=0;i<N;++i){for(int j=0;j<N;++j)a[i][j]=Rat(eq[i].a[j]);a[i][N]=Rat(eq[i].rhs);}
  for(int c=0;c<N;++c){
    int p=c;while(p<N&&a[p][c].n==0)++p;if(p==N)return false;std::swap(a[p],a[c]);
    Rat pivot=a[c][c];for(int j=c;j<=N;++j)a[c][j]=a[c][j]/pivot;
    for(int i=0;i<N;++i)if(i!=c&&a[i][c].n){Rat q=a[i][c];for(int j=c;j<=N;++j)a[i][j]=a[i][j]-q*a[c][j];}
  }
  for(int i=0;i<N;++i)x[i]=a[i][N];return true;
}

static Rat value(const Plane& p,const Point& x){Rat q(-p.rhs);for(int i=0;i<N;++i)if(p.a[i])q=q+x[i];return q;}

static int rank(std::vector<std::array<int,N>> rows){
  std::vector<std::array<Rat,N>> a(rows.size());for(size_t i=0;i<rows.size();++i)for(int j=0;j<N;++j)a[i][j]=Rat(rows[i][j]);
  int r=0;for(int c=0;c<N&&r<(int)a.size();++c){int p=r;while(p<(int)a.size()&&!a[p][c].n)++p;if(p==(int)a.size())continue;std::swap(a[p],a[r]);Rat q=a[r][c];for(int j=c;j<N;++j)a[r][j]=a[r][j]/q;for(size_t i=0;i<a.size();++i)if((int)i!=r&&a[i][c].n){q=a[i][c];for(int j=c;j<N;++j)a[i][j]=a[i][j]-q*a[r][j];}++r;}return r;
}

struct Cell { std::vector<int> candidates; std::vector<int8_t> signs; };

int main(int argc,char** argv){
  // Optional second positional argument refines the global box grid.  The
  // local arrangement is independent of D because every threshold becomes
  // sum(y_i)=integer after x_i=(box_i+y_i)/D.  D must be a multiple of five
  // so that the fixed lower endpoint m=1/5 is itself a box boundary.
  int D=argc>2?std::stoi(argv[2]):5;
  if(D%5||D<5||D>10)throw std::runtime_error("denominator must be 5 or 10");
  auto planes=all_planes();std::set<Point> vertex_set;
  std::array<Plane,N> chosen;
  for(int a=0;a<(int)planes.size();++a)for(int b=a+1;b<(int)planes.size();++b)
   for(int c=b+1;c<(int)planes.size();++c)for(int d=c+1;d<(int)planes.size();++d)
    for(int e=d+1;e<(int)planes.size();++e){chosen={planes[a],planes[b],planes[c],planes[d],planes[e]};Point x;if(!solve(chosen,x))continue;bool ok=true;for(auto q:x)ok=ok&&Rat(0)<=q&&q<=Rat(1);if(ok)vertex_set.insert(x);}
  std::vector<Point> vertices(vertex_set.begin(),vertex_set.end());
  std::cout<<"vertices "<<vertices.size()<<"\n"<<std::flush;
  Cell start;start.candidates.resize(vertices.size());std::iota(start.candidates.begin(),start.candidates.end(),0);std::vector<Cell> cells{std::move(start)};
  auto hs=walls();
  for(size_t h=0;h<hs.size();++h){std::vector<Cell> next;for(auto& cell:cells){Cell neg,pos;neg.signs=pos.signs=cell.signs;neg.signs.push_back(-1);pos.signs.push_back(1);bool hn=false,hp=false;for(int id:cell.candidates){Rat q=value(hs[h],vertices[id]);if(q<=Rat(0))neg.candidates.push_back(id);if(q.n>=0)pos.candidates.push_back(id);hn=hn||(q<Rat(0));hp=hp||(q>Rat(0));}if(hn)next.push_back(std::move(neg));if(hp)next.push_back(std::move(pos));}cells=std::move(next);std::cout<<"wall "<<h+1<<" cells "<<cells.size()<<"\n"<<std::flush;}
  std::map<int,int> dist;int minv=1000,maxv=0;
  for(auto& cell:cells){std::vector<std::array<int,N>> constraints;for(int axis=0;axis<N;++axis){std::array<int,N> lo{},hi{};lo[axis]=-1;hi[axis]=1;constraints.push_back(lo);constraints.push_back(hi);}for(size_t h=0;h<hs.size();++h){auto a=hs[h].a;if(cell.signs[h]>0)for(auto& q:a)q=-q;constraints.push_back(a);}int count=0;for(int id:cell.candidates){std::vector<std::array<int,N>> active;Point x=vertices[id];for(int axis=0;axis<N;++axis){if(x[axis]==Rat(0)){std::array<int,N>a{};a[axis]=-1;active.push_back(a);}if(x[axis]==Rat(1)){std::array<int,N>a{};a[axis]=1;active.push_back(a);}}for(size_t h=0;h<hs.size();++h)if(value(hs[h],x).n==0){auto a=hs[h].a;if(cell.signs[h]>0)for(auto&q:a)q=-q;active.push_back(a);}if(rank(active)==N)++count;}dist[count]++;minv=std::min(minv,count);maxv=std::max(maxv,count);}
  std::cout<<"cells "<<cells.size()<<" vertex_distribution";for(auto [k,v]:dist)std::cout<<' '<<k<<':'<<v;std::cout<<"\n";

  // Per-global-box arrangements: each subset has at most one relevant level.
  // Count all boxes and S4 orbits of (box, chamber sign pattern).
  std::vector<std::array<int,4>> perms;
  std::array<int,4> perm{0,1,2,3};
  do { perms.push_back(perm); } while(std::next_permutation(perm.begin(),perm.end()));
  std::set<std::string> orbit_keys;
  long long total_box_cells=0;int maximum_box_cells=0,minimum_box_cells=100000;
  std::map<int,int> box_cell_distribution;
  for(int b0=D/5;b0<D;++b0)for(int x0=0;x0<D;++x0)for(int x1=0;x1<D;++x1)
   for(int x2=0;x2<D;++x2)for(int x3=0;x3<D;++x3){
    std::array<int,5> base{x0,x1,x2,x3,b0};
    std::vector<Plane> relevant;
    for(int mask=1;mask<(1<<N);++mask){
      int s=popcount(mask);if(s<2)continue;int lower_sum=0;for(int i=0;i<N;++i)if((mask>>i)&1)lower_sum+=base[i];
      int level=D-lower_sum;if(1<=level&&level<s){Plane p;p.mask=mask;p.rhs=level;for(int i=0;i<N;++i)p.a[i]=(mask>>i)&1;relevant.push_back(p);}
    }
    Cell initial;initial.candidates.resize(vertices.size());std::iota(initial.candidates.begin(),initial.candidates.end(),0);std::vector<Cell> boxcells{std::move(initial)};
    for(auto& wall:relevant){std::vector<Cell> next;for(auto& cell:boxcells){Cell neg,pos;neg.signs=pos.signs=cell.signs;neg.signs.push_back(-1);pos.signs.push_back(1);bool hn=false,hp=false;for(int id:cell.candidates){Rat q=value(wall,vertices[id]);if(q<=Rat(0))neg.candidates.push_back(id);if(q.n>=0)pos.candidates.push_back(id);hn=hn||(q<Rat(0));hp=hp||(q>Rat(0));}if(hn)next.push_back(std::move(neg));if(hp)next.push_back(std::move(pos));}boxcells=std::move(next);}
    int count=boxcells.size();total_box_cells+=count;maximum_box_cells=std::max(maximum_box_cells,count);minimum_box_cells=std::min(minimum_box_cells,count);box_cell_distribution[count]++;
    for(auto& cell:boxcells){
      std::array<int8_t,32> sign_by_mask{};for(size_t j=0;j<relevant.size();++j)sign_by_mask[relevant[j].mask]=cell.signs[j];
      std::string best;
      for(auto p:perms){
        std::array<int,4> pb{};for(int old=0;old<4;++old)pb[p[old]]=base[old];
        std::array<int8_t,32> ps{};
        for(int mask=1;mask<32;++mask)if(sign_by_mask[mask]){int pmask=mask&(1<<4);for(int old=0;old<4;++old)if((mask>>old)&1)pmask|=1<<p[old];ps[pmask]=sign_by_mask[mask];}
        std::string key;key.push_back(char('0'+b0));for(int q:pb)key.push_back(char('0'+q));for(int mask=1;mask<32;++mask)key.push_back(ps[mask]<0?'-':ps[mask]>0?'+':'0');
        if(best.empty()||key<best)best=key;
      }
      orbit_keys.insert(best);
    }
   }
  std::cout<<"denominator "<<D<<" global_boxes "<<(D-D/5)*D*D*D*D
           <<" total_relevant_cells "<<total_box_cells
           <<" s4_cell_orbits "<<orbit_keys.size()<<" min_per_box "<<minimum_box_cells
           <<" max_per_box "<<maximum_box_cells<<" distinct_box_counts "
           <<box_cell_distribution.size()<<"\n";
  if(argc>1){
    std::ofstream out(argv[1]);if(!out)throw std::runtime_error("cannot open output");
    out<<"VERTICES "<<vertices.size()<<"\n";
    for(auto& x:vertices){for(int i=0;i<N;++i){if(i)out<<' ';out<<x[i].n;if(x[i].d!=1)out<<'/'<<x[i].d;}out<<"\n";}
    out<<"ORBITS "<<orbit_keys.size()<<"\n";for(auto& key:orbit_keys)out<<key<<"\n";
    out<<"END\n";
    std::cout<<"wrote "<<argv[1]<<"\n";
  }
}
