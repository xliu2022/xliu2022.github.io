// Independent exact check of the finite exclusion in universal_nine_exact_certificate.cpp.
// This implementation directly enumerates all balanced direction matrices and
// derives fork completions from disjoint 9-bit hyperedges, not coordinate formulas.
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
using U=uint64_t;
const unsigned ALL=(1u<<27)-1;
array<pair<int,int>,27> pair_vertices;
array<unsigned,54> exchange_edges;
vector<U> triangles;
array<vector<U>,27> forks;

void initialize() {
  for(int a=0;a<3;a++)for(int b=0;b<3;b++)pair_vertices[3*a+b]={a,3+b};
  for(int b=0;b<3;b++)for(int c=0;c<3;c++)pair_vertices[9+3*b+c]={3+b,6+c};
  for(int c=0;c<3;c++)for(int a=0;a<3;a++)pair_vertices[18+3*c+a]={6+c,a};
  for(int p=0;p<27;p++) {
    auto [u,v]=pair_vertices[p];
    exchange_edges[p]=((7u<<(3*(u/3)))&~(1u<<u))|(1u<<v);
    exchange_edges[p+27]=((7u<<(3*(v/3)))&~(1u<<v))|(1u<<u);
  }
  for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int c=0;c<3;c++) {
    U cycle=(1ULL<<(3*a+b))|(1ULL<<(9+3*b+c))|(1ULL<<(18+3*c+a));
    triangles.push_back(cycle);triangles.push_back(cycle<<27);
  }
  for(int p=0;p<54;p++)for(int q=p+1;q<54;q++) {
    unsigned e=exchange_edges[p],f=exchange_edges[q];
    if(e&f)continue;
    unsigned left=511u&~(e|f);
    if(__builtin_popcount(left&7)!=1||__builtin_popcount((left>>3)&7)!=1
       ||__builtin_popcount((left>>6)&7)!=1)continue;
    int a=__builtin_ctz(left&7),b=__builtin_ctz((left>>3)&7),c=__builtin_ctz((left>>6)&7);
    forks[9*a+3*b+c].push_back((1ULL<<p)|(1ULL<<q));
  }
  for(auto const& list:forks)assert(list.size()==6);
}

void run(unsigned missing,int orbit) {
  auto begin=chrono::steady_clock::now();
  unsigned available=ALL&~missing;
  int n=__builtin_popcount(available);
  array<array<vector<int>,10>,3> matrices;
  array<int,3> sizes;
  for(int i=0;i<3;i++) {
    int allowed=(available>>(9*i))&511;
    sizes[i]=__builtin_popcount(allowed);
    for(int s=allowed;;s=(s-1)&allowed) {
      matrices[i][__builtin_popcount(s)].push_back(s);
      if(!s)break;
    }
  }
  long long balanced=0,compatible=0,parity=0,covered=0,free_passed=0;
  for(int na=0;na<=sizes[0];na++) {
    int nb=na+sizes[2]-sizes[0],nc=na+sizes[2]-sizes[1];
    if(nb<0||nb>sizes[1]||nc<0||nc>sizes[2])continue;
    for(int A:matrices[0][na])for(int B:matrices[1][nb])for(int C:matrices[2][nc]) {
      ++balanced;
      unsigned forward=A|(B<<9)|(C<<18);
      U arcs=forward|((U)(available&~forward)<<27);
      U cycle_arcs=0;
      int cycles=0;
      bool valid=true;
      for(U triangle:triangles)if((triangle&arcs)==triangle) {
        if(triangle&cycle_arcs){valid=false;break;}
        cycle_arcs|=triangle;++cycles;
      }
      if(!valid)continue;
      ++compatible;
      U residual=arcs&~cycle_arcs;
      int count=__builtin_popcountll(residual);
      int needed_pure=13-cycles-count/2;
      if(count%2||needed_pure<0)continue;
      ++parity;
      U covered_arcs=0;
      int free_count=0;
      for(auto const& list:forks) {
        int ways=0;U unique=0;
        for(U fork:list)if((fork&arcs)==fork){++ways;unique=fork;}
        if(!ways)++free_count;
        if(ways==1&&!(unique&cycle_arcs))covered_arcs|=unique;
      }
      if(covered_arcs!=residual)continue;
      ++covered;
      if(free_count>=3*needed_pure)++free_passed;
    }
  }
  double seconds=chrono::duration<double>(chrono::steady_clock::now()-begin).count();
  cout<<"{\"missing_mask\":"<<missing<<",\"orbit_size\":"<<orbit
      <<",\"all_balanced_orientations\":"<<balanced
      <<",\"triangle_compatible\":"<<compatible<<",\"parity_compatible\":"<<parity
      <<",\"residual_covered\":"<<covered<<",\"enough_free_edges\":"<<free_passed
      <<",\"seconds\":"<<seconds<<"}"<<endl;
}

int main(int argc,char**argv) {
  if(argc!=2)return 2;
  initialize();
  ifstream input(argv[1]);if(!input)return 2;
  unsigned missing;int orbit;
  while(input>>missing>>orbit)run(missing,orbit);
}
