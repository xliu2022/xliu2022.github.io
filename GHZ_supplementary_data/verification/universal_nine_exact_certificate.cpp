// Exact integer enumeration excluding a 14-rigid 3-graph on nine vertices.
// Build: g++ -std=c++17 -O3 universal_nine_exact_certificate.cpp -o checker
// Run: checker universal_nine_missing_orbits.txt > certificate.jsonl
// The supplied orbit file covers every missing-pair set of size at most 3.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
using Mask=uint32_t;
const Mask ALL=(1u<<27)-1;
const int Q=14;
vector<Mask> transverse_pms;

void init_pms() {
  array<int,3> b={0,1,2};
  do {
    array<int,3> c={0,1,2};
    do {
      Mask P=0;
      for(int a=0;a<3;a++)P|=1u<<(9*a+3*b[a]+c[a]);
      transverse_pms.push_back(P);
    }while(next_permutation(c.begin(),c.end()));
  }while(next_permutation(b.begin(),b.end()));
  assert(transverse_pms.size()==36);
}

bool no_extra_transverse(Mask selected,int prescribed) {
  int actual=0;
  for(Mask P:transverse_pms)if((P&selected)==P)++actual;
  return actual==prescribed;
}

bool choose_pure(const vector<Mask>& candidates,int start,int need,
                 Mask selected,int chosen) {
  if(!need)return true;
  if((int)candidates.size()-start<need)return false;
  for(int i=start;i<(int)candidates.size();i++) {
    Mask P=candidates[i];
    if(P&selected)continue;
    Mask next=selected|P;
    if(!no_extra_transverse(next,chosen+1))continue;
    if(choose_pure(candidates,i+1,need-1,next,chosen+1))return true;
  }
  return false;
}

struct ForkSearch {
  array<Mask,27> adj{};
  array<array<int,27>,27> tag{};
  vector<Mask> pure;
  int needed_pure=0;
  long long leaves=0;
  bool search(Mask rem,Mask selected) {
    if(!rem) {
      ++leaves;
      return choose_pure(pure,0,needed_pure,selected,0);
    }
    int pick=-1,minimum=100;
    for(Mask z=rem;z;z&=z-1) {
      int p=__builtin_ctz(z), degree=__builtin_popcount(adj[p]&rem);
      if(!degree)return false;
      if(degree<minimum){minimum=degree;pick=p;if(degree==1)break;}
    }
    for(Mask z=adj[pick]&rem;z;z&=z-1) {
      int other=__builtin_ctz(z), t=tag[pick][other];
      Mask next=selected|(1u<<t);
      if(!no_extra_transverse(next,0))continue;
      if(search(rem&~((1u<<pick)|(1u<<other)),next))return true;
    }
    return false;
  }
};

struct Stats {
  long long flow_candidates=0, triangle_compatible=0, parity_compatible=0;
  long long residual_covered=0, enough_free_edges=0, enough_pure_pms=0;
  long long degree_compatible=0, fork_leaves=0, realizations=0;
};

Stats run(Mask missing) {
  Stats S;
  int n=27-__builtin_popcount(missing);
  int avA=511&~missing,avB=511&~(missing>>9),avC=511&~(missing>>18);
  int nAB=__builtin_popcount(avA),nBC=__builtin_popcount(avB),nCA=__builtin_popcount(avC);
  for(int B=avB;;B=(B-1)&avB) {
    int nb=__builtin_popcount(B);
    for(int C=avC;;C=(C-1)&avC) {
      int nc=__builtin_popcount(C);
      // Every vertex has hypergraph degree Q.  Summing over each base
      // block forces its outgoing arc count to equal the number of
      // undirected exchange pairs between the other two blocks.
      int requiredA=nb+nAB-nCA;
      if(requiredA>=0&&requiredA<=nAB&&nb-nc==nBC-nAB) {
        int force1=0,force0=0;
        array<array<int,2>,9> witness;
        for(int ab=0;ab<9;ab++) {
          if(!(avA>>ab&1)){witness[ab]={-1,-1};continue;}
          int a=ab/3,b=ab%3,plus=0,minus=0;
          for(int c=0;c<3;c++) {
            if(!(avB>>(b*3+c)&1)||!(avC>>(c*3+a)&1))continue;
            int bc=B>>(b*3+c)&1,ca=C>>(c*3+a)&1;
            if(bc&&ca)plus|=1<<c;
            if(!bc&&!ca)minus|=1<<c;
          }
          if(__builtin_popcount(plus)>1)force0|=1<<ab;
          if(__builtin_popcount(minus)>1)force1|=1<<ab;
          witness[ab][0]=minus?__builtin_ctz(minus):-1;
          witness[ab][1]=plus?__builtin_ctz(plus):-1;
        }
        if(!(force0&force1)) {
          int freeA=avA&~(force0|force1);
          for(int sub=freeA;;sub=(sub-1)&freeA) {
            int A=force1|sub;
            if(__builtin_popcount(A)==requiredA) {
              ++S.flow_candidates;
              Mask used=0;
              int cycles=0;
              bool overlap=false;
              for(int ab=0;ab<9;ab++) {
                if(!(avA>>ab&1))continue;
                int c=witness[ab][A>>ab&1];
                if(c<0)continue;
                int a=ab/3,b=ab%3;
                Mask triangle=(1u<<ab)|(1u<<(9+b*3+c))|(1u<<(18+c*3+a));
                if(used&triangle){overlap=true;break;}
                used|=triangle;++cycles;
              }
              if(!overlap) {
                ++S.triangle_compatible;
                Mask rem=ALL&~(missing|used);
                int residual=__builtin_popcount(rem);
                int pure_needed=Q-1-cycles-residual/2;
                if(residual%2==0&&pure_needed>=0) {
                  ++S.parity_compatible;
                  ForkSearch F;
                  F.needed_pure=pure_needed;
                  Mask free_tags=0,legal_tags=0;
                  for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int c=0;c<3;c++) {
                    int count=0,fp=-1,fq=-1;
                    bool overlaps_cycle=false;
                    auto add=[&](int p,int q,bool present) {
                      if(present&&!(missing&((1u<<p)|(1u<<q)))) {
                        ++count;fp=p;fq=q;
                        if(used&((1u<<p)|(1u<<q)))overlaps_cycle=true;
                      }
                    };
                    for(int u=0;u<3;u++)if(u!=a) {
                      int v=3-a-u,p=u*3+b,q=18+c*3+v;
                      add(p,q,!(A>>p&1)&&(C>>(c*3+v)&1));
                    }
                    for(int u=0;u<3;u++)if(u!=b) {
                      int v=3-b-u,p=a*3+u,q=9+v*3+c;
                      add(p,q,(A>>p&1)&&!(B>>(v*3+c)&1));
                    }
                    for(int u=0;u<3;u++)if(u!=c) {
                      int v=3-c-u,p=18+u*3+a,q=9+b*3+v;
                      add(p,q,!(C>>(u*3+a)&1)&&(B>>(b*3+v)&1));
                    }
                    int tag=9*a+3*b+c;
                    if(count<=1&&!overlaps_cycle) {
                      legal_tags|=1u<<tag;
                      if(!count)free_tags|=1u<<tag;
                      else {
                        F.adj[fp]|=1u<<fq;F.adj[fq]|=1u<<fp;
                        F.tag[fp][fq]=F.tag[fq][fp]=tag;
                      }
                    }
                  }
                  bool covered=true;
                  for(Mask z=rem;z;z&=z-1)if(!F.adj[__builtin_ctz(z)]){covered=false;break;}
                  if(covered) {
                    ++S.residual_covered;
                    if(__builtin_popcount(free_tags)>=3*pure_needed) {
                      ++S.enough_free_edges;
                      for(Mask P:transverse_pms)if((P&free_tags)==P)F.pure.push_back(P);
                      if((int)F.pure.size()>=pure_needed) {
                        ++S.enough_pure_pms;
                        // Independent necessary regular-degree check.
                        array<int,9> d{},out{},legal_degree{};
                        for(int p=0;p<27;p++)if(!(missing>>p&1)) {
                          int u,v,forward;
                          if(p<9){u=p/3;v=3+p%3;forward=A>>p&1;}
                          else if(p<18){int j=p-9;u=3+j/3;v=6+j%3;forward=B>>j&1;}
                          else {int j=p-18;u=6+j/3;v=j%3;forward=C>>j&1;}
                          ++d[u];++d[v];++out[forward?u:v];
                        }
                        for(Mask z=legal_tags;z;z&=z-1) {
                          int tag=__builtin_ctz(z),a=tag/9,b=tag/3%3,c=tag%3;
                          ++legal_degree[a];++legal_degree[3+b];++legal_degree[6+c];
                        }
                        bool degrees=true;
                        const array<int,3> block_out={nBC,nCA,nAB};
                        for(int v=0;v<9;v++) {
                          int exchange_degree=block_out[v/3]+d[v]-2*out[v];
                          int needed=Q-1-exchange_degree;
                          if(needed<0||needed>legal_degree[v])degrees=false;
                        }
                        if(degrees) {
                          ++S.degree_compatible;
                          if(F.search(rem,0)) {
                            ++S.realizations;
                            cerr<<"REALIZATION CANDIDATE: missing="<<missing<<" A="<<A<<" B="<<B<<" C="<<C<<"\n";
                          }
                          S.fork_leaves+=F.leaves;
                        }
                      }
                    }
                  }
                }
              }
            }
            if(!sub)break;
          }
        }
      }
      if(!C)break;
    }
    if(!B)break;
  }
  return S;
}

int main(int argc,char**argv) {
  if(argc!=2){cerr<<"usage: checker missing_orbits.txt\n";return 2;}
  init_pms();
  ifstream input(argv[1]);
  if(!input){cerr<<"cannot open orbit file\n";return 2;}
  Mask missing;int orbit_size;long long total=0;int cases=0;
  while(input>>missing>>orbit_size) {
    auto begin=chrono::steady_clock::now();
    Stats S=run(missing);total+=S.realizations;++cases;
    double seconds=chrono::duration<double>(chrono::steady_clock::now()-begin).count();
    cout<<"{\"missing_mask\":"<<missing<<",\"missing_count\":"<<__builtin_popcount(missing)
        <<",\"orbit_size\":"<<orbit_size<<",\"flow_candidates\":"<<S.flow_candidates
        <<",\"triangle_compatible\":"<<S.triangle_compatible
        <<",\"parity_compatible\":"<<S.parity_compatible
        <<",\"residual_covered\":"<<S.residual_covered
        <<",\"enough_free_edges\":"<<S.enough_free_edges
        <<",\"enough_pure_pms\":"<<S.enough_pure_pms
        <<",\"degree_compatible\":"<<S.degree_compatible
        <<",\"fork_leaves\":"<<S.fork_leaves
        <<",\"realizations\":"<<S.realizations
        <<",\"seconds\":"<<seconds<<"}"<<endl;
  }
  cout<<"{\"summary\":true,\"cases\":"<<cases<<",\"q\":"<<Q
      <<",\"realizations\":"<<total<<"}"<<endl;
  return total?1:0;
}
