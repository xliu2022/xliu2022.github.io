#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <set>
#include <vector>
using namespace std;
const array<array<int,2>,36> splits={{
{7,4},
{56,0},
{448,0},
{131,128},
{292,32},
{88,80},
{259,0},
{52,0},
{200,0},
{13,12},
{322,0},
{176,0},
{69,4},
{42,32},
{400,0},
{25,8},
{134,0},
{352,256},
{73,8},
{290,256},
{148,128},
{265,8},
{162,32},
{84,80},
{273,256},
{138,0},
{100,0},
{97,0},
{22,4},
{392,128},
{193,128},
{274,0},
{44,8},
{161,0},
{82,0},
{268,8}
}};
vector<uint32_t> edges;
array<vector<int>,18> inc;
vector<vector<int>> solutions;
uint64_t nodes=0;
void visit(uint32_t used, vector<int>& chosen){
 ++nodes;
 if(used==((1u<<18)-1)){solutions.push_back(chosen);return;}
 vector<int> best;
 size_t minsize=100;
 for(int v=0;v<18;v++)if(!(used>>v&1)){
  vector<int> choices;
  for(int i:inc[v])if(!(used&edges[i]))choices.push_back(i);
  if(choices.empty())return;
  if(choices.size()<minsize){best=choices;minsize=choices.size();}
 }
 for(int i:best){chosen.push_back(i);visit(used|edges[i],chosen);chosen.pop_back();}
}
int main(){
 for(int c=0;c<12;c++)for(int j=0;j<2;j++)for(int k=0;k<3;k++){
  auto [e,o]=splits[3*c+k];uint32_t edge=0;
  for(int v=0;v<9;v++)if(e>>v&1)edge|=1u<<(9*(j^((o>>v)&1))+v);
  edges.push_back(edge);
 }
 assert(set<uint32_t>(edges.begin(),edges.end()).size()==72);
 for(int c=0;c<12;c++){
  uint32_t used=0;for(int k=0;k<6;k++){assert(!(used&edges[6*c+k]));used|=edges[6*c+k];}
  assert(used==((1u<<18)-1));
 }
 for(int i=0;i<72;i++)for(int v=0;v<18;v++)if(edges[i]>>v&1)inc[v].push_back(i);
 for(auto&v:inc)assert(v.size()==12);
 vector<int> chosen;visit(0,chosen);
 assert(solutions.size()==12);
 set<vector<int>> actual,expected;
 for(auto s:solutions){sort(s.begin(),s.end());actual.insert(s);}
 for(int c=0;c<12;c++){vector<int>s;for(int k=0;k<6;k++)s.push_back(6*c+k);expected.insert(s);}
 assert(actual==expected);
 cout<<"{\"vertices\":18,\"edges\":72,\"perfect_matchings\":12,\"nodes\":"<<nodes<<",\"complete\":true}\n";
}
