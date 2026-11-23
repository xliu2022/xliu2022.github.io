#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <set>
#include <vector>
using namespace std;
const array<array<int,2>,40> splits={{
{2084,2080},
{522,8},
{208,128},
{1281,256},
{13,0},
{642,640},
{1312,1024},
{2128,2048},
{22,4},
{1096,0},
{161,32},
{2816,512},
{37,32},
{2184,2048},
{1600,0},
{274,256},
{772,512},
{1042,0},
{2113,2112},
{168,0},
{560,32},
{3200,1024},
{321,64},
{14,0},
{1544,0},
{292,0},
{2114,2048},
{145,16},
{392,384},
{112,0},
{2054,2048},
{1537,1536},
{1282,1024},
{580,576},
{2177,0},
{56,0},
{265,8},
{1092,1088},
{162,160},
{2576,0}
}};
vector<uint32_t> edges;
array<vector<int>,24> inc;
vector<vector<int>> solutions;
uint64_t nodes=0;
void visit(uint32_t used, vector<int>& chosen){
 ++nodes;
 if(used==((1u<<24)-1)){solutions.push_back(chosen);return;}
 vector<int> best;
 size_t minsize=100;
 for(int v=0;v<24;v++)if(!(used>>v&1)){
  vector<int> choices;
  for(int i:inc[v])if(!(used&edges[i]))choices.push_back(i);
  if(choices.empty())return;
  if(choices.size()<minsize){best=choices;minsize=choices.size();}
 }
 for(int i:best){chosen.push_back(i);visit(used|edges[i],chosen);chosen.pop_back();}
}
int main(){
 for(int c=0;c<10;c++)for(int j=0;j<2;j++)for(int k=0;k<4;k++){
  auto [e,o]=splits[4*c+k];uint32_t edge=0;
  for(int v=0;v<12;v++)if(e>>v&1)edge|=1u<<(12*(j^((o>>v)&1))+v);
  edges.push_back(edge);
 }
 assert(set<uint32_t>(edges.begin(),edges.end()).size()==80);
 for(int c=0;c<10;c++){
  uint32_t used=0;for(int k=0;k<8;k++){assert(!(used&edges[8*c+k]));used|=edges[8*c+k];}
  assert(used==((1u<<24)-1));
 }
 for(int i=0;i<80;i++)for(int v=0;v<24;v++)if(edges[i]>>v&1)inc[v].push_back(i);
 for(auto&v:inc)assert(v.size()==10);
 vector<int> chosen;visit(0,chosen);
 assert(solutions.size()==10);
 set<vector<int>> actual,expected;
 for(auto s:solutions){sort(s.begin(),s.end());actual.insert(s);}
 for(int c=0;c<10;c++){vector<int>s;for(int k=0;k<8;k++)s.push_back(8*c+k);expected.insert(s);}
 assert(actual==expected);
 cout<<"{\"vertices\":24,\"edges\":80,\"perfect_matchings\":10,\"nodes\":"<<nodes<<",\"complete\":true}\n";
}
