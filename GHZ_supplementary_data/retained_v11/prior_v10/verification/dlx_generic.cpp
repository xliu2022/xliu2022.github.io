#include <vector>
#include <cstdint>
#include <algorithm>
struct DLX {
 int n,cap,total=0; uint64_t nodes=0,limit; bool stopped=false;
 std::vector<int>L,R,U,D,C,S;
 DLX(int nv,int ne,const int*ee,int cc,uint64_t lim):n(nv),cap(cc),limit(lim){
 int sz=3*ne+n+1;L.resize(sz);R.resize(sz);U.resize(sz);D.resize(sz);C.resize(sz);S.assign(n+1,0);
 for(int j=0;j<=n;j++){L[j]=j-1;R[j]=j+1;U[j]=D[j]=j;C[j]=j;}L[0]=n;R[n]=0;
 int z=n+1;
 for(int e=0;e<ne;e++){
  for(int t=0;t<3;t++){int c=1+ee[3*e+t];C[z+t]=c;U[z+t]=U[c];D[z+t]=c;D[U[c]]=z+t;U[c]=z+t;S[c]++;L[z+t]=z+(t+2)%3;R[z+t]=z+(t+1)%3;}z+=3;
 }
 }
 void cover(int c){L[R[c]]=L[c];R[L[c]]=R[c];for(int i=D[c];i!=c;i=D[i])for(int j=R[i];j!=i;j=R[j]){U[D[j]]=U[j];D[U[j]]=D[j];S[C[j]]--;}}
 void uncover(int c){for(int i=U[c];i!=c;i=U[i])for(int j=L[i];j!=i;j=L[j]){S[C[j]]++;U[D[j]]=j;D[U[j]]=j;}L[R[c]]=c;R[L[c]]=c;}
 void dfs(){if(++nodes>limit){stopped=true;return;}if(!R[0]){total++;if(total>=cap)stopped=true;return;}int c=R[0];for(int j=R[c];j;j=R[j])if(S[j]<S[c]){c=j;if(!S[c])break;}if(!S[c])return;cover(c);for(int r=D[c];r!=c&&!stopped;r=D[r]){for(int j=R[r];j!=r;j=R[j])cover(C[j]);dfs();for(int j=L[r];j!=r;j=L[j])uncover(C[j]);}uncover(c);}
};
extern "C" int count_generic(int nv,int ne,const int*ee,int cap,uint64_t limit,uint64_t*out){DLX x(nv,ne,ee,cap,limit);x.dfs();out[0]=x.nodes;out[1]=x.stopped;return x.total;}
