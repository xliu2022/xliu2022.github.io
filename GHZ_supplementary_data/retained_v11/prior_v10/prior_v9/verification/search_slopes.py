import ctypes,time,random,json
from pathlib import Path
lib=ctypes.CDLL(str(Path(__file__).with_name('dlx_slope.so').resolve()));f=lib.solve_slope
f.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_int,ctypes.c_uint64,ctypes.POINTER(ctypes.c_uint64)];f.restype=ctypes.c_int

def count(m,S,q,limit=100000):
 a=(ctypes.c_int*len(S))(*S);o=(ctypes.c_uint64*2)();ct=f(m,len(S),a,q,len(S)+1,limit,o)
 return dict(m=m,S=list(S),q=q,count=ct,nodes=o[0],complete=not o[1])

def slopeclasses(p):
 seen=set();out=[]
 for q in range(1,p-1):
  if q in seen:continue
  orb={q,(-q-1)%p,pow(q,-1,p),(-pow(q,-1,p)-1)%p,pow((-q-1)%p,-1,p),(-pow((-q-1)%p,-1,p)-1)%p}
  seen|=orb;out.append(q)
 return out

if __name__=='__main__':
 random.seed(222)
 R=[];best=[];start=time.time()
 for m in [13,17,19,23,29,31,37,41,43,47,53,59,61,67,71]:
  qs=slopeclasses(m);st=time.time();calls=0;timeouts=0
  while time.time()-st<20:
   q=random.choice(qs);S=[0,1]+sorted(random.sample(range(2,m),4))
   r=count(m,S,q);calls+=1
   if r['complete'] and r['count']==6:
    print('FOUND12',r,flush=True);R.append(r);Path(__file__).with_name('twelve_seeds.json').write_text(json.dumps(R,indent=2));break
   if not r['complete'] and r['count']<7:
    timeouts+=1
    if r['count']<=6:best.append(r)
  print('END',m,qs,'calls',calls,'timeouts',timeouts,'seconds',time.time()-st,flush=True)
  Path(__file__).with_name('slope_incomplete.json').write_text(json.dumps(best,indent=2))
  if R:break
