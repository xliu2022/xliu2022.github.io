#!/usr/bin/env python3
"""Independent complete exact-cover checks of small cyclic lifts.

Compile once with:
  g++ -O3 -std=c++17 -shared -fPIC dlx_generic.cpp -o dlx_generic.so
Then run this file. Large-output rigidity is proved analytically.
"""
from pathlib import Path
import ctypes,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from cyclic_modular_lift import construct
ROOT=Path(__file__).resolve().parent
lib=ctypes.CDLL(str(ROOT/'dlx_generic.so'))
solve=lib.count_generic
solve.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int),
               ctypes.c_int,ctypes.c_uint64,ctypes.POINTER(ctypes.c_uint64)]
solve.restype=ctypes.c_int

def exact_cpp(M):
    m=len(M[0]);t=len(M);E=sorted(set(e for row in M for e in row))
    aa=(ctypes.c_int*(3*len(E)))(*(e[j]+m*j for e in E for j in range(3)))
    out=(ctypes.c_uint64*2)()
    count=solve(3*m,len(E),aa,t+1,10000000,out)
    assert not out[1] and count==t, 'Incomplete search or extra perfect matching'
    return {'m':m,'states':t,'perfect_matchings':count,'search_nodes':out[0],'complete':True}

if __name__=='__main__':
    records=[]
    for m,lam,S,L in [(3,1,[0,1,2],2),(5,1,[0,1,2,3],2),
                      (5,1,[0,1,2,3],3),(5,1,[0,1,2,3],4)]:
        record={'base_m':m,'sheets':L,**exact_cpp(construct(m,lam,S,L))}
        records.append(record);print(json.dumps(record),flush=True)
    (ROOT/'cyclic_lift_cpp_complete_checks.json').write_text(json.dumps(records,indent=2)+'\n')
