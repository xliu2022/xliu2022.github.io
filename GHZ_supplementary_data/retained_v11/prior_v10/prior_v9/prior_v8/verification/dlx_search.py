import ctypes
from pathlib import Path
lib=ctypes.CDLL(str(Path(__file__).with_name('dlx.so').resolve()));f=lib.solve
f.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_int),ctypes.c_int,ctypes.c_uint64,ctypes.POINTER(ctypes.c_uint64)];f.restype=ctypes.c_int

def count_fast(m,S,cap=None,limit=10000000):
 a=(ctypes.c_int*len(S))(*S);o=(ctypes.c_uint64*2)();ct=f(m,len(S),a,cap or len(S)+1,limit,o)
 return dict(m=m,S=list(S),count=ct,nodes=o[0],complete=not o[1])

if __name__=='__main__':
 for m,S in [(23,[0,1,2,5,11]),(29,[0,1,2,5,17]),(31,[0,1,2,5,11]),
             (37,[0,1,2,5,21]),(41,[0,1,2,5,16]),(43,[0,1,2,5,38]),
             (47,[0,1,2,11,23])]:
  result=count_fast(m,S)
  assert result['complete'] and result['count']==5,result
  print(result)
