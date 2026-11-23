from time import monotonic

def count_shifts(m,S,cap=None,seconds=None):
 S=sorted(set(S));n=len(S)*m
 V=(1<<(3*m))-1;E=(1<<n)-1
 edges=[];masks=[];inc=[0]*(3*m)
 for a in range(m):
  for s in S:
   e=(a,m+(a+s)%m,2*m+(2*a+s)%m)
   j=len(edges);edges.append((a,s));mask=sum(1<<v for v in e);masks.append(mask)
   for v in e:inc[v]|=1<<j
 conflict=[0]*n
 for j,mask in enumerate(masks):
  rem=mask
  while rem:
   bit=rem&-rem;rem-=bit;conflict[j]|=inc[bit.bit_length()-1]
 total=0;nodes=0;solutions=[];start=monotonic();timeout=False
 def rec(V,E,chosen):
  nonlocal total,nodes,timeout
  nodes+=1
  if nodes%16384==0 and seconds is not None and monotonic()-start>seconds:timeout=True;return True
  if V==0:
   total+=1
   solutions.append([edges[j]for j in chosen]);return cap is not None and total>=cap
  rem=V;bc=n+1;best=0
  while rem:
   bit=rem&-rem;rem-=bit
   cand=inc[bit.bit_length()-1]&E;ct=cand.bit_count()
   if ct<bc:bc=ct;best=cand
   if bc<2:break
  while best:
   bit=best&-best;best-=bit;j=bit.bit_length()-1
   if rec(V^masks[j],E&~conflict[j],chosen+[j]):return True
  return False
 stopped=rec(V,E,[])
 return dict(m=m,S=S,count=total,complete=not stopped,nodes=nodes,seconds=monotonic()-start,timeout=timeout,solutions=solutions)

if __name__=='__main__':
 import json
 for m,S in [(23,[0,1,2,5,11]),(29,[0,1,2,5,17]),(29,[0,1,3,7,20])]:
  r=count_shifts(m,S,6);print(json.dumps(r),flush=True)
