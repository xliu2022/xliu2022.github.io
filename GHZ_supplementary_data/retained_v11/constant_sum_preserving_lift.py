"""Odd-sheet lift preserving constant sums and the number of perfect matchings.

Input: k disjoint rigid constant-sum transversals on cyclic Z_m,
represented by permutations P_i with edges (a,P_i[a],-a-P_i[a]).
The theorem certifies rigidity. Only the structural conditions are checked here.
"""
import math

def crt(a,x,m,L):
    if L==1:return a
    return a+m*(((x-a)*pow(m,-1,L))%L)

def construct(m,permutations,L=3,mark=0):
    assert m%2 and L%2 and L>=1 and math.gcd(m,L)==1
    assert 0<=mark<m
    assert 1<=len(permutations)<=m
    seen=set()
    for P in permutations:
        assert sorted(P)==list(range(m))
        assert sorted((a+P[a])%m for a in range(m))==list(range(m))
        for a,b in enumerate(P):
            assert(a,b)not in seen;seen.add((a,b))
    Q=[]
    for i,P in enumerate(permutations):
        row=[None]*(L*m)
        for a,b in enumerate(P):
            alpha=int(a==mark and i>0)
            for x in range(L):row[crt(a,x,m,L)]=crt(b,(x+alpha)%L,m,L)
        assert sorted(row)==list(range(L*m))
        assert sorted((a+row[a])%(L*m)for a in range(L*m))==list(range(L*m))
        Q.append(row)
    assert len({(a,b)for row in Q for a,b in enumerate(row)})==len(Q)*L*m
    return Q
