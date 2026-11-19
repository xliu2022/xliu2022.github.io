"""Small exact polynomial implementation; all coefficients are Fractions."""
from fractions import Fraction as Q
from math import prod

class Poly:
    def __init__(self, terms=None, dim=3):
        self.dim = dim
        self.terms = {e: Q(c) for e, c in (terms or {}).items() if c}
    @classmethod
    def constant(cls, c, dim=3):
        return cls({(0,)*dim:Q(c)},dim)
    @classmethod
    def variable(cls, i, dim=3):
        e=[0]*dim; e[i]=1
        return cls({tuple(e):Q(1)},dim)
    def _coerce(self, other):
        return other if isinstance(other,Poly) else Poly.constant(other,self.dim)
    def __add__(self,other):
        other=self._coerce(other); assert self.dim==other.dim
        terms=self.terms.copy()
        for e,c in other.terms.items(): terms[e]=terms.get(e,Q(0))+c
        return Poly(terms,self.dim)
    __radd__=__add__
    def __neg__(self): return Poly({e:-c for e,c in self.terms.items()},self.dim)
    def __sub__(self,other): return self+-self._coerce(other)
    def __rsub__(self,other): return self._coerce(other)+-self
    def __mul__(self,other):
        other=self._coerce(other); assert self.dim==other.dim
        terms={}
        for e,c in self.terms.items():
            for f,d in other.terms.items():
                g=tuple(a+b for a,b in zip(e,f)); terms[g]=terms.get(g,Q(0))+c*d
        return Poly(terms,self.dim)
    __rmul__=__mul__
    def __truediv__(self,other): return self* (Q(1)/Q(other))
    def __pow__(self,n):
        assert isinstance(n,int) and n>=0
        out=Poly.constant(1,self.dim)
        for _ in range(n): out=out*self
        return out
    def derivative(self,i):
        terms={}
        for e,c in self.terms.items():
            if e[i]:
                f=list(e); f[i]-=1; terms[tuple(f)]=c*e[i]
        return Poly(terms,self.dim)
    def at(self,point):
        assert len(point)==self.dim
        return sum((c*prod(Q(x)**k for x,k in zip(point,e)) for e,c in self.terms.items()),Q(0))

def book_polynomial(D,q1,q2):
    a,s,d=[Poly.variable(i) for i in range(3)]
    z=1-a-s; u=(s+d)/2; v=(s-d)/2
    return 3*a*((1-a)**2-(s*s+d*d)/2-z*z/D)+Q(3,2)*(s*s-d*d)*z+q1*u**3+q2*v**3+Q(4,9)*z**3/(D*D)

def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
            -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
            +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
