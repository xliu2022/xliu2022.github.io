from fractions import Fraction as F

if not __debug__:
    raise RuntimeError("Run this verifier without Python's -O flag.")

import sympy as sp

# ---------------------------------------------------------------------------
# I. A 42-row exact heterogeneous-satellite chain from 2/25 to 36/251.
# None denotes s_i = infinity, hence alpha_i = 1.
# ---------------------------------------------------------------------------
ROWS = [
((4,8,5,15,4),((0,3),(0,4),(1,2))),
((3,4,200,20,9),((0,1),(1,4),(2,3))),
((100,None,3,8,5),((0,2),(1,4),(2,3))),
((10,30,9,50,3),((0,3),(1,4),(2,4))),
((30,7,500,3,1000),((0,3),(1,3),(2,4))),
((200,8,3,1000,200),((0,3),(0,4),(1,4),(2,3))),
((1000,None,100,3,12),((0,1),(0,2),(1,3),(2,4))),
((200,3,200,200,30),((0,2),(0,3),(1,3),(2,4))),
((50,500,4,8,None),((0,1),(0,4),(1,3),(2,4))),
((12,1000,50,None,4),((0,1),(1,2),(2,3),(3,4))),
((None,500,4,20,500),((0,1),(0,3),(1,4),(2,4))),
((4,2,None,None),((0,2),(1,3),(2,3))),
((5,50,None,200,1000),((0,3),(1,2),(2,4),(3,4))),
((8,500,15,500,None),((0,3),(1,2),(1,4),(3,4))),
((1000,200,1000,8,100),((0,1),(0,4),(1,2),(2,3))),
((500,12,1000,500,200),((0,2),(0,4),(1,3),(2,3))),
((500,None,500,1000,30),((0,2),(1,3),(1,4),(2,3))),
((2,None,None,None),((0,1),(1,3),(2,3))),
((4,3,9,None),((0,2),(1,3),(2,3))),
((3,15,4,100),((0,3),(1,2),(1,3))),
((50,4,3,500),((0,1),(0,3),(2,3))),
((20,5,None,3),((0,2),(0,3),(1,2))),
((3,6,30,100),((0,3),(1,2),(2,3))),
((200,4,100,4),((0,2),(0,3),(1,2))),
((3,10,1000,1000),((0,3),(1,2),(2,3))),
((None,3,100,30),((0,1),(0,2),(2,3))),
((None,4,None,8),((0,2),(0,3),(1,2))),
((None,15,4,None),((0,2),(0,3),(1,3))),
((1000,6,200,9),((0,1),(0,2),(2,3))),
((1000,1000,20,6),((0,1),(0,2),(1,3))),
((None,100,7,None),((0,2),(0,3),(1,3))),
((500,500,12,500),((0,3),(1,2),(1,3))),
((200,100,500,1000),((0,1),(0,3),(1,2))),
((4,3,8),((0,2),(1,2))),
((12,4,3),((0,1),(0,2))),
((4,30,3),((0,1),(1,2))),
((3,20,5),((0,1),(1,2))),
((3,15,7),((0,1),(1,2))),
((12,5,4),((0,1),(0,2))),
((3,12,50),((0,1),(1,2))),
((9,4,12),((0,2),(1,2))),
((10,5,9),((0,2),(1,2))),
]

def satellite_data(svals, edges):
    alpha = [F(1) if s is None else F(s, s-1) for s in svals]
    S = sum(alpha, F(0))
    K = sum((alpha[i]*alpha[j] for i,j in edges), F(0))
    L = max(
        sum((alpha[j] for i0,j in edges if i0 == i), F(0))
        + sum((alpha[i0] for i0,j in edges if j == i), F(0))
        for i in range(len(alpha))
    )
    rho = 1 - K/S
    T = F(1, 2*S)
    kappa = 4*(1-rho)/(2-rho)**2
    B = T*kappa
    xstar = 4*rho*T
    return S,K,L,rho,T,B,xstar

M = F(2,25)
for row_id,(svals,edges) in enumerate(ROWS,1):
    m = len(svals)
    deg = [0]*m
    eset = {tuple(sorted(e)) for e in edges}
    for i,j in eset:
        deg[i] += 1; deg[j] += 1
    assert min(deg) >= 1 and max(deg) <= 2
    assert all(not ({tuple(sorted((i,j))),tuple(sorted((i,k))),tuple(sorted((j,k)))} <= eset)
               for i in range(m) for j in range(i+1,m) for k in range(j+1,m))
    S,K,L,rho,T,B,xstar = satellite_data(svals,edges)
    assert 0 < rho <= F(1,4)
    assert S <= F(31,5)
    assert L <= 3
    assert B <= M and xstar <= M < T
    M = T
assert M == F(36,251)
print(f'I. verified {len(ROWS)} exact satellite rows; endpoint = {M}')

# The uniform analytic box used for every row.
rho,L,S = sp.symbols('rho L S', nonnegative=True)
A0 = -16*L*rho**2 + 16*L*rho + S*rho**3 - 48*rho**2 + 112*rho - 64
B0 = L*rho + 7*rho - 8
R1 = (-16*L**2*rho**3 + 16*L**2*rho**2 + L*S*rho**4
      -128*L*rho**3 + 320*L*rho**2 - 192*L*rho
      +5*S*rho**4 - 8*S*rho**3 -208*rho**3 +784*rho**2 -960*rho +384)
R2 = (-4*L**2*rho**3 + 4*L**2*rho**2
      -56*L*rho**3 +152*L*rho**2 -96*L*rho
      +S*rho**4 -2*S*rho**3 -100*rho**3 +388*rho**2 -480*rho +192)
# Exact worst-corner values and derivative certificates.
assert sp.factor(R1.subs({L:3,S:sp.Rational(31,5),rho:sp.Rational(1,4)})) == sp.Rational(16947,160)
assert sp.factor(R2.subs({L:3,S:sp.Rational(31,5),rho:sp.Rational(1,4)})) == sp.Rational(64103,1280)
assert sp.simplify(sp.diff(R1,rho).subs({L:3,S:sp.Rational(31,5)}) - sp.Rational(8,5)*(124*rho**3-1473*rho**2+2360*rho-960)) == 0
assert sp.simplify(sp.diff(R2,rho).subs({L:3,S:sp.Rational(31,5)}) - sp.Rational(2,5)*(62*rho**3-2373*rho**2+4400*rho-1920)) == 0
# Rational inequalities used in the analytic-box proof.
assert F(9) + F(31,320) + F(25) - F(64) < 0
assert F(10,4) - 8 < 0
assert F(18) + F(31,320) + F(80) - F(192) < 0
assert F(9,2) + F(38) - F(96) < 0
assert F(124,64) + F(590) - F(960) < 0
assert F(62,64) + F(1100) - F(1920) < 0
print('   verified analytic-box corner identities and monotonicity bounds')

# ---------------------------------------------------------------------------
# II. Analytic P3 tail to 5/32.
# Parameters are (n,6,n), with the s=6 vertex at the center of P3.
# ---------------------------------------------------------------------------
n = sp.symbols('n', integer=True, positive=True)
Tn = sp.Rational(5)*(n-1)/(4*(8*n-3))
Bn = 30*n*(n-1)/(14*n-3)**2
Xn = 5*(n-1)*(2*n-3)/(8*n-3)**2
Tprev = 5*(n-2)/(4*(8*n-11))
assert Bn.subs(n,8) <= sp.Rational(36,251)
assert Xn.subs(n,8) <= sp.Rational(36,251) < Tn.subs(n,8)
assert sp.simplify((Tprev-Xn) - 5*(72*n**2-211*n+114)/(4*(8*n-11)*(8*n-3)**2)) == 0
assert sp.simplify((Tprev-Bn) - 5*(4*n**3-20*n**2-87*n-18)/(4*(8*n-11)*(14*n-3)**2)) == 0
assert (72*n**2-211*n+114).subs(n,8) > 0
assert (4*n**3-20*n**2-87*n-18).subs(n,8) > 0
assert sp.diff(4*n**3-20*n**2-87*n-18,n).subs(n,8) > 0
assert sp.limit(Tn,n,sp.oo) == sp.Rational(5,32)
print('II. verified the P3 reciprocal tail and limit 5/32')

# ---------------------------------------------------------------------------
# III. Exact overlap certificates for the ordinary-3-link compiler.
# p_q(t)=0 is the unique fixed point.  For s in (1/2,1),
# theta_{q'}(0)<s^2 iff q'<(2s-1)/(2s(1-s)).
# ---------------------------------------------------------------------------
t,q = sp.symbols('t q')
p = t**3 + (1+4*q)*t**2 + (3-4*q)*t - 1
OVERLAPS = [
    (F(2,9),   F(8,27),  F(7,12)),
    (F(8,27),  F(3,8),   F(3,5)),
    (F(3,8),   F(12,25), F(8,13)),
    (F(12,25), F(5,9),   F(7,11)),
    (F(5,9),   F(21,32), F(13,20)),
    (F(21,32), F(55,72), F(43,64)),
    (F(55,72), F(552,625), F(20,29)),
]
for q0,q1,s in OVERLAPS:
    threshold = (2*s-1)/(2*s*(1-s))
    test = s*s
    pval = F(test**3) + (1+4*q0)*test**2 + (3-4*q0)*test - 1
    assert q1 < threshold
    assert pval < 0
assert p.subs({q:sp.Rational(2,9),t:sp.Rational(1,3)}) < 0
assert p.subs({q:sp.Rational(552,625),t:sp.Rational(1,2)}) == -sp.Rational(41,5000)
print('III. verified seven exact ordinary-link overlap certificates and the final crossing of 1/2')

# Mixed-pair margin in the ordinary-link compiler.
r = sp.symbols('r', positive=True)
cubic = 2*r**3 - 4*r**2 + 7*r - 2
assert sp.discriminant(sp.diff(cubic,r),r) < 0
assert sp.simplify(cubic.subs(r,sp.sqrt(2)/2)) == -4 + 4*sp.sqrt(2)
print('   verified the uniform mixed-pair margin for all q in [0,1], t <= 1/2')

# ---------------------------------------------------------------------------
# IV. Explicit wrapper endpoints used by the rectangle ladder.
# beta * (b)_2/a^2.
# ---------------------------------------------------------------------------
ENDPOINTS = [
    (F(1,2),8,5,F(5,32)),
    (F(1,2),6,4,F(1,6)),
    (F(2,3),6,4,F(2,9)),
    (F(1,2),12,9,F(1,4)),
    (F(1,2),30,25,F(1,3)),
]
for beta,a,b,target in ENDPOINTS:
    assert beta*F(b*(b-1),a*a) == target
# Ordinary 3-density sequence: complete multipartite values plus 8/27.
assert F((3-1)*(3-2),3*3) == F(2,9)
assert F((4-1)*(4-2),4*4) == F(3,8)
assert F((5-1)*(5-2),5*5) == F(12,25)
assert F((6-1)*(6-2),6*6) == F(5,9)
assert F((8-1)*(8-2),8*8) == F(21,32)
assert F((12-1)*(12-2),12*12) == F(55,72)
assert F((25-1)*(25-2),25*25) == F(552,625)
print('IV. verified all exact transversal-wrapper endpoints and ordinary density values')

print('\nHIGH MODULE EXACT CERTIFICATES PASSED')
