#!/usr/bin/env python3
"""Exact rational checks for the minimum 2-degree interval proofs.

Checks the two four-uniform endpoint transfers, the six-uniform bridge,
the five-uniform Bernstein, continuity-endpoint, and polynomial identities,
and the initial overlap comparisons in higher uniformities.
Uses only the Python standard library.
"""
from fractions import Fraction as F
from math import comb, factorial


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def phi(theta, x, z):
    return (theta + 3*x*z*z)/(1+z)**3


def psi(theta, g, z):
    return (3*g*z+3*theta*z*z)/(1+z)**3


def derivative_core(theta, g, z):
    return g+2*(theta-g)*z-theta*z*z


def bernstein(theta, g, lo, hi):
    power = [3*g, 8*theta-9*g, 3*g-9*theta, 3*theta]
    shift = [sum((power[j]*comb(j,i)*lo**(j-i)*(hi-lo)**i
                  for j in range(i,4)), F(0)) for i in range(4)]
    return tuple(sum((shift[j]*F(comb(i,j),comb(3,j))
                      for j in range(i+1)), F(0)) for i in range(4))


def verify_four():
    for q,t,z,bounds,target in [
        (F(8,27),F(1,3),F(7,6),(F(548,1521),F(61,169)),F(9,25)),
        (F(3,8),F(9,25),F(5,4),(F(31,81),F(769,2025)),F(3,8)),
    ]:
        got=((1+2*q*z)/(1+z)**2,(z*z+t)/(1+z)**2)
        check(got==bounds and min(got)>target, 'four-uniform endpoint')
    check(F(27,20)**2-2*F(8,27)*F(27,20)-1==F(9,400),'root upper bound')
    check(F(729,2209)<F(1,3),'first left overlap')
    check(73<9**2,'second root upper bound')
    check(97<10**2 and F(196,529)<F(3,8),'continuous-family left overlap')
    for beta,b,a,target in [(F(1,2),4,6,F(1,6)),(F(2,3),4,6,F(2,9)),
                             (F(1,2),9,12,F(1,4)),(F(1,2),25,30,F(1,3))]:
        check(beta*b*(b-1)/a**2==target,'exact symmetric fixed point')


def verify_profiles():
    for profile,target in [((2,1,1),F(3,16)),((2,2),F(3,8))]:
        r=sum(profile)
        value=F(factorial(r),r**r)
        for a in profile:
            value*=F(a**a,factorial(a))
        check(value==target,'single-profile density')


def verify_six():
    theta=F(3,16); z=F(10,11)
    f=lambda z:(1-z)*(1+z)**3
    check(f(F(5,6))-1==F(35,1296),'left root bracket')
    check(f(F(6,7))-1==-F(204,2401),'right root bracket')
    check(theta/(3*z*(1-z))==F(121,160)<1,'auxiliary density range')
    p1=8*z*z*theta-9*z*theta+3*z+3*theta-3
    p2=17*z**3*theta-18*z*z*theta+12*z*z+3*z*theta-12*z+3*theta
    check((p1,p2)==(-F(9,1936),-F(6651,21296)),'mixed-pair endpoints')
    # Both derivative polynomials are increasing and positive for z >= 5/6.
    low=F(5,6)
    check(F(3,16)*(16*low+7)==F(61,16)>0,'first derivative lower bound')
    check(F(3,16)*(51*low*low+92*low-61)==F(613,64)>0,'second derivative lower bound')
    value=theta/f(z)
    check(value==F(14641,49392),'six-uniform bridge endpoint')
    check(value-F(5,18)==F(307,16464)>0,'interpolation overlap')


def verify_five():
    theta,g=F(8,27),F(3,32)
    first=(F(9,32),F(1805,1728),F(65,3456),F(235,1152))
    second=(F(235,1152),F(1345,3456),F(3085,1728),F(2129,288))
    check(bernstein(theta,g,F(0),F(3,2))==first,'first Bernstein row')
    check(bernstein(theta,g,F(3,2),F(3))==second,'second Bernstein row')
    check(min(first+second)==F(65,3456)>0,'positive polynomial bound')
    z=F(2,5)
    check(derivative_core(theta,g,z)==F(4501,21600)>0,'first psi monotonicity')
    check(psi(theta,g,z)<phi(theta,0,z)==F(1000,9261)<F(1,9),'first left bound')
    check(min(phi(theta,1,F(3,2)),psi(theta,g,F(3,2)))==F(31,200)>F(3,20),
          'first right bound')
    # Coefficient of theta in P; completing its quadratic factor.
    check([0,F(8),F(-9),F(3)]==[0,3*F(9,4)+F(5,4),-9,3],
          'parameter monotonicity identity')
    check(F(4,9)>F(8,27),'continuous theta range')
    theta=F(4,9);z=F(1,2)
    check(derivative_core(theta,g,z)==F(1,3)>0,'middle psi monotonicity')
    check(psi(theta,g,z)==F(91,648) and phi(theta,0,z)==F(32,243),
          'middle left values')
    check(max(psi(theta,g,z),phi(theta,0,z))<F(3,20),'middle left overlap')
    check(phi(F(1),F(3,4),F(1))==F(13,32),'middle phi endpoint')
    check(psi(F(1),g,F(1))==F(105,256),'middle psi endpoint')
    check(min(F(13,32),F(105,256))>F(3,8),'middle right overlap')
    # Exact coefficients of the sum of squares in ascending order.
    theta,g=F(1),F(3,8)
    actual=[8*3*g,8*(8*theta-9*g),8*(3*g-9*theta),8*3*theta]
    squares=[9*F(17,18)**2+F(35,36),24*F(3,2)**2-18*F(17,18),
             -48*F(3,2)+9,F(24)]
    check(actual==squares==[9,37,-63,24],'last polynomial identity')
    z=F(2,5)
    check(derivative_core(theta,g,z)==F(143,200)>0,'last psi monotonicity')
    check(psi(theta,g,z)==F(465,1372)<phi(theta,0,z)==F(125,343)<F(3,8),
          'last left overlap')
    z=F(8,9)
    check(phi(theta,1,z)==psi(theta,g,z)==F(2457,4913),'last endpoint equality')
    check(phi(theta,1,z)-F(1,2)==F(1,9826)>0,'last right overlap')


def verify_higher():
    p = 6
    ratio = F(2*(p-2)**p*(p+2)**(p-1), p**p*factorial(p+1))
    check(ratio == F(262144,229635) > 1, 'higher-uniformity initial ratio')
    initial = 3*F(3,5)**5
    start = F(factorial(6),2*7**4)
    check(initial == F(729,3125) and start == F(360,2401),
          'seven-uniform interval endpoints')
    check(initial-start == F(625329,7503125) > 0, 'seven-uniform overlap')


if __name__=='__main__':
    for function in (verify_four,verify_profiles,verify_six,verify_five,verify_higher):
        function()
        print(function.__name__+': passed')
