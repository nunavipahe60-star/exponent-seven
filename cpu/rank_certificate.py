"""Directed-interval verification of the E=7 rank and conductor bands.

Inputs: Dirichlet GRH, fundamental discriminant -q, 3.1e20<=q<=4*4339**7,
and elementary 7-class group. The six-step Fourier certificate supplies
the sufficient, smaller upper bound 4*3169**7.
Source of analytic inequalities: Lamzouri--Li--Soundararajan, Section 5.1
(just before the maximization over characters) and inequality (5.1).
Only the Python standard library is required.
"""
from decimal import Decimal as D, getcontext
from functools import cache
from math import isqrt

getcontext().prec=40
class I:
    def __init__(self,x,y=None):
        self.lo=D(x); self.hi=D(y) if y is not None else self.lo
    def __add__(self,r):
        r=iv(r); return widen(self.lo+r.lo,self.hi+r.hi)
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self,r): return self+(-iv(r))
    def __rsub__(self,r): return iv(r)-self
    def __mul__(self,r):
        r=iv(r); v=[a*b for a in (self.lo,self.hi) for b in (r.lo,r.hi)]
        return widen(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,r):
        r=iv(r)
        if r.lo<=0<=r.hi: raise ArithmeticError("interval contains zero")
        return self*widen(1/r.hi,1/r.lo)
    def __rtruediv__(self,r): return iv(r)/self
    def ln(self): return widen(self.lo.ln(),self.hi.ln())
    def exp(self): return widen(self.lo.exp(),self.hi.exp())
    def sqrt(self): return widen(self.lo.sqrt(),self.hi.sqrt())

def widen(a,b): return I(a.next_minus(),b.next_plus())
def iv(x): return x if isinstance(x,I) else I(x)
# Standard rational enclosure of pi, obtainable from Machin's formula.
PI=I('3.141592653589793238462643383279502884197',
     '3.141592653589793238462643383279502884198')
def primes_to(n):
    flags=bytearray([1])*(n+1)
    flags[:2]=bytes(2)
    for k in range(2,isqrt(n)+1):
        if flags[k]: flags[k*k::k]=bytes((n-k*k)//k+1)
    return [k for k in range(2,n+1) if flags[k]]

PRIMES=primes_to(10000)
@cache
def coefficients(x,cut,upper):
    """log L bound = constant + log(q/pi) * coefficient."""
    z=I(x); t=z.ln(); s=z.sqrt(); e=1 if upper else -1
    a=(1/t-e*2/(s*t*t))/((1+e/s)*(1+e/s))
    b=e*(1+(1 if upper else (1+1/s))/t)/((s+e)*t)
    # For a quadratic character prime-power weights are positive. Extremize
    # by chi(p)=-1 throughout for lower, and +1 at unknown p for upper.
    constant=I(0)
    for p in PRIMES:
        if p>x: break
        sign=(-1 if (p<cut or not upper) else 1)
        n=p;k=1;c=sign;lp=I(p).ln()
        while n<=x:
            w=(1-k*lp/t)/(n*k)+a*lp*(1/I(n)-1/z)
            constant+=c*w
            n*=p;k+=1;c*=sign
    return constant,b

def bound(q,x,cut,upper):
    c,b=coefficients(x,cut,upper)
    return (c+b*(I(q)/PI).ln()).exp()

def first_possible(q):
    return next(p for p in PRIMES if 4*p**7>=D(q))

def check():
    lower='3.1e20'; upper=str(4*4339**7)
    a=bound(upper,2100,first_possible(lower),False)
    b=bound(upper,4000,first_possible(lower),True)
    assert (a*I(lower).sqrt()/PI).lo>7**10
    assert (b*I(upper).sqrt()/PI).hi<7**14
    print('Uniform L lower and upper:',a.lo,b.hi)
    # The release also states a four-decimal inequality on the smaller cap
    # supplied by the verified sixth rank-independent Fourier row.
    cap=4*3169**7
    cut=first_possible(lower)
    assert coefficients(2100,cut,False)[1].hi<0
    assert coefficients(4000,cut,True)[1].lo>0
    cap_lo=bound(cap,2100,cut,False).lo
    cap_hi=bound(cap,4000,cut,True).hi
    assert D('0.1084')<cap_lo and cap_hi<D('0.1875')
    print('Six-step cap L lower and upper:',cap_lo,cap_hi)
    # Upper character inequalities give increasing bounds on h in q, since
    # their coefficient b>0. Lower character inequalities give increasing
    # bounds on h as well, since 1/2+b>0.
    steps={
        11:([('1.154e21',3500),('1.254e21',3800),('1.259e21',3800)],'3.044e21',1600),
        12:([('5.557e22',3600),('7.71e22',5800),('7.873e22',5900)],'1.538e23',1800),
        13:([('2.677e24',3800),('4.641e24',8800),('4.792e24',9000),
             ('4.802e24',9000),('4.807e24',9000)],'7.743e24',2000)}
    for r,(sequence,qhi,xhi) in steps.items():
        qlo=lower
        for qnext,x in sequence:
            p=first_possible(qlo)
            assert coefficients(x,p,True)[1].lo>0
            assert (I(qnext).sqrt()*bound(qnext,x,p,True)/PI).hi<7**r
            qlo=qnext
        assert (I('0.5')+coefficients(xhi,0,False)[1]).lo>0
        assert (I(qhi).sqrt()*bound(qhi,xhi,0,False)/PI).lo>7**r
        print('Rank, q interval and first admissible prime:',r,qlo,qhi,first_possible(qlo))

if __name__=='__main__': check()
