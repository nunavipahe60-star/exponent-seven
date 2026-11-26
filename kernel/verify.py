"""Replay the E=7 Fourier-positive certificates. Standard Python >=3.11.
Interval primitives adapted from the earlier Exponent Six verifier.
"""
from decimal import Decimal as D, ROUND_FLOOR, getcontext
from fractions import Fraction
from math import isqrt

getcontext().prec = 48


class Interval:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = D(lo)
        self.hi = D(lo if hi is None else hi)
        assert self.lo <= self.hi

    @classmethod
    def round(cls, lo, hi):
        return cls(lo.next_minus(), hi.next_plus())

    def __add__(self, other):
        return self.round(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        values = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return self.round(min(values), max(values))

    def __truediv__(self, other):
        assert other.lo > 0
        return self * self.round(D(1)/other.hi, D(1)/other.lo)

    def log(self):
        assert self.lo > 0
        return self.round(self.lo.ln(), self.hi.ln())

    def exp(self):
        return self.round(self.lo.exp(), self.hi.exp())

    def sqrt(self):
        assert self.lo >= 0
        return self.round(self.lo.sqrt(), self.hi.sqrt())

    @staticmethod
    def minimum(a,b):
        return Interval(min(a.lo,b.lo),min(a.hi,b.hi))


def rat(q):
    return Interval(q.numerator) / Interval(q.denominator)


def primes_below(x):
    sieve = bytearray([1])*x
    sieve[:2] = bytes(2)
    for p in range(2,isqrt(x-1)+1):
        if sieve[p]:
            sieve[p*p:x:p] = bytes((x-1-p*p)//p + 1)
    return (p for p in range(2,x) if sieve[p])



from collections import defaultdict
from functools import lru_cache
from pathlib import Path
import json

one = Interval(1)

@lru_cache(None)
def ilog(x):
    return Interval(x).log()

@lru_cache(None)
def isqrt_ball(x):
    return Interval(x).sqrt()

def kernel(log_n,width,c):
    m=len(c)-1
    z=log_n/width
    lo=max(0,int(z.lo))
    hi=min(m-1,int(z.hi))
    assert lo<=hi
    parts=[c[i]+(z-Interval(i))*(c[i+1]-c[i]) for i in range(lo,hi+1)]
    return Interval(min(a.lo for a in parts),max(a.hi for a in parts))

def remainder(x,width,c):
    # Left rectangle for decreasing a(u) and exact mean for 1-f(u).
    m=len(c)-1
    integral=one-c[1]
    pieces=32
    step=width/Interval(pieces)
    d1=(Interval('-1.5')*step).exp()
    d2=(Interval(-2)*step).exp()
    a1=(Interval('-1.5')*width).exp()
    a2=(Interval(-2)*width).exp()
    interp=[rat(Fraction(2*j+1,2*pieces)) for j in range(pieces)]
    for i in range(1,m):
        slope=c[i+1]-c[i]
        for j in range(pieces):
            a=Interval(2)*a1/(one-a2)
            integral=integral+step*a*(one-c[i]-slope*interp[j])
            a1=a1*d1
            a2=a2*d2
    return (integral+Interval(4)/(Interval(3*x)*isqrt_ball(x)*(one-one/Interval(x*x)))).hi

def verify(row):
    band,k,low,high,first_upper,T,S,left=row
    v=left+[-w for w in reversed(left)]
    m=len(v)
    den=sum(w*w for w in v)
    c=[rat(Fraction(sum(a*b for a,b in zip(v,v[j:])),den)) for j in range(m+1)]
    width=ilog(T)/Interval(m)
    Jhi=remainder(T,width,c)
    total=Interval(0)
    tail=Interval(0)
    losses={}
    for p in primes_below(T):
        sign=-1
        n=p
        loss=Interval(0)
        A=Interval(0)
        while n<T:
            term=Interval(2*sign)*ilog(p)*kernel(ilog(n),width,c)/isqrt_ball(n)
            A=A+term
            if sign==-1:
                loss=loss+Interval(2)*term
            sign=-sign
            n*=p
        total=total+A
        if p>S:
            tail=tail+Interval(max(D(0),loss.hi))
        elif 4*p**7>=low:
            losses[p]=loss.hi
    if k==1:
        worst=ilog(high)
    else:
        assert first_upper<=S
        costs=[]
        for p,loss in losses.items():
            if p>first_upper:
                continue
            top=sorted([d for ell,d in losses.items() if ell>p and d>0],reverse=True)[:k-2]
            cost=ilog(min(high,4*p**7))+Interval(loss)
            for d in top:
                cost=cost+Interval(d)
            costs.append(cost.hi)
        assert costs
        worst=Interval(max(costs))
    # log(pi)-psi(3/4)=log(pi)+gamma-pi/2+3log2 > 2.23058.
    # pi in [3.14159265,3.14159266]; gamma > H_1000-log(1000)-1/2000.
    gam=rat(sum((Fraction(1,j) for j in range(1,1001)),Fraction()))-ilog(1000)-rat(Fraction(1,2000))
    const=Interval('3.14159265').log()+gam-Interval('3.14159266')/Interval(2)+Interval(3)*ilog(2)
    assert const.lo>D('2.23058')
    margin=total+Interval('2.23058')-Interval(Jhi)-tail-worst
    shown=margin.lo.quantize(D('.000001'),rounding=ROUND_FLOOR)
    assert margin.lo>0,(band,k,T,S,margin.lo)
    return shown

if __name__=='__main__':
    rows=json.loads(Path(__file__).with_name('filters.json').read_text(encoding='utf-8'))
    by_band=defaultdict(dict)
    # The only initial unconditional-in-rank inputs are the ARSX bounds.
    cap=39*10**29
    for row in rows:
        b,k,lo,hi,P,T,S,_=row
        assert 2<=S<T<lo
        if k==1 and b==0:
            assert hi==cap
            cap=4*S**7
        if k>1:
            assert by_band[b][1]<=P
        if b==0 and k>1:
            assert hi==cap
        margin=verify(row)
        print(f'{b:2} {k:3} {hi:32} {T:6} {S:6} {margin:11}',flush=True)
        by_band[b][k]=S
