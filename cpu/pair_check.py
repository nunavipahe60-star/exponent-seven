"""CPU-only two-split-prime divisor cross-check for the E=7 rank bands.
Requires python-flint (tested 0.9.0 with FLINT 3.6.0) and Python 3.11+.
"""
import argparse
from math import gcd,isqrt,prod
from pathlib import Path
from time import perf_counter
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent/'lib'))
from flint import fmpz
from rank_certificate import primes_to

BANDS={
  11:(1259*10**18,3044*10**18,2371,2389),
  12:(7873*10**19,1538*10**20,2713,2729),
  13:(4807*10**21,7743*10**21,3089,3119)
}
PRIMES=primes_to(3300)
def ceilroot(n):
    return 0 if n<=0 else isqrt(n-1)+1

def local_options(p,lo,hi):
    a=4*p**7
    options=[]
    for y in range(1,isqrt(a//lo)+1):
        if gcd(p,y)!=1: continue
        y2=y*y
        # q=3 mod 8 and q integral imply this congruence.
        m=8*y2
        xl=max(1,ceilroot(a-hi*y2)); xh=isqrt(a-lo*y2)
        if xh<xl: continue
        residues=range(xl,xh+1) if xh-xl<m else range(m)
        if not any((r*r-a+3*y2)%m==0 for r in residues): continue
        signs=[]
        for ell in (3,5,7,11,13):
            s=0
            squares={k*k%ell for k in range(ell)}
            for qr in range(1,ell):
                if (-qr)%ell in squares:continue
                if (a-qr*y2)%ell in squares:s|=1<<qr
            signs.append(s)
        if all(signs):options.append((y,xl,xh,signs))
    return options

def divisors(n):
    fs=fmpz(n).factor()
    # Verify the full product with FLINT's proof-based primality function.
    assert prod(int(p)**int(e) for p,e in fs)==n
    assert all(p.is_prime() for p,e in fs)
    ds=[1]
    for p,e in fs:
        pp=int(p); power=pp; pre=len(ds)
        for _ in range(e):
            ds.extend(ds[j]*power for j in range(pre))
            power*=pp
    return ds

def run(rank,limit=None,selection=None,count=False):
    low,high,maxp,maxs=BANDS[rank]
    opts={p:local_options(p,low,high) for p in PRIMES if 4*p**7>=low and p<=maxs}
    factors=0; pairs=0; res=0; hit=0; start=perf_counter()
    ps=[p for p in opts if p<=maxp and (selection is None or p%5==selection)]
    for p in ps:
        for s in opts:
            if s<=p:continue
            for y,xmin,xmax,qr in opts[p]:
                for v,umin,umax,tr in opts[s]:
                    pairs+=1
                    if not all(a&b for a,b in zip(qr,tr)):continue
                    g=gcd(y,v); y0=y//g;v0=v//g
                    V=v0*v0*p**7-y0*y0*s**7; K=abs(V)
                    assert K>0
                    end=v0*xmax+y0*umax
                    begin=v0*xmin+y0*umin
                    mn=(2*K+end-1)//end;mx=min(isqrt(K),2*K//begin)
                    if mn>mx:continue
                    factors+=1
                    if count:continue
                    for a in divisors(K):
                        if a<mn or a>mx:continue
                        b=K//a; d=(a if V>0 else -a)
                        vx=b+d; yu=b-d
                        if vx%v0 or yu%y0:continue
                        x=vx//v0;u=yu//y0
                        if not xmin<=x<=xmax or not umin<=u<=umax:continue
                        qn=4*p**7-x*x
                        if qn%(y*y):continue
                        q=qn//(y*y)
                        assert q*v*v==4*s**7-u*u and low<=q<=high
                        if q%8!=3:continue
                        res+=1
                        if all(pow(-q%ell,(ell-1)//2,ell)==ell-1 for ell in PRIMES if 2<ell<p):
                            hit+=1;print("HIT",rank,p,s,y,v,x,u,q,flush=True)
                    if limit and factors>=limit:return pairs,factors,res,hit,perf_counter()-start
    return pairs,factors,res,hit,perf_counter()-start

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('rank',type=int,nargs='+')
    ap.add_argument('--limit',type=int);ap.add_argument('--count',action='store_true')
    ap.add_argument('--subset',type=int)
    args=ap.parse_args()
    for rank in args.rank:
        print(rank,run(rank,args.limit,args.subset,args.count),flush=True)
