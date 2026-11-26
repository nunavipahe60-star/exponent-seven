"""Verify the 40 A.7 entries by primitive reduced forms and ideal powers.

This checks positive entries. Completeness at small discriminants uses
the published enumeration, and completeness in the gap uses the sieve.
"""
import csv
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path

CYCLIC = [71,151,223,251,284,343,463,467,487,587,604,811,827,859,892,
          1163,1171,1372,1483,1523,1627,1787,1852,1948,1987,2011,2083,
          2179,2251,2467,2707,3019,3067,3187,3907,4603,5107,5923]
RANK_TWO = [63499,118843]
TABLE = Path(__file__).resolve().parent.parent / "known_exponent_seven.csv"
HEADERS = ["discriminant", "is_fundamental", "fundamental_discriminant",
           "conductor_f", "class_group", "class_number"]


def reduced_forms(D):
    for A in range(1,isqrt((-D)//3)+1):
        for B in range(-A,A+1):
            C,rem = divmod(B*B-D,4*A)
            if rem or A>C or gcd(A,gcd(B,C))!=1:
                continue
            if B<0 and (A==abs(B) or A==C):
                continue
            yield A,B,C


def bezout(a,b):
    if b==0:
        return (1 if a>=0 else -1),0
    q,r = divmod(a,b)
    x,y = bezout(b,r)
    return y,x-q*y


def least_value_after_seventh_power(D,A,B):
    # In the order Z[w], w^2=t*w-n. The chosen ideal has a basis
    # A, m+w; changing m by a multiple of A preserves its lattice.
    t = D%2
    n = (t*t-D)//4
    m = (-(B+t)//2)%A
    a,b = 1,0
    for _ in range(7):
        v = [(a*A,0),(a*m,a),(b*A,A),(b*m-n,b+m+t)]
        content = gcd(*(gcd(i,j) for i,j in v))
        v = [(i//content,j//content) for i,j in v]
        index = gcd(*(u*y-vv*x for (u,vv),(x,y) in combinations(v,2)))
        coord,vertical = 0,0
        for x,y in v:
            c,d = bezout(vertical,y)
            coord,vertical = c*coord+d*x,c*vertical+d*y
        assert vertical == 1
        a,b = index,coord%index
        assert all((x-b*y)%a==0 for x,y in v)
    # Binary norm form on the ideal divided by its lattice index.
    first,mid,last = a,2*b+t,(b*b+t*b+n)//a
    assert mid*mid-4*first*last==D
    # Gauss reduction. The leading coefficient is the least value
    # represented by the finished reduced positive definite form.
    while True:
        q=(mid+first)//(2*first)
        mid-=2*q*first
        last=(mid*mid-D)//(4*first)
        if first<=last:
            break
        first,mid,last=last,-mid,first
    return first


def field_and_conductor(discriminant):
    """Compute the fundamental discriminant and f in D = D_K f^2."""
    if discriminant >= 0 or discriminant % 4 not in (0, 1):
        raise ValueError(f"Invalid negative discriminant {discriminant}")
    squarefree = discriminant
    for p in range(2, isqrt(-discriminant) + 1):
        while squarefree % (p*p) == 0:
            squarefree //= p*p
    fundamental = squarefree if squarefree % 4 == 1 else 4*squarefree
    quotient, rem = divmod(discriminant, fundamental)
    conductor = isqrt(quotient)
    assert rem == 0 and conductor**2 == quotient
    return fundamental, conductor


def main():
    with TABLE.open(encoding="ascii", newline="") as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == HEADERS, "Unexpected A.7 CSV columns"
        rows = list(reader)
    indexed = {int(row["discriminant"]): row for row in rows}
    expected = {-d for d in CYCLIC+RANK_TWO}
    assert len(indexed) == len(rows) == 40 and set(indexed) == expected
    total = fields = 0
    for D, row in indexed.items():
        forms = list(reduced_forms(D))
        h = len(forms)
        assert h == (49 if -D in RANK_TWO else 7)
        assert all(least_value_after_seventh_power(D,a,b)==1 for a,b,c in forms)
        fundamental, f = field_and_conductor(D)
        computed = [str(D), "true" if f == 1 else "false", str(fundamental),
                    str(f), "C7" if h == 7 else "C7 x C7", str(h)]
        assert None not in row and [row[key] for key in HEADERS] == computed, D
        fields += f == 1
        total += h
    assert fields == 33 and total == 364
    print(f"A.7 CSV verified: {fields} fields, {len(rows)} orders, "
          f"{len(CYCLIC)} cyclic and {len(RANK_TWO)} rank-two groups; "
          f"all {total} class seventh powers principal.")


if __name__=="__main__":
    main()
