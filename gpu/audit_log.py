"""Audit prime coverage, bounds, retained traces and exact rejection witnesses.

Pass one or more complete --full logs from trace_sieve.py. Distinct wheels and
mask orders must produce exactly the same list of (p,x) pairs.
"""
import argparse
import csv
from math import isqrt, prod
from pathlib import Path
import re

LOW = 310000000000000000000
UPPER_PRIME = 4597


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n)+1))


def read(path):
    raw = Path(path).read_bytes()
    text = raw.decode("utf-16" if raw[:2] in (b"\xff\xfe",b"\xfe\xff") else "utf-8-sig")
    rows = list(csv.reader(line for line in text.splitlines() if re.match(r"^\d+,", line)))
    hits = [tuple(map(int, m)) for m in re.findall(r"^trace p=(\d+) x=(\d+)\r?$", text, re.M)]
    return rows, sorted(hits)


def check(path, upper_prime=UPPER_PRIME):
    rows, hits = read(path)
    expected = [p for p in range(2, upper_prime+1) if prime(p) and 4*p**7>=LOW]
    if [int(row[0]) for row in rows] != expected:
        raise ValueError("Incomplete prime list")
    widths = []
    wheels = [(3,5,7,11,13,17,19,23), (3,5,7,11,13,17,19,29)]
    for p,row in zip(expected, rows):
        n = isqrt(4*p**7 - LOW)
        # Recompute both the number of wheel classes and the number of tiles
        # independently of the lifting routine and the producer's tile count.
        def valid(wheel):
            def classes(ell):
                squares = {j*j % ell for j in range(1,ell)}
                return sum((i*i-4*p**7) % ell not in squares for i in range(ell))
            return (prod(classes(ell) for ell in wheel) == int(row[2])
                    and n//(32*prod(wheel))+1 == int(row[3]))
        wheels = [wheel for wheel in wheels if valid(wheel)]
        if n != int(row[1]) or int(row[3]) != int(row[4]) or not wheels:
            raise ValueError(f"Incomplete trace domain at {p}")
        if int(row[7]) != sum(a==p for a,b in hits):
            raise ValueError(f"Incorrect hit count at {p}")
        widths.append(n)
    if len(hits) != len(set(hits)):
        raise ValueError("Duplicate trace")
    print(f"{len(rows)} primes, {sum(widths)} positive trace positions, {len(hits)} retained")
    print("GPU seconds recorded by launches:", round(sum(float(r[5]) for r in rows),3))
    for p,x in hits:
        f = x*x-4*p**7
        if not 1<=x<=isqrt(4*p**7 - LOW):
            raise ValueError("Trace outside domain")
        def symbol(ell):
            return pow(f%ell, (ell-1)//2, ell)
        if any(symbol(ell)==1 for ell in range(3,252,2) if prime(ell)):
            raise ValueError("Reported survivor fails the GPU presieve")
        w = next((ell for ell in range(253,p,2) if prime(ell) and symbol(ell)==1), None)
        if w is None:
            raise ValueError(f"Unresolved trace: p={p}, x={x}")
        print(f"p={p} x={x} witness={w} F_mod_witness={f%w}")
    return hits


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("logs", nargs="+")
    ap.add_argument("--upper-prime", type=int, choices=(3169,4597), default=UPPER_PRIME)
    args = ap.parse_args()
    found = None
    for name in args.logs:
        result = check(name, args.upper_prime)
        if found is not None and result != found:
            raise ValueError("Different retained traces")
        found = result
