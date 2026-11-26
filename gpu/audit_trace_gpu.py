"""Independent CPU/PTX comparison for large traces and clipped high tiles.

Run beside the research tree with trace_sieve.py in the parent directory or
alongside this file. Only a few small GPU blocks are launched.
"""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.append(str(Path(__file__).resolve().parent.parent))
import trace_sieve as production


def true_values(p, residues, first, length, high, kept):
    # Construct square sets directly; no Legendre power or producer mask code.
    squares = {e: {u*u % e for u in range(1, e)} for e in kept}
    return sorted(x for tile in range(first, first+length)
                  for bit in range(32) for r in residues
                  if 0 < (x := int(r)+(tile*32+bit)*production.M) <= high
                  and all((x*x-4*p**7) % e not in squares[e] for e in kept))


def device_test():
    gpu = production.Device(None)
    modulus = production.M
    assert modulus == 111546435 and 32*modulus < 2**32
    p = 4339
    first, length = 2701, 2
    limit = first*32*modulus + (32+17)*modulus + 90000000
    assert 2**40 < limit < math.isqrt(4*p**7-production.LOW)
    rw = production.wheel(p)
    residues = rw[::max(1,len(rw)//73)].copy()
    keep = set(production.TEST[:7])
    masks = production.masks(p)
    offset = 0
    for ell in production.TEST:
        if ell not in keep:
            masks[offset:offset+ell] = np.uint32(2**32-1)
        offset += ell
    answer = true_values(p, residues, first, length, limit, keep)
    hits, _ = gpu.run(residues, masks, limit, [(first, length)])
    assert answer == hits and len(hits) >= 8
    assert any(x//(32*modulus) == first for x in hits)
    assert any(x//(32*modulus) == first+1 for x in hits)
    print('Genuine prime, seven active test primes, high clipped spans:', len(hits), hits[0], hits[-1])

    # A square test parameter gives F=0 at a known nonzero trace. This
    # exercises every unmodified production mask and every integer-width step.
    # It is a control for the sieve, not a hypothetical exponent-seven field.
    n = 47
    p = n*n
    planted = 2*n**7
    idx = planted//(32*modulus)
    r = planted % modulus
    rw = production.wheel(p)
    assert np.any(rw == r), 'real wheel missed the planted zero-symbol trace'
    residues = np.array([int(r), int(rw[0]), int(rw[min(193, len(rw)-1)])], dtype=np.uint32)
    for cap in (planted-1, planted, planted+1):
        hits, _ = gpu.run(residues, production.masks(p), cap, [(idx, 1)])
        answer = true_values(p, residues, idx, 1, cap, set(production.TEST))
        assert hits == answer and (planted in hits) == (planted <= cap)
    print('All genuine test masks, square control, high-end clipping:', planted, 'tile', idx, 'lane', (planted//modulus)%32)

    # Full natural positives from the exploratory 251-prime-cutoff scan.
    # Recompute each with an independent CPU square-set oracle, along with
    # the entire neighboring 96 lanes for several control residues.
    natural = [(3209, 2044007472677), (3853, 3871258885871),
               (4129, 5310890329408), (4133, 7753344152825),
               (4561, 3028955847366)]
    for p, x in natural:
        cap = math.isqrt(4*p**7-production.LOW)
        assert x <= cap
        rw = production.wheel(p)
        r = x % modulus
        assert np.any(rw == r)
        tile = x//(modulus*32)
        residues = np.unique(np.asarray([r, rw[0], rw[len(rw)//2], rw[-1]], dtype=np.uint32))
        answer = true_values(p, residues, tile-1, 3, cap, set(production.TEST))
        hits, _ = gpu.run(residues, production.masks(p), cap, [(tile-1, 3)])
        assert answer == hits and x in answer
        print('Natural positive, true p and all genuine masks:', p, x, 'tile', tile)


if __name__ == '__main__':
    device_test()
