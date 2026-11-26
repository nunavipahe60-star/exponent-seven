# Independent CPU check by a two-prime norm identity

The input consists of the three rank bands from `rank_certificate.py` and the certified first- and second-split-prime bounds from the Fourier certificates. Both split primes supply norm equations. Write them as

    4 p^7 = x^2 + q y^2,      4 s^7 = u^2 + q v^2,     p < s.

Take x,u,y,v positive. An odd prime fundamental discriminant in the gap must have q = 3 mod 8 because 2 is inert. Thus x and y have the same parity, and so do u and v. Set g = gcd(y,v), y0 = y/g and v0 = v/g. The integers

    a = |v0 x - y0 u| / 2,    b = (v0 x + y0 u) / 2

satisfy ab = K = |v0^2 p^7 - y0^2 s^7|, with a <= b. This constant is nonzero because p and s are distinct and the valuations of p^7/s^7 at p and s are odd. Conversely, every necessary trace pair occurs among the divisors a <= sqrt(K). Write V = v0^2 p^7 - y0^2 s^7 and reconstruct

    v0 x = b + sign(V) a,    y0 u = b - sign(V) a.

For each rank and pair (p,s), use the integer-square-root endpoints for x,u supplied by the q interval. Bound y,v by their integer square-root limits. For each (p,y) keep the possibilities for q modulo 3,5,7,11,13 which are inert and allow the norm equation. Require the two sets of local possibilities to meet at each prime. The two bounds on v0 x + y0 u bound a = 2K/(v0 x + y0 u), occasionally eliminating a tuple before factoring. For every remaining tuple, factor K with FLINT, multiply the prime powers back, and call the proven `is_prime` routine on each base. Generate all divisors, reconstruct both traces, and check every odd prime smaller than p strictly for inertness.

The three bands require (p1,p2) <= (2371,2389), (2713,2729), and (3089,3119), respectively. These imply (y_max,v_max) = (36,37), (7,7), (1,1). The divisor check includes intermediate q that might be composite. An empty output after the last test excludes prime q as well.

On Python 3.11+, install python-flint 0.9.0, place the three Python scripts in the same folder, and run `python replay_pairs.py`. It divides each rank among the four nonzero residues modulo 5 of the first prime. All possible first primes are at least 853. Alternatively, `python pair_check.py 11 12 13` runs sequentially. `python pair_check.py --count 11 12 13` counts tuples and factor calls without performing them. The printed counts denote local tuples, completed factor calls, reconstructed representations passing q = 3 mod 8, and representations meeting all inertness tests. The script also prints elapsed seconds.

A positive control uses the known discriminant -5923 in Table A.7. In Python, use the following commands in the same folder.

    import pair_check as check
    check.BANDS[1] = (5923, 5923, 19, 29)
    print(check.run(1))

This recovers the record `HIT 1 19 29 758 3413 13128 2243 5923`. It tests both norm equations, the parity convention, the sign in the factor identity, and endpoint inclusion.
