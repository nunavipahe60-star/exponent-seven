# Claims ledger

**nunavipahe60-star** coordinates the research and release in collaboration with Claude, Codex and Gemini. The reviewer roles follow the exponent-six project.

| Role | Assignment |
| --- | --- |
| **R1 = Claude** | Independent review of the analysis, algebra and finite exclusion. Review of this exponent-seven repository is **open / pending R1**. |
| **R2 = Codex** | Derivation, implementation, production computation and preparation of the release. |
| **R3 = Gemini** | Independent replay and source audit of the exponent-seven verification package, as attested by the maintainer. |

A claim ID such as **R1** in the last table identifies a result. In the review columns **R1**, **R2** and **R3** identify the reviewers above. The dependency column gives sufficient inputs for each claim.

**R3 verified** records the maintainer's report that R3 ran the original default replay and audited the source line by line. For the GPU and CPU production runs, the ledger distinguishes a checked program or saved record from the person responsible for the actual exhaustive run. **External** identifies a theorem or tabulation accepted from its cited source. The paper, CSV cross-check, release wrapper and sharper presentation bound introduced for this release have separate provenance. A substantive revision to a claim or its inputs reopens the affected review.

## External inputs

| ID | Input | Source and scope | Depends on | Status |
| --- | --- | --- | --- | --- |
| E1 | The field list through absolute fundamental discriminant 3.1e20 contains exactly 33 fields of exponent seven. | EKN, Theorem 1 and Table 1. The low-range classification is unconditional. [EKN](https://arxiv.org/abs/1803.02056) | Published enumeration | External |
| E2 | Under Dirichlet GRH the remaining field search has upper endpoint 3.9e30. | EKN, Theorem 3 and Table 2. ARSX use this endpoint in Theorem 1.3(6). [EKN](https://arxiv.org/abs/1803.02056) · [ARSX](https://arxiv.org/abs/2605.26020) | Dirichlet GRH | External |
| E3 | Completeness of the fundamental exponent-seven list in the gap implies completeness of ARSX Table A.7 for all orders. The table lists the 40 positive instances. | ARSX, Theorem 1.3(6), Sections 5-6 and Table A.7. Their conductor calculation includes field exponents one and seven. [ARSX](https://arxiv.org/abs/2605.26020) | Published conductor argument and table | External |
| E4 | The primitive Dirichlet-character explicit formula and conditional short Euler-product inequalities apply with the odd quadratic character, conductor q and the stated truncation parameters. | Carneiro et al., Lemma 9, and Lamzouri, Li and Soundararajan, corrected version, Section 5. [Guinand-Weil formula](https://arxiv.org/abs/2404.08380) · [LLS](https://arxiv.org/abs/1309.3595v3) | Dirichlet GRH for the conditional applications | External |

## Algebraic layer

| ID | Claim | Evidence | Depends on | R3 verification and R1 status |
| --- | --- | --- | --- | --- |
| B1 | A field of exponent seven has elementary abelian 7-class group. Genus theory in the gap gives fundamental discriminant D = -q for a prime q congruent to 3 modulo 4. | Paper, Section 1. The one-prime genus condition also allows the small exceptional fundamental discriminants -4 and -8 before the lower endpoint is imposed. | Classical genus theory, lower endpoint of E1 | Mathematical argument **open / pending R1** |
| B2 | Each split prime p yields positive x,y with 4p^7 = x^2 + qy^2. Thus q is at most 4p^7 and p is at least 701. The prime 2 is inert, so q is 3 modulo 8. | Paper, Section 1. The full trace sieve includes 701. The separate boundary computation uses y = 1 and excludes 701 by the prime 137. [Endpoint code](cpu/exceptional_701.py) | B1, unique factorization of ideals | **R3 verified** for source audit. Algebraic deduction **open / pending R1** |
| B3 | A nonzero square residue of x^2 - 4p^7 at an odd prime below the least split prime excludes that trace. If this value is zero modulo the testing prime, a factor of y can account for it. | Paper, Sections 1 and 4. The sieve retains both zero and nonsquare residues. [GPU source](gpu/trace_sieve.py) · [Witness audit](gpu/audit_log.py) | B1, B2, quadratic splitting law | **R3 verified** for implementation. Proof review **open / pending R1** |

## Analytic layer

| ID | Claim | Evidence | Depends on | R3 verification and R1 status |
| --- | --- | --- | --- | --- |
| A1 | With T greater than the split-prime threshold S, the autocorrelation has nonnegative transform, the zero sum has the GRH sign, and each unknown prime above S contributes at most its individual worst loss. There is no ramified prime in the support. | Paper, Section 2. The kernel notes give the normalization, odd gamma factor, prime-power signs and tail bound. [Kernel argument](kernel/README.md) | B1, E4, Dirichlet GRH | **R3 verified** for the associated verifier source. Analytic derivation **open / pending R1** |
| A2 | Starting at 3.9e30, the first 256-cell certificate with T = 20550 forces the least split prime to be at most 4597, with certified margin at least 0.059634. Six steps give successive thresholds 4597, 3643, 3313, 3203, 3181 and 3169. The rank-independent higher-index rows give bounds for the subsequent split primes. | Positive margins from directed interval evaluation with integer kernel coefficients. [Witnesses](kernel/filters.json) · [Verifier](kernel/verify.py) · [Certified rows](kernel/certificate.txt) | E2, A1, B2 | **R3 verified** for certificate replay and source. **Open / pending R1** |
| A3 | The class number formula and short Euler products restrict a gap field to ranks 11, 12 or 13. The respective q intervals are [1.259e21, 3.044e21], [7.873e22, 1.538e23] and [4.807e24, 7.743e24]. Inserting these intervals in the Fourier witnesses gives respective upper bounds (2371,2389), (2713,2729) and (3089,3119) for the first two split primes. | Paper, Sections 2-3. The archived broad-cap calculation gave L(1,chi) below 0.18752520. At the smaller cap 4 times 3169 to the seventh power the same interval inequalities give L(1,chi) above 0.10919961 and below 0.18665978, which implies the paper's 0.1084 and 0.1875 bounds. [Rank verifier](cpu/rank_certificate.py) · [Rank-labeled Fourier rows](kernel/certificate.txt) | A2, A1, B1, B2, E4, class number formula | **R3 verified** for the original rank bands and Fourier rows. R2 checked the sharper release bound. **Open / pending R1** |

The paper's two auxiliary observations explain the role of rank in choosing a search. Rank one gives q at most 4 times the cube of the first split prime times the second. The vectors (1,t,t squared,t cubed) over an extension of the field of seven elements obstruct a general degree-six signed-relation argument at high rank. Both observations have independent mathematical review **open / pending R1**. The direct GPU exclusion uses A2, B2 and B3 for every rank.

## GPU trace exclusion

| ID | Claim | Evidence | Depends on | R3 verification and R1 status |
| --- | --- | --- | --- | --- |
| S1 | The positive integer trace interval and CRT wheels include every genuine trace. The last tile is clipped at the exact integer square root. Each host wheel and its device kernel jointly test all odd primes through 251, retaining zero residues. | Paper, Section 4. PTX uses a wheel through 23. The CUDA C++ variant uses 29 in place of 23, direct CRT and a different mask order. [Sieve](gpu/trace_sieve.py) · [High-trace controls](gpu/audit_trace_gpu.py) | A2, B2, B3 | **R3 verified** for source. **Open / pending R1** |
| S2 | The two wide runs cover 497 prime rows and 1,622,545,300,159,089 raw positions each, retaining the same five (p,x) pairs. The six-step short run covers 324 rows and 311,248,928,889,382 positions, retaining zero. | R2 production records. R3's default replay checked the prime lists, endpoints, wheel classes, tile totals, hit counts and agreement of the wide traces. [Primary](logs/trace-primary.log) · [Cross-check](logs/trace-crosscheck.log) · [Short run](logs/trace-short.log) | S1 | **R3 verified** for original record replay and code. **Open / pending R1** |
| S3 | The five wide survivors, in increasing prime order, have exclusion witnesses 263, 271, 257, 269 and 257. In every row the witness prime is smaller than p and the residue is a nonzero square. No trace remains. | Recomputed directly from the saved integers. [Exact witness audit](gpu/audit_log.py) · [Recorded values](logs/coverage-and-witnesses.txt) | B3, S2 | **R3 verified** for replay and source. **Open / pending R1** |

## CPU two-prime cross-check

| ID | Claim | Evidence | Depends on | R3 verification and R1 status |
| --- | --- | --- | --- | --- |
| P1 | Eliminating q from two seventh-power norm equations and dividing the two coefficients by their gcd gives a nonzero signed factorization. The parity conditions make both half-factors integral. Enumerating the allowed divisors recovers every possible pair of positive traces. | Paper, Section 4. The input intervals bound each trace and the smaller divisor. The input prime pair is ordered. [Derivation](cpu/CPU-Divisor-Check.md) · [Implementation](cpu/pair_check.py) | B2, A3 | **R3 verified** for source. Algebraic argument **open / pending R1** |
| P2 | In the three rank bands the maxima for the two coefficients (y,v) are (36,37), (7,7) and (1,1). Local and interval tests retain 533,439 integers to factor across 587,350 allowed local tuples. The largest factorization input has 89 bits. FLINT proves primality of the returned bases and the program multiplies the prime powers back to the input. | Factor calls by rank are 493,357, 39,339 and 743. [Divisor source](cpu/pair_check.py) · [Production totals](logs/cpu-divisors.log) | P1, A3, B3 | **R3 verified** for original source audit. Production run by R2. **Open / pending R1** |
| P3 | The divisor search reconstructs 13,115, 5,802 and 655 in-band representations by rank. All three final counts are zero after the inertness tests. As a positive control, the same reconstruction recovers q = 5923 with split primes 19 and 29. | Independent arithmetic from the GPU trace runs. Run the CPU option to regenerate all 19,572 representations counted with multiplicity. [CPU record](logs/cpu-divisors.log) · [Positive control](cpu/CPU-Divisor-Check.md) | P2, B3 | **R3 verified** for original source audit. Production run by R2. **Open / pending R1** |

## Results

| ID | Claim | Evidence | Depends on | R3 verification and R1 status |
| --- | --- | --- | --- | --- |
| R1 | Under Dirichlet GRH there are exactly 33 imaginary quadratic fields of class-group exponent seven. | EKN below the gap. The first Fourier row and wide GPU exclusion suffice in the gap. The six-step short GPU route and the CPU divisor route also exclude it. | E1, E2, B1-B3, A1-A2 and S1-S3. Alternative finite exclusion through A3 and P1-P3. | Computational evidence **R3 verified** to the scopes above. Theorem review **open / pending R1** |
| R2 | Under Dirichlet GRH there are exactly 40 imaginary quadratic orders of Picard exponent seven. The two rank-two discriminants are -63499 and -118843. The other 38 groups have order seven. The largest absolute discriminant is 118843. | ARSX for conductors. The positive-entry check enumerates 40 reduced-form class groups and verifies all 364 class seventh powers. The release CSV records the 33 maximal orders and seven nonmaximal orders. [Positive-entry verifier](cpu/check_table.py) · [Data](known_exponent_seven.csv) | R1 (field result), E3, positive-entry computation | **R3 verified** for the original 40 entries and 364 powers. R2 checked the release CSV. Theorem review **open / pending R1** |

The all-eight-exponents corollary combines R2, ARSX Theorem 1.3(1)-(5) and the separate result on exponent six. The latter is an external input for this corollary only. [Exponent-six archive](https://doi.org/10.5281/zenodo.24489676)
