# Imaginary quadratic fields and orders of class-group exponent seven

<!-- Zenodo DOI badge: insert the badge for this repository when its DOI has been assigned. -->

Under Dirichlet GRH, exactly **33 imaginary quadratic fields** have ideal class group of exponent seven. The 33 maximal orders and seven nonmaximal orders give **40 imaginary quadratic orders** of Picard exponent seven. The two groups of type $C_7^2$ have discriminants $-63499$ and $-118843$. The other 38 groups are cyclic of order seven. The largest absolute discriminant is $118843$.

The new calculation excludes every missing fundamental discriminant of exponent seven in the published gap $3.1\times10^{20}<|D|\leq3.9\times10^{30}$. The low-range field enumeration and the passage to arbitrary conductors are published inputs. This closes Aiken, Ross, Shvydkoy and Xue (ARSX), Theorem 1.3(6) and Table A.7. Their results for the other exponents, together with the separate exponent-six result closing Theorem 1.3(5), give completeness of Tables A.1 through A.8 under Dirichlet GRH. The exponent-seven proof is independent of the exponent-six proof. [EKN](https://arxiv.org/abs/1803.02056) · [ARSX](https://arxiv.org/abs/2605.26020) · [Exponent-six archive](https://doi.org/10.5281/zenodo.24489676)

## Evidence

| File | Contents |
| --- | --- |
| [Classification_Proof.pdf](Classification_Proof.pdf) and [LaTeX source](Classification_Proof.tex) | The field and order classifications, the split-prime and rank bounds, and two finite exclusions |
| [known_exponent_seven.csv](known_exponent_seven.csv) | All 40 order discriminants, fundamental status, conductor, class-group structure and class number |
| [kernel/](kernel/) | The explicit-formula normalization, 57 integer Fourier witnesses, interval verifier and certified margins |
| [cpu/](cpu/) | Short Euler-product intervals, reduced forms and seventh ideal powers, the endpoint 701, and the two-prime divisor enumeration |
| [gpu/](gpu/) | PTX and CUDA C++ trace sieves, wheel-coverage and witness audit, and high-trace positive controls |
| [logs/](logs/) | Two complete wide GPU records, one complete short GPU record, the CPU divisor totals and timings |
| [replay.py](replay.py) and [verify.sh](verify.sh) | Certificate replay, optional regeneration of the searches, and the one-command release check |
| [CLAIMS.md](CLAIMS.md) | Dependencies, computational evidence and independent-review status for each claim |
| [SHA256SUMS.txt](SHA256SUMS.txt) | SHA-256 digests of all other release files |
| [CITATION.cff](CITATION.cff) and [LICENSE](LICENSE) | Citation metadata and the Apache License, Version 2.0 |

## Three ways to check the gap

| Route | Certified search domain | Recorded result |
| --- | --- | --- |
| Direct GPU | The first 256-cell Fourier certificate gives $p_1\leq4597$ at support parameter $T=20550$. Each of the two sieve implementations covers 497 primes and 1,622,545,300,159,089 positive trace positions. | Both yield the same five traces after the tests through 251. For each trace, a smaller prime gives a nonzero square residue and excludes the proposed least split prime. |
| Short GPU | Six rank-independent Fourier steps give $p_1\leq3169$. The scan covers 324 primes and 311,248,928,889,382 positive trace positions. | Zero traces survive the tests through 251. The recorded trace run took 63.84 seconds on an RTX 2060. |
| Two-prime CPU | The class number formula and short Euler products give ranks 11, 12 and 13. Per-rank Fourier bounds on $p_1,p_2$ permit a signed-divisor search from two norm equations. | 533,439 factorizations of integers up to 89 bits produced 19,572 in-band representations counted with multiplicity. Zero pass the final inertness tests. The six-core i7-10750H run took 877.91 seconds. |

The GPU test allows zero modulo a sieving prime. A nonzero quadratic residue of $x^2-4p^7$ excludes a trace for every admissible $y$ at once. The two GPU implementations use different CRT wheels, constructions of residue tables, mask orders and device kernels. They share a host program and analytic input. The CPU route uses integer factorization, divisor enumeration and proven primality tests on the factors. The known discriminant $-5923$ is a positive control for that reconstruction. [CPU divisor method](cpu/CPU-Divisor-Check.md) · [Coverage and witnesses](logs/coverage-and-witnesses.txt)

## Short independent review

Run this from the repository root in Bash or Git Bash, with Python 3.11 or later and assertions enabled.

```sh
./verify.sh
```

The script checks the release manifest and then runs `python replay.py`. The default replay needs only Python's standard library. It evaluates the 57 positive-kernel inequalities with outward-rounded intervals, checks the LLS rank bounds, validates the CSV, and enumerates the reduced forms for all 40 listed discriminants. All 364 ideal-class seventh powers are principal. It also audits all three GPU records, checking prime coverage, integer trace endpoints, wheel sizes, tile totals and the final modular witnesses. Matching recorded candidates provides a direct comparison of the two wide runs.

The default replay reads the saved sieve outputs. Fresh exhaustive enumeration requires the full commands below. The mathematical argument also uses Dirichlet GRH, the cited theorems, the explicit-formula specialization and the reductions from ideals to search predicates. Independent review of these steps is tracked in the ledger. [Claims and dependencies](CLAIMS.md)

## Regenerate the searches

For the CPU search, the recorded environment used Python 3.11 and python-flint 0.9.0. The command reruns the analytic checks and the three rank bands.

```sh
python -m pip install python-flint==0.9.0
python -B replay.py --cpu
```

For the GPU search, install NumPy and use an NVIDIA Turing-or-newer GPU with a compatible driver. The primary kernel uses driver-assembled PTX. The recorded runs used NumPy 2.4.3, an RTX 2060 and driver 576.52. The following Git Bash commands write new records in a sibling directory.

```sh
python -m pip install numpy==2.4.3
repo=$(pwd)
mkdir -p ../exponent-seven-reruns
(cd ../exponent-seven-reruns && python -B "$repo/replay.py" --gpu)

python -B gpu/trace_sieve.py --short --selftest > ../exponent-seven-reruns/short.log
python -B gpu/audit_log.py --upper-prime 3169 ../exponent-seven-reruns/short.log
```

For the second wide implementation, set `NVRTC_LIBRARY` to the filename of the installed NVRTC shared library, including its directory. The recorded CUDA C++ run used NVRTC 12.9.86. Run these commands after creating the output directory above.

```sh
python -B gpu/trace_sieve.py --full --selftest --crosscheck --nvrtc "$NVRTC_LIBRARY" > ../exponent-seven-reruns/crosscheck.log
python -B gpu/audit_log.py logs/trace-primary.log ../exponent-seven-reruns/crosscheck.log
```

## Credits and citation

**nunavipahe60-star** coordinated the project in collaboration with Claude (R1), Codex (R2) and Gemini (R3). R2 produced the computations and certificates. The maintainer reports that R3 independently replayed and audited the verification package. R1's independent review of this release is pending. The ledger distinguishes that audit from additions made during release preparation. For citation, use the repository metadata. [Claims ledger](CLAIMS.md) · [Citation metadata](CITATION.cff)

The analytic inputs are the Guinand-Weil formula of Carneiro, Milinovich, Quesada-Herrera and Ramos and the conditional character estimates of Lamzouri, Li and Soundararajan. Holmin, Jones, Kurlberg, McLeman and Petersen provide the earlier odd class-group data and the wider conjectural context. [Guinand-Weil formula, Lemma 9](https://arxiv.org/abs/2404.08380) · [LLS, corrected version, Section 5](https://arxiv.org/abs/1309.3595v3) · [Odd class-group data](https://arxiv.org/abs/1510.04387)
