# Certified split primes in the exponent-seven gap

Run `python verify.py` with Python 3.11 or newer. The program uses only the standard library. The file `filters.json` contains integer coefficients for every row, and `certificate.txt` gives the output of a successful verification. Each output row lists the rank label (0 means rank independent), the split-prime index k, the input conductor cap, the support T, the bound S for the k-th split prime, and a rigorous lower margin rounded downwards.

## Hypotheses

Write D=-q for a fundamental discriminant in the ARSX gap, q prime, q >= 3.1*10^20, with ideal class group exponent seven. Assume the Dirichlet GRH for the primitive odd character chi_D. As no q is in the support, every prime has sign +1 (split) or -1 (inert). No ramification penalty is taken. Each split p satisfies q <= 4p^7, by taking a generator of the seventh power of a prime ideal above p. Indeed it is nonrational by unique ideal factorization, and its norm identity is 4p^7=x^2+q y^2 with y nonzero.

The rank-independent rows start at Q=3.9*10^30, from ARSX Theorem 1.3(6). The labeled rows additionally take the following input conductor intervals. The Fourier verifier checks the analytic implications of these intervals. Their independent proof requires the class number formula and the separate rank argument.

| Rank | Lower | Upper |
|---:|---:|---:|
| 11 | 1.259e21 | 3.044e21 |
| 12 | 7.873e22 | 1.538e23 |
| 13 | 4.807e24 | 7.743e24 |

## Guinand–Weil inequality

For each filter, take its stored first half followed by its reversed negation. Let g take these integer values in order on consecutive cells of length (log T)/m, and vanish elsewhere. Write f(u) = integral g(t)g(t+u) dt / integral g(t)^2 dt. The autocorrelation f is even, supported on [-log T,log T], piecewise linear, bounded in absolute value by 1, and satisfies f(0)=1. Its Fourier transform h(t) for the exponential exp(itu) is nonnegative for real t. The transform continues to an entire function decaying as O((1+|Re t|)^(-2)) in each bounded horizontal strip, because the step function transform decays as O((1+|Re t|)^(-1)) there. Under GRH the zero sum in the primitive-character explicit formula is nonnegative. In this normalization the result is

    2 sum_{n<T} Lambda(n) chi_D(n) f(log n)/sqrt(n)
      <= log(q/pi) + psi(3/4) + J,
    J = integral_0^infinity exp(-u/2) (1-f(u))/sinh(u) du.

This follows from the digamma integral for the real part of psi(3/4+it/2), with Fourier inversion and the substitution v=2u. The integrand in the intermediate integral of (1-cos(tv/2)) against h(t) is nonnegative. The source for the explicit formula is Carneiro, Milinovich, Quesada-Herrera and Ramos, Lemma 9, equation (3.1), https://arxiv.org/html/2404.08380v2#S3.Thmlemmatheorem9 .

## Larger Fourier support than the desired prime threshold

For each p<T let A_p be the local left-side contribution if chi(p)=-1. Let B_p be the inert contribution minus the split contribution. They are evaluated including every p^a<T. Thus under the true signs the sum is exactly sum A_p minus sum B_p over the split primes. For p>S the verifier deducts max(0,B_p) regardless of its sign. This permits T to exceed S substantially and is important for the numerical improvements. No independence assumptions on the signs are needed.

For the first split prime, suppose none splits up to and including S. After deducting that tail and log Q, adding log(pi)-psi(3/4), and deducting J, a positive value contradicts the explicit formula. Hence p1<=S and the next conductor cap is 4S^7. For the k-th prime with k>=2, let a be the actual first split prime. Its previously certified cap is P<=S, and 4a^7 >= the lower conductor endpoint. If there are at most k-1 split primes through S, the additional worst deduction is at most

    max over admissible a<=P of
      [log(min(Q,4a^7)) + B_a
       + sum of up to k-2 largest positive B_p over a<p<=S].

In these rows this deduction includes the conductor term. The earlier first-prime bound ensures a maximum over a nonempty admissible set. A positive remaining margin proves p_k<=S.

## Numerical rigor

Autocorrelation values at grid points are rational dot products of integers. The program encloses all needed transcendental and arithmetic evaluations using 48-digit Python Decimal intervals padded outwards. It uses correct rounding of Decimal ln, exp and sqrt. The summands include all integer prime powers below T.

To enclose J, write d=(log T)/m. On the first cell use exp(-u/2)/sinh(u) <= 1/u to bound its integral by 1-f(d). Divide each other cell into 32 equal subcells. Multiply the left endpoint value of the decreasing positive factor by the cell length and the exact mean of the linear function 1-f. The tail is bounded by 4/(3 T^(3/2) (1-T^(-2))). Finally use log(pi)-psi(3/4)=log(pi)+gamma-pi/2+3log(2)>2.23058. The verifier checks this with short bounds on pi and gamma > H_1000-log(1000)-1/2000. The harmonic inequality follows by summing the integral trapezoid error for 1/t.

The filters were discovered by optimizing the piecewise linear function of the autocorrelation coefficients with penalty terms for unknown prime signs. Discovery used floating point. Verification depends only on the integer witnesses. The computed upper bounds make no claim of global optimality.
