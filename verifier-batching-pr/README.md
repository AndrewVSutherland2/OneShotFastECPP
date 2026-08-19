# Batched reference verifiers: vsmallECPP.py and voneshot.py

Copy-ready replacements for the reference verifiers of
[ShortPrimalityProofs](https://github.com/AndrewVSutherland/ShortPrimalityProofs) and
[OneShotPrimalityProofs](https://github.com/AndrewVSutherland/OneShotPrimalityProofs),
making each an O(n^2 (log n)^2)-bit-operation algorithm *as implemented* (FFT
multiplication assumed), with no precomputed tables.  Accept/reject behavior is
unchanged.  Written by Claude (Fable 5) under the direction of A. V. Sutherland;
see reports/ecpp-varieties/ in this repository (Remark 6.4) for the analysis.

## The two defects in the current implementations

1. **Sequential batch-product builds** (both files).  `prime_divisors` builds each
   batch product with `Q = 1; for t in batch: Q *= t`, which is quadratic in the
   batch size (each step re-scans the growing partial product).  With batch_bits
   comparable to the input and ~1.44 n^2 bits of primes, this is Theta(n^3 / log n)
   word operations at the top level -- the true asymptotic order of both verifiers,
   though with so small a constant that it is invisible below n ~ 10^6 bits.

2. **Per-level repetition of the pi(n^2) pass** (vsmallECPP.py only).  Each of the
   Theta(log n) levels independently sweeps all primes <= n^2, so each level does
   Omega(n^2)-bit work regardless of its size (profiled: ~30 ms per level at
   10^200 whether the level has 665 or 45 bits), for Theta(n^2 (log n)^3) total
   even with tree-built products.

## The changes

Both files:

- `balanced_product(xs)`: product by power-of-two pairing, O(M(S) log k) for k
  factors of total size S (sequential fast path for k <= 8, where sequential is
  already O(M(S))).  Used for every product of primes: the trial-division batch
  products, the radical R in the exact-order step, and (voneshot) `is_smooth`.

vsmallECPP.py additionally:

- `batch_prime_divisors(os, primes)` replaces the per-level `prime_divisors`
  calls: a prime can divide some o_i only if it divides D = prod(os) -- an
  O(n)-bit integer for a valid chain -- so ONE sweep of remainder trees against
  D finds the at most log2(D) primes dividing any level, which are then
  distributed to the levels by direct division.  One pass instead of one per
  level.
- Two cheap pre-screens before the batched pass, both implied by validity (so
  rejecting on them is sound) and needed so the pass runs on a
  certificate-independent budget against adversarial inputs: at most n levels
  (a valid chain has at most 1 + log2 n), and every o_i below 2n+2 bits (a
  valid o_i is below the Hasse bound of its modulus, hence below p_0^2).
- The CRT/lcm counterexample from the repository README
  (p_0 = 2098153 * 2102167) added to the `_INVALID` self-test vectors.
- Docstring cost paragraph updated (the previous per-level claim of
  O~(n log o_i) undercounted: a level cannot touch every prime below n^2 in
  o(n^2) bit operations).

With these, the whole-verifier costs are O(n^2 (log n)^2) bit operations and
O(n^2) bits of memory, worst case over certificates, with the (implicit) prime
table construction included -- matching the repositories' stated targets.

## Behavior preservation and tests

- Self-tests of both files pass (`--test`), including the fold-attack vectors
  (voneshot) and the new CRT counterexample (vsmallECPP).
- All 20 published short chains verify True in old and new (10^10 ... 10^200).
- All 51 published one-shot certificates plus this repository's 18 local
  certificates (including the n^2-smooth restricted ones, empty q-list) verify
  True in old and new.
- Tamper sweep: every single-field +1 mutation of the 10^10 and 10^200 chains
  is rejected by both, with identical verdicts.
- Randomized differential: `batch_prime_divisors` agrees with per-level
  `prime_divisors` on 120 random multi-level inputs.
- Adversarial inputs: a 5000-bit o field and a 300-level garbage chain are
  rejected by both (the new pre-screens make the rejection immediate instead of
  after an n^2-scale scan).

## Measurements (Ryzen AI Max+ 395, CPython 3.12)

| workload                              | old      | new      |
|---------------------------------------|----------|----------|
| all 20 short chains                   | 1.12 s   | 0.50 s   |
| the 10^200 chain alone (665 bits)     | 0.143 s  | 0.063 s  |
| all 69 one-shot certificates          | 0.60 s   | 0.62 s   |
| reject a 5000-bit forged o field      | 57 ms    | < 0.1 ms |

The short-chain speedup is the pass count (levels -> 1); at challenge sizes the
sequential-product defect is invisible (its constant is ~1/(178 log n) words),
so voneshot timings are unchanged within noise -- the point of that fix is the
asymptotic claim, not present-day speed.
