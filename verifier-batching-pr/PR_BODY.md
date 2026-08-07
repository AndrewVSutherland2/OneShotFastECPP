Two implementation changes to the reference verifier so that its stated O(n^2 (log n)^2) bit-operation bound (FFT multiplication assumed) holds as implemented, worst case over certificates and with no precomputed tables. Accept/reject behavior is unchanged.

**1. One trial-division pass for the whole chain.** Previously each of the Θ(log n) levels swept all π(n²) primes independently, so every level cost Ω(n²) bit operations regardless of its size (profiled: ~30 ms per level on the 10^200 chain whether the level has 665 or 45 bits). A prime can divide some o_i only if it divides D = ∏ o_i — an O(n)-bit integer for a valid chain — so `batch_prime_divisors` runs a single sweep of remainder trees against D and distributes the ≤ log₂ D hit primes to the levels by direct division. Two cheap pre-screens, both implied by validity (hence sound to reject on), keep the pass on a certificate-independent budget against adversarial input: at most n levels (a valid chain has ≤ 1 + log₂ n), and every o_i below 2n+2 bits (a valid o_i is below the Hasse bound of its modulus, hence below p₀²).

**2. Balanced products of primes.** Batch products (and the radical R in the exact-order step) were built by sequential accumulation, which is quadratic in the product size — Θ̃(n³) in total, with a constant small enough to be invisible below n ~ 10⁶ bits but the true asymptotic order. `balanced_product` builds them by power-of-two pairing, O(M(S) log k), with a sequential fast path for ≤ 8 factors.

Also: the CRT/lcm counterexample from the README (p₀ = 2098153·2102167) is added to the `_INVALID` self-test vectors, and the docstring's cost paragraph is updated (the per-level Õ(n·log o_i) claim undercounted — a level cannot touch every prime below n² in o(n²) bit operations).

**Testing** (verdicts identical to the current verifier in every case):
- `--test` self-suite, including the new CRT counterexample;
- all 20 published chains in certs.csv verify True (2.2× faster: 0.50 s vs 1.12 s total; the 10^200 chain 63 ms vs 143 ms);
- tamper sweep: every single-field +1 mutation of the 10^10 and 10^200 chains rejected by both old and new;
- 120 randomized differentials: `batch_prime_divisors` agrees exactly with the previous per-level `prime_divisors`;
- adversarial inputs (a 5000-bit forged o field; a 300-level garbage chain) rejected by both, now in <0.1 ms instead of after an n²-scale scan.

A companion patch for OneShotPrimalityProofs' voneshot.py (the balanced-product change; its single pass is already chain-global) is in [AndrewVSutherland2/OneShotFastECPP, `verifier-batching-pr/`](https://github.com/AndrewVSutherland2/OneShotFastECPP/tree/short-ecpp/verifier-batching-pr), alongside the complexity analysis (`reports/ecpp-varieties/`, Remark 6.4).

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01QGLZscZDf6ymeo8TopbHhk
