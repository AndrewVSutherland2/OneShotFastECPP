# OneShotFastECPP — design notes

Working notes for a CM-method ("fast ECPP" of Atkin–Morain + Shallit; see Morain
arXiv:math/0502097 and Enge arXiv:2404.05506, `fastecpp.pdf`) approach to the
one-shot elliptic-curve primality-proving challenge
(https://github.com/AndrewVSutherland/OneShotPrimalityProofs).

## Goal & algorithm

Given a probable prime `p` (target sizes here are **small**, ~128–384 bits, and
eventually larger), produce a *one-shot* primality certificate: a single CM
discriminant `D<0` and curve `E/F_p` whose order `p+1∓t` contains a **sufficiently
large smooth divisor**, so primality follows in one step (no downrun).

For fixed `p`, iterate over fundamental imaginary-quadratic discriminants `D`:
1. **D-search / Cornacchia** — find `D` with `4p = t² − v²D = t² + |D|v²` solvable. *(component a — DONE)*
2. **Smoothness** — test whether the n⁴-smooth part of `p+1−t` or `p+1+t` exceeds `L≈√p`. *(component b — DONE)*
3. **Class polynomial** — compute `H_D` (best class invariant) mod `p` with `classpoly`. *(component d)*
4. **Root of `H_D` over `F_p`** — big-`F_p` root-finding (`ff_poly` is word-size only). *(component c — DONE)*
5. Build `E` with trace `±t` from the root, find a point of order `m`, emit the certificate. *(component d — DONE)*

Time-critical pieces: (a) Cornacchia over *many* `D`, (b) smoothness of few-hundred-bit
integers, (c) root-finding of a large `H_D` over `F_p`.

**Key difference from textbook fastECPP:** their `p` is thousands of bits and `D` is
tiny relative to `p`; ours is the opposite (`p` small, `D` up to `10^10`+), which moves
where the time goes (see below).

## Repo layout

- `classpoly_v1.0.3/` — Drew's class-polynomial library (Hilbert CRT method + many class
  invariants). Builds into `./local`. **Modified by us** (see "classpoly changes"); also holds
  Drew's `class_inv_mpz.c` (big-p invariant→j, linked by the ecpp tools).
- `ff_poly_v2.0.0/` — word-size `F_p` / `F_p[x]` library; classpoly's dependency.
- `zp_poly/` — Drew's large-p `F_p[x]` library (Harvey–Sutherland); needed by `class_inv_mpz.c`.
- `phi_files/` — modular polynomials. A 46 MB subset is committed (`phi_files_manifest.txt`);
  the full 2.2 GB database stays external (see INSTALL).
- `local/` — project-local install prefix for ff_poly + zp_poly (gitignored).
- `work/` — scratch, class-poly output, and `pcache/` prime-product segments (gitignored).
- `tests/` — classpoly correctness suite vs PARI/GP (Tests 1/2/3).
- `ecpp/` — **the ECPP code**: dscan (a), smooth (b), fproot (c), cminv/cm_method + curve +
  oneshot (d), plus validation drivers (roottest, invjtest, cmjtest, asmtest, zpbench).
- `certs/` — voneshot-verified example certificates (2²⁵⁵−19, NIST P-256, secp256k1, 10⁸⁰/⁹⁰/¹⁰⁰+ε).
- `Makefile` (top level) — `make -j`: ff_poly → classpoly → zp_poly → ecpp; `make test`.
- `setenv.sh` — env vars to run everything in-tree (phi dir, PATH, `work/pcache`).
- `fastecpp.pdf` — Enge's 2024 fastECPP-over-MPI paper (reference, not committed).

## Build & test

```sh
make -j            # ff_poly -> classpoly -> zp_poly -> ecpp tools (~15 s)
make test          # classpoly suite vs PARI over |D|<=1000 (Tests 1/2/3); ~2.5 min
make test MAXD=200 # smaller/faster
. ./setenv.sh && ./ecpp/oneshot p=<prime>      # the headline tool
python3 ecpp/test_dscan.py     # validate the discriminant scan vs brute force + PARI
python3 ecpp/test_smooth.py    # validate the smoothness engine vs PARI
```

Environment: `. ./setenv.sh` points classpoly at `phi_files/`, sets `work/` output dirs and
`work/pcache`, and puts the tools on PATH. Toolchain: gcc 13, GMP 6, PARI/GP 2.18.
Box: **16 physical cores** (32 vCPUs, HT).

## Component (a): CM-discriminant search — DONE

Code: `ecpp/cornacchia.{c,h}` (solver) + `ecpp/dscan.c` (the scan) + `ecpp/cmsearch.c`
(a simpler baseline that enumerates *all* fundamental D and Legendre-filters).

### Cornacchia (`cornacchia.{c,h}`)
Per-`p` context precomputes `4p`, `L=⌊2√p⌋`, and Tonelli-Shanks data. `cornacchia_solve(d)`
does Legendre filter → modular sqrt → Euclidean descent → square test, for `4p=t²+d·v²`,
`d=|D|` (`d ≡ 0 or 3 mod 4`). `cornacchia_solve_x0(d, x0)` takes a **precomputed** `√(−d) mod p`
(so the scan supplies the multiplicative root and never runs a per-D Tonelli). Follows
classpoly's `mpz_cornacchia4` but is thread-safe (all scratch in the ctx, one ctx/thread).

**Descent is mpn + a small-quotient fast path**: `(2p, x0)` reduced by Euclid until `b≤L`;
most quotients are 1 or 2, done by `mpn_sub_n`, only q≥3 uses `mpn_tdiv_qr`; pointer-rotated,
no allocation. Root recovered from `(4p−b²)/d` being a perfect square (divisibility-by-`d`
is the cheap common-case filter). ~1.7–1.9× over the mpz version.

### The scan (`dscan.c`) — fastECPP substeps 1+2, residue-factor variant
Given `p` and bound `B`, enumerates **every negative fundamental `D` with `|D|<B` all of
whose prime-discriminant factors `q*` are QR mod p** and reports solvable traces `t`.

- **Factor base** = QR prime discriminants: `q* = q` if `q≡1 (4)`, `−q` if `q≡3 (4)`, plus
  `−4, ±8`, keeping those with `(q*/p)=1`. Crucial identity: `(q*/p) = (p mod q / q)`
  (small Legendre) — determination uses **Legendre symbols, never Tonelli**.
- **Multiplicative sqrt**: precompute `√(q*) mod p` for each factor-base element (the *only*
  Tonelli-Shanks). A DFS over products of distinct factor-base primes maintains
  `√(D mod p) = ∏ √(q*)` incrementally (one modmul per D); Cornacchia runs with that root.
- **Why residue-factor restriction (not F_p² for all D):** genus theory — if all `(qᵢ*/p)=1`
  then `p` splits in the genus field and `Pr[4p=t²+|D|v² solvable] = g/h = 2^(t−1)/h` (t =
  #prime-disc factors), a `2^t` boost over the generic `1/(2h)`. So all-QR-factor D are the
  *highest-yield* ones AND keep `√` in `F_p`. Observed yield ~10% at small B (vs ~0.6% for
  arbitrary D).

Validated by `test_dscan.py`: DFS enumeration == an independent linear brute force (0
missing/extra) and solvability == PARI's principal-form oracle (`qfbsolve`), both p mod 4.

### Parallelism (both substeps)
- **Substep 1 (factor-base build)**: sieve → parallel Legendre determination → parallel
  Tonelli roots (per-thread `cornacchia_ctx`). ~16× on the roots.
- **Substep 2 (scan)**: parallelized by **seed-based load balancing**. A naive parallel-for
  over the first prime plateaus at 3.6× (the smallest-prime subtree is ~27% of the work). Fix:
  a cheap serial `gen_seeds` recurses only into big subtrees and emits `~SEED_D`-sized
  independent "seed" chunks (batching the millions of tiny large-prime subtrees), then a flat
  `omp parallel for` over seeds (each a pure serial `dfs_serial`, no shared state, no task
  suspension). **13.5× @32 threads** (~85% of the 16-core ceiling; `SEED_D=2000` default, `sd=` to tune).

### Benchmarks (32 threads, B=1e8, seed=1)
| p bits | build (substep 1) | scan (substep 2) | scan µs/D | solvable rate |
|---|---|---|---|---|
| 128 | 0.61 s | 0.26 s | 0.049 (20 M D/s) | 1/601 |
| 256 | 2.10 s | 0.63 s | 0.102 (9.8 M D/s) | 1/523 |
| 384 | 3.10 s | 0.75 s | 0.176 (5.7 M D/s) | 1/708 |

Scan is ~24× faster than the original mpz-serial baseline (mpn descent ×1.7 · parallel ×13.5).

### Findings & open levers (for component a)
- **Substep 1 now dominates** (Tonelli roots ~π(B)/2). Cause: enumerating *all* D<B pushes
  the factor base up to B. **Bounding the factor base** (small primes + products, à la fastECPP)
  would slash it and focus the scan on high-yield D — **deferred by user until obviously needed.**
- Yield falls ~`1/√B` (small D far higher-yield; `Pr[solvable]≈2^(t−1)/√|D|`).
- **Lehmer/half-gcd does NOT help our descent**: it's short (avg steps `≈0.29·bits`, over a
  1–2 limb span), so there's no span to batch. Half-gcd is a big win only at their 1000s-of-bits
  `p`. Per-step small-operand division is our cost; mpn+fast-path is the right tool. A hand-rolled
  1-limb-quotient division could add ~2× per descent (fiddly, deferred).
- Scan scaling likely capped by memory bandwidth on the `mpz` sqrt array (~800 MB @B=1e8);
  **packing sqrts as contiguous limbs** would help — deferred.
- The "à la Bach" the user first mentioned is actually **Shallit** (the L⁴ fastECPP idea);
  the multiplicative-sqrt is credited to Crandall–Pomerance §8.4.3.
- Remainder-tree for `p mod q`: measured negligible at our small `p` (`p mod q` is 0.014 µs/prime
  @256-bit, ~9% of the already-cheap determine phase). Worth it only at much larger `p`. Deferred.

## Component (b): n⁴-smoothness testing — DONE (v1, validated)

Code: `ecpp/smooth.{c,h}` (engine) + `ecpp/smoothtest.c` (driver) + `ecpp/test_smooth.py`.

### The requirement changed: n⁴-smooth one-shot certificates
The challenge (github.com/AndrewVSutherland/OneShotPrimalityProofs, verifier `voneshot.py`)
now accepts an **n⁴-smooth** certificate `(p, A, x₀, m, q₁…q_k)`: `m = ord(P)` on a Montgomery
curve `E/F_p`, with
- **`m` is n⁴-smooth**, `n = bitlength(p)` (every prime factor of `m` is `≤ n⁴`);
- **`L < m < L·r`**, `L = ⌊√p⌋+1+⌊√(4⌊√p⌋)⌋ = (p^{1/4}+1)²`, `r` = least prime of `m`;
- `m ≤ p+1+⌊2√p⌋` (Hasse); `q₁<…<q_k` = the primes of `m` in `(n², n⁴)` (verifier trial-divides
  only to `n²` and needs the larger primes listed).

Since `m | #E = p+1∓t`, a valid `m` exists **iff the n⁴-smooth part `S` of `p+1−t` or `p+1+t`
exceeds `L`** (then `m` = smallest divisor of `S` above `L` lands in `(L,L·r)`; proof: `m/r ≤ L`
by minimality). **So component (b) = for each solvable `(t,v)`, compute the n⁴-smooth part of the
two candidate orders and test `S > L`.** This is a big relaxation vs. full smoothness (only ~half
the bits of `#E` need be smooth), which is why usable `D` are now far more common — but we run the
test on *every* candidate, so it must be fast.

### Sizing (n = bitlength(p), y = n⁴ = smoothness bound)
| p | n²=trial bd | y=n⁴ | ‖P‖=∏_{q≤y}q | need S > | u=logN/logy | est. yield/candidate |
|---|---|---|---|---|---|---|
| 128b | 2¹⁴ | 2²⁸ | 48 MiB | 2⁶⁴ | 4.6 | ~6% (measured 6.3%) |
| 256b | 2¹⁶ | 2³² | 738 MiB | 2¹²⁸ | 8 | **~1.3×10⁻³ (measured)** |
| 384b | 2¹⁷·² | 2³⁴·³ | 3.7 GiB | 2¹⁹² | 11.2 | **< 3×10⁻⁶ (measured, 0/326222)** |

Yield = Pr[n⁴-smooth part of a ~random `#E` exceeds `√p`]. Harder as `p` grows (`y=n⁴` is only
polynomial in `n`). So we need ~`1/yield` candidates → 128b: tens; 256b: ~10³–10⁴; 384b: **~2×10⁵**.
2 candidates per solvable `D`. At 128/256-bit the `dscan` scan is cheap (≪1 s) so **smoothness
testing is the cost**, as the user flagged. **At 384-bit the picture shifts**: the number of solvable
`D` up to `B` grows only like `~0.6·√B` (yield ∝ 1/√B — small `|D|` dominate), so ~10⁵ candidates
needs `B≈10¹⁰` and ~2×10⁵ needs `B≈2.5×10¹⁰` — a `dscan` whose factor-base build (π(B) Tonelli
roots, tens of GB of √-table) becomes the bottleneck. So at the frontier, *producing* candidates
(component a) costs as much as testing them; the deferred "bound the factor base / favor small |D|"
lever starts to matter. Smoothness testing itself stays cheap (measured 262 µs/candidate @384-bit).

### Algorithm: Bernstein batched remainder tree (chosen over Pollard–Strassen)
Per-candidate single-number smoothness (Pollard–Strassen with a wheel) costs `Õ(√y)` modmuls
≈ `√(2³²)=2¹⁶` × log ≈ 10⁶–10⁷ modmuls ≈ 0.1–1 s **each** — a non-starter at 256/384-bit scale.
Instead, **batch** with Bernstein ("How to find smooth parts of integers"):
1. Precompute the single big integer `P = ∏_{prime q≤y} q` (‖P‖ ≈ 1.44·y bits). Depends only on
   `n`, so **built once per bit-length and cached to disk** (`smooth_base_save/load`).
2. For a batch `N₀…N_{k-1}`, a **remainder tree** gives `rᵢ = P mod Nᵢ` for all `i`: build the
   product tree of the `Nᵢ` (root `X=∏Nᵢ`), reduce `P mod X` **once**, then descend.
3. The y-smooth part is `Sᵢ = gcd(Nᵢ, rᵢ^{2^s} mod Nᵢ)` with `2^s ≥ bitlen(Nᵢ)` (`s≈9`; ~9
   squarings + a gcd per candidate). **Exact & deterministic**: for a prime `q|Nᵢ`, `q|rᵢ ⇔ q|P
   ⇔ q≤y`, so the gcd captures exactly the primes `≤y` with full multiplicity.

Cost per candidate ≈ (one big `P mod X`, `Õ(‖P‖)`, amortized over the batch) + `Õ(bitlen·log k)`
descent + ~9 modmuls. The only step scaling with `y` is `P` (built once) and the top reduction
(per batch); everything else is quasi-linear in the number size. Crossover vs Pollard–Strassen is
~400 candidates; at 10³–10⁵ batches this is ~25–300× faster. **This is the same primitive the
verifier uses** (remainder tree + repeated-squaring `is_smooth`), so the generator mirrors it.

Why not a QS/NFS-style sieve over candidates: our solvable `t` are sparse and scattered over a
`~2^{129}`-wide range, so there's no contiguous array to step primes through. Small primes *could*
be stripped by the `t ≡ p+1 (mod q)` test, but the cost is dominated by the many primes near `y`,
which the remainder tree handles anyway — so a sieve pre-pass buys little. (Deferred.)

### Implementation (`smooth.{c,h}`)
- `smooth_base_build(y, nth)` — parallel segmented sieve (banded, chunked, `2` + odd wheel) →
  `uint64` prime array → parallel product tree (OpenMP tasks). `smooth_base_save/load` cache `P`.
- `smooth_parts(sb, N[], k, S[], nth)` — the remainder tree (parallel level-by-level build +
  descent) then parallel extraction. Reduces `P mod X` once at the root.
- `cert_bounds(p, …)` → `L`, Hasse, `n`, `n²`, `n⁴`. `factor_smooth` (trial + Brent-rho) and
  `build_m` (assemble `m∈(L,L·r)` largest-prime-first, list `q_i∈(n²,n⁴)`) for the rare winners.
- **Cache gotcha (fixed):** the `P` disk cache uses native limbs + a 64-bit count, NOT
  `mpz_out_raw` — the latter's **4-byte size field silently overflows for P > 2³¹ bytes** (i.e.
  p ≳ 384-bit), yielding a corrupt `P` that's missing primes. 256-bit P (0.8 GiB) is under the cap;
  384-bit P (3.7 GiB) is over it, so this bit the first 384-bit cache. Native-limb format has no cap.
- Driver `smoothtest`: `pbuild` (build/cache `P`), `parts` (stdin ints → smooth parts; unit test),
  `gate` (reads `dscan … dump`, forms `N=p+1∓t`, batch-gates, prints winners, self-checks each).

### Validation (`make -C ecpp && python3 ecpp/test_smooth.py`)
- **Test A**: batched smooth-part vs PARI factor-based smooth part over mixed random integers
  (y ∈ {10³,10⁶,2²⁸}, batches ≤ 400) — exact.
- **Test B (end-to-end, 128-bit)**: 872 candidates from `dscan`; every reported winner
  **independently verified with PARI** (`m|N`, n⁴-smooth, `L<m<L·r≤`Hasse, `q_i` exact) and the
  gate decision (`S>L`) matches PARI on the full candidate set. 55/55 winners OK, 0 mismatches.
- Disk-cache round-trip (fresh `P` vs loaded `P`) bit-identical. Every winner is also self-checked
  in C (`check_winner`, PARI-independent). Algorithm is size-independent (only `n²,n⁴` are `uint64`
  — fine to ≫1000-bit), so correctness carries to 256/384-bit.

### Measured performance (32 vCPU box, **heavily contended** during these runs — lower bounds)
- Build+cache `P` (y=2³², 203 280 221 primes → 738 MiB), 16 threads: **62 s one-time** (reused
  across all `p` of that bit-length; `mpz_out_raw` cache loads in 0.17 s). π(2³²) exact.
- 256-bit gate, cached `P`, `smooth_parts` on a batch:
  - 10 514 candidates: **3.5 s** (338 µs/candidate) — was 17.7 s before the parallel reduction.
  - 28 984 candidates: **4.6 s** (160 µs/candidate) — the top reduction is ~batch-independent, so
    per-candidate cost falls with batch size (amortization confirmed).
- 384-bit gate, cached `P` (3.7 GiB, loads 0.66 s), validated vs PARI. Amortization scales with batch:
  60 196 cand → 385 µs, 141 114 → 262 µs, **326 222 → 181 µs/candidate** (59 s). Over 326 222 real
  candidates: **0 winners ⇒ yield < 3×10⁻⁶**. So smoothness testing at 384-bit is cheap and correct,
  but **no certificate was found** — the frontier bottleneck is producing enough candidates: the
  B=2.5×10¹⁰ `dscan` needed to make 3×10⁵ of them peaked at **58 GiB** (its √-table), and can't push
  much past B≈4×10¹⁰ without OOM. Finding a 384-bit certificate needs the factor-base-bounding lever
  (or a lower-RAM scan), not faster smoothness testing.
- The single `P mod X` top reduction dominates a batch; it's now parallel (`reduce_big_mod`):
  limb-aligned zero-copy chunks reduced mod `X` concurrently, recombined with `2^{g·chunk} mod X`.
  ~5× here (contended); ~16× on an idle box. Extraction (≈9 squarings+gcd/candidate) is ~free.
- Practical upshot: finding a 256-bit certificate's `(D,t,m)` needs a batch of ~few×10³ candidates
  (yield 1.3×10⁻³) → **~1–2 s** of smoothness work after the one-time `P` build.

### Open levers (component b)
- **Prefer large batches** (amortize the ~batch-independent top reduction): one batch of all
  candidates for a given `p`; per-candidate cost then → the quasi-linear descent + squarings.
- **384-bit**: `P` = 3.7 GiB (955 576 801 primes; measured build 6:53 at 8 threads, contended, 28 GiB
  peak from the top FFT multiply). Fits in RAM, cached to disk; yield ~10⁻⁵ ⇒ ~10⁵-candidate batches
  where the parallel `reduce_big_mod` matters most. (Build peak RAM = the top product-tree multiply;
  a streaming "binary-counter" merge would cut it — deferred, fine to 384-bit.)
- Tune `y` below `n⁴` if the yield/‖P‖ trade favors it (smaller `P`, fewer qualifying `N`).
- Extra candidates from `D=−3` (6 twists) / `D=−4` (4 twists) — more orders per solvable `t`.
- Strip tiny primes with the `t ≡ p+1 (mod q)` sieve test only if profiling shows it helps.
- Squeeze the last of the 16× on the top reduction (serial `pw`-chain + recombine are small but
  nonzero); revisit on an idle box.

### Scaling ceiling (for "eventually larger" p)
`‖P‖ ≈ 1.44·n⁴` bits, so `P` in RAM caps out around **~900-bit p** (n⁴≈2³⁸, P≈40 GiB); at 1024-bit
n⁴≈2⁴⁰ ⇒ P≈200 GiB > 128 GiB. Past that, `reduce_big_mod` already shows the fix: it consumes `P` as
independent limb-chunks, so a **disk-resident `P` streamed in chunks** (accumulating `Σ Bⱼ·(2^{jw}
mod X)`) gives the same result with O(‖X‖) RAM. Also, as p grows the yield falls (n⁴ only polynomial
in n), so we need more candidates — which stresses `dscan` too (its factor-base build ~π(B)); the
deferred "bound the factor base / prefer small |D|" lever then starts to matter. 128–384-bit (the
current target) is comfortable: P ≤ ~4 GiB, fully in RAM.

## Component (c): root of H_D over F_p — DONE (core, validated)

Code: `ecpp/fproot.{c,h}` (engine) + `ecpp/roottest.c` (validation) + `ecpp/hdroot.c` (end-to-end).

`H_D` is degree `d = h(D)` over `F_p`, `p` few-hundred bits (beyond `ff_poly`'s word size). We need
**one** root `j₀` (a valid j-invariant with CM by D). Because `D` was chosen with `4p = t²+|D|v²`
solvable, `p` splits in the ring class field and **`H_D` splits completely over `F_p`** — so there is
no distinct-degree step; we go straight to **equal-degree splitting** for linear factors.

### The inner loop (per the user's recipe)
Pick random `a`, form `b = (x+a)^{(p−1)/2} mod h`, split by `g = gcd(h, b−1)`, keep the
smaller-degree side, recurse to degree 1 → root `= −h₀`.
- **`(x+a)·f mod h`** (h = xᵈ + h′ monic): a shift, a scalar multiply, and subtracting `lead·h` —
  **no polynomial multiplication**. `O(d)` per step; negligible.
- **`f² mod h`** dominates. `f²` is formed by **Kronecker substitution** (pack coefficients into one
  big integer at limb-aligned slots, square with GMP `mpn_sqr`, unpack), then each coefficient is
  reduced mod p by **Montgomery REDC**, then the degree-(2d−2) result is reduced mod h.
- **All scalar F_p arithmetic is Montgomery** (multi-limb REDC), coefficients kept in Montgomery
  form throughout. The Montgomery limb count `s` is padded (`64s ≥ n+34`) so a Kronecker product
  coefficient (a sum of ≤ d products, `< d·p²`) still lies below `p·R` and REDC returns it directly.
- **mod-h reduction is Barrett/Newton** (not schoolbook): precompute `rev_d(h)⁻¹ mod x^{d−1}` once
  per h (Newton), then reduce with two short Kronecker products. This makes the squaring the
  dominant cost (as intended); schoolbook is kept for `d < 48` (faster there) and as a reference.

### Getting H_D mod p (classpoly)
`compute_classpoly(D, inv=0, p, file)` (the validated Test-3 ECRT path) writes monic `H_D mod p`
in the **j-invariant** (`inv=0`), degree `h(D)`, coeffs in `[0,p)`, low-degree-first. `hdroot`
invokes the `classpoly` binary, parses the file, loads into the Montgomery poly, and finds a root.
(For smaller class invariants — Weber, double-eta — `H_D` is smaller but a big-`F_p` invariant→j
map would be needed; j-invariant first keeps the root == j₀ directly. Deferred.)

### Validation & performance (`make -C ecpp && ecpp/roottest`)
- **`roottest`**: builds `h = ∏(x−rᵢ)` from known roots, finds a root, checks membership + an
  independent Horner root check; a `pari` mode cross-checks the *entire* root set vs `polrootsmod`.
  Green over p ∈ {128,160,256,384}-bit, d ∈ {5…640}, the schoolbook/Barrett boundary, 20-trial
  sweeps. Montgomery layer separately checked vs GMP (200k random mul/add/sub/inv).
- **Barrett speedup**: d=640 root-find dropped 5.1→0.81 s (256-bit) / 12.6→1.85 s (384-bit), ~6–7×;
  scaling is now near-linear in d.
- **End-to-end (`hdroot`, 256-bit p, real class polynomials)**: D=−24447 (h=92) root in 0.069 s,
  D=−1274827 (h=102) in 0.063 s, D=−103044 (h=104) in 0.072 s — every root verified `H_D(j₀)≡0`.
  (`H_D mod p` itself via classpoly is 0.05–0.7 s here, depending on |D|.)

### Full pipeline verified end-to-end (certificate assembly prototyped in PARI)
From `j₀`: build `E/F_p` with that j-invariant, take the twist of order `N=p+1∓t`, put it in
Montgomery form (2-torsion root `α` with `3α²+a₄` a QR ⇒ `s=1/√(3α²+a₄)`, `A=3αs`), find a point of
order `m` (random `Q`, `[N/m]Q`, check order), and emit `(p, A, x₀=s(x_P−α), m, q_i)`. Prototyped in
PARI over three 256-bit smoothness winners, feeding **j₀ straight from the C root-finder** →
**`voneshot.py` returns True** (p proven prime) for D=−15607 (h=39) and D=−2390772 (**h=528**).

**Two compatibility filters the assembly needs** (both number-theoretic, characterized on the 13
winners): the Montgomery model exists iff **N ≡ 0 mod 4** (6/13 winners; the N≡2-mod-4 ones have no
2-torsion QR); and a point of order `m` exists iff **m | exponent(E)** (fails when `E` is
`Z/2 × Z/(N/2)` and `m`'s 2-part = v₂(N)). So a Montgomery one-shot cert wants `N≡0 mod 4` and an `m`
whose 2-part ≤ v₂(N)−1 — cheap extra filters for component (b)/the assembler. This validates the whole
chain **a→b→c→certificate**; the remaining work is a C implementation of the assembly (Montgomery
x-only ladder over big `F_p`, reusing `fp_ctx`) for speed — the math is confirmed correct.

## classpoly integration & validation

Builds/runs entirely in-tree (no `/usr/local`). `make test` runs three suites (green over
|D|≤1000):
- **Test 1** (`compare_pari.py`): classpoly vs PARI `polclass` over Z for every shared invariant
  — 3037/3037 exact.
- **Test 2** (`verify_jinv.py`): PARI-unsupported invariants (Atkin A₄₁/₄₇/₅₉/₇₁, Ramanujan t,
  single-eta, w₃,₁₃) — map class-poly roots mod a split prime to j via `invtoj`, check they cover
  the Hilbert roots — 1450/1450.
- **Test 3** (`compare_modp.py`): `H_D mod (2²⁵⁵−19)` == `(H_D over Z) mod (2²⁵⁵−19)` — 3792/3792.
  This validates the large-modulus ECRT path that component (c) will consume.

### classpoly changes we made
- Env-var directory hooks: `CLASSPOLY_PHI_DIR`, `CLASSPOLY_H_DIR`, and `CLASSPOLY_TMPDIR` (now a
  *base* dir; classpoly makes a private per-process `mkdtemp` subdir under it, default `/tmp`).
- **Parallel-safety fix**: CRT scratch was deterministic per-computation → concurrent runs in a
  shared dir clobbered each other (SIGABRT). Now self-isolating (per-process /tmp dir, auto-cleaned).
- Added `invtoj` binary + `ff_all_j_from_inv()` (invariant→j over F_p, power-aware, non-aborting).
- Code-review fixes: an 864-byte-per-mod-p-run ECRT leak (`ecrt_clear` never called + missing
  `free(ecrt->m)`); a 5×`mpz_t` leak + a stray cwd-writing "hack" in `classpoly_inv_setup`;
  `bach_gcd` `abs`→`labs` (long truncation); a missing `break` in `classpoly_load`.
- Invariant codes match PARI's `polclass` except w₅,₇ (PARI 35 / cp 535) and w₃,₁₃ (39 / 539).

## Gotchas / conventions

- **PARI/GP via subprocess**: pipe the script on **stdin** (`gp -q < script`), NOT `gp file.gp`
  (that drops into interactive mode and hangs). One statement per line (gp reads line-by-line).
  `polclass` needs a big stack: `default(parisizemax,"8G")`.
- **classpoly not parallel-safe across a shared cwd/TMPDIR** unless self-isolating (now fixed);
  any harness running many instances still runs each in its own dir out of habit.
- **Build with `-Wall -Wextra`** (the vendored makefiles didn't; that's how we found the bugs).
  `-march=native` is fine (fixed box).
- Fundamental discriminant `d=|D|`: `d≡3 mod4` squarefree, OR `d≡0 mod4` with `d/4≡1,2 mod4`
  squarefree. Prime-discriminant 2-part `∈ {1,−4,8,−8}`.

## Component (d): best invariant, invariant→j, curve assembly, and `oneshot` — DONE

Code: `ecpp/cminv.{c,h}` (invariant→j) + `ecpp/cm_method.c` + `ecpp/curve.{c,h}` (assembly) +
`ecpp/oneshot.c`; classpoly's `class_inv_mpz.c` (Sutherland) is used for the hardwired formulas.

**Best invariant, not j.** Forcing the Hilbert (j) class polynomial is wasteful: the best class
invariant (Weber, Atkin, double-η) is **25–50× faster** to compute (measured: h=528, 1.53→0.06 s;
h=609, 4.03→0.08 s) and often halves the degree. classpoly picks it with `inv=-1`.

**Invariant→j over big F_p** (`cm_j_from_inv`): the value `f₀` maps to `j` = a root of `Φ_inv(f₀^e, ·)`
with a per-invariant power `e`. **Now that the `zp_poly` library is in the tree (see below), the
production path links classpoly's authoritative `mpz_j_from_inv` (`class_inv_mpz.c`) directly**; a
self-contained `fproot`-based implementation (`cm_j_from_inv_ref`, Sutherland's formulas re-ported +
`Φ_inv` root-find via `invj.{c,h}`) is retained as a cross-validation reference and covers the generic
single-η range (400–499) that `class_inv_mpz.c` does not dispatch. **Three-way validated** — for 22
invariant families × random `f₀` at word size, `mpz_j_from_inv` ∈ `invtoj`'s root set and the
reference's root set == `invtoj`'s (81/81). Found & fixed in `class_inv_mpz.c`: `mpz_j_from_u8`
returned `J` unreduced (≈7·|p| bits) and never `mpz_init`'d `T3` — flagged for upstream.

**zp_poly** (Harvey–Sutherland large-`p` `F_p[x]` library, added to the tree with an in-tree makefile
staging into `./local` like ff_poly): provides the `zp_poly_find_root`/`bipoly_eval_mod_mpz` machinery
`class_inv_mpz.c` needs. Its `zp_poly_find_split_root` is the same EDS algorithm as `fproot`
(plus radical shortcuts for deg ≤ 3). Head-to-head on identical split polys (`ecpp/zpbench`):
zp_poly ~5–25% faster at 128/256-bit, parity at 384-bit (fproot slightly ahead at deg 800) — so
**`fproot` is kept for the `H_D` root-find** (equal-or-better where it matters, Montgomery-form
integration with `curve.c`), and zp_poly serves the invariant→j path.

**`cm_method D p`** → the j-invariant of `E/F_p` with CM by `D` and trace ±t: Cornacchia → best-inv
`H_D^inv` → `fproot` root f₀ → `cm_j_from_inv` → verify a candidate's curve has trace ±t (disambiguates
the ≤2 j-roots and self-checks). Every output verified to be a genuine Hilbert `H_D` root.

**`mont_assemble` (`curve.{c,h}`)**: j₀ → Montgomery `A` (solve `256(A²−3)³=j₀(A²−4)` for `A²`,
`A=√u`); Montgomery x-only ladder; determine which of `E_A`/twist has order `N`; find `x₀` = x-coord
of `[N/m]·`(random point on that curve) with order exactly `m`. Produces the `(A, x₀)` voneshot wants.

**`oneshot p`** ties it together (all C; shells to `dscan`, `cm_method`, `classpoly`): scan → keep
`N≡0 mod 4` orders → batch n⁴-smoothness → for a winner, `cm_method` gives j₀ and `mont_assemble`
builds `(A, x₀)`; emit `(p, A, x₀, m, q_i)`. Verified end-to-end by **`voneshot.py`** for 2²⁵⁵−19
(~6 s with cached P) and fresh random primes. Montgomery needs `N≡0 mod 4` and `m|exponent(E)`;
`oneshot` filters the former and skips winners failing the latter.

## Roadmap
1. ~~Cornacchia + discriminant scan (a)~~ — DONE, parallel, validated.
2. ~~n⁴-smoothness testing (b)~~ — DONE: Bernstein batched remainder tree (`smooth.{c,h}`),
   parallel, cached `P`, validated vs PARI (128/256-bit end-to-end, all winners checked).
3. ~~Root-finding of large `H_D` over `F_p` (c)~~ — DONE: Montgomery + Kronecker + Barrett EDS
   (`fproot.{c,h}`), validated vs PARI, wired to classpoly end-to-end.
4. ~~Wire it together (d)~~ — DONE: best invariant + invariant→j (`cm_method`) + Montgomery curve
   assembly (`mont_assemble`) + `oneshot`; certificates verified by `voneshot.py`. `git clone && make
   && . ./setenv.sh && ./ecpp/oneshot p=<prime>`.
5. Deferred optimizations (do when obviously needed): bound the factor base; pack factor-base
   sqrts as contiguous limbs; hand-rolled 1-limb-quotient division; SIMD the per-D updates.

### Prime-product cache ladder + near-miss top-up (2026-07-01, post-merge)
`oneshot` no longer builds one exact-n⁴ `P` per bit-length. Policy: use any cached power-of-2
product `2^j ≥ n⁴` (oversize primes are skipped by `build_m` when assembling `m`), else build and
cache `y' = 2^⌊log₂ n⁴⌋ ≤ n⁴` — at most the exact cost, one file per octave of `n` (e.g. `P_2³²`
serves every 227–304-bit prime). The gap `(y', n⁴]` is recovered by **near-miss top-up**: candidates
with `S ∈ (L/n⁴², L]` get bounded Pollard rho on the cofactor `N/S` (a missing prime is ≤ n⁴ ≈ 2³⁴
⇒ ~2¹⁷ mulmods), folding gap primes into `S` — no second prime product needed. Measured (10⁹⁰+289,
loaded box): P 199 s → 0.3 s (reused `P_2³²`), 1236 near-misses rechecked, cert = the SAME winner as
the exact-y run, voneshot True. Caches live in `$ONESHOT_PCACHE_DIR` (setenv.sh → `work/pcache`,
survives reboots); sizes ~0.172·y bytes: 2³² → 739 MiB, 2³³ → 1.4 GiB, 2³⁴ → 2.9 GiB.
Motivation: P build was 89%/74%/22% of the fresh-bit-length 10⁸⁰/10⁹⁰/10¹⁰⁰ runs — far above the
"never spend more building P than on H_D + root-finding" balance point; with the ladder the
amortized cost → 0. (10¹⁰⁰ breakdown, original run: P 247.5 s, dscan 597 s [~75% factor-base
build/Tonelli, ~25% DFS scan], gate 26 s, H_D^f deg 35084 ≈ 80 s, root-find ≈ 170 s, assembly ~10 s.)

### Adaptive smoothness ladder (2026-07-01, replaces the fixed-P flow in oneshot)
Per Drew's proposal: instead of building `P = ∏_{q≤n⁴} q` up front, `oneshot` climbs a power-of-2
ladder of prime-product **segments** `(y_{j-1}, y_j]` starting just above `n²` (final rung capped at
`n⁴` exactly), keeping a running smooth part `S = ∏` per-segment parts for every candidate. The
ladder deepens only when cumulative testing work ≥ `c`·(P bits so far + next segment), `c=1` default
(work-normalized version of his bits-tested trigger); otherwise the pool widens via incremental
`dscan Bmin=` chunks (both growths geometric ⇒ within a small constant of hindsight-optimal). A
winner at any rung stops the run — and certificates rarely need primes near `n⁴` (measured on a
256-bit pool: y=2²⁸→3, 2³⁰→6, 2³²→13 winners), so runs stop octaves early. Segments are
bit-length independent and cached individually (`oneshot_Pseg_<lo>_<hi>.bin`); legacy full-P caches
are used as the ladder base when present. The near-miss rho top-up is no longer needed in oneshot
(the ladder reaches n⁴ exactly); `smooth_topup` remains available in smooth.{c,h}.
**Measured cold runs** (no caches, loaded box): 256-bit 2.7 s (was ~65 s), 2²⁵⁵−19 quick-start 3.5 s,
10⁸⁰+129 4.5 s (was 160 s), 10⁹⁰+289 23.6 s (was 269 s) — the 10⁹⁰ run stopped at y=2²⁸ with a
14k pool instead of ever touching the 2³¹⁺ rungs. Warm ≈ cold now ("first run as fast as the nth").

### Parallel classpoly (multi-job ECRT) + sub-quadratic gcd (2026-07-01)
Two improvements motivated by the 10¹⁰⁰+267 winner (D=−2557415807, **h=35085**):

**Parallel H_D across CRT primes.** classpoly is single-threaded and ff_poly is not thread-safe, but
the computation is embarrassingly parallel across the CRT primes, and Drew's ECRT layer already had
distributed-job machinery (`ecrt_init(jobs, jobid)`, `ecrt_dump`, `ecrt_merge` — the accumulators are
pure sums). We wired it up process-parallel: `CLASSPOLY_JOBS=W` + `CLASSPOLY_JOBID=j` makes a classpoly
run a worker (prime indices ≡ j−1 mod W, dumps partial coefficient sums to `$CLASSPOLY_ECRT_DIR`);
`CLASSPOLY_JOBID=0` merges the dumps and writes the polynomial. `cm_method jobs=N` orchestrates
workers+merge (gated on |D| ≥ 10⁸ — for small h the repeated per-worker setup outweighs the win);
`oneshot` passes its thread count through. **h=35085: 158 s → 28.9 s (16 workers), byte-identical
output**; full `cm_method` on that winner 634 s → 188 s (both under external load). The residual gap
to W× is the per-worker setup (class group, phi loads, ECRT precomputation), which is repeated.

**Sub-quadratic gcd in the root-finder.** `fproot`'s schoolbook gcd/division is O(d²) and dominates
EDS at large degree. zp_poly's `zp_poly_xgcd` is a genuine half-gcd (quad-matrix recursion), but its
modmul engine is slower than fproot's Montgomery/Kronecker/Barrett at our sizes — full delegation to
`zp_poly_find_split_root` measured ~parity. The hybrid wins: keep fproot's exponentiation, delegate
gcd + exact division to zp_poly above degree 1024 (`fpoly_gcd_fast`/`fpoly_divexact_fast`).
Measured: d=8000 31.3→25.7 s, d=16000 81→66 s, d=32000 152→127 s (~16–19%; the powmod dominates
what remains). Validated vs PARI `polrootsmod` across the threshold; suites green.
