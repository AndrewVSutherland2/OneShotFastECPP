# CLAUDE.md — OneShotFastECPP

CM-method ("fast ECPP") approach to one-shot elliptic-curve primality proofs.
**Read `design.md` for the full technical picture, decisions, and roadmap.**

> **CURRENT MAIN TASK (2026-08-19): the ecpp-varieties write-up.**  Drew has a
> detailed review from GPT 5.6 Sol to work through.  Read
> **`reports/ecpp-varieties/HANDOFF.md`** first — it has the full context:
> the deliverable (`ecpp-varieties.tex`, sole maintained source; PDF committed
> alongside), house conventions, the corrections already made (do not
> regress), where every number is reproduced from, sources, and pending items.

## Status (2026-07-02) — COMPLETE, end to end
- **`oneshotECPP p=<prime>` → n⁴-smooth one-shot ECPP certificate**, verified by the challenge's
  `voneshot.py`. Cold runs (fresh pcache): 10⁶⁰+7 0.8 s, 10⁸⁰+129 5.5 s, 10⁹⁰+289 23.6 s,
  10¹⁰⁰+267 7.9 min (B=4×10⁹ scan + a degree-35085 H_D dominate). Verified certs in `certs/`
  (2²⁵⁵−19, NIST P-256, secp256k1, 10⁶⁰…10¹⁰⁰+ε).
- **Adaptive smoothness ladder** (Drew's design): no upfront `P=∏_{q≤n⁴}q` build — prime-product
  *segments* `(y_{j−1},y_j]` climb from just above n² (final rung capped at n⁴), candidates keep a
  running smooth part across rungs, pool widens via incremental `dscan Bmin=` chunks; deepen-vs-widen
  decided by a work trigger (`c=1` default). Winners usually appear octaves below n⁴. Segments are
  bit-length independent, cached in `$ONESHOT_PCACHE_DIR` (setenv.sh → `work/pcache`).
- **Component (a) discriminant search**: `dscan` — Cornacchia + QR-factor-base scan, parallel,
  validated vs PARI + brute force. ~10 M D/s @256-bit. `Bmin=` for incremental chunks.
- **Component (b) smoothness**: `smooth.{c,h}` — Bernstein batched remainder tree, ranged segment
  products, parallel top reduction; winner → `m∈(L,L·r)` + `q_i∈(n²,n⁴)` via `build_m`.
  Validated vs PARI end-to-end (every winner independently checked).
- **Component (c) root of `H_D` over big `F_p`**: `fproot.{c,h}` — H_D splits completely, so straight
  to equal-degree splitting; Montgomery arithmetic, Kronecker squaring, Barrett/Newton reduction,
  cheap `(x+a)·f`. Validated vs PARI `polrootsmod` (deg ≤ 640, up to 384-bit; deg-35085 in production).
- **Component (d) CM engine + assembly**: `cm_method D p` → j-invariant with CM by D (best class
  invariant via classpoly `inv=-1`, 25–50× faster than j; invariant→j via classpoly's authoritative
  `mpz_j_from_inv` = `class_inv_mpz.c` linked against **zp_poly** (both in-tree now); fproot-based
  `cm_j_from_inv_ref` kept as cross-check + for the generic single-η range; 3-way validated 81/81
  vs `invtoj`). `mont_assemble` (`curve.{c,h}`) builds Montgomery `(A,x₀)` with a point of order `m`
  (needs `N≡0 mod4` + `m|exponent`).
- **Isogeny-volcano descent** (`cm_method ells=2,7,...`): for each prime ℓ | n₁ (≤ 97), walk the
  ℓ-volcano to its floor (3 non-backtracking walkers on classical `Φ_ℓ` roots, cap v_ℓ(v)+2) —
  the floor curve has cyclic ℓ-Sylow, so it is always Montgomery-representable (kills the
  p≡3 mod 4, N≡4 mod 8 obstruction — no intake filter) and m may use the full ℓ-power of N
  (oneshotECPP reduces S only by n₁'s >97 part). Validated vs PARI (floor/trace/group/representability).
- **phi bundle**: 46 MB subset committed via git add -f (`phi_files_manifest.txt`): every family's
  `Φ_{ℓ,f}` for ℓ≤71 + all `Ψ_f` maps. Classical `Φ_ℓ` are needed only at the chosen invariant's
  LEVEL primes (e.g. `phi_j_59` for A₅₉) + `Φ_2`, never for enumeration. Full 2.2 GB DB external;
  misses degrade gracefully (skip winner). Fresh-clone validated (make + tests + oneshotECPP).
- Fixed upstream-worthy bugs: `class_inv_mpz.c` `mpz_j_from_u8` (unreduced J + uninit T3);
  vendored-header warnings (qform/prime/mpzutil/iqclass).
- **Parallel classpoly**: `CLASSPOLY_JOBS`/`CLASSPOLY_JOBID`/`CLASSPOLY_ECRT_DIR` env → multi-job
  ECRT (workers split CRT primes mod W, `ecrt_dump` partial sums; JOBID=0 merges). `cm_method jobs=N`
  orchestrates (auto-gated |D|≥1e8); oneshotECPP passes threads. h=35085: 158→29 s, byte-identical.
- **Hybrid root-finder**: fproot keeps its Montgomery/Kronecker/Barrett powmod but delegates gcd +
  exact division to zp_poly's half-gcd above degree 1024 (~16–19% at d≥8000; full zp delegation
  measured parity — its modmul is slower, its gcd faster).
- **Parallel root-finding (OpenMP)**: the EDS squarings are sequentially dependent, so the
  parallelism lives inside each product — task-parallel Karatsuba over polynomial halves
  (`fpoly_mul_par`, 3^depth leaf Kronecker products, depth ≤ 3 → up to 27 tasks) for the
  squarings and Barrett's two big multiplications, above degree ~4096.  ff_poly is not involved
  in this step, so OpenMP is safe (per Drew).
- **Big-machine scaling (2026-07-03)**: oneshotECPP's threads default is `omp_get_max_threads()`
  and is passed through to the `dscan`/`cm_method jobs=` shell-outs (was: hardcoded 16 —
  capped a 256-core production box at 16 cores for the scan phases). dscan's factor-base
  sieve reuses smooth.c's parallel `sieve_primes_range` (now exported; dscan links smooth.o)
  and the QR survivor collection fills slots via a parallel prefix sum; the prime array is
  64-bit — the old `uint32_t` silently truncated once B > 2³² ≈ 4.3e9 (reachable: Bmax
  defaults to 2e10). Validated: dump/dumpscan sorted-identical pre/post (5 configs, B ≤ 1e8),
  test_dscan.py + test_smooth.py green, (2³², 2³²+6e4) window exact vs PARI (enumeration,
  qfbsolve solvability, 4p = t²+dv² arithmetic), 256-bit oneshotECPP cert re-verified in PARI.
  Known remaining ceilings at 256 cores: classpoly jobs setup overhead, fproot's ~27-task cap.
- **Smoothness at scale (2026-07-04)**: `smooth_parts_multi` fuses a batch against ALL ladder
  rungs — one product tree, one shared chunk-power table (parallel prefix doubling; the old
  `reduce_big_mod` recombination ladder was G *serial* |X|-sized mulmods ≈ 90 s single-core per
  big batch at 256 threads, growing with thread count), one flat (segment×chunk) job list, one
  descent (gcd vs ∏P_s = ∏ per-seg parts; segments disjoint). oneshotECPP overlaps norm-equation
  scanning with smoothness testing: dscan is a forked child + drain pthread (parses dump,
  precomputes n1/n1h off-thread), next window spawned speculatively before test_batch, killed
  if a cert lands first. Phase instrumentation on stderr ([smooth segs], [scan join ... stall],
  [cm_method], dscan one-line summary in dump mode). Validated: fused vs per-candidate
  reference exact, gates green, 256-bit cert byte-identical (19.7→14.9 s), 10^110+7 B=8e9:
  identical pool (39878) + outcome, scan+smooth critical path 855→390 s @32 vCPU.
  Cold-start ceiling left: segment-build product-tree top muls (355 s/9 rungs at 10^110 here,
  once per rung then pcached) — multi-part segment forests if it matters.
- **Bmax default is memory-aware (2026-07-04)**: largest doubling of 2e10 whose dscan factor
  base fits ~30% of physical RAM (~67 B/prime at 383 bits: 128 GB → 2e10, 1 TB → 8e10,
  1.5 TB → 1.6e11); `B=` still overrides. The y=n⁴ exhaust message now suggests the next B
  and the ~GB-per-1e9-of-B factor-base cost. (10^115+79 exhausted the old fixed 2e10 cap.)
- **Descent-chain prototype (2026-07-22)**: `ecpp/chain_prove.py` + `ecpp/vchain.py` — Drew's
  b-descent variant: chain of GK certificates, bits(q) ≤ b·bits(p) per step (b=0.85 default),
  cofactor peeled by trial division + gmp-ecm rounds over the dscan pool; heuristically
  L_p(1/3, (2/3)(18(1-b))^(1/3)) vs the one-shot's p^(1/4+o(1)), verification
  O(n² log n)·1/(1-b²) (Pomerance-order, log log better than one-shot). vchain: stdlib-only,
  gcd-guarded affine ladder (sound for composite p), MR-12-bases base case < 2^64. Certs in
  `certs/chain/` (256-bit ~3 s vs one-shot ~15 s; 10^100+267 3.7 s vs 7.9 min; 1024-bit 128 s /
  9 steps / verify 0.4 s — far beyond one-shot reach). All steps PARI-cross-checked; tamper
  tests rejected. Analysis + results: `reports/descent-chain/index.html`. Needs a gmp-ecm
  binary (`ecm=` arg or CHAIN_ECM env; built from release tarball, no autotools on dev box).
- **SEA supply + crossover race (2026-07-29)**: `ecpp/sea_supply.py` — random Montgomery curves
  + PARI ellap as the candidate-order supply (heuristically p^(1/8+o(1)): flat poly cost per
  candidate, no quadratic scan, no winner classpoly), same `smoothtest gate` n⁴ window; winners
  assembled to voneshot certs with verify() as acceptance oracle (needs `voneshot.py` at repo
  root). Raced vs oneshotECPP on a 192-core spot box, same prime each size: 256: CM 4.1 s / SEA
  68 s; 384: 2306/615; 416: 1774/2805; 448: 7182 / stopped 20035 s no winner. Per-candidate:
  SEA flat 7–14 core-s, CM avg 2.2–4.5 rising — expected-cost crossover extrapolates to
  ~480–512 bits (walls flip on density luck from 384 up; observed fudge over ρ(u/2): 1–8×).
  No SEA-side accelerations implemented (stage-1 gate/batch tail need ellsea internals).
  Report: `reports/sea-crossover/`; certs in `certs/sea/`; raw cands in `work/race-archive/`.

- **Short ECPP records (2026-08-18, frontier 2026-08-28)**: `ecpp/shortECPP.py` — CM prover
  for the ShortPrimalityProofs *short ECPP* format; produced the table entries for
  nextprime(10^c), c = 210..310 (~3,500 core-h, one day) and the frontier extension
  c = 320..400 (largest: 10^400+69, 1329 bits; ~28,000 core-h, ~$420 of spot across six
  preemptions; 10^390 was a ×6 outlier at ~13,700 core-h over five ladder passes).  Frontier
  fixes: exact-order cofactor sampling needs the ℓ-volcano floor for EVERY ℓ | gcd(m, v) —
  build_level now retries via cm_method ells=<those ℓ> (was: ells=2 only at p≡3 mod 4);
  get_j retries rootless runs and logs cm_method failures; opt-in SHORTECPP_DEEP=1 adds a
  B1=3e7 tail tier (unused insurance so far); ECM-only idle padding (P-1 is deterministic).  Per level: dscan (new `maxfb=` bounds the factor base —
  essential for 1e5-candidate pools at 1000+ bits) → `smoothtest parts` batch strip (2^32
  primorial, built on demand) → staged gmp-ecm under a Bayesian index scheduler (coverage
  ~ln s, failure discounts, round-robin sub-sweeps, economic widen trigger) → classpoly/
  cm_method winner build (2-volcano floor when p ≡ 3 mod 4 needs it) → short.gp tail below
  135 bits; winner journal + `resume=1` survive preemptions; dual verification before write.
  Campaign report: `reports/short-ecpp-records/`.  Certs: `certs/short/`.  Gotcha fixed en
  route: subprocesses need classpoly on PATH (ECM_ENV handles it).  SEA/CM crossover for
  this criterion extrapolates to ~2000 bits (see report).

- **Short ECPP — the revised (August 2026) format**: revised with Drew over 2026-08-21/22:
  smoothness bound B = ⌈n²/log₂n⌉; radical cap log₂ rad(m_i) ≤ ⌈n/log₂n⌉ (exact as a
  bit-length test); p₁..p_{k+1} B-rough with recursion floor B² < p_{i+1} < √p_i and a
  self-certifying terminal prime p_{k+1} < B² (rough below B² ⟹ 1 or prime); payload
  unchanged.  Verification O(n² log n) worst case, self-contained (primorial-gcd recovery —
  gcd(P_B mod o, o) = rad(m)), ~5× Pomerance.  Tools (all on the settled convention):
  `ecpp/vshort2.py` (working verifier), `short2.gp`, `repair_short2.py`, `shortECPP.py
  v2=1` (gated on vshort2), `assemble_short2.py` (zero-arg = validate the tracked table).
  Full 31-chain table migrated: repairs dual-verified (revised + original), then forced
  terminal-prime truncation — the shipped `certs/short2/certs.csv` passes the revised
  verifier only (21/31 end in a terminal prime, intentionally original-invalid).  Existence swept for all p ≤ 4096.  Radical-capped stage merged
  upstream (ShortPrimalityProofs cc2caf1); terminal-prime + settled conventions on the
  AndrewVSutherland2 fork main (verifier renamed vshortECPP.py), upstream PR pending (Drew).
  Spec note: `reports/short2-spec/`; paper: Def 6.2 + Prop 6.4 in `reports/ecpp-varieties/`.

- **Short ECPP prover revisit (2026-11-14, branch `shortecpp-improve`)**: supply and round
  handling reworked, scheduler untouched.  dscan's factor base is capped from the first scan
  (`fb_cap(B)` = B/1000 in [1e6, 1e8]; `fbmin=` lets the cap grow across widenings without
  repeats; `roots=` caches the square roots per modulus, self-checking) — 1030 bits, B = 2e10:
  61,315 D for 1.9 core-h vs 85,954 for ~210 uncapped.  Strip bound 2^32 -> 2^26 (10x faster,
  same n²-smooth parts); `smoothtest parts n2= prp=1` returns the n²-smooth part and the tail's
  Fermat verdict, so the Python classify pool and `n2_parts` are gone; ECM rounds handle
  results as they complete and stop at the first winner (`EcmRunner`), several tails per
  gmp-ecm process at shallow tiers.  `ecpp/shortsim.py` = parameter-free model of one level
  (Poisson medium primes + Buchstab cofactor + measured/Dickman ECM success, exact
  compound-Poisson DP); it reproduces the 1030-bit campaign log tier by tier and says a
  per-candidate stopping policy (most orders get no ECM, the best few go to B1 = 3e6+) cuts
  ECM cost per winner several-fold.  Tested with real gmp-ecm on 249,648 orders of
  nextprime(2^768): 17 winners for 26.9 core-h of ECM (model: 24 for 32; 20.9 with the
  measured ECM success table `reports/shortecpp-improve/data/pi_calib.json`) — 1.6 core-h per
  winner; the prover with its old scheduler, same prime and box, needed 59 core-h of search
  for its first winner (148 for the 765-bit campaign prime).  (In the prover since
  2026-11-15, see the stopping-rule entry below.)  Consequence: construction is now the bottleneck, and it is root-finding —
  `cm_method` at 769 bits, 192 cores: h = 11,364 (inv 159) class poly 0.5 s + root 110 s;
  degree 58,520: 3.4 s + 364 s; h = 149,392 (inv 26): 525 s + 816 s.  The class-group
  decomposition of arXiv:1009.1082 is not in the tree (classpoly has only the word-size
  consumer `ff_classpoly_root` / `tower_split`); it is the next thing to build.  Gotcha found
  on the way: the old "classify" phase was 24-round Miller-Rabin in Python on every prime
  tail, not the process pool — now only tails below sqrt p get it.
  Engine benchmark (c8a.48xlarge): Zen-targeted GMP + tuned gmp-ecm gains <= 3% with all
  cores busy (GMP 6.3 configure falls back to k8 on Zen 5); standalone avx-ecm's IFMA path
  gives wrong results, its non-IFMA build is correct and 1.8x faster in stage 1 (1.35-1.5x per
  curve with gmp-ecm stage 2); yafu's IFMA ECM not built.  Report: `reports/shortecpp-improve/`.

- **Class-group decomposition (2026-11-14, same branch; arXiv:1009.1082 Algorithm 1)**: classpoly
  can now compute, in place of H_D^inv, its decomposition through a subgroup G of the class group
  (V of degree m = #cosets, interpolants W_k, n = |G|), and `cm_method` uses it by default
  (`tower=0` off, `tn=<n>` forces |G|): root of V, then of U = X^n + yX^(n-1) + sum W_k(y)/V'(y) X^k.
  Env: `CLASSPOLY_TOWER=1` (decompose, G chosen by a CRT + root-finding cost model), `=2` (also
  sign-free, below), `CLASSPOLY_TOWER_N`, `CLASSPOLY_TOWER_C` (exact subgroup, tests),
  `CLASSPOLY_TOWER_SLACK` (percent added to the height bound).  Output = the F/G format that
  `classpoly_load` / `ff_classpoly_root` already read (see classpoly_crt.c).  At 768 bits on 192
  cores the 17 winners of the Nov experiment build in a median of 11.7 s (359 s total, 7.3
  core-h); the four timed with the plain polynomial went from 2,072 s / 49 core-h to 145 s / 2.6.
  Things that are NOT in the paper and matter:
  (1) subgroups: G = <a_i^c_i> with c_i | n_i, valid iff for every i with c_i < n_i the power
  relation a_i^n_i = prod a_j^r_ij has c_j | r_ij (orbits are then the positions with fixed
  e_i mod c_i, whatever direction the enumeration walked) -- the paper's prefix form is not
  enough, e.g. relative orders (2,1894,3) need c = (2,2,3) for 12 x 947;
  (2) sign-free variant for the sqrt invariants (Weber f, ...): odd coefficients stored times the
  trace T plus S = T^2 (header line `S=`), which skips classpoly's trace CRT (23-58 s PER WORKER
  at h = 117,668 in multi-job mode);
  (3) the height bound is HEURISTIC for inv != j and can be too small (measured 3.0% over the
  unscaled bound, ramified w_{5,7}): 5% slack by default, and cm_method retries once with 20%,
  then falls back to the plain polynomial; a wrong decomposition is only caught by V/U not
  splitting or the trace check, so keep those checks;
  (4) `fp_find_root_bounded` (gives up on non-splitting input by evidence in bits) and a parallel
  first splitting step (lg(deg) shifts at once, 1.7x) for degrees 48..4096.
  Validated: 908 decompositions over Z vs PARI (21 invariants, sqrt cases, both subgroup forms),
  `ecpp/test_tower.py` (also mod P, multi-job, sign-free, retry/fallback via `fault=`).  Bugs
  fixed on the way: `qform_table_free` dangling pointer (crash on re-entering compute_classpoly),
  `fpoly_rem` with deg f < deg g.  Known and unchanged: the PLAIN height bound is ~2x too large
  when the level is ramified (sqrt of the class polynomial).  Root-finding is what is left: the
  slow builds are class numbers with a large prime factor (16 x 9337: 110 s), and fproot's
  parallel multiplication above degree 4096 gains only 1.1-3x.
- **yafu IFMA ECM (2026-11-14)**: yafu 3.1.9 built with `USE_AVX512IFMA=1` on c8a (Zen 5) is
  correct and 4-5x cheaper per curve than gmp-ecm 7.0.5 on 1030-bit inputs with all cores busy
  (0.059 / 0.298 / 0.946 vs 0.274 / 1.206 / 4.938 core-s at B1 = 5e4 / 2.5e5 / 1e6, same B2 via
  `-avxecm_stg2 1`), finding planted factors at least as often (388 vs 333 of 864).  Its ECM
  (twisted Edwards, FFT stage 2) is NOT the standalone avx-ecm, whose IFMA build is wrong.  Runs
  8 curves per call (vector width); expression goes on stdin when not a tty.  Not wired into the
  prover.  Report sections 7-8 of `reports/shortecpp-improve/`.
- **Stopping-rule scheduler (2026-11-15, same branch)**: levels of >= 500 bits
  (`POLICY_MIN_BITS`) are run by `ecpp/shortsched.py`; smaller ones, and `sched=old`, by the old
  tiered rounds.  The level is treated as a bandit problem.  `shortsim.IndexTable` gives every
  candidate state (coverage class, sweeps done, tail bits) its break-even reward r, the value of
  a win in core-seconds at which the next sweep pays; the table does not depend on the price of
  fresh orders and takes ~2 s per level.  `Sched` works candidates in increasing r, and only
  while r <= R, the reward at which a fresh order is worth the scan's marginal price (from the
  same table and the histogram of orders seen); R rises as the scan widens and parked
  candidates come back.  `level_loop` runs it against a driver: `RuleDriver` in shortECPP.py
  (dscan chunks doubling B on a share of the cores, gmp-ecm units) or `SimDriver` (virtual time,
  tails drawn from the model: `shortsched.py sim <pool.pkl> threads=..`).  Things that matter:
  (1) a candidate's class is just round(log2 s): the divisors of the smooth part s cover every
  delta in (0, log2 s), so no divisor enumeration is needed (s is the B-smooth part in the
  revised format, so the strip is called with n2 = V2_BCAP there);
  (2) sweeps are cut into units of a few seconds and the sweeps of one B1 may be in flight
  together, otherwise the winner's deep sweeps run one after the other while the other cores
  do work the sequential rule would never have done; a larger B1 waits for the smaller;
  (3) the supply price depends strongly on the prime (2,422 vs 6,303 solvable D below 6e7 for
  two 768-bit primes), so nothing about it is assumed: it is measured chunk by chunk;
  (4) gmp-ecm's cost is calibrated per run (`ecm_machine_factor`: one B1=5e4 benchmark with all
  threads busy, scaling `shortsim.SEC1030_LOADED`); measured/modelled cost came out 0.96-1.15
  per B1; success probabilities are `ecpp/ecm_pi.json` (= the report's pi_calib.json);
  (5) the log prints the wins the model expected from the work done -- at a level's first win
  that number is Exp(1) if the model is right, which is the calibration check.
  Measured, 13 fresh 768-bit primes on the 16-core box (32 threads): whole chain in a mean of
  202 s (median 112 s), search 1.5 thread-hours on average, expected wins at the win 0.65 +- 0.29
  (so the model is not optimistic); the old scheduler needed 59 and 148 core-h of search on a
  c8a for two primes of this size.  Simulation (768 bits, c8a costs): machine cost per level
  1.8 core-h at 4 threads, 1.6 at 32, 1.6-1.9 at 192 (about 30 s wall), 2.7 at 768; 1030 bits on
  192 cores: 37 core-h, 685 s wall (the August campaign: 812 core-h).  Real, 1030 bits:
  nextprime(10^310) proved on the 16-core box in 5,081 s, 43.9 thread-hours in all (search:
  181,428 orders, scan 13.4 + ECM 29.8 thread-h = 32 c8a core-h; expected wins at the win 0.66;
  top-level D = -155,888,345,711, A_47 class polynomial of degree 251,808 built in 105 s, nearly
  all CRT at 32 threads).  The last scan chunk arrived 5 min before the win (14% of the search
  hardly used): chunks double B (simulated growth of 1.5 or 1.25 is no better, so that stays).
  **On a c8a.48xlarge (192 cores, 2026-11-15, `data/rule_c8a/`)**: 26 primes of 768 bits in
  44 s per chain on average (median 36 s; search 28 s, 1.46 core-h of work, cores 97% busy;
  1.56 core-h of CPU per prime in all), 3 primes of 1030 bits in 167 / 1,157 / 690 s (search
  6.8 / 42.9 / 24.6 core-h; 10^310+733 is the 167 s one), all 29 verified twice in 53 minutes.
  Simulation and reality agree (768: 32 s and 1.6-1.9 core-h simulated), and over all 38 real
  768-bit runs the expected wins at the first win average 0.80 (exact model: 1 +- 0.16); first
  wins come from the B1 = 2.5e5..3e6 sweeps, as the model says.  What is left at 1030 bits is
  the build: 305 s and 195 s for two of the three, i.e. root-finding at prime degree 19,583 and
  10,891 on ~50 of the 192 cores.  Build side, found on the
  way: `build_level` skips hits with an odd prime dividing gcd(m, v) on the H_D curve (they can
  only work from the volcano floor; 14.8 s -> 2.9 s when it happens) and, when the curves with
  the computed j have the wrong order (`j_order_ok`), recomputes j from the plain class
  polynomial instead of dropping the winner; classpoly's decomposition also takes subgroups
  given by sign-invariant lattices (see compute_classpoly.c; 26 s -> 12 s for the one real
  winner in 29 that needed it).  Height bound: per Drew, optimise for speed and rely on the
  downstream checks; the measured true/bound ratio for non-j invariants has median 0.77 and a
  tail to 1.03, so the 5% slack + retry is left as is (slack tuning moves <= 5% of CRT time).
  Tests: `ecpp/test_sched.py`.  Report sections 9-10 of `reports/shortecpp-improve/`.
  Two things the repeated runs showed: (a) for a given prime the SAME discriminant tends to win
  on every machine and engine (9 of 13 primes at least twice in three sets of runs), so reruns
  of a prime are not independent tests of the model -- count primes, not runs; (b) 8 of 26
  primes won before the model expected 0.12 wins (an exponential law gives 3): the model may
  underrate easy winners.  The strip does not scale beyond 8 threads (4,000 orders: 0.41 s on 8
  or 16 threads, 0.53 s on 32), so the rule gives it `STRIP_THREADS = 8` and frees the scan's
  other threads for ECM meanwhile (event "scanthreads").
- **yafu as the ECM engine (2026-11-15, same branch)**: `shortECPP.py engine=yafu yafu=<binary>`
  (or `SHORTECPP_YAFU`; `engine=auto` uses yafu when a binary is given and passes a planted-
  factor self-test) runs the stopping-rule levels with yafu 3.1.9's AVX-512 ECM; the tiered
  rounds (levels below 280 bits since 2026-11-16, 500 before) and P-1 stay with gmp-ecm.  Per engine, `shortsim.ENGINES` holds the cost of a curve
  at 1030 bits (`SEC1030_YAFU`: 0.0036 / 0.0148 / 0.0599 / 0.302 / 0.95 / 2.69 / 9.9 / 36 core-s
  at B1 = 2e3 .. 3e7, i.e. 3.3-5.2x cheaper than gmp-ecm), its growth with size (yafu:
  bits^1.36 up to 1040 bits = 20 words of 52 bits, 14% above that law beyond), B2 per B1, the
  success table (`ecpp/ecm_pi_yafu.json`, 182 cells from `data/yafu_pi_calib.py`, 2.5 h on the
  16-core box) and the curve quantum (8: `Sched(quantum=8)`, sweeps in multiples of 8,
  `default_tiers("yafu")`).  Calibration result: a yafu curve succeeds ~30% more often than a
  gmp-ecm curve at B1 = 2e3, 15-20% more from 1.1e4 to 2.5e5, equally from 1e6; yafu's rates match the
  two-stage Dickman curve with a factor-24 gain within 15% everywhere.  Calling conventions that
  matter: `-avxecm_stg2 1` (B2 a function of B1; mode 2 takes B2 from a memory limit),
  `-aprcl_p 1` (otherwise every prime cofactor gets an APR-CL proof, 2 s at 1000 bits),
  `-logfile /dev/null -session /dev/null`, an empty working directory (it reads yafu.ini from
  cwd), expressions `ecm(n, curves)` one per line on stdin, `-v -v -v` only when per-lane hits
  are wanted (calibration).  Curves are seeded from clock + pid + system randomness, `-seed`
  does not control them.  A thread of the 16-core box does as much yafu work as a c8a core up to
  B1 = 2.5e5 (gmp-ecm: 1/1.35).  yafu is NOT in the repo: clone bbuhrow/yafu (commit 4d88b0b),
  ytools, ysieve; gmp-ecm 7.0.6 from its tarball; `make yafu USE_AVX512IFMA=1 ECM=1
  GMP_INCDIR=/usr/include/x86_64-linux-gnu GMP_LIBDIR=/usr/lib/x86_64-linux-gnu
  ECM_INCDIR=<ecm> ECM_LIBDIR=<ecm>/.libs`.
  Simulated (c8a costs, each engine with its own tables): a 768-bit level 0.77 core-h on 192
  threads (gmp-ecm 1.81) and 0.52 on 32 (1.45); a 1030-bit level 12.6 core-h on 192 threads
  (gmp-ecm 41.3).  Real, 16-core box, the 13 primes of the 2026-11-15 entry: chain 41 s on
  average (median 34 s) against 202 s (112 s) with gmp-ecm, search work 0.24 thread-h against
  1.52.  nextprime(10^310) on that box with yafu: 3,581 s for the chain, 25.5
  thread-h of CPU (gmp-ecm: 5,081 s, 43.9), on a slow draw (1.26 expected wins at the win; the
  gmp-ecm run had 0.66), i.e. 22 against 66 thread-h of search per expected win; its winner is
  D = -83,632,771, the August campaign's.  In that run 27% of the ECM time was "above the
  rule" (cores waiting for a scan chunk; simulation: 1%) -- the scan's share of the cores
  should be sized from the observed backlog, not done yet (tried 2026-11-16 and dropped, see
  the next entry).  Report section 11.
- **The whole ShortPrimalityProofs table on a c8a with yafu (2026-11-16, same branch)**: every
  entry nextprime(10^c), c = 10..320, is proved faster than the README records, by 5x to 750x
  in CPU time; c = 160..320 take 17 minutes and 42 core-hours in all on 192 cores (recorded:
  4,180 core-hours), c = 10..150 take 19 s; 10^310+733 takes 268 s.  The thirteen 768-bit
  primes: chain 21 s mean / 14 s median, 0.52 core-h of CPU (gmp-ecm, same box: 34 s / 27 s,
  1.11).  Data `reports/shortecpp-improve/data/yafu_c8a_{1,2,3}/` (three passes, 95 min of
  aws-spot-7), `bench_table.py`, report section 12.  c = 330..400 NOT rerun (estimate: 5-10
  hours of one box).  What the runs turned up, all fixed:
  (1) yafu's `ecm(n, curves)` drops an input after the batch of 8 curves in which a factor
  appears, even if the cofactor is composite (gmp-ecm goes on with the cofactor).  The model's
  state "j sweeps done, tail of t bits" means a tail that survived the sweeps IN FULL, so
  crediting the curves run and re-offering the rest is wrong too; `EcmRunner._run_yafu` reruns
  the remaining curves on the cofactor until the sweep is complete (`parse_yafu(...).ran` =
  curves run, from yafu's progress lines "ecm: K/TOTAL curves").
  (2) Small levels are sized to their work: `level_params(nb, threads, bound)` gives the strip
  bound (2^16 below 200 bits, 2^20 below 320, else 2^26; `set_strip_bound`), the first scan
  (1e5 below 200 bits ... 1e7 below 500) and the thread cap (2^((nb-180)/40)); the gp tail
  takes over only below `cm_min_bits` = 2 * bits(recursion floor) + 24 (64-94 bits in the
  revised format, was 135).  `POLICY_MIN_BITS = 280`: for the rule below 500 bits the model's
  bound D follows the size (`shortsim.medium_bits` = bits/3, capped at 160; `set_grid` sets
  it), the ladder is cut where its curves would reach D (`default_tiers(engine, nbits)`), the
  reward grid scales as (bits/500)^8, the machine benchmark is skipped below 400 bits, and a
  shared ECM process stops growing at 0.05 expected wins (`Sched(unit_pw=)`).  Model check at
  299 bits by brute force (`data/small_model_check.py`: every order of a scan through the
  whole ladder): 94 winners against 97 modelled, sweep by sweep 9/33/37/15 against
  14.1/38.5/31.2/15.0.  The "expected wins at the first win" check is NOT valid at small
  sizes: the first handful of units carry about three expected wins, so the count at the win
  averages 3 whatever the model.
  (3) classpoly aborted silently when a modular polynomial was missing: the repo carried
  Phi_l for l <= 71 only, and about one winner in ten at 768-1000 bits needs more (Phi_101
  and Phi_79 for Weber f, Phi_73 and Phi_83 for w2w5e1 were the cases seen).  Since
  2026-11-16 the bundle goes to l = 109 for f, f2 (Drew's call) and for eight double eta
  families with small files (w2w3e1, w2w3e2, w2w5e1, w2w5e2, w2w7e1, w2w7e2, w2w13e1,
  w3w3e1: my addition, a separate commit): of 3,000 random fundamental D in [1e7, 1e10], 99
  asked for a missing file before and 31 do now (20 of them Atkin functions).  `git archive` to a spot box ships only that subset: COPY
  THE FULL `phi_files/` (2.2 GB, 3 min) to any box that runs the prover.  The message now goes
  to stderr whatever the verbosity and `get_j` gives up on such a D at once
  (`_CM_UNBUILDABLE`).
  (4) `shortsim.IndexTable` was computed in threads to no effect (GIL); it now runs its
  classes in child interpreters fed over pipes (`_solve_in_processes`; identical tables,
  4.0 s -> 0.6 s per level on the c8a).  Not `multiprocessing` with a fork server: that
  re-imports the caller's main script, and an unguarded script reruns itself in every worker.
  Measured and NOT changed: on 192 cores with yafu 40-100% of the ECM core-seconds are handed
  out "above the rule" (cores without work below R while a scan chunk runs); giving the scan
  more cores instead costs more per win in simulation (0.78 against 0.65 core-h at 768 bits),
  because it buys the next, dearer chunk sooner while the work above the rule is mostly work
  brought forward.  Next lever is the build: root-finding seconds are
  0.0145 x^0.93 (bits/768)^2.74 on the c8a within a factor 1.25, x = the larger part of the
  most balanced factorisation of h(D) (`data/build_fit.py`, 31 builds); builds are a third of
  search + build at 768 bits with yafu (median 2.7 s, mean 5.7 s, tail to 90 s; at 1030 bits
  median 13 s, mean 53 s, tail to 305 s), and x grows with |D| (mean 150 near 1e7, 6,500 near
  1e11).  PARI gives 20,000 class numbers in 3-7 s, so a candidate's build cost can be
  charged at intake; the scheduler does not do that yet.  At the win it does:
  `RuleDriver.should_build` (called by `level_loop` when a winner is in hand) keeps the
  winner in reserve and lets the search go on while its estimated build (`build_estimate`:
  the law above, scaled by this run's own builds) exceeds R/threads + the median estimate
  over a sample of the pool, for at most the time that build would have taken; `try_build`
  takes the cheapest winner first.  On 192 cores at 768 bits the threshold is ~10 s; on 32
  threads it is a minute, so nothing is passed there.  Unit-tested and simulated, NOT yet
  exercised on a c8a.

- **The frontier entries and 10^410 (2026-11-17, same branch)**: on one c8a with yafu,
  nextprime(10^c) for c = 330..400 take 7 min, 44 min, 2.0 h, 43 min, 41 min, 9 min, 5 min and
  87 min: 6.0 h and 1,020 core-h in all (recorded: 51 h, 8,290 core-h); every table entry from
  10^10 to 10^400 is now faster than recorded (`data/yafu_c8a_5`, `bench_table.py`).
  **nextprime(10^410) = 10^410 + 67 (1362 bits) is proved**: `certs/short2/nextprime1e410.csv`,
  5 levels, 28,153 s on 192 cores, 953 core-h; search 4.8 h (362,752 orders), build 3.0 h, of
  which 2.6 h for a root at prime degree 470,279 (h = 4 x 470,279; `build_estimate` said
  3.65 h).  It is not in the upstream table (Drew decides about a PR there).  What that prime
  forced:
  (1) the supply has TWO arms.  A first run stalled: 10^410 + 67 has 5,692 solvable D below
  1e9 (10^400 + 69: 20,888), and the factor-base cap stopped at 1e8, so past B = 1e11 the
  only way to more orders was a wider B.  Now `RuleDriver.next_chunk` offers "widen" (B ->
  2B at the present cap) and "raise" (cap -> 4 cap over the range already scanned: the
  discriminants passed over for a prime above the cap, Drew's suggestion) and takes the
  cheaper by predicted core-seconds per order; the prediction of a raise is corrected by the
  last one's outcome (`raise_fit`).  The log prints "next orders: widen|raise (prices)".  Cap
  ceiling = a quarter of RAM for the square roots (`fb_cap_max`: 2e10 on a c8a, 6e9 on the
  16-core box, at 1362 bits); default scan bound 1e15.  dscan: `Bmin=B fbmin=<old cap>` is a
  cap-only scan and skips old-prime nodes with no room for a new prime (exact: earlier scan
  + raise = full scan at the new cap, test_dscan).  Neutral at 1000 bits (same price per
  order as before), for 10^410 the cap went up six times (1e6 to 4.1e9) and the
  362,752 orders cost 13% less each than the first run's 356,128 (next orders 5.9 against ~7.7 core-s).
  (2) levels of 1300 bits and more get 192 curves at B1 = 3e7 on offer instead of 64 (the
  stalled run had 2,829 candidates through the whole ladder).  B1 = 1e8 was timed and
  rejected: 4.4x the cost of 3e7 for ~2x the success on 120-160-bit factors.
  (3) the model above 1100 bits: expected wins at the first win over the nine primes
  10^330..10^410 are 0.61, 3.58, 4.41, 1.79, 0.92, 0.20, 0.05, 0.66, 0.73 (mean 1.44; exact
  1 +- 0.33): two slow draws in a row looked like a break and were not.  Checked on the way:
  yafu finds planted factors equally often in numbers of 780..1330 bits and labels prime
  cofactors correctly; "splits found/expected by B1" (now in every search line) is 1.05 at
  B1 = 2e3 rising to ~1.5 at 1e6+ (the expectation ignores that deep candidates are
  selected).  A deeper calibration of yafu at B1 = 3e6, 1e7, 3e7 (598,344 curves,
  `data/yafu_pi_calib_deep.json`) gives 0.87, 0.78, 0.95 of the hits the table in use
  expects: NOT yet folded into `ecpp/ecm_pi_yafu.json` (do it, then rerun test_sched).
  (4) `RuleDriver.should_build` has still never fired in a real run (25 runs of 768-bit
  primes drew no costly winner); at 1362 bits it rightly accepted a 3-hour build because
  another winner was 6 hours away.  `get_j`'s timeout follows the build estimate.
  (5) since 2026-11-18: levels of >= 900 bits checkpoint their search every 10 min
  (`work/shortcamp/<tag>.search.pkl`, `RuleDriver.checkpoint` / `_load_ckpt`, `level_loop`'s
  `restore`; env `SHORTECPP_CKPT` seconds, `SHORTECPP_CKPT_BITS`), and `resume=1` takes the
  search up from it (lost: units in flight and the roots cache); a build is NOT checkpointed.
  `ecm_pi_yafu.json` = first + deep calibration pooled (`data/fold_yafu_calib.py`), and
  `EcmModel.pi` pools a cell's ratio to theory with its neighbours up to 100 hits (the old
  "cells with >= 8 hits" cut selected on the outcome: 1.44x / 1.49x theory beyond 110 bits at
  B1 = 5e4 / 2.5e5); the gmp-ecm table is read the same way (model vs campaign log as good
  or better: 0.67 winners per reference pass, 20.0 modelled against 17 found at 768 bits).
  Reward grid to 3e8 core-s from 1300 bits.  nextprime(10^c), c = 420..450 (offsets 1261,
  1587, 709, 2269) launched 2026-11-18 00:23 UTC on four spot boxes (7, 16: c8a; 15, 12:
  c7a = Zen 4, yafu's IFMA build runs there); run/pull scripts in the session scratchpad
  (`big/`), off-box copies in `ckpt_big/`.  A fresh box needs: `git archive` + full
  `phi_files` + `~/.local/ecmpkg` + `~/yb/yafu/yafu` + `apt install python3-numpy`.
  Spot notes: a build uses ~20 of 192 cores for hours: re-claim with `--idle-ok --force` or
  the watchdog stops the box; `spot claim` again with new hours works on one's own box.
  Cost of the two days: about $62 (c8a at ~$3.06/h in us-west-2, the only region with
  capacity).

## Build & test
```sh
make -j                      # ff_poly -> classpoly (incl. zp_poly) -> ecpp (all in-tree, ~15 s)
make test [MAXD=1000]        # classpoly vs PARI (Tests 1/2/3); ~2.5 min
. ./setenv.sh                # env: phi dir, PATH, work/pcache; needed by oneshotECPP/cm_method
./ecpp/oneshotECPP p=<decimal>   # THE tool: prime -> certificate (or pbits=<n> seed=<s>)
./ecpp/oneshotECPP pbits=256 seed=1 threads=16 c=1 B0=4000000 B=20000000000  # all knobs
python3 ecpp/test_dscan.py   # validate the discriminant scan vs PARI + brute force
python3 ecpp/test_smooth.py  # validate smoothness engine vs PARI (Test A + end-to-end gate)
python3 ecpp/test_tower.py   # class-group decomposition: classpoly CLASSPOLY_TOWER + cm_method vs PARI (~2 min)
python3 ecpp/test_sched.py   # stopping-rule scheduler: break-even table vs reference DP, queue mechanics, simulated levels (SHORTECPP_YAFU=<binary> adds the live yafu checks)
make -C ecpp voneshot && python3 ecpp/test_voneshot.py  # C one-shot verifier vs voneshot.py (differential)
ecpp/voneshot $(cat certs/25519.txt)   # C verifier, same CLI as voneshot.py; --bench <certfile> times it
make -C ecpp vshort2 && python3 ecpp/test_vshort2.py  # C short-ECPP verifier (faithful + --fast) vs vshort2.py
ecpp/vshort2 --bench certs/short2/certs.csv 1329   # faithful vs one-ladder-per-level, in MR rounds
./ecpp/roottest pari 256 200 # validate F_p root-finder vs PARI polrootsmod
./ecpp/cm_method D=-24447 pbits=256 seed=1   # j-invariant with CM by D (needs env)
./ecpp/dscan pbits=256 B=1e8 threads=32 dump # scan alone (Bmin=|verify|dumpscan|sd=|seed=)
./ecpp/smoothtest gate pbits=256 seed=1 y=4294967296 < dump.txt  # fixed-y gate (tests/analysis)
make -C ecpp tests           # invjtest cmjtest asmtest zpbench (validation drivers)
```
Toolchain: gcc 13, GMP 6, PARI/GP 2.18. Box: **16 physical cores** (32 vCPUs, HT) — so ~16× is the
parallel ceiling. Verify certs: `python3 voneshot.py $(cat certs/<f>.txt)` (voneshot from the
challenge repo, github.com/AndrewVSutherland/OneShotPrimalityProofs).

## Must-know gotchas
- **PARI/GP**: pipe on **stdin** (`gp -q < f`), never `gp f.gp` (hangs interactive). One statement
  per line. `polclass` needs `default(parisizemax,"8G")`.
- **classpoly** self-isolates its CRT scratch under `/tmp`. Dir hooks are env vars:
  `CLASSPOLY_PHI_DIR`, `CLASSPOLY_H_DIR`, `CLASSPOLY_TMPDIR` (a *base* dir).
- `class_inv_mpz.c` lives in `classpoly_v1.0.3/` so its quote-includes bypass `-isystem`; keep the
  vendored headers warning-clean or the ecpp build gets noisy.
- **P caches**: native-limb format (`mpz_out_raw`'s 4-byte size field overflows > 2³¹ bytes).
  Segment files carry `lo` (v2 header); v1 full-product files still load. Self-check on load.
- Compile new C with **`-Wall -Wextra`** (found real bugs in the vendored code this way).
- Deferred optimizations (don't do until obviously needed, per user): packing factor-base sqrts
  as contiguous limbs, hand-rolled 1-limb-quotient division.  (Bounding the factor base is done
  for shortECPP — capped from the first scan, see the 2026-11-14 entry; `oneshotECPP` still
  scans uncapped.)
- Half-gcd/Lehmer does **not** help our (short, small-`p`) Cornacchia descent — don't reach for it.

## Conventions
- New ECPP code lives in `ecpp/`; keep `test_dscan.py` + `test_smooth.py` green.
- Every solver/scan result is self-checked (`4p == t²+|D|v²`), every winner independently verified
  (PARI cross-checks in tests; `voneshot.py` for certificates).
- Detailed findings are also in the session memory files (recalled automatically).

When addressing codex-review findings, follow the Round protocol in AGENTS.md -> Code Review Rules.
