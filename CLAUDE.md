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
- **`oneshot p=<prime>` → n⁴-smooth one-shot ECPP certificate**, verified by the challenge's
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
  (oneshot reduces S only by n₁'s >97 part). Validated vs PARI (floor/trace/group/representability).
- **phi bundle**: 46 MB subset committed via git add -f (`phi_files_manifest.txt`): every family's
  `Φ_{ℓ,f}` for ℓ≤71 + all `Ψ_f` maps. Classical `Φ_ℓ` are needed only at the chosen invariant's
  LEVEL primes (e.g. `phi_j_59` for A₅₉) + `Φ_2`, never for enumeration. Full 2.2 GB DB external;
  misses degrade gracefully (skip winner). Fresh-clone validated (make + tests + oneshot).
- Fixed upstream-worthy bugs: `class_inv_mpz.c` `mpz_j_from_u8` (unreduced J + uninit T3);
  vendored-header warnings (qform/prime/mpzutil/iqclass).
- **Parallel classpoly**: `CLASSPOLY_JOBS`/`CLASSPOLY_JOBID`/`CLASSPOLY_ECRT_DIR` env → multi-job
  ECRT (workers split CRT primes mod W, `ecrt_dump` partial sums; JOBID=0 merges). `cm_method jobs=N`
  orchestrates (auto-gated |D|≥1e8); oneshot passes threads. h=35085: 158→29 s, byte-identical.
- **Hybrid root-finder**: fproot keeps its Montgomery/Kronecker/Barrett powmod but delegates gcd +
  exact division to zp_poly's half-gcd above degree 1024 (~16–19% at d≥8000; full zp delegation
  measured parity — its modmul is slower, its gcd faster).
- **Parallel root-finding (OpenMP)**: the EDS squarings are sequentially dependent, so the
  parallelism lives inside each product — task-parallel Karatsuba over polynomial halves
  (`fpoly_mul_par`, 3^depth leaf Kronecker products, depth ≤ 3 → up to 27 tasks) for the
  squarings and Barrett's two big multiplications, above degree ~4096.  ff_poly is not involved
  in this step, so OpenMP is safe (per Drew).
- **Big-machine scaling (2026-07-03)**: oneshot's threads default is `omp_get_max_threads()`
  and is passed through to the `dscan`/`cm_method jobs=` shell-outs (was: hardcoded 16 —
  capped a 256-core production box at 16 cores for the scan phases). dscan's factor-base
  sieve reuses smooth.c's parallel `sieve_primes_range` (now exported; dscan links smooth.o)
  and the QR survivor collection fills slots via a parallel prefix sum; the prime array is
  64-bit — the old `uint32_t` silently truncated once B > 2³² ≈ 4.3e9 (reachable: Bmax
  defaults to 2e10). Validated: dump/dumpscan sorted-identical pre/post (5 configs, B ≤ 1e8),
  test_dscan.py + test_smooth.py green, (2³², 2³²+6e4) window exact vs PARI (enumeration,
  qfbsolve solvability, 4p = t²+dv² arithmetic), 256-bit oneshot cert re-verified in PARI.
  Known remaining ceilings at 256 cores: classpoly jobs setup overhead, fproot's ~27-task cap.
- **Smoothness at scale (2026-07-04)**: `smooth_parts_multi` fuses a batch against ALL ladder
  rungs — one product tree, one shared chunk-power table (parallel prefix doubling; the old
  `reduce_big_mod` recombination ladder was G *serial* |X|-sized mulmods ≈ 90 s single-core per
  big batch at 256 threads, growing with thread count), one flat (segment×chunk) job list, one
  descent (gcd vs ∏P_s = ∏ per-seg parts; segments disjoint). oneshot overlaps norm-equation
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
  root). Raced vs oneshot on a 192-core spot box, same prime each size: 256: CM 4.1 s / SEA
  68 s; 384: 2306/615; 416: 1774/2805; 448: 7182 / stopped 20035 s no winner. Per-candidate:
  SEA flat 7–14 core-s, CM avg 2.2–4.5 rising — expected-cost crossover extrapolates to
  ~480–512 bits (walls flip on density luck from 384 up; observed fudge over ρ(u/2): 1–8×).
  No SEA-side accelerations implemented (stage-1 gate/batch tail need ellsea internals).
  Report: `reports/sea-crossover/`; certs in `certs/sea/`; raw cands in `work/race-archive/`.

- **Short ECPP records (2026-08-18)**: `ecpp/short_prove.py` — CM prover for the
  ShortPrimalityProofs *short ECPP* format; produced the table entries for nextprime(10^c),
  c = 210..310 (largest: 1030 bits, beating Bernstein's 1025-bit AKS example) in one day on
  spot instances (~3,500 core-h).  Per level: dscan (new `maxfb=` bounds the factor base —
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
  `ecpp/vshort2.py` (working verifier), `short2.gp`, `repair_short2.py`, `short_prove.py
  v2=1` (gated on vshort2), `assemble_short2.py` (zero-arg = validate the tracked table).
  Full 31-chain table migrated (repairs + forced truncation) and dual-verified:
  `certs/short2/certs.csv`.  Existence swept for all p ≤ 4096.  Radical-capped stage merged
  upstream (ShortPrimalityProofs cc2caf1); terminal-prime + settled conventions on the
  AndrewVSutherland2 fork main (verifier renamed vshortECPP.py), upstream PR pending (Drew).
  Spec note: `reports/short2-spec/`; paper: Def 6.2 + Prop 6.4 in `reports/ecpp-varieties/`.

## Build & test
```sh
make -j                      # ff_poly -> classpoly (incl. zp_poly) -> ecpp (all in-tree, ~15 s)
make test [MAXD=1000]        # classpoly vs PARI (Tests 1/2/3); ~2.5 min
. ./setenv.sh                # env: phi dir, PATH, work/pcache; needed by oneshot/cm_method
./ecpp/oneshot p=<decimal>   # THE tool: prime -> certificate (or pbits=<n> seed=<s>)
./ecpp/oneshot pbits=256 seed=1 threads=16 c=1 B0=4000000 B=20000000000  # all knobs
python3 ecpp/test_dscan.py   # validate the discriminant scan vs PARI + brute force
python3 ecpp/test_smooth.py  # validate smoothness engine vs PARI (Test A + end-to-end gate)
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
- Deferred optimizations (don't do until obviously needed, per user): **bounding the factor base**
  (dscan's Tonelli wall — ~75% of large-B scans), packing factor-base sqrts as contiguous limbs,
  hand-rolled 1-limb-quotient division.
- Half-gcd/Lehmer does **not** help our (short, small-`p`) Cornacchia descent — don't reach for it.

## Conventions
- New ECPP code lives in `ecpp/`; keep `test_dscan.py` + `test_smooth.py` green.
- Every solver/scan result is self-checked (`4p == t²+|D|v²`), every winner independently verified
  (PARI cross-checks in tests; `voneshot.py` for certificates).
- Detailed findings are also in the session memory files (recalled automatically).
