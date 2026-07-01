# CLAUDE.md — OneShotFastECPP

CM-method ("fast ECPP") approach to one-shot elliptic-curve primality proofs.
**Read `design.md` for the full technical picture, decisions, and roadmap.**

## Status (2026-07-01) — COMPLETE, end to end
- **`oneshot p=<prime>` → n⁴-smooth one-shot ECPP certificate**, verified by the challenge's
  `voneshot.py`. Cold runs (no caches): 256-bit ~3 s, 10⁸⁰+129 4.5 s, 10⁹⁰+289 24 s,
  10¹⁰⁰+267 ~19 min (B=10¹⁰ scan + a degree-35085 H_D). Verified certs in `certs/`
  (2²⁵⁵−19, NIST P-256, secp256k1, 10⁸⁰/10⁹⁰/10¹⁰⁰+ε).
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
  (needs `N≡0 mod4` + `m|exponent`; oneshot filters/skips accordingly).
- **phi bundle**: 46 MB subset committed via git add -f (`phi_files_manifest.txt`): every family's
  `Φ_{ℓ,f}` for ℓ≤71 + all `Ψ_f` maps. Classical `Φ_ℓ` are needed only at the chosen invariant's
  LEVEL primes (e.g. `phi_j_59` for A₅₉) + `Φ_2`, never for enumeration. Full 2.2 GB DB external;
  misses degrade gracefully (skip winner). Fresh-clone validated (make + tests + oneshot).
- Fixed upstream-worthy bugs: `class_inv_mpz.c` `mpz_j_from_u8` (unreduced J + uninit T3);
  vendored-header warnings (qform/prime/mpzutil/iqclass). Root-finding stays on `fproot`
  (zp_poly parity per `ecpp/zpbench`; zp_poly's fast gcd would win at deg ≫ 10³ — future lever).

## Build & test
```sh
make -j                      # ff_poly -> classpoly -> zp_poly -> ecpp (all in-tree, ~15 s)
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
  hand-rolled 1-limb-quotient division, zp_poly half-gcd for deg≫10³ root-finds.
- Half-gcd/Lehmer does **not** help our (short, small-`p`) Cornacchia descent — don't reach for it.

## Conventions
- New ECPP code lives in `ecpp/`; keep `test_dscan.py` + `test_smooth.py` green.
- Every solver/scan result is self-checked (`4p == t²+|D|v²`), every winner independently verified
  (PARI cross-checks in tests; `voneshot.py` for certificates).
- Detailed findings are also in the session memory files (recalled automatically).
