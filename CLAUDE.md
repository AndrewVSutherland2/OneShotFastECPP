# CLAUDE.md — OneShotFastECPP

CM-method ("fast ECPP") approach to one-shot elliptic-curve primality proofs.
**Read `design.md` for the full technical picture, decisions, and roadmap.**

## Status (2026-06-30)
- **classpoly + ff_poly**: build in-tree, validated vs PARI (`make test`, 3 suites green over |D|≤1000).
- **Component (a) CM-discriminant search**: DONE in `ecpp/` — Cornacchia + factor-base scan
  (`dscan`), parallel & validated. Scan ~10 M D/s @256-bit (32 threads).
- **Component (b) n⁴-smoothness testing**: DONE in `ecpp/smooth.{c,h}` + `smoothtest`. Challenge now
  accepts **n⁴-smooth** certs (n=bitlen p): need the n⁴-smooth part of `p+1∓t` to exceed `L≈√p`.
  Bernstein batched remainder tree (cached prime-product `P=∏_{q≤n⁴}q`), parallel top reduction.
  Validated vs PARI end-to-end (128/256-bit, every winner checked). 256-bit yield ~1.3×10⁻³.
- **Component (c) root of `H_D` over `F_p`**: DONE in `ecpp/fproot.{c,h}` + `roottest` + `hdroot`.
  H_D splits completely (CM D chosen so), so straight to equal-degree splitting: `(x+a)^((p-1)/2) mod h`,
  gcd, keep smaller side, to degree 1. Big-`F_p` **Montgomery** arithmetic; `f² mod h` via **Kronecker**
  squaring + **Barrett/Newton** reduction; cheap `(x+a)·f` (shift+scalar+adds). Validated vs PARI;
  end-to-end via classpoly (`inv=0` = Hilbert/j). ~0.07 s/root at degree 100, 256-bit.
- **Component (d) CM engine + assembly + oneshot**: DONE. `cm_method D p` → j-invariant of `E/F_p`
  with CM by D (picks the **best class invariant**, `inv=-1`, 25–50× faster than j; converts back to j
  via `cminv.{c,h}` = Sutherland's `class_inv_mpz.c` formulas + `fproot` for the Atkin/η modular-poly
  root-find; validated 120/120 vs word-size `invtoj`). `mont_assemble` (`curve.{c,h}`) builds the
  Montgomery `(A,x₀)` with a point of order `m` (needs `N≡0 mod4` + `m|exponent`). **`oneshot p`** ties
  it all together → certificate verified by `voneshot.py` (2²⁵⁵−19 in ~6 s w/ cached P).
- **Full pipeline WORKS**: `git clone && make && . ./setenv.sh && ./ecpp/oneshot p=<prime>`. README+INSTALL.
- Target `p`: small (128–384 bit), `D` up to `10^10`+.
- **phi bundle**: only a 28 MB subset (`phi_files_manifest.txt`) is committed (git add -f); full DB (2.2 GB)
  external.
- **zp_poly landed** (Drew): in-tree at `zp_poly/`, staged into `./local` by the top-level make.
  `cm_j_from_inv` now calls classpoly's **authoritative `mpz_j_from_inv`** (`class_inv_mpz.c`, linked
  vs zp_poly + classpoly objects); the fproot-based `cm_j_from_inv_ref` is kept as cross-check + for
  the generic single-η range class_inv_mpz.c doesn't dispatch. 3-way validated (81/81 vs `invtoj`).
  Fixed in class_inv_mpz.c: `mpz_j_from_u8` unreduced J + uninit T3 (report upstream). Root-finding
  stays on `fproot`: `ecpp/zpbench` shows zp_poly +5–25% at ≤256-bit, parity at 384-bit — not worth
  switching (root-find is a minor pipeline cost; Montgomery form is integrated with curve.c).

## Build & test
```sh
make                         # ff_poly -> ./local, then classpoly (+ invtoj)
make test [MAXD=1000]        # classpoly vs PARI (Tests 1/2/3); ~2.5 min
cd ecpp && make              # cmsearch, dscan, smoothtest, roottest, hdroot
python3 ecpp/test_dscan.py   # validate the discriminant scan
python3 ecpp/test_smooth.py  # validate n⁴-smoothness engine vs PARI (Test A + end-to-end gate)
./ecpp/roottest pari 256 200 # validate F_p root-finder (EDS+Kronecker+Barrett) vs PARI polrootsmod
. ./setenv.sh; ./ecpp/hdroot D=-24447 pbits=256 seed=1   # end-to-end: H_D mod p -> root j0 (needs env)
./ecpp/dscan pbits=256 B=100000000 threads=32   # run the scan (dump|verify|dumpscan|sd=|seed=)
# smoothness: build/cache P once per bit-length, then gate solvable D (form p+1∓t, test smooth part>L):
./ecpp/smoothtest pbuild y=4294967296 threads=16 save=P_2e32.bin            # P for 256-bit (n⁴=2³²)
./ecpp/dscan pbits=256 seed=1 B=2e7 dump | ./ecpp/smoothtest gate pbits=256 seed=1 load=P_2e32.bin
```
Toolchain: gcc 13, GMP 6, PARI/GP 2.18. Box: **16 physical cores** (32 vCPUs, HT) — so ~16× is the
parallel ceiling. `. ./setenv.sh` to run classpoly interactively in-tree.

## Must-know gotchas
- **PARI/GP**: pipe on **stdin** (`gp -q < f`), never `gp f.gp` (hangs interactive). One statement
  per line. `polclass` needs `default(parisizemax,"8G")`.
- **classpoly** self-isolates its CRT scratch under `/tmp` now (was a parallel-crash hazard). Dir
  hooks are env vars: `CLASSPOLY_PHI_DIR`, `CLASSPOLY_H_DIR`, `CLASSPOLY_TMPDIR` (a *base* dir).
- Compile new C with **`-Wall -Wextra`** (found real bugs in the vendored code this way).
- Deferred optimizations (don't do until obviously needed, per user): **bounding the factor base**,
  packing factor-base sqrts as contiguous limbs, hand-rolled 1-limb-quotient division.
- Half-gcd/Lehmer does **not** help our (short, small-`p`) Cornacchia descent — don't reach for it.

## Conventions
- New ECPP code lives in `ecpp/`; keep `test_dscan.py` green (brute-force + PARI oracle, both p mod 4).
- Every solver/scan result is self-checked (`4p == t²+|D|v²`) and cross-checked vs PARI in tests.
- Detailed findings are also in the session memory files (recalled automatically).
