# OneShotFastECPP

Generate **n⁴-smooth one-shot ECPP certificates** by the CM method ("fast ECPP").
Given a prime `p`, it produces a certificate `(p, A, x₀, m, q₁…q_k)` that proves `p`
prime and is verifiable in quasi-quadratic time by
[`voneshot.py`](https://github.com/AndrewVSutherland/OneShotPrimalityProofs).

## Quick start

```sh
git clone <this repo> && cd OneShotFastECPP
make -j                           # ff_poly -> classpoly -> zp_poly -> ecpp tools (~1 min)
. ./setenv.sh                     # point classpoly at the bundled modular polynomials
./ecpp/oneshot p=$(python3 -c 'print(2**255-19)')
```

Output is a single line `p A x₀ m q₁ … q_k`:

```
57896044618658097711785492504343953926634992332820282019728792003956564819949 \
34272590291611932410078115265478194150993293036704524951758556452491882420894 \
39019257195751330428207690201832939607319582892511312056714375592481312290863 \
260621764559582591965379563027288448623 342527 864107 1396061 15802043 177075901
```

(Certificates are not unique — your run may print a different valid one.)
Verify it with the challenge repo's verifier:

```sh
python3 voneshot.py 57896044618658097711785492504343953926634992332820282019728792003956564819949 34272590...
# True
```

`oneshot` accepts `p=<decimal>` or `pbits=<n> [seed=<s>]` (a random n-bit prime),
plus `threads=<t>`, `c=<work ratio>`, `B0=`/`B=` (initial/max discriminant-scan
bounds), and `pcache=<file>`.  The smoothness bound climbs a power-of-2 ladder
starting just above n² — prime-product *segments* are built on demand (and cached
in `work/pcache/`, shared across all prime sizes), the candidate pool widens
geometrically in between, and the run stops at the first rung that yields a
certificate.  Certificates rarely need primes anywhere near n⁴, so even a fully
cold run takes seconds at 256 bits.

## What it does

For a fixed prime `p` it runs the CM method end to end:

1. **Discriminant search** — find a CM discriminant `D<0` with `4p = t² + |D|v²`
   solvable (Cornacchia over a factor base).
2. **Smoothness** — keep curve orders `N = p+1∓t` with `N ≡ 0 mod 4` whose
   **n⁴-smooth part exceeds `L = (p^{1/4}+1)²`** (`n = ⌈log₂ p⌉`); this gives a
   smooth `m | N`.  Batched with a remainder tree.
3. **Class polynomial + root** — compute `H_D mod p` in the **best class invariant**
   (via `classpoly`), find a root over `F_p`, and convert it to a `j`-invariant.
4. **Curve + point** — build the Montgomery curve `E_A/F_p` with that `j` and order
   `N`, find a point of order `m`, and emit `(p, A, x₀, m, q_i)`.

See **[`design.md`](design.md)** for the full technical writeup and performance.

## Programs (`ecpp/`)

| program | purpose |
|---|---|
| **`oneshot`** | prime → one-shot ECPP certificate (**the main tool**) |
| `cm_method` | `D`, `p` → the `j`-invariant of a curve `E/F_p` with CM by `D`, trace ±t (picks the best class invariant, converts back to `j`) |
| `dscan` | CM-discriminant search (Cornacchia + factor base), parallel |
| `smoothtest` | batched n⁴-smoothness testing (Bernstein remainder tree) |
| `roottest` | validate the big-`F_p` root finder against PARI |

## Example certificates (`certs/`, all verified by `voneshot.py`)

| prime | file | notes |
|---|---|---|
| 2²⁵⁵ − 19 | `certs/25519.txt` | Curve25519 field prime |
| 2²⁵⁶ − 2²²⁴ + 2¹⁹² + 2⁹⁶ − 1 | `certs/p256.txt` | NIST P-256 field prime |
| 2²⁵⁶ − 2³² − 977 | `certs/k256.txt` | secp256k1 (Bitcoin) field prime |
| 10⁸⁰ + 129 | `certs/1e80p129.txt` | |
| 10⁹⁰ + 289 | `certs/1e90p289.txt` | |
| 10¹⁰⁰ + 267 | `certs/1e100p267.txt` | |

Verify any of them with `python3 voneshot.py $(cat certs/<file>)`.

## Requirements

- **gcc** 13+ and **GMP** 6+ (the build and `oneshot` itself).
- **PARI/GP** 2.x — only for the correctness test suites, *not* for `oneshot`.
- Modular polynomials: a **46 MB subset** is bundled (`phi_files/`, see
  [`INSTALL`](INSTALL)); it covers the class invariants and levels the CM method
  uses over the 128–384-bit range.

## Verifying the build

```sh
make test                         # classpoly vs PARI (Tests 1/2/3)
python3 ecpp/test_smooth.py       # n⁴-smoothness engine vs PARI
python3 ecpp/test_dscan.py        # discriminant scan vs PARI + brute force
./ecpp/roottest pari 256 200      # F_p root finder vs PARI polrootsmod
```

## Layout

```
Makefile        setenv.sh          design.md
ff_poly_v2.0.0/   classpoly_v1.0.3/    (vendored: word-size F_p / class polynomials)
zp_poly/        (vendored: large-p F_p[x] library; classpoly's invariant->j needs it)
phi_files/      (28 MB subset of modular polynomials; see INSTALL)
ecpp/           (this project: discriminant search, smoothness, root-finding,
                 invariant->j, curve assembly, oneshot)
tests/          (classpoly correctness suites vs PARI)
```
