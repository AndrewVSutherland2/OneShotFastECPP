# Tests for the classpoly build

Run everything with `make test` from the project root, or directly:

```sh
./run_tests.sh [MAXD]      # default MAXD = 1000
```

## Test 1 — `compare_pari.py`: classpoly vs PARI/GP `polclass`

For every fundamental imaginary quadratic discriminant `D` with `|D| <= MAXD`,
and every class invariant supported by **both** classpoly and PARI/GP, we
compute the class polynomial `H^inv_D(X)` over `Z` with each tool and compare
the coefficient vectors **exactly**.

The invariant codes coincide between the two tools (both come from the
Enge–Sutherland framework) except for the double-eta quotients `w_{5,7}` and
`w_{3,13}`, where PARI uses `35`/`39` and classpoly uses `535`/`539` (the
`500 + p1*p2` encoding). The `INVARIANTS` table in the script records the
`(pari_code, classpoly_code)` pairs.

A `(D, inv)` pair is only **compared** when both tools produce a polynomial.
When exactly one tool accepts it we record a *validity difference* — the two
libraries impose slightly different domain conditions on some invariants; this
is expected and reported separately (see below), not as a failure.

### Result over `|D| <= 1000`

```
compared (both tools produced a polynomial): 3037
  exact matches : 3037
  MISMATCHES    : 0
validity differences: 17 PARI-only, 42 classpoly-only
```

All 26 shared invariant families agree on every comparable discriminant
(`j`, Weber `f`/`f^2`, `gamma_2`, the single/double-eta quotients including
`w_{5,7}` and `w_{3,13}`, and the Atkin invariants `A_3 … A_31`).

The validity differences are all benign domain-heuristic differences:

* **42 classpoly-only**: all are `w_{3,13}`. classpoly computes the minimal
  polynomial of this double-eta even when PARI's "generates the full ring class
  field" condition fails. The degree is `h` when `w_{3,13}` generates the field
  (e.g. `D=-87`, `h=6`, degree 6) and `h/2` in the Fricke square-root cases
  where both 3 and 13 ramify (e.g. `D=-195=-3·5·13`, `h=4`, degree 2). These are
  independently checked by Test 2 (PARI cannot validate them).
* **17 PARI-only**: classpoly's `inv_good_discriminant` is stricter than PARI
  for a handful of Weber-`f` and double-eta cases, plus `D=-4` `gamma_2`
  (classpoly only computes `j` for `D >= -4`).

## Test 2 — `verify_jinv.py`: invariant roots → j-invariants

PARI cannot validate the invariants classpoly supports but it does not (Atkin
`A_41/A_47/A_59/A_71`, the Ramanujan `t`/`t^2`, the single-eta quotients, and the
`w_{3,13}` subfield/Fricke cases). For those we check classpoly against the
Hilbert class polynomial directly. For each fundamental `D` (class number `h`):

1. pick a word-size prime `p` that splits completely in the ring class field of
   `D` (i.e. `H_D` splits into `h` distinct linear factors mod `p`);
2. `JH` := the `h` roots of `H_D` mod `p` (the CM j-invariants), via PARI;
3. for each invariant, compute `H^inv_D(X)` mod `p` with classpoly, find its
   roots `{x_i}` mod `p`, and map each `x_i` to the j-invariant(s) `J` with
   `Phi_inv(x_i, J) = 0` using classpoly's own f→j modular polynomials (via the
   `invtoj` helper, which calls the new `ff_all_j_from_inv`);
4. check every `x_i` maps to at least one element of `JH`, and that the mapped
   j-invariants cover `JH` exactly.

This pins down `H^inv_D` as a genuine class polynomial for `D`: its roots are the
invariant values at the CM points, mapping onto precisely the CM j-invariants.

### Result over `|D| <= 1000`

```
checks run : 1450
  passed   : 1450
  FAILED   : 0
discriminants skipped (no split prime in window): 0
```

All invariants pass, including the controls (`gamma_2`, `w_{2,3}`, `A_3`, `A_31`,
`w_{5,7}`) and the PARI-unsupported `t`, `t^2`, single-eta `w_3^12/w_5^6/w_7^4/
w_13^2`, `w_{3,13}` (77/77 — including the 42 cases PARI refused in Test 1), and
Atkin `A_41` (120), `A_47` (113), `A_59` (115), `A_71` (109).

The `invtoj` helper and `ff_all_j_from_inv` (added to `class_inv.c`) are also the
building blocks for going from a class-polynomial root to a curve in the ECPP
pipeline.

## Test 3 — `compare_modp.py`: mod-p output vs the over-Z reduction

classpoly can emit a class polynomial over `Z` or reduced modulo a given `P` (the
explicit-CRT / ECRT path). ECPP uses the mod-`P` path with `P` the large prime
being certified, so it must be exact. Test 1 anchors the over-Z output to PARI;
this test anchors the **mod-P output** to the over-Z output. For

```
p = 2^255 - 19   (255 bits -- well beyond word size, so it exercises the
                  multi-word ECRT reduction)
```

and every fundamental `D` with `|D| <= MAXD` and every supported invariant, it
checks `(over-Z coefficients) mod p == (mod-p coefficients)` term by term.

### Result over `|D| <= 1000`

```
compared : 3792
  matches : 3792
  MISMATCHES : 0
over-Z-only : 0    mod-p-only : 0
```

`D = -3, -4` are skipped (classpoly emits the linear Hilbert polynomial over Z
only for `D >= -4`). Everything else agrees exactly, validating the large-modulus
path against the (PARI-anchored) over-Z path.

## Note on parallelism

classpoly used to derive its CRT scratch filenames (`CRT_*.crt`) deterministically
in a single shared directory, so two processes sharing `CLASSPOLY_TMPDIR` would
clobber each other's scratch files and abort. This is now fixed: `crt_dir()`
creates a **private per-process scratch directory** with `mkdtemp()` under
`$CLASSPOLY_TMPDIR` → `$TMPDIR` → `/tmp` (removed at exit), so concurrent classpoly
processes never collide. The harness relies on this and launches all jobs from a
shared working directory with no per-task setup.

## Environment

The harness sets the three directory hooks itself; for interactive use,
`source ../setenv.sh` from the project root points them into the project tree.
