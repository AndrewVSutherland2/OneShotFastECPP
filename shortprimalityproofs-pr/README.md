# Contribution prepared for AndrewVSutherland/ShortPrimalityProofs

The PR could not be opened from this session: the GitHub token here is a fine-grained PAT
scoped to existing repos, so both `POST /repos/AndrewVSutherland/ShortPrimalityProofs/forks`
and `POST /user/repos` return 403, and it has no push access to the upstream repo. Everything
is ready to land as two commits against `main` (tested to apply cleanly at upstream d541190).

## To land it

```sh
git clone git@github.com:AndrewVSutherland/ShortPrimalityProofs.git
cd ShortPrimalityProofs
git checkout -b short-ecpp-tools
git am /path/to/000*.patch          # or: git pull /path/to/short-ecpp-tools.bundle short-ecpp-tools
git push -u origin short-ecpp-tools # then open the PR from the branch
```

Alternatively, granting `AndrewVSutherland2` push access (or fork permission for the token)
lets the PR be opened directly from a future session.

## What the two commits contain

**0001 — Add vsmallECPP.py, short.gp, short8all.txt and certs.csv**

| file | contents |
|---|---|
| `vsmallECPP.py` | stdlib-only short-ECPP verifier, quasi-quadratic; sound for composite input |
| `short8all.txt` | all 201,072 short ECPPs with p ≤ 2⁸ (every prime 5 ≤ p ≤ 251 admits one) |
| `short.gp` | toy GP prover: SEA on random curves, in the style of `oneshot.gp` |
| `certs.csv` | short ECPPs for the least prime above 10^c, c = 10, 20, …, 100 |

plus a README section with the primality argument, a note on why the order is taken modulo
every prime divisor (with the composite counterexample for the weaker reading), the resources
list and the challenge table.

**0002 — README typo and hygiene fixes** (separable; drop it if unwanted): `thie`→`this`,
`m_i p_{i+i}`→`m_i p_{i+1}`, `a_i`→`A_i`, `n` defined once from `p_0`, and
`A_i ≠ ±2` strengthened to `gcd(A_i²−4, p_i) = gcd(B_i, p_i) = 1` so the reduction is an
elliptic curve modulo every prime divisor of `p_i`.

## Verification performed

* `python3 vsmallECPP.py --test` — valid vectors (single-level and a 48-bit two-level chain)
  plus six tamper classes, including the composite CRT/lcm pseudo-certificate.
* All 201,072 entries of `short8all.txt` and all 10 entries of `certs.csv` re-verified from a
  clean checkout after applying the patches.
* Every `certs.csv` level independently cross-checked in PARI (`ellisoncurve`, `ellorder`).
* `short.gp` re-run from the patched tree on a fresh prime (10²⁵+13), output verified.

Provenance of these files in this repo: `ecpp/vsmallECPP.py`, `ecpp/short.gp`,
`ecpp/gen_short.py` (the p ≤ 2⁸ enumerator), `ecpp/collect_short.py`, `certs/short/`.
