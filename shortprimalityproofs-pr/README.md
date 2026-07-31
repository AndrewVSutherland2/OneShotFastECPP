# Files to copy into AndrewVSutherland/ShortPrimalityProofs

Everything in `files/` goes at the repo root. Four are new; `README.md` replaces the existing one.
Built and tested against upstream `d541190`.

| file | new? | what it is |
|---|---|---|
| `files/short.gp` | new | the GP prover (SEA on random curves), in the style of `oneshot.gp` |
| `files/vsmallECPP.py` | new | short-ECPP verifier, quasi-quadratic, Python stdlib only |
| `files/short8all.txt` | new | all 201,072 short ECPPs with p ≤ 2⁸, one per line, space-separated |
| `files/certs.csv` | new | the ten certificates for the least prime above 10^c, comma-separated |
| `files/README.md` | replaces | adds the resources list, challenge table, and the primality argument |

## The README changes

Everything the existing README already said is preserved; the additions are the primality
argument, a note on why the order is taken modulo every prime divisor (with a composite
counterexample for the weaker reading), the resources list, and the challenge table.

Five edits touch the existing definition block:

1. `thie` → `this`
2. `n:=\lceil \log_2 p_i\rceil` → `p_0` (the first bullet said `p_i`, the third said `p_0`)
3. `m_ip_{i+i}` → `m_ip_{i+1}` (and the now-redundant restatement of `n` dropped from that bullet)
4. `a_i\ne \pm 2\bmod p_i` → `\gcd(A_i^2-4,p_i)=1` (also fixes the lowercase `a_i`)
5. `B_i,y_i\in [0,p_i-1]` → the same with `\gcd(B_i,p_i)=1`

Edits 4 and 5 are the only substantive ones, and are easy to drop if unwanted: they make the
reduction an elliptic curve modulo *every* prime divisor of `p_i`, which is where the order
condition and the Hasse bound get applied — with `A_i ≡ ±2 mod ℓ` for some `ℓ | p_i` the
reduction mod `ℓ` is a singular cubic, and "order" is not defined if the point reduces to the
singular point. For prime `p_i` the two forms agree (`gcd(A²−4,p)=1` ⟺ `A ≢ ±2 mod p`, and
`gcd(B,p)=1` ⟺ `B ≢ 0`), so no valid certificate is affected either way; `vsmallECPP.py`
checks the gcd form regardless, as `voneshot.py` does for `A`.

## Verification performed

* `python3 vsmallECPP.py --test` — valid vectors (single-level, and a 48-bit two-level chain)
  plus six tamper classes, including the composite CRT/lcm pseudo-certificate.
* All 201,072 lines of `short8all.txt` and all 10 lines of `certs.csv` verify.
* Every level of every `certs.csv` entry independently cross-checked in PARI
  (`ellisoncurve`, `ellorder`).
* `short.gp` re-run from the assembled tree on a fresh prime (10²⁵+13); output verified.

Sources in this repo: `ecpp/short.gp`, `ecpp/vsmallECPP.py`, `ecpp/gen_short.py` (the p ≤ 2⁸
enumerator and chain builder), `ecpp/collect_short.py`, `certs/short/`.
