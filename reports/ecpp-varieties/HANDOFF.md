# HANDOFF — the ecpp-varieties write-up (read this first in a new session)

**Current main task (2026-08-19): Drew has a detailed review of this write-up
from GPT 5.6 Sol and wants to work through it.  That review discussion is the
next unit of work.**  Everything below is the context a fresh session needs.

## Sol-review progress (2026-08-19 session)

Sol's review is committed here as `fable_ecpp_review.{tex,pdf}`.  Triage went to
Drew as three lists; **list 1 (16 agreed items) is applied and committed**:

- CM winner-classpoly cost fixed in §4.3.2 (Õ(|D|) = p^(1/4), only degree/root
  at p^(1/8)); open problem 3 rewritten — output-sensitive enumeration alone no
  longer closes the 1/8-vs-1/4 gap, a sub-Õ(|D|) commit is also needed.
- 10^27 Pomerance estimate: 170 GPU-h (was "ten"); restricted-256: 5–10
  core-years = days on 3000 cores (was "months"); both also in §8.
- "Always arrangeable" now correctly via prime-power fillers (k_i = 1 trivial
  exactness tree), not poly-size fillers (table dagger, Prop 6.2, intro item 4).
- Remark 6.3 repairs: trial division Ω(n³/log n); Pollard–Strassen with search
  radius B = n² → Õ(n)/level for EVERY fixed c (stronger than Sol's c ≤ 8 fix).
- b→1 limit demoted to qualitative; AH92 = genus-2 Jacobians (not this format);
  L(p) = upper bound (not attained max); "exactly two known ways" and other
  absolutes hedged; codim-(n/2) → probability 2^(-n/2+o(n)).
- Sizes: one-shot ≈2.5n–3n (Σlog q_i ≤ log m); short ≲5n — measured 3.2n–4.7n
  over the 21 local chains (c=10..100 + 210..310).
- Hierarchy paragraph → partial order (n⁴ one-shot vs short incomparable;
  classical not a superset; existence scoped to Pomerance-containing formats).
- AKS: table cell labeled "det. decision"; find/verify-inversion sentence added;
  quartic "floor" softened to known-verifier floor (3 spots).
- 40 bib entries got Crossref-VERIFIED DOIs; challenge repos pinned (DANGER3
  3cd7be6, OneShot 188f53e, Short 0408b2b, accessed 2026-08-19); OSF cited by
  branch; \texorpdfstring in 6 math titles (0 bookmark warnings); date → Aug 19.

**Round 2 (Drew's verdicts applied)**: A ✓ Table-1 short cell now reads
"O(n²(log n)²); O(n² log n)†" (worst case first, dagger reworded).  B ✓
rebut-and-clarify sanctioned and done: Prop 6.4's statement now names "the
budgeted-ECM search above (channels (1)–(3))", and a paragraph after (H1)–(H3)
proves the exposure equivalence (S = N/q y-smooth ⟺ q exposable; discarded
cofactor arbitrary within y-smoothness; the q ≤ y large-filler channel needs an
n²-smooth part > √p/y, probability p^(-1/4+o(1)) — NB round-1 note said
p^(-1/8), wrong — negligible vs 1/L(1/3); alarm-budget = soft threshold,
constants only).  C ✓ minimal split: §7.2 typical-output-vs-worst-case
sentences (strict 2q: Θ(n²) size, O(n³ log n) verification) + Table-1 caption
clause; no separate rows (Drew).  **F ✓ AKS Appendix A REMOVED entirely
(Drew's call)** — his rationale, now IN Remark 1.1: finding is cheaper than
verifying, so a Las Vegas verifier regenerates certificate data from the empty
string in the same total time (MR witnesses dispatch composites, Mil76/Rab80
now cited); the certificate buys determinism only, size is beside the point.
Remark 1.1 is the sole self-contained AKS treatment (tuple shape, S={1} with
d,e params, concrete 2^1024+643 instance, quartic-floor argument absorbed from
A.4, Ber98 kept for the perfect-power test).  §8 item 8 reframed three-way →
find/verify two-way (size ≤ verification time always; small certs can be
worthless).  Uncited-since-ever bibitems Bach90 + BLS12 removed.  **Round 3 (D/E resolved)**: E = prose only, confirmed (the partial-order text
from round 2 stands, no diagram).  D = keep all of §6 EXCEPT the 1.31/1.39
archaeology: Remark 6.5 rewritten as "how much the constant means" (convention
warnings kept, anonymized; old constants and "earlier revisions" language
gone), §6.5's pointer to the prototype's coarser convention dropped.  Also:
PR #4 merge reflected in Remark 6.3 (past tense; merged Aug 2026); and a REAL
overclaim in Remark 6.3 fixed — the factored-fillers paragraph promised
worst-case O(n² log n), but factoring alone does not tame the exactness tree
(a property of the group-theoretic checks, not of how factors were learned):
now states cap-dependence, and the three-regimes sentence says "capped AND
factored ... outright".

**Spec question (Drew, 2026-08-19): what format change makes Remark 6.3
irrelevant for short proofs?**  Answer delivered in conversation: (1) list the
filler data in the payload (kills the factor-recovery source: no primorial, no
precomputation; roughness of p_{i+1} needn't be checked — descent induction);
(2) bound the filler's distinct-prime count — cleanest k_i = 1, prime-power
fillers m_i = q_i^{e_i}: payload delta is ONE ≤2·log₂n-bit integer per level
(A_i, x_i, o_i, q_i), parse forced via e_i = v_{q_i}(o_i), window keeps
r_i = q_i.  Worst-case O(n² log n), O(n) memory, ≈5× Pomerance constant.
KEY: a size cap m_i ≤ n^c alone would BREAK unconditional existence (the
Pomerance 2^κ fallback violates any polynomial cap) — prime powers keep it
verbatim.  Finder: unaffected at L(1/3) scale (tiny prime-power fillers, esp.
2^a, complete the window; O(1)-factor candidate penalty).  Migration cost
measured on the 21 local chains (79 levels): 35% of levels are already
prime-power, 63% have k ≤ 2, but only 1/21 chains conforms wholesale — a k=1
spec means re-proving the table.  Restricted one-shot: (log n)² is intrinsic
(unbounded-k fillers are its content).  List 3 = rebuttals, no action.

**Spec thread, round 2 (2026-08-20)**: Drew refines k_i = 1 to "squarefree
part of m_i bounded by n/log n" — read as: the RADICAL rad(m_i) = prod of
distinct primes has ≤ n/log₂n BITS.  Agreed, and it is the sharp knee:
peel p' first, [o/R]P ≤ b_i/2 bits, tree = rad_bits·⌈lg k⌉ ≤
(n/log₂n)·log₂n = n ladder bits/level ⟹ O(n·M(n)) = O(n² log n) total.
Subsumes k=1 for n ≳ 128; keeps 2^κ existence (rad = 2).  Caveats: (a)
must be the radical, not m/□ (square fillers evade); (b) fixes exactness
only — listed-primes payload still wanted for the recovery side (listing ≤
rad_bits ≤ n/log n bits/level, O(n) total); (c) small-n floor needed: for
n ≲ 128, n/log₂n < 2log₂n = one legal prime's size (the n=34 "violation"
is a legal k=1 level).  COMPLIANCE over the 31 upstream chains (120
levels, certs.csv @0408b2b): 9 levels in 8 chains violate rad_bits ≤
n/log₂n, ALL at n ≤ 333 (c ≤ 100 era); worst = 10^50+151 level 0 (46-bit
radical vs 22.6, needed c = 2.03; also the ONLY level violating the exact
functional rad_bits·(1+⌈lg k⌉) ≤ n: 184 > 167); next-worst c = 1.53.
From n = 366 up all levels sit at c ≤ 0.83; the c ≥ 210 records have ≥ 2×
margin.  n/ln n instead: 3 violations (n = 34, 100, 167).  A slack
constant retains everything: 2n/log₂n loses only 10^50 (by < 1 bit),
2.1n/log₂n (or 3n/log₂n) retains all 31.  Re-proving the 8 offenders
(all seconds-to-minutes-scale finds) is trivial anyway.

## The deliverable

`reports/ecpp-varieties/ecpp-varieties.tex` — a survey, written for Drew
(A. V. Sutherland) for a general mathematical audience, of algorithms and
complexity for the elliptic-curve primality-certificate families:

1. Pomerance proofs (github.com/AndrewVSutherland/DANGER3),
2. one-shot ECPP (github.com/AndrewVSutherland/OneShotPrimalityProofs),
3. the restricted n^2-smooth one-shot (the k = 0 stratum; Drew's variant),
4. short ECPP (github.com/AndrewVSutherland/ShortPrimalityProofs),
5. classical ECPP,

each with: problem statement, finding algorithms + heuristic complexity (SEA
and CM supplies treated in parallel throughout), verification algorithm +
rigorous bit-complexity, best results to date.  Single self-contained amsart
file, no BibTeX; build with two `pdflatex` passes; the compiled PDF is
committed alongside and should be recompiled and re-committed with any edit.
The old HTML version is retired (index.html here is a redirect stub); **the
.tex is the sole maintained document**.  Everything lives on branch
`short-ecpp`, pushed to github.com/AndrewVSutherland2/OneShotFastECPP.

## House conventions for the document (Drew-directed)

- All complexity in **bit operations**; multiplication counts only where the
  exact constant is the point (Pomerance's (5/2+o(1)) log2 p).  M(t) = t log t.
- n = ceil(log2 p); L(p) = q+1+floor(2 sqrt q), q = floor(sqrt p).
- Summary Table 1 lists finding costs SEA-first; no descent-chain row; Pratt
  finding = L_p(1/3,(64/9)^{1/3}) via NFS; AKS = the AKS–Bernstein certificate
  row (Remark 1.1 ONLY — Appendix A removed 2026-08-19 per Drew, Round 2 above); ECPP row: (log p)^{6+o(1)} SEA (Drew agreed
  with the strict-downrun accounting), (log p)^{4+eps} CM.

## Review-hardened corrections already made (do NOT regress these)

- LP19's c_0 is effectively computable but was never computed; sharpest
  citable deterministic bound is (log p)^6 (2+log log p)^{c_0} (quoted in
  Remark 1.1).  AKS–Bernstein certificate: S = {1} always suffices (Ber07
  Thms 5.1–5.3), size (log p)^{1+o(1)}, verification (log p)^{4+o(1)} (a
  structural quartic floor: e = Omega~((log p)^2) forced by the counting
  argument), randomized finding (log p)^{2+o(1)}.  [Appendix A, which held the
  definition/existence/7-step verifier/trade-off discussion, was REMOVED
  2026-08-19 per Drew — the essentials live in Remark 1.1; see Round 2 above.]
- Short-ECPP verification: Prop 6.2 (peel the big prime first; cost
  parametrized by filler size/prime count), Remark 6.3 (tightened specs:
  factored fillers OR cap m_i <= n^c; empirical compliance of the 20 published
  chains — largest filler 50 bits at 10^90 level 0, max exponent n^{6.11} at
  10^50 level 0, c = 7 retains all; per-level n_i convention breaks 23/70
  levels across 17/20 chains — data generated by tools/checkcap.gp), and the
  all-inclusive worst case O(n^2 (log n)^2) for the uncapped format (batched
  across the chain; the reference verifier's per-level organization + its
  sequential batch products did NOT meet this as implemented — patched
  verifiers in ../../verifier-batching-pr/; upstream PR
  AndrewVSutherland/ShortPrimalityProofs#4 MERGED as of 2026-08-19: upstream
  main@0408b2b (= the commit pinned in the paper's bibliography) carries our
  vsmallECPP.py byte-identical; paper text updated to past tense).
- **Prop 6.4 + Remark 6.5 (latest, 2026-08-19)**: the L(1/3) finding constant
  for short ECPP / ratio-b chains is (3 sigma)^{1/3}, sigma = 1-b — i.e.
  (3/2)^{1/3} ~ 1.14 at b = 1/2 — for BOTH supplies: with ln ln y =
  (2/3) ln ln p tracked correctly the CM scan coefficient (8 sigma/9)^{1/3} is
  2/3 of the total and never binds.  The earlier constants 1.31 (SEA) and 1.39
  (CM), from the descent-chain session report, used ln ln y ~ ln ln p (a
  (3/2)^{1/3} inflation) and an imposed, non-binding scan/peeling balance.
  Full derivation with heuristics (H1)-(H3) is the proof of Prop 6.4.
- Corrected en route (watch for stale echoes): "largest filler 35 bits" was
  wrong (only the 10^200 chain had been sampled); GK/SEA total is
  (log p)^{6+o(1)} not 5+eps; one-shot verification is the worst case over
  certificates (k-parametrized refinement in 6.2).

## Reproducibility of every number in the paper

- Dickman rho values: `tools/dickman_rho.py` (log-space RK; beyond u ~ 9 use
  published tables — see its docstring).
- 20-chain compliance data (Remark 6.3): `tools/checkcap.gp` — run with a
  clone of ShortPrimalityProofs; **pipe on stdin** (`gp -q < checkcap.gp`),
  never `gp file.gp` (hangs interactive).
- Restricted one-shot certificates (§5.3, first generic ones): finder
  `ecpp/restricted.gp` (oneshot.gp with bound n^2), certs + search stats in
  `certs/restricted/nextprime10c.txt`, parallel driver `tools/drive.sh`.
  10^60+7 was never completed (predicted 10^2–10^3 core-h; stopped at 4.8e3
  curves).  All verified by voneshot.py with empty q-list.
- Race/chain measurements: quoted from `reports/sea-crossover/` and
  `reports/descent-chain/` (historical session reports, superseded on the
  L(1/3) constants by Prop 6.4 — do not re-import their constants).
- Verifier benchmarks/differential tests (verifier-batching-pr/README.md).

## External sources used (re-fetch as needed; scratchpad clones are gone)

- The three challenge repos (definitions + reference verifiers + leaderboards)
  — clone fresh; upstream vsmallECPP.py was last seen unchanged 2026-08-17.
- Pomerance 1987: math.dartmouth.edu/~carlp/PDF/paper62.pdf (scan is pp.
  315–320 only; the 7/2→5/2 remarks on pp. 321–322 are missing from it).
- Bernstein quartic: cr.yp.to/primetests/quartic-20060914-ams.pdf.
- Lenstra–Pomerance Gaussian periods: math.dartmouth.edu/~carlp/aks041411.pdf.
- ECPP records: t5k.org/top20/page.php?id=27 (R(109297), May 2025).

## Open/pending items beyond the review

- PR #4 (batched vsmallECPP.py) MERGED upstream (confirmed 2026-08-19); the
  companion voneshot.py patch (verifier-batching-pr/) still awaits a fork of
  OneShotPrimalityProofs if Drew wants it as a PR.
- OneShotFastECPP PR #12 (this whole branch) was closed unmerged; landing the
  branch on main is undecided.
- Open problems list in §8 of the paper is the research agenda (10^27+103
  Pomerance triple; the p^{1/8} vs p^{1/4} scan gap; 256-bit restricted cert;
  the three-way trade-off).

## Working style expected (from this session's experience)

Verify every review claim against sources or computation before accepting or
rebutting — several confident claims in earlier drafts fell to exactly this
discipline (the 35-bit filler, the 1.31/1.39 constants).  Drew reads on his
phone via the refserve links; send updated PDFs with SendUserFile and keep
commits on short-ecpp pushed.
