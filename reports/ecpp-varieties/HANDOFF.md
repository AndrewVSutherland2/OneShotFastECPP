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

**Spec thread, rounds 3-4 (2026-08-20/21) — CONVERGED DESIGN**: Drew round 3:
tie BOTH knobs to one value: smoothness bound Y = 2n²/log₂n AND rad(m_i) ≤ Y
as an INTEGER, payload unchanged.  Verified: everything O(n² log n)
self-contained (primorial of Y builds in O(n² log n) — the log² died with the
smoothness shrink; gcd(P_Y mod o, o) = rad(m); radical ≤ Y factors by trial
division to √Y; k ≤ log₂Y ≈ 2log₂n so the tree is noise) — but compliance
catastrophic: 2/31 chains (radicals up to 2^57 vs caps 2^15-2^18; needed c up
to 1.6e12).  Taxonomy: 21/120 levels fail via a filler prime in (Y, n²] (16
chains — fail ANY Y≈n²/log n design); 30 more via the integer radical cap.
Key dichotomy established: {unchanged payload, keep published certs, fully
self-contained O(n²log n)} — pick two; self-containment with unchanged
payload FORCES Y = O(n²/log n) (primorial build is Θ(Y·log²)-ish otherwise;
on-the-fly chunking is also log²).  Drew round 4: cap log₂rad(m) instead of
rad(m) — i.e. rad(m_i) < 2^{c·n/log₂n} with Y = c·n²/log₂n.  Verified: still
O(n²log n) (tree ≤ (cn/log₂n)·log₂n = cn ladder bits/level; the big radical
is factored by ONE batched remainder-tree pass across the chain vs the ≤Y
prime tree, O(n²log n) total — or per-level Pollard–Strassen radius Y,
Õ(n²√log n), tree-free).  Bits-cap is non-binding on ALL published data at
c ≥ 3 (max radical 57 bits vs cap ≥ 90); the binding constraint is now the
smoothness bound (filler primes > Y).  COMPLIANCE (both knobs at c):
c=1: 10/31, c=2: 14/31, c=3: 20/31, **c=4: 24/31 with records 9/11 passing
and the 2 failures (10^240 lev2 ~192b, 10^280 lev3 ~114b) trivial**, c=6:
28/31.  At c=4 the full re-find list: 10^90 lev0 (299b, orig ~350 s — the
only level-0), 10^20 lev1 + 10^60 lev1 (toys), 10^130/10^140 lev2, 10^240
lev2, 10^280 lev3 — ALL local-box, no spot campaign.  RECOMMENDED FINAL:
Y = ⌊4n²/log₂n⌋, rad(m_i) < 2^{⌈4n/log₂n⌉}, payload/floor/window unchanged.
Next steps if Drew confirms: patch vsmallECPP.py to the new design, re-find
the 7 broken chains locally, update paper §6/Remark 6.3 wholesale.

**Spec thread, round 5 (2026-08-21) — DECISION + EXECUTION (Drew AFK 8h,
directed "proceed and update the specification to match")**: Drew rejected
the ln normalization ("funny to take natural log of a binary log") and chose
**c = 1: B = ⌊n²/log₂n⌋, radical cap ⌊log₂ rad(m_i)⌋ < n/log₂n**, floor
p' > n² explicit, payload unchanged.  Compliance at c=1: 10/31 chains pass
(c = 70,150,160,180,200,250,260,270,290,310), 21 repair, incidental per-level
74% (⟹ finder penalty ≈1.35×).  BUILT + VALIDATED (commit 7f01659):
ecpp/vshort2.py (v2 reference verifier, primorial-gcd recovery; self-tests;
differential = exactly the predicted 10/21 split on all 31 chains),
ecpp/short2.gp (adapted finder, shortcert2from(p,ntop) for repairs),
ecpp/repair_short2.py (prefix-preserving migration driver, dual-validates
v2+v1 — note v2 ⊂ v1: old verifier accepts all v2 certs),
short_prove.py v2=1 flag (CM prover with v2 filler filter + short2.gp tails;
smoke-tested at 200 bits, 46 s, dual-verified).  Paper updated (Def 6.2
def:short2, Prop prop:short2verif, Remark 6.3 compressed to resolved-design
note, Table 1 short row = O(n²log n)† revision note, hierarchy adjusted
(restricted ⊄ v2-short), §6.4 gains the RECORDS TABLE c=210..310 (the paper
had been stale — records merged upstream 2026-08-18!) + migration paragraph
[placeholders SMALLREPAIRS/BIGREPAIR to fill]).  reports/short2-spec/
index.html = upstream adoption package (README clause text, checklist,
[MIGRATION-RESULTS-TABLE] to fill).  REPAIRS RUNNING on the local box:
driver-2 = 20 chains ex-698 (work/short2repair2/, sequential big-first);
698 = short_prove v2 (work/short2repair/short_prove_698.log) after the SEA
attempt was killed as redundant (also: ~21 vCPU of ANOTHER session's jobs
(smallcensus, backtrack_rat_d) share the box — do not kill).  **MIGRATION COMPLETE (2026-08-21 ~16:15)**: all 31 chains in
certs/short2/certs.csv, every one verified under BOTH v2 and v1
(assemble_short2.py).  19 repairs by short2.gp: 650 s aggregate (max 254 s);
the two level-0 rebuilds by short_prove v2=1 on the shared box: 10^170 =
15,106 s / 6 threads (22,654 lev-0 cands), 10^210 = 47,373 s / 10 threads
(33,876 cands, ~130 ecm core-h, winner D=-501101320 h~22385 — its v1 find
had been a 20 s luck fluke; ordinary-luck price paid this time).  Paper and
spec placeholders filled, compiled clean.  UPSTREAM PACKAGE PREPARED (2026-08-21, evening): branch `radical-cap` on the
AndrewVSutherland2/ShortPrimalityProofs fork (local clone
/home/claude/ShortPrimalityProofs), two commits on top of upstream 0408b2b:
fd5f613 = pure revert of short.gp to Drew's pristine 15a87c8 version (the
IslayResearch optimizations from ca3ccd1 — twists, backtracking, factorint
flags, sctryorder refactor — removed per Drew: short.gp is a demonstration,
simple not fast; ONE deliberate deviation: a parisizemax 2^32 line, without
which the usage examples crash at 10^30 absent a gprc); 1cdbf11 = the v2
revision: README (definition clause + revision paragraph + resources + all
21 repaired <details> entries with updated summaries/level counts),
vsmallECPP.py (ported v2 verifier, upstream-style header, new self-tests),
short.gp (+v2 bounds on the simple base; smoke-tested to 10^60, chains
verify), certs.csv (migrated 31), short8all.txt (v2 exhaustive: 55,056
tuples, generator validated by reproducing the v1 201,072 EXACTLY under v1
rules; all 55,056 + all 31 chains verified by the new verifier),
parallel_short.py (dead knobs dropped; smoke-tested end-to-end).  Drew opens
the PR: https://github.com/AndrewVSutherland2/ShortPrimalityProofs/pull/new/radical-cap
— arm watch-codex at PR-open.  After upstream merge: bump the paper's
[Short] bib pin from 0408b2b to the new commit.

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

## Spec thread, round 6 (2026-08-22): the terminal-prime revision (v2.1)

Drew (post-v2-merge, PR #6 = upstream cc2caf1): raise the recursion floor
n² → B² and allow a terminal prime p_{k+1} < B² (self-certifying: B-rough
below B² is 1 or prime).  Analysis confirmed + measured: 21/31 chains
shorten by exactly one level by PURE TRUNCATION (129→108 levels; absorbed
primes 15-32 bits; 10^10 becomes single-level); truncation is FORCED (old
long forms invalid); short8all provably unchanged (v2.1 enumeration below
2^8 = v2's, equal count + containment); finder simplifies (channel (a'),
ispseudoprime deleted).  Drew's spec catch: p_{k+1}'s B-roughness must be
EXPLICIT (consequence of other conditions for recursive p_i, not forced for
p_{k+1}; without it, junk like p_{k+1}=9 and radical-cap dodges pass an
exists-reading).  Final wording (Drew's): "the p_i are odd integers and
p_1,...,p_{k+1} have no prime factors ≤ B" — scoped to i ≥ 1 because p_0=5
has B=5 (n=3) and would be excluded!  Fork branch `terminal-prime`
(85a2215 + 6146fd1 pycache cleanup — NB the v2 PR accidentally committed
__pycache__, now removed + .gitignore): README clauses, verifier (new
terminal/floor logic + tests incl. forced-truncation rejections), certs.csv
truncated, short.gp channel (a') [bug found+fixed: (a') needed the m>1
guard — factor(1)[1,1] errors when R > L at small moduli], parallel_short
untouched.  Validated: 31/31 truncated verify; 55,056 short8all verify;
21/31 old forms reject; fresh chains 10^10/30/60 verify (1/2/3 levels).
Local side: vshort2.py ported to v2.1, certs/short2/certs.csv truncated,
paper Def 6.2 rewritten to final form + Prop sketch + migration truncation
sentence + [Short] pin → cc2caf1 (RE-BUMP after the terminal-prime PR
merges).  PR link: github.com/AndrewVSutherland2/ShortPrimalityProofs/pull/new/terminal-prime

## Spec thread, round 7 (2026-08-22 afternoon): conventions settled, upstream-ready

Drew merged terminal-prime as FORK PR #1 (fork-first review flow; upstream PR
comes after fork review).  His README edits: B := ⌈n²/log₂n⌉ (ceiling), and —
after a slip that made the cap "log₂rad ≤ B" (VACUOUS: log₂rad < n/2 + log₂B
≪ B always; caught+flagged, short8all ballooned to 162,042 as evidence) — the
FINAL radical cap: **log₂ rad(m_i) ≤ ⌈n/log₂n⌉**.  Key fact: for squarefree
rad this is EXACTLY bit_length(rad) ≤ ⌈n/log₂n⌉ = the old floored test's
accept-set (squarefree 2-power boundary can't arise), so NOTHING migrates:
certs.csv unchanged (31/31, the 0.046-bit 10^20 marginal now clears by ~1
bit), short8all = the same 55,056 file (ceil-B's new smooth prime 11 at n=5
can't enter any filler without busting the cap — proven + regenerated
byte-identical).  Verifier renamed vsmallECPP.py → **vshortECPP.py** (Drew;
parallel_short + README refs updated).  Fork main = 6d100e5 all-green (tests,
31/31, 55,056/55,056, fresh chains).  Local: vshort2.py ported (ceil-B,
radlim = ⌈n/lg⌉, bit-length test); paper Def 6.2/Prop sketch updated to the
settled forms, compiles clean.  REMAINING: Drew reviews fork main → opens the
upstream PR → after merge, bump the paper's [Short] pin (currently cc2caf1)
and note the verifier rename if the paper ever names it (it doesn't).
