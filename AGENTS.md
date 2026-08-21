<!-- Drop-in "## Code Review Rules" section for AGENTS.md. Scope: ONLY what
     a strong reviewer cannot infer from the code -- this repo's conventions
     and this loop's protocol. The generic contract (severity gate, cite
     file:line, demonstrate errors, LGTM sentinel) lives in the workflow
     prompt; do not duplicate it here. Grow this file from observed review
     misses, like a test suite -- not prophylactically. -->

## Code Review Rules

### Claims
- Scope claims to exactly what this repo's code and data establish:
  explicit construction vs existence, unconditional vs conditional (name
  the hypothesis). Flag overclaim in comments, docs, and paper text.

### Generated artifacts
- Shipped tables (.tex/.txt) are generated: flag hand edits. Verifiers
  rebuild from structured sources and compare exactly, totals asserted as
  n of n. Every claiming entry has a machine-checked witness.

### Certificates
- Certificates are checkable by an included verifier independent of the
  pipeline that produced them; format changes re-verify everything shipped.

### Searches
- Completeness claims state bounds and every pruning condition; each
  condition is proved or independently checked and matches the claim.

### C kernels
- Check arithmetic (widths, intermediate products, overflow) at the STATED
  parameter ranges, not typical inputs. Performance claims need a recorded
  measurement.

### Reproducibility
- Checkout-relative paths; zero-argument out-of-tree runs work or the
  README says why; baselines reproducible via documented env switches.
  In paper/, numeric values trace to data files; no prose nits.

### Round protocol
- Read git log BASE..HEAD first: "Codex round N:" commits record
  resolutions. Do not re-raise an addressed finding unless the recorded
  resolution is defective -- then answer that justification specifically.
  Disputes not settled by computation: state both positions once, defer
  to the maintainer.

<!-- Repo-specific: add 1-3 lines here, e.g. "Cassels conditions cited and
     re-verified whenever sieve code changes." -->
