#!/usr/bin/env bash
# Unit-test driver for the classpoly build.  Invoked by `make test`.
#
#   ./run_tests.sh [MAXD]      # default MAXD=1000
#
# Test 1: classpoly vs PARI/GP polclass over Z (shared invariants).
# Test 2: invariant class poly roots mod a word-size prime -> j-invariants.
# Test 3: classpoly mod 2^255-19 vs the over-Z polynomial reduced mod 2^255-19.
# All over fundamental imaginary quadratic discriminants D with |D| <= MAXD.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
MAXD="${1:-1000}"
JOBS="${JOBS:-24}"

if [[ ! -x "$ROOT/classpoly_v1.0.3/classpoly" ]]; then
    echo "classpoly binary missing; building..." >&2
    make -C "$ROOT"
fi

rc=0
echo "### Test 1: classpoly vs PARI/GP polclass  (|D| <= $MAXD) ###"
python3 "$HERE/compare_pari.py" --maxd "$MAXD" --jobs "$JOBS" || rc=1

echo
echo "### Test 2: invariant class poly roots -> j-invariants  (|D| <= $MAXD) ###"
python3 "$HERE/verify_jinv.py" --maxd "$MAXD" --jobs "$JOBS" || rc=1

echo
echo "### Test 3: classpoly mod 2^255-19 vs over-Z reduction  (|D| <= $MAXD) ###"
python3 "$HERE/compare_modp.py" --maxd "$MAXD" --jobs "$JOBS" || rc=1

echo
if [[ $rc -eq 0 ]]; then echo "ALL TESTS PASSED"; else echo "SOME TESTS FAILED"; fi
exit $rc
