#!/bin/bash
# For each top-level winner of nextprime(2^768): presentation, candidate subgroups with height bounds, and the choice made for 192 jobs.
cd /home/claude/OneShotFastECPP && . ./setenv.sh >/dev/null
P=$(cat reports/shortecpp-improve/data/p768.txt)
out=$1; T=$(mktemp -d); : > $out
for D in $(python3 -c "
import json
P=open('reports/shortecpp-improve/data/p768.txt').read().strip()
[print(r['d']) for r in map(json.loads, open('reports/shortecpp-improve/data/exp768_winners.jsonl')) if r['p']==P]" | sort -un); do
  echo "=== D=-$D" >> $out
  ( cd $T; CLASSPOLY_TOWER=2 CLASSPOLY_JOBS=192 CLASSPOLY_JOBID=1 CLASSPOLY_ECRT_DIR=$T classpoly -$D -1 $P $T/o.txt 1 2>&1 | grep -E "Using inv|square-root|relative order|Decompos|  n=|bits =|Selected|Evaluating|Plain|Sign-free" >> $out; rm -f $T/*.ecrt )
done
rm -rf $T; echo LIST_DONE >> $out
