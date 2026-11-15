#!/bin/bash
# Time cm_method on the 17 journalled winners of nextprime(2^768), alone on the box.
# usage: cmtime.sh <outlog> [extra cm_method args]
cd ~/OneShotFastECPP && . ./setenv.sh >/dev/null
P=$(cat reports/shortecpp-improve/data/p768.txt)
out=$1; shift
for D in $(python3 -c "import json;[print(json.loads(l)['d']) for l in open('reports/shortecpp-improve/data/exp768_winners.jsonl')]" | sort -un); do
  /usr/bin/time -f "D=$D wall=%e user=%U sys=%S" timeout 3000 ./ecpp/cm_method D=-$D p=$P jobs=192 "$@" >> $out 2>&1
done
echo CM_DONE >> $out
