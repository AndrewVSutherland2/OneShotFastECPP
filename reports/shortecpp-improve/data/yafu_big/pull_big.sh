#!/bin/bash
# pull_big.sh [nopkl]: off-box copy of the four runs (log, progress, certificate, journal, levels
# done, search checkpoint) and one status line per run.  A copy replaces the previous one only
# when the transfer completed.
S=/tmp/claude-1000/-home-claude-OneShotFastECPP/f1db9bbb-8f84-4630-ad02-28f9e40f5dc3/scratchpad
MAP=${MAP:-"7:450 15:430"}
set -o pipefail
one() {
  N=$1; C=$2; D=$S/ckpt_big/c$C; T=$D/.tmp
  mkdir -p $D; [ -d "$T" ] && find "$T" -type f -delete; mkdir -p $T
  if timeout 240 spot ssh $N "bash ~/rule/pack.sh $C $3" 2>/dev/null | tar x -C $T 2>/dev/null; then
    cp -a $T/. $D/ && echo "$(date -u +%T) ok" > $D/last_pull
    st=ok
  else
    st="PULL FAILED (kept the copy of $(cat $D/last_pull 2>/dev/null))"
  fi
  f=$D/big/c$C.log
  pk=$(ls -la --time-style=+%H:%M $D/shortcamp/c$C.search.pkl 2>/dev/null | awk '{printf "%.0fMB@%s", $5/1e6, $6}')
  echo "[$N c$C] $st; winners $(grep -c WINNER $f 2>/dev/null); ckpt ${pk:-none}; $(tail -n 1 $D/big/progress 2>/dev/null | cut -c1-90)"
  grep -v WINNER $f 2>/dev/null | grep "s: \|search:\|built\|reserve\|resumed\|chain complete\|verif" | tail -n 1 | cut -c1-250
  grep "dscan\[" $f 2>/dev/null | tail -n 1 | sed 's/.*dscan/  dscan/' | cut -c1-175
}
for m in $MAP; do one ${m%%:*} ${m##*:} $1 > $S/ckpt_big/.st_${m##*:} 2>&1 & done
wait
for m in $MAP; do cat $S/ckpt_big/.st_${m##*:}; done
