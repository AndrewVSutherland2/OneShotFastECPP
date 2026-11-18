#!/bin/bash
# restart_after_ckpt.sh C [now]: restart the prover of run C from its next checkpoint, so that it
# runs the code now on disk.  Without "now" it first waits for the scan chunk in flight to arrive
# (a new top-level dscan line), so that no scan work is thrown away.  Does nothing if a winner is
# in hand by then.  Run in its own tmux session on the box.
C=$1; MODE=$2
LOG=~/rule/big/c$C.log
BITS=$(grep -m1 -o 'proving p0 ([0-9]* bits' $LOG | grep -o '[0-9]* bits' | grep -o '[0-9]*')
note() { echo "$(date -u +%F_%T) restart helper: $*" >> ~/rule/big/progress; }
PID=$(pgrep -f "^python3 ecpp/shortECPP.py p=[0-9]+ threads=192 v2=1 tag=c$C " | head -n 1)
[ -z "$PID" ] && { note "no prover found for c$C"; exit 1; }
note "watching c$C (pid $PID, top level $BITS bits, mode ${MODE:-after-chunk})"
n0=$(grep -c "level $BITS bits: dscan\[" $LOG)
if [ "$MODE" != now ]; then
  while [ "$(grep -c "level $BITS bits: dscan\[" $LOG)" -le "$n0" ]; do
    kill -0 $PID 2>/dev/null || { note "prover gone while waiting for the chunk"; exit 0; }
    sleep 20
  done
fi
k0=$(grep -c "checkpoint:" $LOG)
while [ "$(grep -c "checkpoint:" $LOG)" -le "$k0" ]; do
  kill -0 $PID 2>/dev/null || { note "prover gone while waiting for a checkpoint"; exit 0; }
  sleep 5
done
if grep -q "WINNER\|level $BITS bits: built" $LOG; then note "a winner is in hand: no restart"; exit 0; fi
note "restarting c$C after: $(grep 'checkpoint:' $LOG | tail -n 1 | cut -c1-100)"
kill -KILL $PID
sleep 3
pkill -x dscan; pkill -x yafu; pkill -x smoothtest
for i in $(seq 1 60); do tmux has-session -t big 2>/dev/null || break; sleep 2; done
pkill -x yafu; pkill -x dscan
tmux new -d -s big "bash ~/rule/run_big.sh $C"
note "relaunched c$C"
