#!/bin/bash
# relaunch.sh N C: bring aws-spot-N (claimed and started) up as the new home of run C and resume it
# from the off-box copy in ckpt_big/cC (search checkpoint, winners, levels done, log so far)
N=$1; C=$2
S=/tmp/claude-1000/-home-claude-OneShotFastECPP/f1db9bbb-8f84-4630-ad02-28f9e40f5dc3/scratchpad
cd $S/big
./setup_box.sh $N
spot ssh $N 'python3 -c "import numpy" 2>/dev/null || sudo apt-get install -y python3-numpy > ~/rule/apt.log 2>&1; python3 -c "import numpy; print(\"numpy\", numpy.__version__)"'
./ship.sh $N
tar c pack.sh | spot ssh $N 'tar x -C ~/rule'
echo "[$N] smoke test"
spot ssh $N 'cd ~/rule/repo && . ./setenv.sh > /dev/null && P=$(python3 -c "print(10**160+303)") && python3 ecpp/shortECPP.py p=$P threads=192 v2=1 tag=smoke out=$HOME/rule/smoke.cert.txt sched=rule engine=yafu yafu=$HOME/yb/yafu/yafu > ~/rule/smoke.log 2>&1; echo "smoke exit $? verified=$(grep -c "verified twice" ~/rule/smoke.log)"'
echo "[$N] state of c$C"
tar c -C $S/ckpt_big/c$C shortcamp | spot ssh $N 'mkdir -p ~/rule/repo/work && tar x -C ~/rule/repo/work && ls -la ~/rule/repo/work/shortcamp | grep -v smoke'
tar c -C $S/ckpt_big/c$C big | spot ssh $N 'tar x -C ~/rule && echo "$(date -u +%F_%T) moved to $(hostname) after an interruption" >> ~/rule/big/progress'
spot ssh $N "tmux new -d -s big 'bash ~/rule/run_big.sh $C' && sleep 20 && tmux ls && tail -n 3 ~/rule/big/progress && grep -n 'resumed\|Traceback\|Error' ~/rule/big/c$C.log | tail -n 3 | cut -c1-330"
echo "[$N] $(date +%T) RELAUNCH_DONE"
