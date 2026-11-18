#!/bin/bash
# setup_box.sh N: repo at HEAD + full phi_files + gmp-ecm package + yafu on aws-spot-N, then build
N=$1
R=/home/claude/OneShotFastECPP
S=/tmp/claude-1000/-home-claude-OneShotFastECPP/f1db9bbb-8f84-4630-ad02-28f9e40f5dc3/scratchpad
set -o pipefail
cd $R
echo "[$N] $(date +%T) repo"
git archive HEAD | spot ssh $N 'mkdir -p ~/rule/repo && tar x -C ~/rule/repo' || echo "[$N] repo FAILED"
echo "[$N] $(date +%T) phi_files"
tar c phi_files | spot ssh $N 'tar x -C ~/rule/repo' || echo "[$N] phi FAILED"
echo "[$N] $(date +%T) ecm + yafu from box 7"
spot ssh 7 'tar c -C ~ .local/ecmpkg yb/yafu/yafu' | spot ssh $N 'tar x -C ~' || echo "[$N] ecm/yafu FAILED"
echo "[$N] $(date +%T) build"
spot ssh $N 'cd ~/rule/repo && make -j > ~/rule/make.log 2>&1; echo "make exit $?"; ls -la ecpp/dscan ecpp/cm_method ecpp/smoothtest classpoly_v1.0.3/classpoly 2>&1 | cut -c1-120; du -sh phi_files; ls phi_files | wc -l; ~/.local/ecmpkg/usr/bin/ecm --version 2>&1 | head -1; ls -la ~/yb/yafu/yafu'
echo "[$N] $(date +%T) SETUP_DONE"
