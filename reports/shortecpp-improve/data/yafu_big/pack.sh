#!/bin/bash
# pack.sh C [nopkl]: tar stream of this run's state (see pull_big.sh on the driver box)
C=$1
cd ~/rule/repo/work
F=$(ls -d shortcamp/c$C.* 2>/dev/null | grep -v '\.tmp$')
[ "$2" = nopkl ] && F=$(echo "$F" | grep -v '\.search\.pkl$')
cd ~/rule
if [ -n "$F" ]; then tar c big -C ~/rule/repo/work $F; else tar c big; fi
