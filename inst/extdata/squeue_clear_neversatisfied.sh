#!/bin/bash
u=$(basename ~)
tc=$(squeue -u $u | grep DependencyNeverSatisfied  | awk '{if (NR > 0){print $1}}')
for c in $tc; do scancel $c; done
