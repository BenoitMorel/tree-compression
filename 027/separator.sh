#!/bin/bash
cat own_example.nex.run1.t | sed 's/tree gen.[0-9]* = \[&U\] //g' | sed 's/ //g'| cat > test.txt
COUNTER = 0
while read p; do
  COUNTER=$((COUNTER+1))
  cat > tree_$COUNTER.nwk << EOF
$p
EOF
done < test.txt
