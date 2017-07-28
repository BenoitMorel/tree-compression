#!/bin/bash

for i in `seq 1 200`;
do
   echo
   echo ----------------------- tree_$i \& tree_$(($i+1)) -----------------------
   ./rf-distance 027/tree_$i.nwk 027/tree_$(($i+1)).nwk
done   
