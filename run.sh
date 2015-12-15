#!/bin/bash

#export GOMP_CPU_AFFINITY=0-3

filename=test.txt
nsamples=500000000
#nsamples=100

#rm -f $filename
echo "Subject: sample=$nsamples" > $filename
echo "" >> $filename

for t in {1..4}; do
  #printf "# threads = $t, " 2>&1 | tee -a $filename ; ( exit ${PIPESTATUS} )
  export OMP_NUM_THREADS=$t
  foo=$( { ./bench $nsamples; } 2>&1 )
  printf "$foo, " 2>&1 | tee -a $filename ; ( exit ${PIPESTATUS} )
  echo "" 2>&1 | tee -a $filename ; ( exit ${PIPESTATUS} )
  echo "============================================" 2>&1 | tee -a $filename ; ( exit ${PIPESTATUS} )
  echo "" 2>&1 | tee -a $filename ; ( exit ${PIPESTATUS} )
done;


sendmail felipe.mt87@gmail.com < $filename
