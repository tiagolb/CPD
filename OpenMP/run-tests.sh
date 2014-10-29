#!/bin/bash

for num in 8 4 2 1
do
   echo "$num Threads:"
   export OMP_NUM_THREADS=$num
   for i in $(ls ../public-instances/*.ini)
   do
      file_var=${i:20}
      echo "./lcs-omp $i > ../result/Threads$num/test-$file_var.out"
      ./lcs-omp $i > ../result/Threads$num/test-$file_var.out
      echo "DONE Exec"
   done
   echo "DONE THREAD"
done
echo "DONE TEST"
