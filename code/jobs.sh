#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -q serial.q
#$ -N _c0.1_p1
#$ -pe parSMP 16
export OMP_NUM_THREADS=16
module load gcc
module load gsl
module load boost
module load mpi
time ./flat_simulation -X 14.28 -Y 14.28 -Z 3 -p 1 -n -1 -c 0.1 -d 0.714 -g 0.9 -e 80.1 -l 20 -r 20 -s 5000 -S 40000000 -P 200000 -F 10 -y 1000000 -x 1000000 -w 1000000 -m 500000 -k 0.002 -f 100 -M 1 -T 0.001 -q 1
