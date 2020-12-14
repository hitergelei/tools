#!/bin/sh

#export OMP_NUM_THREADS=20
#cvd=`get_avail_gpu.py`
export CUDA_VISIBLE_DEVICES=

lmp_mpi -in input-md.in > out 
