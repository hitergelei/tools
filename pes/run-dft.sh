#!/bin/sh

#export OMP_NUM_THREADS=20
#cvd=`get_avail_gpu.py`
export CUDA_VISIBLE_DEVICES=

mpiexec.hydra -np $OMP_NUM_THREADS vasp_std > out
