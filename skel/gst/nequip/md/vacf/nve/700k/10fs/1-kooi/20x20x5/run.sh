#!/bin/sh

#export OMP_NUM_THREADS=20
cvd=`get_avail_gpu.py`
export CUDA_VISIBLE_DEVICES=$cvd
export KMP_AFFINITY=none
#export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:100

lmp_mpi -in input.in > out 
#python ase-md.py > out

#rm -rf lmp-result.traj
lmp2traj.py &
sleep 10
