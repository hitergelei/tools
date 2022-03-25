#!/bin/sh

#export OMP_NUM_THREADS=20
cvd=`get_avail_gpu.py`
export CUDA_VISIBLE_DEVICES=$cvd

lmp_mpi -in input.in > out 
#python ase-md.py > out

#rm -rf lmp-result.traj
lmp2traj.py &
sleep 10
