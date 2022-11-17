#!/bin/sh

# GPU Number
cvd=`get_avail_gpu.py`
export CUDA_VISIBLE_DEVICES=$cvd
export KMP_AFFINITY=none
export KMP_BLOCKTIME=0
#export CUDA_VISIBLE_DEVICES=""

#export NSLOTS=4
export NSLOTS=$OMP_NUM_THREADS
#export OUTPUT=out

#export EXEC=vasp_std
#export INPUT=
#time mpiexec.hydra -np $NSLOTS $EXEC $INPUT > $OUTPUT

export EXEC=python
export INPUT=pho-phonopy.py
time $EXEC $INPUT

#export EXEC=lmp_mpi
#export INPUT=input-rlx.in
#time $EXEC -in $INPUT > $OUTPUT
#lmp2traj.py


# time sleep 100

