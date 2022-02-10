#!/bin/bash

# queue request
#$ -N bn_band
#$ -pe mpich 280
#$ -q Twice
##$ -l h=!'pcs7-18'
##$ -l h=!'pcs7-19'

#$ -S /bin/bash
#$ -V

#$ -cwd

#$ -o mpi.log
#$ -e mpi.log
#$ -j y

echo "Got $NSLOTS slots."
cat $TMPDIR/machines

#######################################################
### INTEL MPI (w/ INTEL COMPILER)
#######################################################

 #export I_MPI_FABRICS=shm:dapl
 cd $SGE_O_WORKDIR
#cat out > outold
#cat OUTCAR > OUTCARold
#cat POSCAR > POSCARold
#cat CONTCAR > CONTCARold

#export EXEC=python
#export INPUTFILE=pho-phonopy.py
#export OUTPUTFILE=job-out
#time $EXEC $INPUTFILE > $OUTPUTFILE
export EXEC=vasp_std
export INPUTFILE=
export OUTPUTFILE=out
time mpirun -machinefile $TMPDIR/machines -np $NSLOTS $EXEC $INPUTFILE > $OUTPUTFILE

#cat CONTCAR > POSCAR

 exit 0

