#!/bin/bash 

#$ -N ptse2
#$ -pe mpi_16 16
#$ -q cnpl

#$ -S /bin/bash
#$ -V

#$ -cwd

#$ -o vasp.log
#$ -e vasp.log
#$ -j y

echo "Got $NSLOTS slots."
cat $TMPDIR/machines

#######################################################
### openmpi (w/ Intel compiler)
#######################################################

cat out > outold
cat OUTCAR > OUTCARold
cat POSCAR > POSCARold
cat CONTCAR > CONTCARold

 if [ "$MPI" = mvapich2 ]; then
   MPIEXEC=mpirun_rsh
 elif [ "$MPI" = openmpi ]; then
   MPIEXEC=mpirun
 elif [ "$MPI" = impi ]; then
   MPIEXEC=mpiexec.hydra
 fi

 if [ "$1" = ncl ]; then
   EXEC=vasp_ncl
 elif [ "$1" = std ]; then
   EXEC=vasp_std
 else
   exit 0
 fi
 
 for i in `seq 1 1 20`
     do
 cat > INCAR <<!
 SYSTEM=ptse2
#
   PREC    =         a      ! determines the energy cutoff ENCUT, |L|M|N|A|H|
   ISTART  =         1      ! job   : 0-new 1-cont 2-samecut
   ICHARG  =         1      ! charge: 0-wave 1-file 2-atom 10-const
   ISPIN   =         1      ! | 1-non spin polarized | 2-spin polarized |
#  MAGMOM  =  0 0 0  0 0 0 90*0.0
   ENCUT   =        400.    ! cut-off energy for plane wave basis set in eV
   NELM    =        200     ! maximum number of electronic SC (selconsistency)
#  NELMDL  =        -6      ! number of non-selconsistency SC steps
   EDIFF   =       1E-07    ! specifies the global break condition for the electronic
   LREAL   =      .FALSE.   ! real space projection .FALSE. or Auto
   IALGO   =         38     ! select algorithm (8=CG for small, 48=RMM for big systems)
   NSW     =         4      ! maximum number of ionic steps
   IBRION  =         1      ! how to ions are updated and moved
   EDIFFG  =       -.5E-03  ! break conditions for the ionic relaxation loop
   ISIF    =         7      ! controls whether the stress tensor is alculated
#  ISYM    =         0      ! switch symmetry on (1,2,3) or off (-1,0)
#  ADDGRID =      .TRUE.
 DOS related values:
   ISMEAR  =         0      ! for semiconductor or insulators ISMEAR=-5, SIGMA=0.05
   SIGMA   =         0.02   ! for metal ISMEAR=1 or =2 SIGMA=0.2
  EMIN   = -10
  EMAX   =  10
  NEDOS  = 20000
  LORBIT = 11
  LORBMOM = .TRUE.
 Write flags
#  LWAVE   =        .FALSE. ! These tags determine whether the orbitals (file WAVECAR),
#  LCHARG  =        .FALSE. ! the charge densities (files CHGCAR and CHG) are written
   LVTOT   =        .FALSE. !
 Non-colinear calculations and spin orbit coupling
   LSORBIT =        .TRUE.  ! switch on SOC and automatically set LNONCOLLINEAR=.TRUE.
   SAXIS   =         0 0 1  ! quantisation axis for spin
   GGA_COMPAT =     .FALSE. ! apply spherical cutoff on gradient field
!



 $MPIEXEC -np $NSLOTS $EXEC  > out

cat out > outold
cat OUTCAR > OUTCARold
cat POSCAR > POSCARold
cat CONTCAR > CONTCARold
cat CONTCAR > POSCAR
cat > INCAR <<!
 SYSTEM=ptse2
#
   PREC    =         a      ! determines the energy cutoff ENCUT, |L|M|N|A|H|
   ISTART  =         1      ! job   : 0-new 1-cont 2-samecut
   ICHARG  =         1      ! charge: 0-wave 1-file 2-atom 10-const
   ISPIN   =         1      ! | 1-non spin polarized | 2-spin polarized |
#  MAGMOM  =  0 0 0  0 0 0 90*0.0
   ENCUT   =        400.    ! cut-off energy for plane wave basis set in eV
   NELM    =        200     ! maximum number of electronic SC (selconsistency)
#  NELMDL  =        -6      ! number of non-selconsistency SC steps
   EDIFF   =       1E-07    ! specifies the global break condition for the electronic
   LREAL   =      .FALSE.   ! real space projection .FALSE. or Auto
   IALGO   =         38     ! select algorithm (8=CG for small, 48=RMM for big systems)
   NSW     =         4      ! maximum number of ionic steps
   IBRION  =         1      ! how to ions are updated and moved
   EDIFFG  =       -.5E-03  ! break conditions for the ionic relaxation loop
   ISIF    =         2      ! controls whether the stress tensor is alculated
#  ISYM    =         0      ! switch symmetry on (1,2,3) or off (-1,0)
#  ADDGRID =      .TRUE.
 DOS related values:
   ISMEAR  =         0      ! for semiconductor or insulators ISMEAR=-5, SIGMA=0.05
   SIGMA   =         0.02   ! for metal ISMEAR=1 or =2 SIGMA=0.2
  EMIN   = -10
  EMAX   =  10
  NEDOS  = 20000
  LORBIT = 11
  LORBMOM = .TRUE.
 Write flags
#  LWAVE   =        .FALSE. ! These tags determine whether the orbitals (file WAVECAR),
#  LCHARG  =        .FALSE. ! the charge densities (files CHGCAR and CHG) are written
   LVTOT   =        .FALSE. !
 Non-colinear calculations and spin orbit coupling
   LSORBIT =        .TRUE.  ! switch on SOC and automatically set LNONCOLLINEAR=.TRUE.
   SAXIS   =         0 0 1  ! quantisation axis for spin
   GGA_COMPAT =     .FALSE. ! apply spherical cutoff on gradient field
!

 $MPIEXEC -np $NSLOTS $EXEC  > out
 cat CONTCAR > POSCAR
 lat=`cat CONTCAR | head -3 | tail -1 | awk {'print $1'}`
 date=`date | awk '{print $1 " " $2 " "  $3 " " $4}'`
 ene=`grep 'E0= ' out | tail -1 | awk {'printf $5'}`
 echo $lat $ene $date  >> lattice.out

 done


 exit 0
