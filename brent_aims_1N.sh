#!/bin/bash

# brent_aims.sh ver1.2 2015-11-22

# This script is written by Yoongu Kang.

# 1. Calculate energy in la, la+h, la+2h.
# 2. From the given three points, find the lowest point of the parabola.


# Change geometry.in for your system.

energy ()
{
a1=`echo "$la" | bc -l`
latt_prist=4.04625  
latt_Hterm=4.07177     
tau1=`echo "$latt_Hterm*9.*sqrt(3.)/2." | bc -l`
tau2=`echo "$tau1+$latt_prist*$a1*1./(2*sqrt(3.))" | bc -l`
a2=`echo "$tau2+$latt_Hterm*1./sqrt(3.)" | bc -l`



cat >./geometry.in <<!
lattice_vector   $latt_Hterm    0.0000000000    0.0000000000
lattice_vector  0.0000000000        $a2         0.0000000000
lattice_vector  0.0000000000    0.0000000000    30.0000000000

atom       1.94619000        $tau1        -0.29164114 Si
  constrain_relaxation x
atom       0.00000000        $tau2         0.29161101 Si
  constrain_relaxation x
atom       0.00000000      0.00681335      0.68564069 Si
  constrain_relaxation x
atom       0.00000000      4.44058970     -0.36845125 Si
  constrain_relaxation x
atom       0.00000000      6.71051483      0.30027624 Si
  constrain_relaxation x
atom       0.00000000     11.22803351     -0.37699572 Si
  constrain_relaxation x
atom       0.00000000     13.47771380      0.35579810 Si
  constrain_relaxation x
atom       0.00000000     18.00404835     -0.32477166 Si
  constrain_relaxation x
atom       0.00000000     20.25007747      0.41917509 Si
  constrain_relaxation x
atom       0.00000000     24.75139728     -0.39574124 Si
  constrain_relaxation x
atom       0.00000000     27.05147852      0.16325829 Si
  constrain_relaxation x
atom       1.94619000      1.03943791     -0.16315952 Si
  constrain_relaxation x
atom       1.94619000      3.33955506      0.39575368 Si
  constrain_relaxation x
atom       1.94619000      7.84070868     -0.41918249 Si
  constrain_relaxation x
atom       1.94619000     10.08674556      0.32475129 Si
  constrain_relaxation x
atom       1.94619000     14.61314242     -0.35580233 Si
  constrain_relaxation x
atom       1.94619000     16.86282161      0.37701860 Si
  constrain_relaxation x
atom       1.94619000     21.38029691     -0.30024946 Si
  constrain_relaxation x
atom       1.94619000     23.65024268      0.36840361 Si
  constrain_relaxation x
atom       1.94619000     28.08392387     -0.68571244 Si
  constrain_relaxation x
atom       0.00000000      0.22659719      2.17260439 H
  constrain_relaxation x
atom       0.00000000      4.39757732     -1.86962222 H
  constrain_relaxation x
atom       0.00000000      6.73136781      1.80222743 H
  constrain_relaxation x
atom       0.00000000     11.24294721     -1.87887566 H
  constrain_relaxation x
atom       0.00000000     13.47120492      1.85778098 H
  constrain_relaxation x
atom       0.00000000     18.02796464     -1.82641442 H
  constrain_relaxation x
atom       0.00000000     20.23933500      1.92101046 H
  constrain_relaxation x
atom       0.00000000     24.66398437     -1.89547967 H
  constrain_relaxation x
atom       0.00000000     27.18747930      1.65685046 H
  constrain_relaxation x
atom       1.94619000      0.90330298     -1.65680132 H
  constrain_relaxation x
atom       1.94619000      3.42678586      1.89547628 H
  constrain_relaxation x
atom       1.94619000      7.85151840     -1.92101541 H
  constrain_relaxation x
atom       1.94619000     10.06288253      1.82639985 H
  constrain_relaxation x
atom       1.94619000     14.61961582     -1.85778362 H
  constrain_relaxation x
atom       1.94619000     16.84787559      1.87889158 H
  constrain_relaxation x
atom       1.94619000     21.35947394     -1.80220950 H
  constrain_relaxation x
atom       1.94619000     23.69325620      1.86958802 H
  constrain_relaxation x
atom       1.94619000     27.86432926     -2.17262716 H
  constrain_relaxation x




!


mpirun -np 4 /opt/FHI-aims/bin/aims.081213.scalapack.mpi.x > 0$a1-out
E=`grep "Total energy," 0$a1-out | tail -1 | awk '{print $10}'`
echo 0$a1  $a2  $tau2  $E >> ./E_latt.dat
}


# Calculate new lattice constant : x1,x2,x3 => la

new_lattice_constant ()
{
c1=`echo "(($x2)-($x1))*(($y2)-($y3))" | bc -l`
c2=`echo "(($x2)-($x3))*(($y2)-($y1))" | bc -l`
c3=`echo "(($x2)-($x1))*($c1)-(($x2)-($x3))*($c2)" | bc -l`
c4=`echo "($c1)-($c2)" | bc -l`
nx=`echo "($x2)-0.5*($c3)/($c4)" | bc -l`
la=`echo "(2.*sqrt(3.)/$latt_prist)*($nx-($latt_Hterm*29./6.*sqrt(3.)))" | bc -l`
}


# main function

la=0.827  
h=0.001 
tol=1    


while [ $tol -ne 0 ]
do

echo la- $la  a2- $nx  h- $h >> ./E_latt.dat

la=`echo "($la)-($h)" | bc -l`
energy $la
x1=`echo "($a2)" | bc -l`
y1=`echo "($E)" | bc -l`

la=`echo "($la)+($h)" | bc -l`
energy $la
x2=`echo "($a2)" | bc -l`
y2=`echo "($E)" | bc -l`

la=`echo "($la)+($h)" | bc -l`
energy $la
x3=`echo "($a2)" | bc -l`
y3=`echo "($E)" | bc -l`

new_lattice_constant

if [ `echo "(($y2) < ($y3)) && (($y2) < ($y1)) " | bc` = "1" ]
then
        h=`echo "($h)/10." | bc -l`
	tol=$(($tol-1))
fi
done

exit 0

