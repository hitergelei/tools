#!/bin/bash

# brent_aims.sh ver1.2 2015-11-22

# This script is written by Yoongu Kang.

# 1. Calculate energy in la, la+h, la+2h.
# 2. From the given three points, find the lowest point of the parabola.


# Change geometry.in for your system.

energy ()
{
a11=`echo "$1" | bc -l`
a21=`echo "$1/2." | bc -l`
a22=`echo "sqrt(3.)*$1/2." | bc -l`
cat >./geometry.in <<!
lattice_vector        $a11          0.00000000      0.00000000
lattice_vector        $a21              $a22        0.00000000
lattice_vector      0.00000000      0.00000000     30.00000000

atom_frac       0.00000000      0.00000000      0.01200558 Si
  constrain_relaxation x
  constrain_relaxation y
atom_frac       0.33333333      0.33333333     -0.01200558 Si
  constrain_relaxation x
  constrain_relaxation y
atom_frac       0.00000000      0.00000000      0.06207127 H
  constrain_relaxation x
  constrain_relaxation y
atom_frac       0.33333333      0.33333333     -0.06207127 H
  constrain_relaxation x
  constrain_relaxation y



!

mpirun -np 4 /opt/FHI-aims/bin/aims.081213.scalapack.mpi.x > $a11-out
E=`grep "Total energy," $a11-out | tail -1 | awk '{print $10}'`
echo $1  $E >> ./E_latt.dat
}


# Calculate new lattice constant : x1,x2,x3 => la

new_lattice_constant ()
{
c1=`echo "(($x2)-($x1))*(($y2)-($y3))" | bc -l`
c2=`echo "(($x2)-($x3))*(($y2)-($y1))" | bc -l`
c3=`echo "(($x2)-($x1))*($c1)-(($x2)-($x3))*($c2)" | bc -l`
c4=`echo "($c1)-($c2)" | bc -l`
la=`echo "($x2)-0.5*($c3)/($c4)" | bc -l`
}


# main function

la=3.892  
h=0.001 
tol=2    


while [ $tol -ne 0 ]
do

echo la $la  h $h >> ./E_latt.dat

la=`echo "($la)-($h)" | bc -l`
energy $la
x1=`echo "($la)" | bc -l`
y1=`echo "($E)" | bc -l`

la=`echo "($la)+($h)" | bc -l`
energy $la
x2=`echo "($la)" | bc -l`
y2=`echo "($E)" | bc -l`

la=`echo "($la)+($h)" | bc -l`
energy $la
x3=`echo "($la)" | bc -l`
y3=`echo "($E)" | bc -l`

new_lattice_constant

if [ `echo "(($y2) < ($y3)) && (($y2) < ($y1)) " | bc` = "1" ]
then
        h=`echo "($h)/10." | bc -l`
	tol=$(($tol-1))
fi
done

exit 0

