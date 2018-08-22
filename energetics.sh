

#enter npoint value to final seq below!!!!!!!!!!!!!!!!!!
for i in `seq 1.0  1.0 20.0`
do

i=$i
#enter vol/atom of system you want to start
ivol=5
#enter vol/atom of system you want to end up
fvol=9
#enter how many points you want between start point & end point
npoint=20
#enter your systems lattice vectors
ia1=7.8202837469019588
ia2=0
ia3=0
ib1=0
ib2=4.8838707152591994
ib3=0
ic1=0
ic2=0
ic3=3.2407019257781196
#enter your systems total atom #
natom=16


vol=`echo "(($ic1*(($ia2*$ib3)-($ia3*$ib2)))+($ic2*(($ia3*$ib1)-($ia1*$ib3)))+($ic3*(($ia1*$ib2)-($ia2*$ib1))))/$natom" | bc -l`
interval=`echo "($fvol-$ivol)/($npoint-1)" | bc -l`
newvol=`echo "($ivol+($interval*($i-1)))" | bc -l`
ff=`echo "$newvol/$vol" | bc -l`
echo 'print('"$ff"'**0.333333333)' > ff.py | bc -l
f=`python ff.py`

a1=`echo "$f*$ia1" | bc -l`
a2=`echo "$f*$ia2" | bc -l`
a3=`echo "$f*$ia3" | bc -l`
b1=`echo "$f*$ib1" | bc -l`
b2=`echo "$f*$ib2" | bc -l`
b3=`echo "$f*$ib3" | bc -l`
c1=`echo "$f*$ic1" | bc -l`
c2=`echo "$f*$ic2" | bc -l`
c3=`echo "$f*$ic3" | bc -l`

#adjust POSCAR input below except lattice vectors!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


cat > POSCAR <<!
C16_vol_$newvol
   1.00000000000000
            $a1                   $a2                   $a3
            $b1                   $b2                   $b3
            $c1                   $c2                   $c3
   C
    16
Selective dynamics
Direct
  0.0884570498509460  0.0000000000000000  0.5000000000000000   T   T   T
  0.0884570498509460  0.5000000000000000  0.5000000000000000   T   T   T
  0.1768824158169276  0.2500000000000000  0.3736290823649631   T   T   T
  0.1768824158169276  0.7500000000000000  0.6263709176350369   T   T   T
  0.3231175841830725  0.2500000000000000  0.1263709176350372   T   T   T
  0.3231175841830725  0.7500000000000000  0.8736290823649631   T   T   T
  0.4115429501490541  0.0000000000000000  0.0000000000000000   T   T   T
  0.4115429501490541  0.5000000000000000  0.0000000000000000   T   T   T
  0.5884570498509459  0.0000000000000000  0.0000000000000000   T   T   T
  0.5884570498509459  0.5000000000000000  0.0000000000000000   T   T   T
  0.6768824158169340  0.2500000000000000  0.1263709176350372   T   T   T
  0.6768824158169340  0.7500000000000000  0.8736290823649631   T   T   T
  0.8231175841830660  0.2500000000000000  0.3736290823649631   T   T   T
  0.8231175841830660  0.7500000000000000  0.6263709176350369   T   T   T
  0.9115429501490541  0.0000000000000000  0.5000000000000000   T   T   T
  0.9115429501490541  0.5000000000000000  0.5000000000000000   T   T   T
!



mpirun -np 2  vasp.5.3.5_nonSOC > out 

cat POSCAR > POSCAR_$newvol
cat CONTCAR > CONTCAR_$newvol
cat OUTCAR > OUTCAR_$newvol
cat out > out_$newvol
grep "d E =" out >> situation.dat

TOT=`grep TOTEN OUTCAR | tail -1 | awk '{print $5}'`
date=`date | awk '{print $1 " " $2 " "  $3 " " $4}'`

echo $newvol $TOT $date >> Results.OUT


done



exit 0
