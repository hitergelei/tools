#!/bin/bash

for la in `seq 0.980 0.002 0.980`        
do
latt_prist=4.68001  
latt_Hterm=4.72611     
tau1=`echo "$latt_Hterm*9.*sqrt(3.)/2." | bc -l`
tau2=`echo "$tau1+$latt_prist*$la*1./(2*sqrt(3.))" | bc -l`
a2=`echo "$tau2+$latt_Hterm*1./sqrt(3.)" | bc -l`

cat >./geometry.in <<!

  lattice_vector         4.72611000        0.00000000        0.00000000
  lattice_vector         0.00000000           $a2            0.00000000
  lattice_vector         0.00000000        0.00000000       30.00000000

            atom         2.36305500          $tau1           0.43126216  Sn
  constrain_relaxation x
            atom         0.00000000          $tau2          -0.43153911  Sn
  constrain_relaxation x
            atom         0.00000000       -0.00894785        0.41807920  Sn
  constrain_relaxation x
            atom         0.00000000        5.45957444       -0.40784635  Sn
  constrain_relaxation x
            atom         0.00000000        8.18831417        0.40700616  Sn
  constrain_relaxation x
            atom         0.00000000       13.64613182       -0.40879126  Sn
  constrain_relaxation x
            atom         0.00000000       16.37256345        0.40858429  Sn
  constrain_relaxation x
            atom         0.00000000       21.82768854       -0.40853699  Sn
  constrain_relaxation x
            atom         0.00000000       24.55359407        0.40896175  Sn
  constrain_relaxation x
            atom         0.00000000       30.01565926       -0.40702866  Sn
  constrain_relaxation x
            atom         0.00000000       32.74835554        0.39939966  Sn
  constrain_relaxation x
            atom         2.36305500        1.35944074       -0.39954290  Sn
  constrain_relaxation x
            atom         2.36305500        4.09206962        0.40712230  Sn
  constrain_relaxation x
            atom         2.36305500        9.55415951       -0.40887705  Sn
  constrain_relaxation x
            atom         2.36305500       12.28007454        0.40859775  Sn
  constrain_relaxation x
            atom         2.36305500       17.73520475       -0.40858451  Sn
  constrain_relaxation x
            atom         2.36305500       20.46162208        0.40883527  Sn
  constrain_relaxation x
            atom         2.36305500       25.91945390       -0.40689638  Sn
  constrain_relaxation x
            atom         2.36305500       28.64818455        0.40799146  Sn
  constrain_relaxation x
            atom         2.36305500       34.11669512       -0.41830543  Sn
  constrain_relaxation x
            atom         0.00000000        0.03630028        2.16091329  H
  constrain_relaxation x
            atom         0.00000000        5.46588506       -2.14813763  H
  constrain_relaxation x
            atom         0.00000000        8.19422700        2.14712236  H
  constrain_relaxation x
            atom         0.00000000       13.64698851       -2.14865703  H
  constrain_relaxation x
            atom         0.00000000       16.37275326        2.14859595  H
  constrain_relaxation x
            atom         0.00000000       21.82355317       -2.14840959  H
  constrain_relaxation x
            atom         0.00000000       24.54912791        2.14900966  H
  constrain_relaxation x
            atom         0.00000000       30.00711635       -2.14721451  H
  constrain_relaxation x
            atom         0.00000000       32.74763907        2.14079007  H
  constrain_relaxation x
            atom         2.36305500        1.36031906       -2.14093295  H
  constrain_relaxation x
            atom         2.36305500        4.10048147        2.14730796  H
  constrain_relaxation x
            atom         2.36305500        9.55860600       -2.14892520  H
  constrain_relaxation x
            atom         2.36305500       12.28421856        2.14847028  H
  constrain_relaxation x
            atom         2.36305500       17.73503264       -2.14859600  H
  constrain_relaxation x
            atom         2.36305500       20.46073036        2.14870102  H
  constrain_relaxation x
            atom         2.36305500       25.91358144       -2.14701284  H
  constrain_relaxation x
            atom         2.36305500       28.64187420        2.14828225  H
  constrain_relaxation x
            atom         2.36305500       34.07140260       -2.16113763  H
  constrain_relaxation x

!

mpirun -np 4 /opt/FHI-aims/bin/aims.081213.scalapack.mpi.x > $la-out
E=`grep "Total energy," $la-out | tail -1 | awk '{print $10}'`
echo $la   $a2   $tau2   $E >> ./E_latt.dat
done

exit 0
