#!/usr/bin/sh

export e_fermi=`grep fermi ../OUTCAR | tail -n 1 | awk '{print $3}'`
echo "e_fermi=$e_fermi"
eigen2bp -10 20 GKMG $e_fermi
gnuplot -p Eigenbp.in
evince band.eps
