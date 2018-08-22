set term post portrait  enhanced color "Helvetica,20"
set output 'Energetics_C.eps'

################################
set title "Energetics of Carbons" font "Curier Bold,18,"
 set size  ratio 1;
 set xrange[4.50000:9.50000]
 set yrange[-9.2:-7.8]
#set lmargin 1.2;set bmargin 1.5
set xtics 5,1,9
set ytics -9.2,0.2,-7.8 
set ylabel "Energy (eV)" font "Curier Bold,30";
 set tics nomirror


plot "oC8.dat" using 1:2 w lp lw 3 lt 1  pt 7 ps .3  lc rgb "red"   title 'oC8',\
 "H18.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "blue"   title 'H18',\
 "C16.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "green"   title 'C16',\
 "gra.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "black"   title 'grapite',\
 "dia.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "#008888"   title  'diamond',\
 "gra2.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "#888800"   title  'grapite2',\
 "sp3.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "#880088"   title  'C16_sp3'
