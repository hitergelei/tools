#set term post portrait  enhanced color "Helvetica,20"
#set output 'epot.eps'

################################
set title "Phonon Dispersion" font "Curier Bold,18,"
set size ratio 0.6
#set xrange[4.50000:9.50000]
set yrange[-6:8]
#set logscale y
set lmargin at screen 0.13
set bmargin 1.0 
#set xtics 0,0.5,2
#set xtics 0.5
set xtics ("GM" 0.00000000, "H" 0.31283850, "P" 0.58376460, "GM" 0.85469060, "N" 1.07590090, "H" 1.29711110)
#set ytics -80,10,-30
#set ytics 100
#set ytics ("10^1" 10, "10^2" 100, "10^3" 1000, "10^4" 10000, "10^5" 100000, "10^6" 1000000)
set ytics 0,1,0
set xtics font "Curier Bold, 16"
set ytics font "Curier Bold, 16"
#set xlabel "steps(x10^3)" font "Curier Bold,16" offset 0,-0.7,0
set ylabel "frequency" font "Curier Bold,16" offset -3,0,0
set grid
#set tics nomirror
#set arrow from 1,-6 to 1,8 nohead lt rgb "gray"
#set object 1 rectangle from 0,10000 to 2,1200000 behind fc rgb "gray" fs solid 0.4 noborder

plot "band.in" using 1:2 w l lw 1.5 lt 1  lc rgb "purple" title 'beta-Ti'
# "H18.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "blue"   title 'H18',\
# "C16.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "green"   title 'C16',\
# "gra.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "black"   title 'grapite',\
# "dia.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "#008888"   title  'diamond',\
# "gra2.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "#888800"   title  'grapite2',\
# "sp3.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "#880088"   title  'C16_sp3'
