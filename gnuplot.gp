#set term post portrait  enhanced color "Helvetica,20"
#set output 'epot.eps'

################################
set title "Temperature" font "Curier Bold,18,"
set size ratio 0.6
#set xrange[4.50000:9.50000]
set yrange[500:1200000]
set logscale y
set lmargin at screen 0.13
set bmargin 1.0 
#set xtics 0,0.5,2
set xtics 0.5
#set ytics -80,10,-30
#set ytics 100
set ytics ("10^1" 10, "10^2" 100, "10^3" 1000, "10^4" 10000, "10^5" 100000, "10^6" 1000000)
set xtics font "Curier Bold, 16"
set ytics font "Curier Bold, 16"
set xlabel "steps(x10^3)" font "Curier Bold,16" offset 0,-0.7,0
set ylabel "T(K)" font "Curier Bold,16" offset -3,0,0
#set grid
#set tics nomirror
set arrow from 0,10000 to 2,10000 nohead lt rgb "gray"
set object 1 rectangle from 0,10000 to 2,1200000 behind fc rgb "gray" fs solid 0.4 noborder

plot "dyn_log.txt" using 1:5 w l lw 2 lt 1  lc rgb "purple" notitle
# "H18.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "blue"   title 'H18',\
# "C16.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "green"   title 'C16',\
# "gra.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "black"   title 'grapite',\
# "dia.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "#008888"   title  'diamond',\
# "gra2.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "#888800"   title  'grapite2',\
# "sp3.dat" using 1:2 w lp lw 3  lt 1  pt 7 ps .3  lc rgb "#880088"   title  'C16_sp3'
