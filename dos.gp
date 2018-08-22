set term post landscape enhanced color "Helvetica,18"
set output 'tri-Pb-X-1.eps'
set yrange[-4:2]
set xrange [0:2.4]
set size  ratio 5;
#set xtics ("-2" -2,"-1" -1,"E_f" 0,"1" 1,"2" 2)
#set xtics ("-10" -10, "-8" -8,"-6" -6, "-4" -4, "-2" -2,"E_f" 0,"2" 2)
#set ytics ("-2" -6, "-4" -4, "-2" -2,"E_f" 0,"2" 2)
set ytics -4,1,2 font ",18"
set xtics 0,1.2,2.4 font ",18"
#set ytics  0,2,10 font ",24"
set ylabel "Energy (eV)" font "Helvetica,20"
set xlabel "PDOS (a.u)" font "Helvetica,20"
set arrow from 0,0 to 2.4,0 nohead lt 3 lw 0.5 lc rgb "black"
#set arrow from 0,-2 to 0,2 nohead lt 3 lw 0.5 lc rgb "black"

set key  at  6.0, 1.5
set key spacing 1

fermi=-4.2574
plot "DOSCAR-Pb-X-1" u ($2):($1-fermi) w l lt 1 lw 2 lc rgb "black" title 's',\
     "DOSCAR-Pb-X-1" u ($6):($1-fermi) w l lt 1 lw 2 lc rgb "blue" title 'p_y',\
     "DOSCAR-Pb-X-1" u ($10):($1-fermi) w l lt 1 lw 2 lc rgb "red" title 'p_z',\
     "DOSCAR-Pb-X-1" u ($14):($1-fermi) w l lt 1 lw 2 lc rgb "yellow" title 'p_x',\
     "DOSCAR-Pb-X-1" u ($18):($1-fermi) w l lt 1 lw 2 lc rgb "green" title 'd_x_y',\
     "DOSCAR-Pb-X-1" u ($22):($1-fermi) w l lt 1 lw 2 lc rgb "#008888" title 'd_y_z',\
     "DOSCAR-Pb-X-1" u ($26):($1-fermi) w l lt 1 lw 2 lc rgb "#880088" title 'd_{z^2-r^2}',\
     "DOSCAR-Pb-X-1" u ($30):($1-fermi) w l lt 1 lw 2 lc rgb "#888800" title 'd_x_z',\
     "DOSCAR-Pb-X-1" u ($34):($1-fermi) w l lt 1 lw 2 lc rgb "#888888" title 'd_{x^2-y^2}'
#    "KS_DOS_total.dat" u (-1*$3):($1-fermi) w l lt 1 lw 2 lc rgb "black" notitle ,\
#    "atom_proj_dos_spin_upSi0001.dat" u   2:($1-fermi) w l lt 3 lw 1 lc rgb "blue" title 'Si_1',\
#    "atom_proj_dos_spin_dnSi0001.dat" u (-1*$2):($1-fermi) w l lt 3 lw 1 lc rgb "blue" notitle,\
#    "atom_proj_dos_spin_upSi0002.dat" u   2:($1-fermi) w l lt 3 lw 1 lc rgb "red" title 'Si_2',\
#    "atom_proj_dos_spin_dnSi0002.dat" u (-1*$2):($1-fermi) w l lt 3 lw 1 lc rgb "red" notitle

