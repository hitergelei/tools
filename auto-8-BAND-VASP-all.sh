#!/bin/sh
#       for VASP energy band ploting tool
#       written by Hyun-Jung Kim,  2008.06.10 summer ; basic program is written
#       modified by HJ Kim.. 2012.May.29 ;small bug fix
#       modified by HJ Kim.. 2013.Jun.19 ;program structure is entirely changed
#       modified by HJ Kim.. 2013.Jul.24 ;automated more
#       modified by HJ Kim.. 2013.Jul.25 ;addup-band-decomposition fortran..
#       modified by HJ Kim.. 2013.Jul.30 ;spin-polarized case implemented ISPIN=2

# THE Things that you should modify..
# SYS, L_Lines(important!!), Y-axis range, The used symmetry points(important)!, PROCAR related tags
# GNUPLOT tags;x-label, arrows ..etc.

SYS="test"
N_Lines=7           # Number of symmetry lines !!!! Check!!! IMPORTANT !!!!
E_Fermi=8.9516      # grep fermi OUTCAR (total energy calculation)
ISPIN=`grep ISPIN OUTCAR | awk '{print $3}'`  # Spin-polarized=2 , No-spin=1 
Yinit=-5;Yend=5   # Y-axis range

### How many symmetry points? : number of lines = LINE#   !!!! Check !!! IMPORTANT !!!!
                       KP_1="  0.0000000000     0.0000000000     0.0000000000 ";KP_1_name="T";KPS_init='T'
LINE1=" output band " ;KP_2=" -0.5000000000     0.0000000000     0.0000000000 ";KP_2_name="W";KPS_2='W'
LINE2=" output band " ;KP_3=" -0.2777777778    -0.4411764706     0.0000000000 ";KP_3_name="R";KPS_3='R'
LINE3=" output band " ;KP_4="  0.0000000000     0.0000000000     0.0000000000 ";KP_4_name="G";KPS_4='{/Symbol \G}'
LINE4=" output band " ;KP_5="  0.2777777778    -0.5588235294     0.0000000000 ";KP_5_name="X";KPS_5='X'
LINE5=" output band " ;KP_6="  0.2777777778    -0.5588235294     0.5000000000 ";KP_6_name="S";KPS_6='S'
LINE6=" output band " ;KP_7="  0.0000000000     0.0000000000     0.5000000000 ";KP_7_name="G";KPS_7='{/Symbol \G}'
LINE7=" output band " ;KP_8="  0.0000000000     0.0000000000     0.0000000000 ";KP_8_name="G";KPS_end='{/Symbol \G}'

LTU="lp  lt 1 lw 2 pt 7 ps .0 lc rgb 'red' "   # Line type
LTD="lp  lt 2 lw 2 pt 6 ps .0 lc rgb 'blue' " 
atom1=3    # atom to be resolved ;Ir1 (-2)
atom2=4    # atom to be resolved ;Ir2 (-2)
PTa1="p pt 6 ps variable  notitle"
PTa2="p pt 8 ps variable  notitle"
spd1=12    # orbital to be resolved ;d
#PT1="p pt 6 ps variable lc rgb 'red' notitle"
#PT1="p pt 6 ps variable lc rgb 'red' notitle"
DOSNAME_atom="DOS_atom_projected.dat"
DOSNAME_spd="DOS_spd_projected.dat"
EIGENVAL="EIGENVAL"
PROCAR="PROCAR"         # PROCAR result from LORBIT=11

rm vasp_band.f vasp_band  
cat >> vasp_band.f <<EOF
        character*80 string(5),filenm,foname
        dimension band(6000,6000,2),reci(6000,1000)
        write(6,*) "Enter file name to be read"
        read(5,*)filenm

        open(unit=8,file=filenm,status='old')
        write(6,*) ' Enter Fermi Energy $ total # of symmetry lines'
        read(5,*)fermi,n_line
        read(8,*) natom,natom,dummy,ispin
        do j=1,4
        read(8,*) string(j)
        enddo
        read(8,*) nelect,nkp,neig

        do n=1, nkp
        read(8,*) (reci(n,j),j=1,4)
        do i=1, neig
        read(8,*) n_eig, (band(n,i,isp),isp=1,ispin)
        enddo
        enddo

        write(6,*)'Number of K-points=',nkp,'& eigenvalues=',neig

        do isp=1,ispin
        nk=0
        do line_n=1,n_line  
        write(foname,'(A,I1,A,I1,A)')'band',isp,'00',line_n,'.out'
        open(unit=9,file=foname,status='unknown')

        do k=1,int(nkp/n_line)   ! write k-point
        nk=nk+1
        write(9,'(I4,2X,3F15.7)',ADVANCE='NO')k,(reci(nk,j),j=1,3)

         do i=1,neig   ! write eigenvalue
          if( band(nk,i,isp) .le. fermi )then
         write(9,'(F12.5,F15.5)',ADVANCE='NO')2.,band(nk,i,isp)-fermi
          else
         write(9,'(F12.5,F15.5)',ADVANCE='NO')0.,band(nk,i,isp)-fermi
          endif
         enddo

         write(9,'()')
        enddo
         close(9)
        enddo   
        enddo
        close(8);
        end
EOF

ifort -o vasp_band vasp_band.f

#for LINE in 1 2 3
#do
cat >>input_band.d<<EOF
$EIGENVAL
$E_Fermi  $N_Lines
EOF
./vasp_band < input_band.d
#mv band.out band100$LINE.out
rm input_band.d
#done

if [ $ISPIN -eq 2 ]; then
cat >> cont.in << EOF
spin   collinear
EOF
fi

NDIVK=`head -2 KPOINTS | tail -1 | awk '{print $1}'`;# let z=($N_lines*3-1)
cat >>cont.in<<EOF
$LINE1 $KP_1  $KP_2  $NDIVK  $KP_1_name $KP_2_name 
$LINE2 $KP_2  $KP_3  $NDIVK  $KP_2_name $KP_3_name
$LINE3 $KP_3  $KP_4  $NDIVK  $KP_3_name $KP_4_name
$LINE4 $KP_4  $KP_5  $NDIVK  $KP_4_name $KP_5_name
$LINE5 $KP_5  $KP_6  $NDIVK  $KP_5_name $KP_6_name
$LINE6 $KP_6  $KP_7  $NDIVK  $KP_6_name $KP_7_name
$LINE7 $KP_7  $KP_8  $NDIVK  $KP_7_name $KP_8_name
$LINE8 $KP_8  $KP_9  $NDIVK  $KP_8_name $KP_9_name
EOF

lattice_vector11=`head -3 POSCAR | tail -1 | awk '{print $1}'`
lattice_vector12=`head -3 POSCAR | tail -1 | awk '{print $2}'`
lattice_vector13=`head -3 POSCAR | tail -1 | awk '{print $3}'`
lattice_vector21=`head -4 POSCAR | tail -1 | awk '{print $1}'`
lattice_vector22=`head -4 POSCAR | tail -1 | awk '{print $2}'`
lattice_vector23=`head -4 POSCAR | tail -1 | awk '{print $3}'`
lattice_vector31=`head -5 POSCAR | tail -1 | awk '{print $1}'`
lattice_vector32=`head -5 POSCAR | tail -1 | awk '{print $2}'`
lattice_vector33=`head -5 POSCAR | tail -1 | awk '{print $3}'`
cat >>geom.in<<EOF
  lattice_vector $lattice_vector11 $lattice_vector12 $lattice_vector13
  lattice_vector $lattice_vector21 $lattice_vector22 $lattice_vector23
  lattice_vector $lattice_vector31 $lattice_vector32 $lattice_vector33
EOF

vasp_band_plotting.pl
rm cont.in geom.in

############GNUPLOT ###############
if [ $ISPIN -eq 1 ]; then
FNAMEU="band_structure.dat"
#LTU="lp  lt 1  pt 7 ps .3"
LTU=$LTU
fi
if [ $ISPIN -eq 2 ]; then
FNAMEU="band_structure_spin_up.dat"
FNAMED="band_structure_spin_down.dat"
LTU=$LTU
LTD=$LTD
fi

K_init_name=`grep "Starting point for band" $FNAMEU | head -1 | awk '{print $7}'`
K_end_name=`grep "Ending   point for band" $FNAMEU | tail -1 | awk '{print $7}'`
K_init_point=`grep "Starting point for band" $FNAMEU | head -1 | tail -1| awk '{print $19}'`
K_end_point=`grep "Ending   point for band" $FNAMEU | tail -1 | head -1| awk '{print $19}'`
K_2_name=`grep "Starting point for band" $FNAMEU  | head -2 | tail -1| awk '{print $7}'`
K_2_point=`grep "Starting point for band" $FNAMEU | head -2 | tail -1| awk '{print $19}'`
K_3_name=`grep "Starting point for band" $FNAMEU  | head -3 | tail -1| awk '{print $7}'`
K_3_point=`grep "Starting point for band" $FNAMEU | head -3 | tail -1| awk '{print $19}'`
K_4_name=`grep "Starting point for band" $FNAMEU  | head -4 | tail -1| awk '{print $7}'`
K_4_point=`grep "Starting point for band" $FNAMEU | head -4 | tail -1| awk '{print $19}'`
K_5_name=`grep "Starting point for band" $FNAMEU  | head -5 | tail -1| awk '{print $7}'`
K_5_point=`grep "Starting point for band" $FNAMEU | head -5 | tail -1| awk '{print $19}'`
K_6_name=`grep "Starting point for band" $FNAMEU  | head -6 | tail -1| awk '{print $7}'`
K_6_point=`grep "Starting point for band" $FNAMEU | head -6 | tail -1| awk '{print $19}'`
K_7_name=`grep "Starting point for band" $FNAMEU  | head -7 | tail -1| awk '{print $7}'`
K_7_point=`grep "Starting point for band" $FNAMEU | head -7 | tail -1| awk '{print $19}'`
K_8_name=`grep "Starting point for band" $FNAMEU  | head -8 | tail -1| awk '{print $7}'`
K_8_point=`grep "Starting point for band" $FNAMEU | head -8 | tail -1| awk '{print $19}'`

Xinit=$K_init_point;Xend=$K_end_point
#Yinit=-10.5;Yend=.5
rm gnuBAND.gp gnuBAND-decomp.gp

cat >> gnuBAND.gp <<EOF
set term post portrait  enhanced color "Helvetica,20"
set output 'BAND-$SYS.eps'

#  $SYS  ###############################
set title "Energy band $SYS" font "Curier Bold,18,"
 set size  ratio 1;
 set xrange[$Xinit:$Xend+0.01]
 set yrange[$Yinit:$Yend]
#set lmargin 1.2;set bmargin 1.5
#set xtics 1,2,31 
#set xtics ("{/Symbol \G}" $Xinit, "X" 0.69772 ,"M" 1.32016,"{/Symbol \G}" $Xend)
set xtics ("$KPS_init" $Xinit, "$KPS_2" $K_2_point, "$KPS_3" $K_3_point, "$KPS_4" $K_4_point, "$KPS_5" $K_5_point, "$KPS_6" $K_6_point, "$KPS_7" $K_7_point, "$KPS_8" $K_8_point, "$KPS_end" $Xend) font "Helvetical,24" offset 0,0.2
set ylabel "Energy (eV)" font "Curier Bold,30";
 x1=$Xinit; x3=$Xend; y3=0.0;unset arrow
 set tics nomirror
 set arrow from x1,y3 to x3,y3 nohead lt 3 lw 1 
 set arrow from $K_2_point,$Yinit to $K_2_point,$Yend nohead lt 3 lw 1 
 set arrow from $K_3_point,$Yinit to $K_3_point,$Yend nohead lt 3 lw 1 
 set arrow from $K_4_point,$Yinit to $K_4_point,$Yend nohead lt 3 lw 1 
 set arrow from $K_5_point,$Yinit to $K_5_point,$Yend nohead lt 3 lw 1 
 set arrow from $K_6_point,$Yinit to $K_6_point,$Yend nohead lt 3 lw 1 
 set arrow from $K_7_point,$Yinit to $K_7_point,$Yend nohead lt 3 lw 1 
 set arrow from $K_8_point,$Yinit to $K_8_point,$Yend nohead lt 3 lw 1 

plot "$FNAMEU" using 1:2 w $LTU notitle,\\
EOF

NEIG=`grep NBANDS OUTCAR | awk '{print $15}'`
for k in `seq 3 $NEIG`
do
cat >> gnuBAND.gp <<EOF
 "$FNAMEU" using 1:$k w $LTU notitle,\\
EOF
done
if [ $ISPIN -eq 1 ]; then
cat >> gnuBAND.gp <<EOF
 "$FNAMEU" using 1:$k w $LTU notitle
EOF
fi
if [ $ISPIN -eq 2 ]; then
cat >> gnuBAND.gp <<EOF
 "$FNAMEU" using 1:$k w $LTU notitle,\\
EOF
for k in  `seq 2 $NEIG`
do
cat >> gnuBAND.gp <<EOF
"$FNAMED" using 1:$k w $LTD notitle,\\
EOF
done
cat >> gnuBAND.gp <<EOF
"$FNAMED" using 1:$k w $LTD notitle
EOF
fi


gnuplot gnuBAND.gp
eps2eps BAND-$SYS.eps a.eps
mv a.eps BAND-$SYS.eps
rm band_decomposed.f 
cat >> band_decomposed.f <<EOF
        character*80 string(5),filenm,str
        dimension eig(300,400),band_dist(300)
        dimension dos_total(300,400,150)
        dimension dos_atom_tot(300,400,150)
        dimension dos_atom_s(300,400,150),dos_atom_p(300,400,150,3)
        dimension dos_atom_d(300,400,150,5)
        dimension dos_total_s(300,400),dos_total_p(300,400)
        dimension dos_total_d(300,400)

        open(unit=11,file="$EIGENVAL",status='old')
        read(11,*) natom, natom, niter, ispin
        do j=1,4
        read(11,*) string(j)
        enddo
        read(11,*) nelect,nkp,neig
        close(11)
        open(unit=12,file="$FNAMEU",status='old')
        do j=1,6
        read(12,*) string(1)
        enddo
        do j=1,$N_Lines*3 !N_Lines*3
        read(12,*) string(1)
        enddo
        do k=1, nkp
        read(12,*) band_dist(k),(eig(k,n),n=1,neig)
        enddo
        close(12)

        open(unit=13,file="$PROCAR",status='old')
        read(13,*) string(1);read(13,*) string(2)
        do k=1,nkp
        read(13,*) string(1)
         do n=1,neig
         read(13,*) string(1)
         read(13,*) string(1)
          do i=1,natom
           read(13,*)nn,dos_atom_s(k,n,i),(dos_atom_p(k,n,i,l),l=1,3),
     &(dos_atom_d(k,n,i,l),l=1,5),dos_atom_tot(k,n,i)
          enddo
           read(13,'(A3,10F7.3)')str,(dos_total(k,n,l),l=2,11)
           dos_total_s(k,n)=dos_total(k,n,2)
           dos_total_p(k,n)=dos_total(k,n,3)+dos_total(k,n,4)+
     &dos_total(k,n,5)
           dos_total_d(k,n)=dos_total(k,n,6)+dos_total(k,n,7)+
     &dos_total(k,n,8)+dos_total(k,n,9)+dos_total(k,n,10)
         enddo
        enddo
        close(13)
!!      atom projected dos
       open(unit=14,file="$DOSNAME_atom",status='unknown')
        do n=1,neig
        write(14,'(A)',ADVANCE="NO")"## k_dist      eigen  atom1"
        write(14,'(A)',ADVANCE="NO")"  atom2  atom3  atom4  atom5"
        write(14,'(A)',ADVANCE="NO")"  atom6  atom7  atom8  atom9"
        write(14,'(A)')" atom10 atom11 atom12"
         do k=1,nkp
         write(14,'(F9.5,F11.5)',ADVANCE='NO')band_dist(k),eig(k,n)
          do i=1,natom
           write(14,'(F7.3)',ADVANCE='NO')dos_atom_tot(k,n,i)
          enddo
          write(14,'()')
         enddo
         write(14,'()')
        enddo
        write(6,*)"Atom projected eigenvalue file generated:
     $ DOS_atom_projected.dat"
!!      orbital projected dos
        open(unit=15,file="$DOSNAME_spd",status='unknown')
        do n=1,neig
         write(15,'(A)',ADVANCE="NO")"# k_dist      eigen     s      p"
         write(15,'(A)',ADVANCE="NO")"      d       s       py     pz"
         write(15,'(A)',ADVANCE="NO")"     px      dxy    dyz    dz2"
         write(15,'(A)')"    dxz    dx2     tot"
         do k=1,nkp
        write(15,'(F9.5,F11.5,3F7.3,1X,F7.3,1X,3F7.3,1X,5F7.3,1X,F7.3)')
     &band_dist(k),eig(k,n),dos_total_s(k,n),dos_total_p(k,n),
     &dos_total_d(k,n),(dos_total(k,n,l),l=2,11)
         enddo
         write(15,'()')
        enddo
        write(6,*)"s,px,py,pz,dxy.. projected eigenvalue file generated:
     $ DOS_spd_projected.dat"
        end
EOF
ifort -o band_decomposed band_decomposed.f
#./band_decomposed 


cat >> gnuBAND-decomp.gp <<EOF
set term post portrait  enhanced color "Helvetica,18"
set output 'BAND-$SYS-decomp.eps'

#  $SYS  ###############################
set title "Energy band $SYS" font "Curier Bold,18,"
 set size  ratio 1;
 set xrange[$Xinit:$Xend]
 set yrange[$Yinit:$Yend]
#set lmargin 1.2;set bmargin 1.5
#set xtics 1,2,31 
#set xtics ("{/Symbol \G}" $Xinit, "X" 0.69772 ,"M" 1.32016,"{/Symbol \G}" $Xend)
set xtics ("$K_init_name" $Xinit, "$K_2_name" $K_2_point, "$K_3_name" $K_3_point, "$K_4_name" $K_4_point,"$K_end_name" $Xend)
set ylabel "Energy (eV)" font "Curier Bold,24";
 x1=$Xinit; x3=$Xend; y3=0.0;unset arrow
 set arrow from x1,y3 to x3,y3 nohead lt 3 lw 1 
 set arrow from $K_2_point,$Yinit to $K_2_point,$Yend nohead lt 3 lw 1 
 set arrow from $K_3_point,$Yinit to $K_3_point,$Yend nohead lt 3 lw 1 
 set arrow from $K_4_point,$Yinit to $K_4_point,$Yend nohead lt 3 lw 1 
# Define palette
h1=0/360.
h2=360/360.
set palette model HSV functions (1-gray)*(h2-h1)+h1,1,0.68
# Define RGB color variable
 rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b)
set palette rgbformulae 33,13,10
#s p d s py pz px dxy dyz dz2 dxz dx2 tot 
#3 4 5 6 7  8  9  10  11  12  13  14  15  

###  "$DOSNAME_atom" every :$Degeneracy using 1:2:((`A="$""$atom1" ;echo $A`+`A="$""$atom2" ;echo $A`)*40) w $PT
#plot "$DOSNAME_atom" using 1:2 w $LT notitle,\\
#     "$DOSNAME_atom" using 1:2:((`A="$""$atom1";echo $A`)*3) w $PTa1,\\
#     "$DOSNAME_atom" using 1:2:((`A="$""$atom2";echo $A`)*3) w $PTa2
plot "$DOSNAME_spd" using 1:2:((`A="$""$spd1";echo $A`)*1) w l lw 7 lc variable notitle

EOF
#gnuplot gnuBAND-decomp.gp
#eps2eps BAND-$SYS-decomp.eps a.eps;mv a.eps BAND-$SYS-decomp.eps
