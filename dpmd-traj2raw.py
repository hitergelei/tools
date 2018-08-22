#!/usr/bin/env python2
from ase.io.trajectory import Trajectory
from ase import Atom, Atoms
import datetime
import sys
from ase.io import read, write
from commands import getstatusoutput
from ase.calculators.calculator import PropertyNotImplementedError

now = datetime.datetime.now()
time = now.strftime('%Y-%m-%d %H:%M:%S')
print("\n      ***** Code by Youngjae Choi @ POSTECH *****")
print("      Moment of code start : "+time)
print("\n##################################################\n")
print("        useage ==> ./dpmd_traj2raw.py 'trajactory file'\n")
print("##################################################")
if len(sys.argv) is 2:
	print(" ")
else:
	print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
	sys.exit(1)

trajfile = sys.argv[1]

print("I'll extract informations from '"+trajfile+"' file.")
print("Raw files for DeePMD code will be writen in 'raw-"+trajfile+".d' directory.")
print(" i.e.) `box.raw`, `coord.raw`, `force.raw`, `energy.raw` and `virial.raw`")
print(" %%%%%%%%%%%%%%%''type.raw'' file need to be generated manually %%%%%%%%%%%%%%%")

traj = Trajectory(trajfile,"r")
natom = len(traj[0])
print("\n#####################################################################")
print("\n*Number of frames :: "+str(len(traj)))
print("  *First frame contains :: "+str(natom)+" atoms.")
print("  *First frame's chemical symbols ::")
print("     => "+str(traj[0].get_chemical_symbols()))
print("     Caution)) Every frame mush have identical number and type of atoms.")
print("  *First frame's pbc :: "+str(traj[0].get_pbc()))
print("  *First frame's lattice vectors ::")
print(""+str(traj[0].get_cell()))
print("\n")
getstatusoutput("rm -rf raw-"+trajfile+".d.old")
getstatusoutput("mv raw-"+trajfile+".d raw-"+trajfile+".d.old")
getstatusoutput("mkdir raw-"+trajfile+".d")

################# box.raw ####################
box_raw = open("raw-"+trajfile+".d/box.raw", "w")
n=0
for atoms in traj:
	n+=1
	print("box.raw :: writing "+str(n)+" th frame.")
	cell_array = atoms.get_cell()
	for i in range(3):
		for j in range(3):
			box_raw.write(str(cell_array[i][j])+" ")
		box_raw.write("    ")
	box_raw.write("\n")

################# type.raw ###################
type_raw = open("raw-"+trajfile+".d/type.raw", "w")
print("type.raw :: writing type.raw file")
for atom in traj[0]:
	type_raw.write(atom.symbol+" ")

################# energy.raw #################
E_raw = open("raw-"+trajfile+".d/energy.raw", "w")
n=0
for atoms in traj:
	n+=1
	print("energy.raw :: writing "+str(n)+" th frame.")
	E_raw.write(str(atoms._calc.results['energy'])+"\n")
E_raw.close()

################# force.raw ##################
F_raw = open("raw-"+trajfile+".d/force.raw", "w")
n=0
for atoms in traj:
	n+=1
	print("force.raw :: writing "+str(n)+" th frame.")
	force_array = atoms._calc.results['forces']
	for i in range(natom):
		for j in range(3):
			F_raw.write(str(force_array[i][j])+" ")
		F_raw.write("    ")
	F_raw.write("\n")
F_raw.close()

################# virial.raw #################
try:
	traj[0].get_stress()
except PropertyNotImplementedError:
	print("\n#####################################################################")
	print("\n        No virial tensor information in trajectory you gave")
	print("\n#####################################################################")
except:
	print("\n#####################################################################")
	print("\n        Something wrong with stress information")
	print("\n#####################################################################")
else:
	V_raw = open("raw-"+trajfile+".d/virial.raw", "w")
	n=0
	for atoms in traj:
		n+=1
		print("virial.raw :: writing "+str(n)+" th frame.")
		stress_array = atoms._calc.results['stress']
		for i in range(3):
			for j in range(3):
				V_raw.write(str(stress_array[i][j])+" ")
			V_raw.write("    ")
		V_raw.write("\n")

################# coord.raw #################
coord_raw = open("raw-"+trajfile+".d/coord.raw", "w")
n=0
for atoms in traj:
	n+=1
	print("coord.raw :: writing "+str(n)+" th frame.")
	#print atoms.get_positions()
	atoms.wrap(eps=0)
	#print atoms.get_positions()
	posi_array = atoms.arrays['positions']
	for i in range(natom):
		for j in range(3):
			coord_raw.write(str(posi_array[i][j])+" ")
		coord_raw.write("    ")
	coord_raw.write("\n")

print("\n#####################################################################")
print("\n#####################################################################")
print("\n      NOTE :: type.raw FILE MUST BE REVISED BY YOU !!!!!!!!!!!")
print("\n                   (symbols to numbers from zero")
print("\n#####################################################################")
print("\n#####################################################################")
