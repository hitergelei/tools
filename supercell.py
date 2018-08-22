#!/usr/bin/env python2
##### CODE BY YOUNGJAE CHOI #####

from ase import Atoms, Atom
from ase.io.trajectory import Trajectory
from ase.build import make_supercell
import sys

print("\n\n#######################################################################################\n")
print("useage ==> ./supercell.py 'trajactory file of single image' 'supercell integer matrix'\n")
print("           EXAMPLE) ./supercell.py amp.traj [200,020,001]\n")
print("#######################################################################################")
if len(sys.argv) is 3:
	print(" ")
else:
	print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
	sys.exit(1)

trajfile = sys.argv[1]
traj = Trajectory(trajfile, "r")

raw = list(sys.argv[2])
matrix = [[int(raw[1]),int(raw[2]),int(raw[3])],
		  [int(raw[5]),int(raw[6]),int(raw[7])],
		  [int(raw[9]),int(raw[10]),int(raw[11])]]

if len(traj) == 1:
	print " "
else:
	print "\nTrajectory file you gave has wrong number of images.(1 required)\nPlease check and try again\n"
	sys.exit(1)

print("I'll make supercell of '"+trajfile+"' file.")
print("Supercell trajectory file will be writen as 'supercell_"+trajfile+"' file.")

atoms = traj[0]
print("\n#############  Primitive cell  #############")
print("Total number of atom = "+str(len(atoms)))
print("Cell =")
print(atoms.get_cell())

supercell = make_supercell(atoms, matrix)
print("\n#############  Supercell cell  #############")
print("Transformation matrix, P ="),;print(matrix)
print ""
print("Total number of atom = "+str(len(supercell)))
print("Cell =")
print(supercell.get_cell())
print ""
print("pbc ="),;print(supercell.get_pbc())
print ""

supertraj = Trajectory("supercell_"+trajfile, "w")
supertraj.write(supercell)
supertraj.close()
