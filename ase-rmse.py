#!/usr/bin/env python

import numpy as np
from ase.io.trajectory import Trajectory as Traj
import sys

print("\n\n#######################################################################################")
print("      %%%%%%%%%%% This code will give you RMSE of two traj files %%%%%%%%%")
print("         useage ==> ./ase-rmse.py 'traj file1' 'traj file2'")
print("              e.g.) ./ase-rmse.py conv-test.traj vasp-result.traj")
print("#######################################################################################")
if len(sys.argv) is 3:
    print("                The Number of arguments is correct.\n\n")
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
    sys.exit(1)

traj1_f = sys.argv[1]
traj2_f = sys.argv[2]
traj1 = Traj(traj1_f)
traj2 = Traj(traj2_f)

if len(traj1) != len(traj2):
	print(" ***** ERROR ***** The # of images is different between two traj files.")
	sys.exit(1)
if len(traj1[0]) != len(traj2[0]):
	print(" ***** ERROR ***** The # of atoms is different between two traj files.")
	sys.exit(1)
	
nimage = len(traj1)
natom = len(traj1[0])
ermse = 0
frmse = 0
for i in range(nimage):
	print("    "+str(i)+"th image / total "+str(nimage)+" images")
	ermse += (traj1[i]._calc.results['energy']-traj2[i]._calc.results['energy'])**2
	for j in range(natom):
		for k in range(3):
			frmse += (traj1[i]._calc.results['forces'][j][k]-traj2[i]._calc.results['forces'][j][k])**2

ermse = ermse / nimage
ermse_per_atom = ermse / natom
frmse = frmse / ( nimage * natom * 3 )

print("****** RESULTS ******")
print("Energy RMSE / atom = "+str(ermse_per_atom*1000)+" meV/atom")
print("Force RMSE         = "+str(frmse*1000)+" meV\n\n")

