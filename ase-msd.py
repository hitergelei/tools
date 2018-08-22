#!/usr/bin/env python

import numpy as np
from ase.io.trajectory import Trajectory
import sys


print("\n\n#######################################################################################")
print("      %%%%%%%%%%% This code will give you Mean Square Displacement %%%%%%%%%%%\n")
print("useage ==> ./ase-msd.py 'init structure traj file(1 image)' 'traj file(N images)")
print("    if init structure traj file has more than one image, it will use traj[0] image.")
print("           EXAMPLE) ./ase-msd.py si-init.traj si.traj")
print("#######################################################################################")
if len(sys.argv) is 3:
	print("          The Number of arguments is correct.\n")
else:
	print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
	sys.exit(1)

tf1 = sys.argv[1]
tf2 = sys.argv[2]
t1 = Trajectory(tf1, "r")
t2 = Trajectory(tf2, "r")
ref = t1[0]
natom = len(ref)
if len(ref) != len(t2[0]):
	print("*****ERROR***** Two traj files have different number of atoms\n")
	sys.exit(1)

add = 0
num = 0
for image in t2:
	for i in range(natom):
		disp = np.linalg.norm(ref.arrays["positions"][i]-image.arrays["positions"][i])
		add += (disp ** 2)
		num += 1
		print num
		print disp
msd = add / num
print("********** RESULT :: MSD = <|r(t)-r(0)|^2> = "+str(msd)+" *************\n\n")

