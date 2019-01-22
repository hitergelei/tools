#!/usr/bin/env python

import numpy as np
from ase.io.trajectory import Trajectory as Traj
import sys

print("\n")
print("#######################################################################################".center(100))
print("%%%%%%%%%%% This code will give you RMSE of two traj files %%%%%%%%%".center(100))
print("useage ==> ./ase-rmse.py 'traj file1' 'traj file2'".center(100))
print("e.g.) ./ase-rmse.py conv-test.traj vasp-result.traj".center(100))
print("#######################################################################################".center(100))
if len(sys.argv) is 3:
    print("The Number of arguments is correct.".center(100))
    print("")
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(100))
    print("\n")
    sys.exit(1)

traj1_f = sys.argv[1]
traj2_f = sys.argv[2]
traj1 = Traj(traj1_f)
traj2 = Traj(traj2_f)

if len(traj1) != len(traj2):
    print(" ***** ERROR ***** The # of images is different between two traj files.".center(100))
    sys.exit(1)
if len(traj1[0]) != len(traj2[0]):
    print(" ***** ERROR ***** The # of atoms is different between two traj files.".center(100))
    sys.exit(1)
    
nimage = len(traj1)
natom = len(traj1[0])
try:
    traj1[0].get_potential_energy()
    traj2[0].get_potential_energy()
except:
    print('No potential energy info.'.center(100))
    e_ok = False
else:
    e_ok = True
try:
    traj1[0].get_forces()
    traj2[0].get_forces()
except:
    print('No force info.'.center(100))
    f_ok = False
else:
    f_ok = True

e1_arr = np.array([])
f1_arr = np.array([])
e2_arr = np.array([])
f2_arr = np.array([])
if e_ok:
    for i in range(nimage):
        e1_arr = np.append(e1_arr, traj1[i].get_potential_energy())
        e2_arr = np.append(e2_arr, traj2[i].get_potential_energy())
    print("calculating energy RMSE...".center(100))
    ermse = np.sqrt(np.mean((e1_arr - e2_arr)**2))
    adjust = np.mean(e1_arr) - np.mean(e2_arr)
    edrmse = np.sqrt(np.mean((e1_arr - adjust - e2_arr)**2))
    ermse_per_atom = ermse / natom
    edrmse_per_atom = edrmse / natom
if f_ok:
    for i in range(nimage):
        f1_arr = np.append(f1_arr, traj1[i].get_forces())
        f2_arr = np.append(f2_arr, traj2[i].get_forces())
    print("calculating force RMSE...".center(100))
    frmse = np.sqrt(np.mean((f1_arr - f2_arr)**2))



print("")
print("       ****** RESULTS ******")
print("")
if e_ok:
    print(("    Energy RMSE / atom            = "+str(ermse_per_atom*1000)+" meV/atom"))
    print(("    Energy difference RMSE / atom = "+str(edrmse_per_atom*1000)+" meV/atom"))
if f_ok:
    print(("    Force RMSE                    = "+str(frmse*1000)+" meV/Ang\n\n"))

