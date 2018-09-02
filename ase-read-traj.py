#!/usr/bin/env python
import sys

print("\n\n")
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(80))
print("___________________________".center(80))
print(" __________|  C o d e  b y  Y.J. Choi  |_________".center(80))
print("|______________ ssrokyz@gmail.com _______________|".center(80))
print("")
print("*******   This code will read traj file   *******".center(80))
print("useage ==> ./ase-read-traj.py 'trajectory file'".center(80))
print("EXAMPLE) ./ase-read-traj.py si.traj".center(80))
print("")
print("     traj  : Traj file you gave as an object")
print("     ntraj : len(traj)")
print("     atoms : traj[0]")
print("     calc  : atoms._calc")
print("")
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(80))
print("")
if len(sys.argv) is 2:
    print("The Number of arguments is correct.".center(80))
    print("\n")
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(80))
    print("\n")
    sys.exit(1)

from ase.io import Trajectory as Traj

traj = Traj(sys.argv[1], "r")
atoms = traj[0]
ntraj = len(traj)
calc = atoms._calc
