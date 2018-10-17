#!/usr/bin/env python

from ase.io.trajectory import Trajectory as Traj
from ase.io import read
import sys
from numpy import array

print("\n")
print("#######################################################################################".center(100))
print("%%%%%%%%%%% This code will covert lammps results to traj file %%%%%%%%%".center(100))
print("useage ==> ./lmp2traj.py 'dump file' 'energy logfile'".center(100))
print("e.g.) ./lmp2traj.py out.dump log.lammps".center(100))
print("The result file name will be lmp-result.traj".center(100))
print("*****NOTE***** There is some issue when ase.io.lammpsrun import dump file. *****NOTE*****".center(100))
print("Make sure that you revised it. (velocity-> vel /1000/units.fs, symbol issue)".center(100))
print("#######################################################################################".center(100))
if len(sys.argv) is 3:
    print("The Number of arguments is correct (n=2).".center(100))
    print("\n")
    dump_inpf = sys.argv[1]
    log_f = sys.argv[2]
elif len(sys.argv) is 1:
    print("No argument is provided. Load 'out.dump' & 'log.lammps' as arguments".center(100))
    print("\n")
    dump_inpf = 'out.dump'
    log_f = 'log.lammps'
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(100))
    print("\n")
    sys.exit(1)

dump_inp = read(dump_inpf,':')
log = open(log_f, "r")
epot = []

j = 0
while True:
    line = log.readline()
    llist = line.split()
    j += 1
    print("Reading "+str(j)+"th line")

    seq = {}
    for (i, x) in enumerate(llist):
        seq[x] = i

    for word in llist:
        if word == 'PotEng':
            epot.append(llist[seq['PotEng']+2])

    if not line: break

if len(epot) != len(dump_inp):
    print(" ***** ERROR ***** The # of images in log file and dump file do not match\n")
    print("   i.e.) len(epot) != len(dump_inp)\n")
    print("   "+str(len(epot))+" != "+str(len(dump_inp))+"\n")
    sys.exit(1)

traj = Traj("lmp-result.traj", "w")
for i in range(len(dump_inp)):
    atoms = dump_inp[i]
    atoms._pbc = array([True, True, True])
    atoms._calc.atoms.pbc = array([True, True, True])
    atoms._calc.results['energy'] = float(epot[i])
    traj.write(atoms)
    print("Writing "+str(i)+"th image")
traj.close()


print("\n\n#######################################################################################")
print("      %%%%%%%%%%% This code will covert lammps results to traj file %%%%%%%%%")
print("         useage ==> ./lmp2traj.py 'dump file' 'energy logfile'")
print("              e.g.) ./lmp2traj.py out.dump log.lammps")
print("                The result file name will be lmp-result.traj")
print("  *****NOTE***** There is some issue when ase.io.lammpsrun import dump file. *****NOTE*****")
print("       Make sure that you revised it. (velocity-> vel /1000/units.fs, symbol issue)")
print("#######################################################################################\n\n")

