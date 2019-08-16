#!/usr/bin/env python

from ase.io.trajectory import Trajectory as Traj
from ase.io import read
import sys
from numpy import array

print("\n")
print("#######################################################################################".center(120))
print("%%%%%%%%%%% This code will covert lammps results to traj file %%%%%%%%%".center(120))
print("useage ==> ./lmp2traj.py 'dump file' 'energy logfile'".center(120))
print("e.g.) ./lmp2traj.py out.dump log.lammps".center(120))
print("The result file name will be lmp-result.traj".center(120))
print("*****NOTE***** There is some issue when ase.io.lammpsrun import dump file. *****NOTE*****".center(120))
print("Make sure that you revised it. (velocity-> vel /1000/units.fs, symbol issue)".center(120))
print("#######################################################################################".center(120))
if len(sys.argv) is 3:
    print("The Number of arguments is correct (n=2).".center(120))
    print("\n")
    dump_inpf = sys.argv[1]
    log_f = sys.argv[2]
elif len(sys.argv) is 1:
    print("No argument is provided. Load 'out.dump' & 'log.lammps' as arguments".center(120))
    print("\n")
    dump_inpf = 'out.dump'
    log_f = 'log.lammps'
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(120))
    print("\n")
    sys.exit(1)

print('Reading dump file...'.center(120))
dump_inp = read(dump_inpf,':',format='lammps-dump')
print('Successively read dump file!'.center(120))
log = open(log_f, "r")
epot = []

while True:
    line = log.readline()
    llist = line.split()
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
    if i % 1000 == 0:
        print("Writing "+str(i)+"th image")
    atoms = dump_inp[i]
    atoms._pbc = array([True, True, True])
    atoms._calc.atoms.pbc = array([True, True, True])
    atoms._calc.results['energy'] = float(epot[i])
    traj.write(atoms)
traj.close()


print("\n\n#######################################################################################")
print("      %%%%%%%%%%% This code will covert lammps results to traj file %%%%%%%%%")
print("         useage ==> ./lmp2traj.py 'dump file' 'energy logfile'")
print("              e.g.) ./lmp2traj.py out.dump log.lammps")
print("                The result file name will be lmp-result.traj")
print("  *****NOTE***** There is some issue when ase.io.lammpsrun import dump file. *****NOTE*****")
print("       Make sure that you revised it. (velocity-> vel /1000/units.fs, symbol issue)")
print("#######################################################################################\n\n")

