#!/usr/bin/env python

from ase.io.trajectory import Trajectory as Traj
from ase.io import read
import sys
from numpy import array

print("\n\n#######################################################################################")
print("      %%%%%%%%%%% This code will covert lammps results to traj file %%%%%%%%%")
print("         useage ==> ./lmp2traj.py 'thermo logfile'")
print("              e.g.) ./lmp2traj.py log.lammps")
print("                The result file name will be termo.log")
print("#######################################################################################")
if len(sys.argv) is 2:
    print("                The Number of arguments is correct.\n\n")
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
    sys.exit(1)

log_f = sys.argv[1]
log = open(log_f, "r")

etot = []
epot = []
ekin = []
temp = []
pres = []

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
		if word == 'TotEng':
			etot.append(llist[seq['TotEng']+2])
		elif word == 'PotEng':
			epot.append(llist[seq['PotEng']+2])
		elif word == 'KinEng':
			ekin.append(llist[seq['KinEng']+2])
		elif word == 'Temp':
			temp.append(llist[seq['Temp']+2])
		elif word == 'Press':
			pres.append(llist[seq['Press']+2])
	if not line: break

log.close()
log = open("thermo.log", "w")
log.write("i     etot        epot     ekin  temperature   pressure")
for i in range(len(etot)):
	log.write("\n"+str(i)+"  "+etot[i]+"  "+epot[i]+"  "+ekin[i]+"  "+temp[i]+"  "+pres[i])
	print("Writing "+str(i)+"th log data")


print("\n\n#######################################################################################")
print("      %%%%%%%%%%% This code will covert lammps results to traj file %%%%%%%%%")
print("         useage ==> ./lmp2traj.py 'thermo logfile'")
print("              e.g.) ./lmp2traj.py log.lammps")
print("                The result file name will be termo.log")
print("#######################################################################################")

