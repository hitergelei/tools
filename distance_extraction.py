#!/usr/bin/env python

from ase.io import read, write
from ase.io.trajectory import Trajectory
import sys
import datetime

print("\n\n##################################################\n")
print("useage ==> ./distance_extraction.py 'trajactory file'\n")
print("##################################################\n")
if len(sys.argv) is 2:
    print(" ")
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
    sys.exit(1)

trajfile = sys.argv[1]

print("\nI'll extract informations from '"+trajfile+"' file.")
print("\nDistance information will be writen in '"+trajfile+"_dist.log' file.")

traj = Trajectory(trajfile,"r")
print("\nTotally, "+str(len(traj))+" images in trajectory file.\n")

logfile = open("distance_"+trajfile+".log", "w")
now = datetime.datetime.now()
time = now.strftime('%Y-%m-%d %H:%M:%S')
atomnum = len(traj[0])
logfile.write("\n***** code by Youngjae Choi *****")
logfile.write("\n\nlog file of "+trajfile+"\nexcuted at "+time)
logfile.write("\nTotal number of images : "+str(len(traj)))
logfile.write("\nTotal number of atoms in each image : "+str(atomnum))
logfile.write("\n\nAtom sequence :\n")

i = 0
for atom in traj[0]:
    i += 1
    logfile.write("\n"+str(i)+". "+atom.symbol+"(AN = "+str(atom.number)+")(AM = "+str(atom.mass)+")")

a=0
b=0
c=0
trajlist = []
for atoms in traj:
    a += 1
    logfile.write("\n\n########################################################")
    logfile.write("\nDistance matrix of image "+str(a)+"/"+str(len(traj)))
    logfile.write("\n########################################################")
    atoms.wrap(0)
    logfile.write("\n\n"+str(atoms.get_all_distances()))
    atomlist = []
    for i in range(atomnum):
        for j in range(atomnum):
            if i < j:
                atomlist.append(atoms.get_all_distances()[i][j])
                trajlist.append(atoms.get_all_distances()[i][j])
    atomlist.sort()
    logfile.write("\n\n sorted distances for image "+str(a)+"/"+str(len(traj))+" :\n")
    logfile.write(str(atomlist))

trajlist.sort()
logfile.write("\n\n\n\n\n\n\n\n\n#######################################################")
logfile.write("\n sorted distances for total images of trajectory :")
logfile.write("\n#######################################################\n\n")
logfile.write(str(trajlist))
logfile.write("\n\nminimum distance : "+str(min(trajlist))+" Angstrom\n")
logfile.write("maximum distance : "+str(max(trajlist))+" Angstrom\n")

logfile.close()
