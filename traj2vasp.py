#!/usr/bin/env python2
from ase.io.trajectory import Trajectory
from ase.io import read, write
from ase.calculators import vasp
import sys
import commands

print("\n\n##################################################\n")
print("useage ==> ./traj2vasp.py 'trajactory file'\n")
print("##################################################\n")
if len(sys.argv) is 2:
	print(" ")
else:
	print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
	sys.exit(1)

trajfile = sys.argv[1]

print("\nI'll extract informations from '"+trajfile+"' file.")
print("Informations will be writen in folder, 'vasp_images_"+trajfile+"'")

traj = Trajectory(trajfile, "r")
steps = len(traj)
print("\nTotally, "+str(steps)+" images of each step will be writen\n")


commands.getstatusoutput("mkdir vasp_images_"+trajfile)
imagenum = 0
for image in traj:
	write("./vasp_images_"+trajfile+"/POSCAR", image)
 	raw = open("./vasp_images_"+trajfile+"/POSCAR", "r")
 	lines = raw.readlines()
	raw.close()
	poscar = open("./vasp_images_"+trajfile+"/POSCAR_%.5d"%imagenum, "w")
	poscar.write("image "+str(imagenum)+"/"+str(steps-1)+"\n")
	poscar.writelines(lines[1:5])
	poscar.writelines(lines[5:])
	poscar.close()
	imagenum += 1
	
	
commands.getstatusoutput("rm -rf ./vasp_images_"+trajfile+"/POSCAR")
