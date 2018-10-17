#!/usr/bin/env python

import numpy as np
import sys


print("\n\n#######################################################################################\n")
print("useage ==> ./vasp-volume.py 'POSCAR or CONTCAR file'\n")
print("           EXAMPLE) ./vasp-volume.py POSCAR-unitcell.vasp\n")
print("#######################################################################################")
if len(sys.argv) is 2:
    print("          The Number of arguments is correct.\n\n")
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
    sys.exit(1)

inp_fname = sys.argv[1]
    
########### cell and # of atoms ###########

inp = open(inp_fname, "r")
lines = inp.readlines()
inp.close()
cell = []
i = 0
natom = 0
for line in lines:
    l = line.split()
    if len(l) != 0:
        if len(l) != 3 and len(l) != 6:
            continue
        if len(l) == 3 or len(l) == 6:
            if i < 3:
                cell.append(l)
                i += 1
            elif i == 3:
                natom += 1
    if len(l) == 0:
        break

cell_array = np.asarray(cell).astype(np.float)

##########  cell volume  ###############

c1 = cell_array[0,:]
c2 = cell_array[1,:]
c3 = cell_array[2,:]

vol = np.dot(c1, np.cross(c2, c3))

########## print ############

print " <<< cell array >>>"
print cell_array
print "\n <<< cell volume = "+str(vol)+" >>>"
print "\n <<< # of atoms = "+str(natom)+" >>>\n"

