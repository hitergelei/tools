#!/usr/bin/env python

import numpy as np
import sys

print("\n")
print("#######################################################################################".center(120))
print("")
print("useage ==> ./vasp-volume.py 'POSCAR or CONTCAR file'".center(120))
print("           EXAMPLE) ./vasp-volume.py POSCAR-unitcell.vasp".center(120))
print("")
print("#######################################################################################".center(120))
print("")
if len(sys.argv) is 2:
    print("The Number of arguments is correct.".center(120))
    print("\n")
else:
    print(">>>>> ERROR <<<<<     The number of arguments is not correct     >>>>> ERROR <<<<<".center(120))
    print("\n")
    sys.exit(1)

inp_fname = sys.argv[1]
    
## Open file
from ase.io import read
atoms = read(inp_fname, -1)

## Data
cell = atoms.get_cell()
vol = np.linalg.det(cell)
num = len(atoms)
vpa = vol / num

## Print
print(" <<< cell array >>>")
print(cell)
print("\n <<< cell volume = "+str(vol)+" >>>")
print(" <<< # of atoms = "+str(num)+" >>>\n")
print(" <<< volume/atom (Ang^3) = "+str(vpa)+" >>>\n")

