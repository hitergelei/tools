import sys
import numpy as np

print("\n\n")
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(80))
print("            ___________________________           ".center(80))
print(" __________|  C o d e  b y  Y.J. Choi  |_________ ".center(80))
print("|______________ ssrokyz@gmail.com _______________|".center(80))
print("")
print("*******   This code will read ASE files   *******".center(80))
print("useage ==> ./ase-read-traj.py 'file' ('format')".center(80))
print("EXAMPLE) ./ase-read-traj.py vasprun.xml (vasp-xml)".center(80))
print("")
print("     alist   : Atoms list you gave as an object")
print("     nlist   : len(alist)")
print("     atoms   : alist[0]")
print("     calc    : atoms._calc")
print("     results : atoms._calc.results")
print("     Traj    : ase.io.trajectory:Trajectiory == Traj")
print("    coord    : atoms.get_positions()")
print("     frac    : atoms.get_scaled_positions()")
print("    forces   : atoms.get_forces()")
print("    energy   : atoms.get_potential_energy()")
print("    stress   : atoms.get_stress()")
print("")
print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(80))
print("")
if len(sys.argv) is 2 or len(sys.argv) is 3:
    print(("The Number of arguments(= %d) is correct." %(len(sys.argv)-1)).center(80))
    print("\n")
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(80))
    print("\n")
    sys.exit(1)

from ase.io import read
from ase.io.trajectory import Trajectory
Traj = Trajectory

alist = read(
    sys.argv[1],
    index  = ':',
    format = None if len(sys.argv) is 2 else sys.argv[2],
    )
try:
    nlist = len(alist)
except:
    pass
try:
    atoms = alist[0]
except:
    pass
try:
    calc = atoms._calc
except:
    pass
try:
    results = calc.results
except:
    pass
try:
    forces = atoms.get_forces()
except:
    pass
try:
    energy = atoms.get_potential_energy()
except:
    pass
try:
    stress = atoms.get_stress()
except:
    pass
try:
    coord = atoms.get_positions()
except:
    pass
try:
    frac = atoms.get_scaled_positions()
except:
    pass
