#!/usr/bin/env python
from ase.io.trajectory import Trajectory
from ase import Atom, Atoms
import datetime
import sys
from ase.io import read, write
from ase.calculators.calculator import PropertyNotImplementedError
import subprocess as sp

now = datetime.datetime.now()
time = now.strftime('%Y-%m-%d %H:%M:%S')
print("")
print("***** Code by Youngjae Choi @ POSTECH *****".center(100))
print(("code started time: "+time).center(100))
print("")
print("############################################################".center(100))
print("useage ==> ./dpmd_traj2raw.py 'trajactory file'".center(100))
print("############################################################".center(100))
if len(sys.argv) is 2:
    print(" ")
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(100))
    print("")
    print("")
    sys.exit(1)

trajfile = sys.argv[1]

print(("I'll extract informations from '"+trajfile+"' file.").center(100))
print(("Raw files for DeePMD code will be writen in 'raw-"+trajfile+".d' directory.").center(100))
print(" i.e.) `box.raw`, `coord.raw`, `force.raw`, `energy.raw` and `virial.raw`".center(100))

traj = Trajectory(trajfile,"r")
natom = len(traj[0])
print("\n#####################################################################")
print("\n*Number of frames :: "+str(len(traj)))
print("  *First frame contains :: "+str(natom)+" atoms.")
print("  *First frame's chemical symbols ::")
print(str(traj[0].get_chemical_symbols()))
print("     Caution)) Every frame mush have identical number and type of atoms.")
print("  *First frame's pbc :: "+str(traj[0].get_pbc()))
print("  *First frame's lattice vectors ::")
print(""+str(traj[0].get_cell()))
print("\n")
sp.call(["rm -rf raw-"+trajfile+".d.old"], shell=True)
sp.call(["mv raw-"+trajfile+".d raw-"+trajfile+".d.old"], shell=True)
sp.call(["mkdir raw-"+trajfile+".d"], shell=True)

################# box.raw ####################
box_raw = open("raw-"+trajfile+".d/box.raw", "w")
n=0
print("box.raw :: writing")
for atoms in traj:
    n+=1
    # print("box.raw :: writing "+str(n)+" th frame.")
    cell_array = atoms.get_cell()
    for i in range(3):
        for j in range(3):
            box_raw.write(str(cell_array[i][j])+" ")
        box_raw.write("    ")
    box_raw.write("\n")

################# type.raw ###################
type_raw = open("raw-"+trajfile+".d/type.raw", "w")
type_ref = open("raw-"+trajfile+".d/ref_type_symbols.txt", "w")
print("type.raw :: writing")
from ss_util import list2numlist as l2nl
symbols = traj[0].get_chemical_symbols()
type_ref.write(str(symbols))
if traj[-1].get_chemical_symbols() != symbols:
    raise ValueError("Chemical symbols seem to be not consistent. Please check")
symbols_num = l2nl(list(symbols))
for nums in symbols_num:
    type_raw.write(str(nums)+" ")

################# energy.raw #################
try:
    traj[0].get_potential_energy()
except PropertyNotImplementedError:
    print("#####################################################################")
    print("        No energy matrix information in trajectory you gave")
    print("#####################################################################")
except:
    print("\n#####################################################################")
    print("\n        Something wrong with energy information")
    print("\n#####################################################################")
else:
    E_raw = open("raw-"+trajfile+".d/energy.raw", "w")
    n=0
    print("energy.raw :: writing")
    for atoms in traj:
        n+=1
        # print("energy.raw :: writing "+str(n)+" th frame.")
        E_raw.write(str(atoms._calc.results['energy'])+"\n")
    E_raw.close()

################# force.raw ##################
try:
    traj[0].get_forces()
except PropertyNotImplementedError:
    print("#####################################################################")
    print("        No forces matrix information in trajectory you gave")
    print("#####################################################################")
except:
    print("\n#####################################################################")
    print("\n        Something wrong with forces information")
    print("\n#####################################################################")
else:
    F_raw = open("raw-"+trajfile+".d/force.raw", "w")
    n=0
    print("force.raw :: writing")
    for atoms in traj:
        n+=1
        # print("force.raw :: writing "+str(n)+" th frame.")
        force_array = atoms._calc.results['forces']
        for i in range(natom):
            for j in range(3):
                F_raw.write(str(force_array[i][j])+" ")
            F_raw.write("    ")
        F_raw.write("\n")
    F_raw.close()

################# virial.raw #################
try:
    traj[0].get_stress()
except PropertyNotImplementedError:
    print("#####################################################################")
    print("        No virial tensor information in trajectory you gave")
    print("#####################################################################")
except:
    print("\n#####################################################################")
    print("\n        Something wrong with stress information")
    print("\n#####################################################################")
else:
    V_raw = open("raw-"+trajfile+".d/virial.raw", "w")
    n=0
    print("virial.raw :: writing")
    for atoms in traj:
        n+=1
        # print("virial.raw :: writing "+str(n)+" th frame.")
        stress_array = atoms._calc.results['stress']
        for i in [[0,5,4],[5,1,3],[4,3,2]]:
            for j in i:
                V_raw.write(str(stress_array[j])+" ")
            V_raw.write("    ")
        V_raw.write("\n")

################# coord.raw #################
coord_raw = open("raw-"+trajfile+".d/coord.raw", "w")
n=0
print("coord.raw :: writing")
for atoms in traj:
    n+=1
    # print("coord.raw :: writing "+str(n)+" th frame.")
    #print atoms.get_positions()
    atoms.wrap(eps=0)
    #print atoms.get_positions()
    posi_array = atoms.arrays['positions']
    for i in range(natom):
        for j in range(3):
            coord_raw.write(str(posi_array[i][j])+" ")
        coord_raw.write("    ")
    coord_raw.write("\n")
print("")
print("")
