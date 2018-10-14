#!/usr/bin/env python
##### CODE BY YOUNGJAE CHOI #####

from ase.io.trajectory import Trajectory
from ase import Atom, Atoms
from ase.calculators.singlepoint import SinglePointCalculator as spc
from ase.calculators.calculator import PropertyNotImplementedError as pnie
from ase.build import make_supercell
import sys
from copy import deepcopy
import numpy as np

print("\n\n#######################################################################################\n")
print("useage ==> ./ase-supercell.py 'trajactory file' 'supercell integer matrix'\n")
print("           EXAMPLE) ./supercell.py traj 2 2 1\n")
print("#######################################################################################")
if len(sys.argv) is 5:
    print(" ")
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
    sys.exit(1)

trajfile = sys.argv[1]
traj = Trajectory(trajfile, "r")
frames = len(traj)

raw1 = int(sys.argv[2])
raw2 = int(sys.argv[3])
raw3 = int(sys.argv[4])
ncell = raw1*raw2*raw3
matrix = [[raw1,0,0],
          [0,raw2,0],
          [0,0,raw3]]

supertraj = Trajectory("supercell_"+str(raw1)+"x"+str(raw2)+"x"+str(raw3)+"_"+trajfile, "w")

print("I'll make supercell of '"+trajfile+"' file with "+str(frames)+" frames.")


for i in range(frames):
    print("Making supercell of "+str(i+1)+"th frame")
    atoms = traj[i]
    super = make_supercell(atoms, matrix)
    calc = atoms._calc
    if calc != None:
        scalc = deepcopy(calc)
        super.set_calculator(scalc)
        scalc.atoms = super
    
        try:
            energy=calc.results['energy']
        except KeyError:
            print("No energy info")
            energy=None
        except:
            print("Something goes wrong with energy info")
            sys.exit()
        else:
            scalc.results['energy']=energy*ncell
            
        try:
            forces=calc.results['forces']
        except KeyError:
            print("No forces info")
            forces=None
        except:
            print("Something goes wrong with forces info")
            sys.exit()
        else:
            newf = forces
            for j in range(ncell-1):
                newf=np.concatenate((newf, forces))
            scalc.results['forces']=newf
        
        #try:
        #    stress=calc.results['stress']
        #except KeyError:
        #    print("No stress info")
        #    stress=None
        #except:
        #    print("Something goes wrong with stress info")
        #    sys.exit()
        #else:
        #    pass
        
    
    
    
    
        #print(calc.__dict__)
        #print(scalc.__dict__)
    
    supertraj.write(super)

print("\n#############  Primitive cell  #############")
print("Total number of atom = "+str(len(atoms)))
print("Cell =")
print(atoms.get_cell())

print("\n#############  Supercell cell  #############")
print("Transformation matrix, P ="),;print(matrix)
print("")
print("Total number of atom = "+str(len(super)))
print("Cell =")
print(super.get_cell())
print("")
print("pbc ="),;print(super.get_pbc())
print("")

supertraj.close()
print("Supercell trajectory file :: 'supercell_"+str(raw1)+"x"+str(raw2)+"x"+str(raw3)+"_"+trajfile+"'\n")

