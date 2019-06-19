#!/usr/bin/env python
from ase import Atoms, Atom
import ase.io
from ase.io import read, write, trajectory
from ase.io.trajectory import Trajectory as Traj
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as Max
from ase.md.velocitydistribution import Stationary
from ase import units
from ase.calculators.vasp import Vasp
from ase.build import make_supercell
import datetime
from ss_util import random_position_generator as RPG
from ase.calculators.lammpsrun import LAMMPS
import os

## Global params
label    = "GeTe-liquid"
atoms    = read('init.vasp')
temp     = 1100 *units.kB
timestep = 10 *units.fs
ttime    = 1000 *units.fs
ptime    = 1000 *units.fs

############# calculator ############
calc = Vasp()
atoms.set_calculator(calc)
calc.read_incar()
from kpoints_gen import get_grid_num, write_KPOINTS
write_KPOINTS(get_grid_num(atoms.get_cell(), 45))
calc.read_kpoints()
    
########### dynamics ###############
Max(atoms, temp* 1.0)
Stationary(atoms)

# from ase.md.verlet import VelocityVerlet
# dyn = VelocityVerlet(atoms, 10*units.fs , trajectory = 'atoms_training_set_'+str(i)+'.traj')  #'dyn' is defined as an object of VV class 
                                                                #'superrhombo_training_set.traj' will be writen
# from ase.md import Langevin
# dyn = Langevin(
    # atoms       = atoms,
    # timestep    = 10 *units.fs,
    # temperature = 1100 *units.kB,
    # friction    = 1e-02,
    # trajectory  = label+'.traj',
    # logfile     = 'log_'+label+'.txt',
    # ) 
from ase.md.npt import NPT
dyn = NPT(
    atoms = atoms,
    timestep = timestep,
    temperature = temp,
    externalstress = 0.,
    ttime = ttime,
    pfactor = (ptime)**2 * 100. *units.GPa,
    trajectory  = label+'.traj',
    logfile     = 'log_'+label+'.txt',
    )
## relax option
dyn.set_fraction_traceless(0.) # 0 --> no shape change but yes volume change
    
############ run ###############
#for i in range(2):        #repeat loop 10 times
#    pot = superrhombo.get_potential_energy()/numofa    
#    kin = superrhombo.get_kinetic_energy()/numofa
#    print('Energy per atom: pot. E. = %f (eV) kin. E. = %f (eV) (Tot. E. = %f (eV))'% (pot, kin, pot + kin))
dyn.run(steps=200000)     #MD simulation of object 'dyn' is performed by 'run' method of VV class
# traj.write(atoms)

