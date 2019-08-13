#!/usr/bin/env python
from ase import Atoms, Atom
import ase.io
from ase.io import read, write, trajectory
from ase.io.trajectory import Trajectory as Traj
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as Max
from ase.md.velocitydistribution import Stationary
from ase import units
from ase.build import make_supercell
import datetime
from ss_util import random_position_generator as RPG
from ase.calculators.lammpsrun import LAMMPS
import os

#### Global params
## cell
label = "lj-Ar-200K-10fs-NPT"
atoms = read('rlx-lj-Ar-2.51K.traj', -1)
for i in range(3):
    for j in range(3):
        if i != j:
            atoms.cell[i,j] = 0.
##
temp   = 5 *units.kB
d_t    = 10 *units.fs
t_step = 5000
############# calculator ############
# from ase.calculators.vasp import Vasp
# calc = Vasp()
# atoms.set_calculator(calc)
# calc.read_incar()
# calc.read_kpoints()
from ase.calculators.lj import LennardJones as LJ
calc = LJ(epsilon=120 *units.kB, sigma=0.34 *units.nm)
atoms.set_calculator(calc)

########### Relaxation #############
# from ase.constraints import StrainFilter as SF
# sf = SF(atoms)
# from ase.optimize.bfgslinesearch import BFGSLineSearch as BFGSLS
# opti_posi = BFGSLS(
    # sf,
    # logfile='rlx-log'+label+'.txt',
    # trajectory = 'rlx-'+label+'.traj',
    # )
# opti_posi.run(1e-3)
    
########### dynamics ###############
Max(atoms, temp * 2)
Stationary(atoms)
# from ase.md.verlet import VelocityVerlet
# dyn = VelocityVerlet(
    # atoms,
    # d_t,
    # trajectory = label+'.traj',
    # logfile = 'log_'+label+'.txt',
    # )
# from ase.md import Langevin
# dyn = Langevin(
    # atoms       = atoms,
    # timestep    = 10 *units.fs,
    # temperature = 100 *units.kB,
    # friction    = 1e-02,
    # trajectory  = label+'.traj',
    # logfile     = 'log_'+label+'.txt',
    # ) 
from ase.md.npt import NPT
dyn = NPT(
    atoms = atoms,
    timestep = d_t,
    temperature = temp,
    externalstress = 0.,
    ttime = 75 * units.fs,
    pfactor = (500. *units.fs)**2 * 100. *units.GPa,
    trajectory  = label+'.traj',
    logfile     = 'log_'+label+'.txt',
    )
# # relax option
dyn.set_fraction_traceless(0) # 0 --> no shape change but yes volume change
dyn.run(steps=t_step)     #MD simulation of object 'dyn' is performed by 'run' method of VV class

