#!/usr/bin/env python

from ase import Atoms, Atom
from ase.calculators.lammpsrun import LAMMPS
from ase.io import read, write
from ase.io.trajectory import Trajectory as Traj

## atoms
atoms = read('POSCAR_1_Kooi')

## calc
calc = LAMMPS(
    specorder=['Ge','Sb','Te'],
    parameters={
        'units':'metal',
        'boundary':'p p p',
        'box':'tilt large',
        'pair_style':'deepmd frozen_model.pb',
        'pair_coeff':' ',
        'mass':['1 72.64', '2 121.76', '3 127.60'],
        },
    files=['frozen_model.pb'],
    # keep_tmp_files=True,
    )
atoms.set_calculator(calc)

## rlx
# from ase.optimize.bfgslinesearch import BFGSLineSearch as BFGSL
# from ase.optimize.gpmin.gpmin import GPMin
from ase.optimize.mdmin import MDMin
dyn = MDMin(
    atoms=atoms,
    trajectory='rlx-tmp.traj',
    # force_consistent=False,
    dt=0.1,
    )
dyn.run(fmax=1e-2)
