#!/usr/bin/env python
from ase import Atoms, Atom
from ase.calculators.lammpsrun import LAMMPS
from ase import Atoms, Atom
from ase.io import read, write
from ase.io.trajectory import Trajectory as Traj
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as Max 
from ase.md.velocitydistribution import Stationary
from ase.md.verlet import VelocityVerlet
from ase.md import Langevin
from ase import units
from ase.build import make_supercell, bulk
import subprocess as sp


########### cell ##########
label = "sim-to-trn-set"
traj = Traj("GST225-hexa-0-1.traj")
			 
traj_w = Traj(label+"-lmp.traj","w")

for atoms in traj:
    ############## calculator ##################

    calc = LAMMPS(parameters={'units':'metal',
                              'boundary':'p p p',
                              'box':'tilt large',
                              'pair_style':'deepmd frozen_model.pb',
                              'pair_coeff':' ',
                              'mass':['1 121.760', '2 72.64', '3 127.60']},
                  files=['frozen_model.pb'],
                  keep_tmp_files=True,
                  )
    atoms.set_calculator(calc)
    print(calc.__dict__)
    atoms.get_potential_energy()
    atoms.get_forces()
    traj_w.write(atoms)
    sp.call(["kill -9 $(nvidia-smi | sed -n 's/|\s*[0-9]*\s*\([0-9]*\)\s*.*/\\1/p' | sort | uniq | sed '/^$/d')"], shell=True)

