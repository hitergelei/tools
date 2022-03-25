#!/usr/bin/env python

import numpy as np
from ase.io import read, write
from RAG import random_atoms_gen as RAG
from subprocess import call
import datetime

## Hyperparams
backbone = read('wrapped-POSCAR_1_Kooi-16x18x4.vasp')
Iter = 1
# strain_max = 1.15
# strain_min = 0.90
strain_mean  = 1.03
strain_sigma = 0.00
vacuum_max = 0.0

## Iteration
for i in range(Iter):
    # Random ratio
    # strain_ratio = np.random.rand(3) * (strain_max - strain_min) + strain_min
    strain_ratio = np.random.normal(strain_mean, strain_sigma, 3)
    vacuum = [0.,0.,np.random.rand(1) * vacuum_max]

    # Random system generation
    atoms = RAG(
        backbone,
        num_spec_dict = {'Ge':2304, 'Sb':2304, 'Te':5760},
        # fix_ind_dict  = {'Te': [ 0, 2, 5, 7, 8,10,13,15,
                                # 17,19,20,22,25,27,28,30,
                                # 32,34,37,39,40,42,45,47,
                                # 49,51,52,54,57,59,60,62]},
        # cutoff_radi   = None,
        cutoff_frac   = 0.9,
        random_frac   = 0.0,
        # strain        = None,
        strain_ratio  = strain_ratio,
        # vacuum        = None,
        # vacuum_ratio  = None,
        # max_trial_sec = 5,
        # log           = True,
        )
    write('POSCAR{}'.format(i), atoms)

