#!/usr/bin/env python

import numpy as np

from ase.io import read, write
backbone = read('cubic-backbone-64.vasp-4x4x4.vasp')

# posi = backbone.get_positions()
# fix_ind = []
# for i in range(len(backbone)):
    # if posi[i][1] < 28.7:
        # fix_ind.append(i)

# from ase.visualize import view
# view(backbone[fix_ind])

# np.save('fix_inds.npy', fix_ind)

from ss_util import random_atoms_gen as RAG
atoms = RAG(
    backbone,
    num_spec_dict = {'Ge': 819, 'Sb': 819, 'Te': 2048, 'V':410},
    fix_ind_dict  = {'Te':list(range(2048, 4096))},
    cutoff_frac   = 0.80,
    random_radi   = 0.0,
    pin_the_fixed = True,
    )

write('cubic-w-rag.vasp', atoms)
