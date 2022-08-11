#!/usr/bin/env python

import numpy as np

from ase.io import read, write
backbone = read('cubic-backbone-64.vasp-7x6x6.vasp')

# posi = backbone.get_positions()
# fix_ind = []
# for i in range(len(backbone)):
    # if posi[i][1] < 28.7:
        # fix_ind.append(i)

# from ase.visualize import view
# view(backbone[fix_ind])

# np.save('fix_inds.npy', fix_ind)

from RAG import random_atoms_gen as RAG
atoms = RAG(
    backbone,
    num_spec_dict = {'Ge': 3225, 'Sb': 3225, 'Te': 8064, 'V': 1614},
    fix_ind_dict  = {'Te':list(range(8064, 16128))},
    cutoff_frac   = 0.80,
    random_radi   = 0.0,
    pin_the_fixed = True,
    strain_ratio  = [1.002355, 1.002355, 1.002355],
    )

write('cubic-w-rag.vasp', atoms)
