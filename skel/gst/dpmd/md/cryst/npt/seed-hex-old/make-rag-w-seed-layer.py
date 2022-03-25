#!/usr/bin/env python

import numpy as np

from ase.io import read, write
backbone = read('wrapped-POSCAR_1_Kooi-16x16x4.vasp')

posi = backbone.get_positions()
fix_ind = []
for i in range(len(backbone)):
    if posi[i][0] < 30.00:
        fix_ind.append(i)

# from ase.visualize import view
# view(backbone[fix_ind])

# np.save('fix_inds.npy', fix_ind)

from ss_util import random_atoms_gen as RAG
atoms = RAG(
    backbone,
    fix_ind_dict  = fix_ind,
    cutoff_frac   = 0.80,
    random_radi   = 0.6,
    pin_the_fixed = True,
    )

write('rag_w_kooi_seed_layer-16x16x4.traj', atoms)
