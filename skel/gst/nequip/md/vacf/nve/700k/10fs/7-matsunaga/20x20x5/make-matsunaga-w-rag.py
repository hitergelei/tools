#!/usr/bin/env python

import numpy as np

from ase.io import read, write
backbone = read('POSCAR_1_Kooi-20x20x5.vasp')

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
    fix_ind_dict  = {'Te':list(range(20*20*5*4, 20*20*5*9))},
    cutoff_frac   = 0.80,
    random_radi   = 0.0,
    pin_the_fixed = True,
    )

write('matsunaga-w-rag-kooi-layer-20x20x5.vasp', atoms)
