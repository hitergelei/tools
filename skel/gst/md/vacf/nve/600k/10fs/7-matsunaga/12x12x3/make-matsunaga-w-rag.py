#!/usr/bin/env python

import numpy as np

from ase.io import read, write
backbone = read('wrapped-POSCAR_1_Kooi-wo-vdw-12x12x3.vasp')

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
    fix_ind_dict  = {'Te':list(range(12*12*3*4, 12*12*3*9))},
    cutoff_frac   = 0.80,
    random_radi   = 0.0,
    pin_the_fixed = True,
    )

write('matsunaga-w-rag-kooi-layer-12x12x3.vasp', atoms)
