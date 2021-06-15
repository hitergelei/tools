#!/usr/bin/env python
import numpy as np

from ase.io import read, write
alist = read('../lmp-results-20000:40000-0.traj', '::10')
posi = np.zeros(alist[0].get_positions().shape)
for atoms in alist:
    posi += atoms.get_positions()
posi /= len(alist)
alist[0].set_positions(posi, apply_constraint=False)
write('init.vasp', alist[0])
