#!/usr/bin/env python
import numpy

from ase.io import read, write
alist = read('lmp-results-:-0.traj', ':')
new_alist = []

posi = alist[0].get_positions()
mask = posi[:,1] > 23.2
mask *= posi[:,1] < 48.7
for i in range(len(alist)):
    new_alist.append(alist[i][mask])

write('only-half.traj', new_alist)
