#!/usr/bin/env python
import numpy

from ase.io import read, write
alist = read('last-10.traj', ':')
new_alist = []

posi = alist[0].get_positions()
mask = posi[:,1] > 23.5
mask *= posi[:,1] < 52.6
for i in range(len(alist)):
    new_alist.append(alist[i][mask])

write('only-half.traj', new_alist)
