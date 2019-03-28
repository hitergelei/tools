#!/usr/bin/env python

import numpy as np

#### Hyper params
cutoff_force = 5

## Read file
from ase.io import read
alist = read('rlxed-ran.traj', ':')

## Gather max forces of each system
max_list = []
over_index = []
over_force = []
for i in range(len(alist)):
    max_tmp = np.amax(np.linalg.norm(alist[i].get_forces(), axis=-1))
    max_list.append(max_tmp)
    if max_tmp > 5:
        over_index.append(i)
        over_force.append(max_tmp)
for i in range(len(over_index)):
    print(over_index[i], over_force[i])



