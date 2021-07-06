#!/usr/bin/env python
import numpy as np

uc_file = 'gete-alpha-prim.vasp'

def rot_mat(mat, new_z):
    # Rotate mat's z-axis to (111)
    new_x = np.array([1.,0.,0.])
    new_y = np.cross(new_z, new_x)
    new_y /= np.linalg.norm(new_y)
    new_x = np.cross(new_y, new_z)
    # R.shape = (1, 3, 3)
    R = np.expand_dims(np.array([new_x, new_y, new_z]).T, axis=0)
    mat = np.matmul(np.transpose(R, [0,2,1]), mat)
    mat = np.matmul(mat, R)
    return mat

from subprocess import call, check_output
out = str(check_output('phonopy-vasp-born', shell=True)).split('\\n')[1:-1]
for i in range(len(out)):
    out[i] = out[i].split()
out = np.array(out, dtype=float).reshape(-1,3,3)

#
from ase.io import read
uc = read(uc_file)
z = np.sum(uc.get_cell(), axis=0)
z /= np.linalg.norm(z)

#
out = rot_mat(out, z)

eps = out[0]
born = out[1:]

print('eps =\n{}\nZ* /e =\n{}'.format(eps, born))
