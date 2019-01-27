#!/usr/bin/env python

from ss_phonopy import *
import pickle as pckl

with open('pickle-x441_d0.050_sym.p', 'rb') as f:
    pho1 = pckl.load(f)
with open('../7.7%/pickle-x441_d0.050_sym.p', 'rb') as f:
    pho2 = pckl.load(f)

mesh = [2, 2, 2]
for pho in [pho1, pho2]:
    pho.set_mesh(
        mesh,
        # is_mesh_symmetry=False,
        is_eigenvectors=True,
        )
    pho.set_total_DOS()

two_dos_plot(
    {'pho1':pho1, 'pho2':pho2},
    # color=['b', 'o'],
    ).show()
