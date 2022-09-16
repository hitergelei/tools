#!/usr/bin/env python
import numpy as np

atoms_file = '10ps.vasp'
# fix_pho_log = 'dynmat.log'

#
from vibration_solver import VibSolver
vs = VibSolver(
    atoms_file   = atoms_file,
    calc         = 'lmp',
    displacement = 0.01,
    plus_minus   = True,
    # cp_files     = [],
    )
# vs._produce_structures()
# vs._get_forces()
# vs.get_force_constants()
# print(vs.get_dynamical_matrix())
w2, eps = vs.get_eigen_sets()
vs.plot_VDOS()

#
from voronoi import Voronoi
voro = Voronoi(
    atoms_file,
    )
# voro.read_fix_phonon_log(fix_pho_log)
voro.get_A_mat()
voro.set_eigen_sets(w2, eps)
voro.plot_LT_VDOS()

# voro.plot_VDOS()
