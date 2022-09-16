#!/usr/bin/env python
import numpy as np

atoms_file = 'init.dump'
# fix_pho_log = 'dynmat.log'

#
from vibration_solver import VibSolver
vs = VibSolver(
    atoms_file   = atoms_file,
    calc         = 'lmp',
    displacement = 0.01,
    # plus_minus   = True,
    # cp_files     = [],
    )
# vs._produce_structures()
# vs._get_forces()
# vs.get_force_constants()
# print(vs.get_dynamical_matrix())
w2, eps = vs.get_eigen_sets()
# vs.plot_VDOS(
    # # nbins=200,
    # # nan_to_zero=True,
    # # plot_imaginary=True,
    # gsmear_std=0.1,
    # )

#
from voronoi import Voronoi
voro = Voronoi(
    atoms_file,
    )
# voro.read_fix_phonon_log(fix_pho_log)
voro.get_A_mat()
voro.set_eigen_sets(w2, eps)
voro.plot_LT_VDOS(
    # nbins=200,
    # nan_to_zero=True,
    # plot_imaginary=True,
    gsmear_std=0.1,
    xlim=(1e-5, 0.09),
    xlim_bp=(1e-5, 0.0115),
    ylim=(0., 23.),
    )

# voro.plot_VDOS()
