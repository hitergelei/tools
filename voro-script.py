#!/usr/bin/env python
import numpy as np

atoms_file = '10ps.vasp'
# fix_pho_log = 'dynmat.log'

#
from vibration_solver import VibSolver
vs = VibSolver(
    atoms_file    = atoms_file,
    calc          = 'lmp',
    displacement  = 0.01,
    # force_cutoff  = 7e-4,
    # plus_minus    = True,
    set_hermitian = False,
    # cp_files    = [],
    )
w2, eps = vs.get_eigen_sets()
vs.plot_VDOS(
    # nbins=200,
    # nan_to_zero=True,
    # plot_imaginary=True,
    gsmear_std=0.02,
    )

# #
# from voronoi import Voronoi
# voro = Voronoi(
    # atoms_file,
    # )
# # voro.read_fix_phonon_log(fix_pho_log)
# voro.get_A_mat()
# voro.set_eigen_sets(w2, eps)
# voro.plot_LT_VDOS(
    # # nbins=200,
    # # nan_to_zero=True,
    # # plot_imaginary=True,
    # gsmear_std=0.02,
    # xlim=(1e-5, 0.37),
    # xlim_bp=(1e-5, 0.55),
    # ylim=(0., 5.5),
    # )
