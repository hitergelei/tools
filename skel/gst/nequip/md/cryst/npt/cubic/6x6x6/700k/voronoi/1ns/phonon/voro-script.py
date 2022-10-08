#!/usr/bin/env python
import numpy as np

atoms_file = '1ns.vasp'
# fix_pho_log = 'dynmat.log'

#
from vibration_solver import VibSolver
vs = VibSolver(
    atoms_file    = atoms_file,
    calc          = 'lmp',
    displacement  = 0.01,
    # force_cutoff  = 7e-4,
    # plus_minus    = True,
    # set_hermitian = True,
    # cp_files    = [],
    )
# dm = vs.get_dynamical_matrix()
# np.save('_dm.npy', dm)

w2, eps = vs.get_eigen_sets()
# vs.plot_VDOS(
    # nbins=200,
    # # nan_to_zero=True,
    # # plot_imaginary=True,
    # gsmear_std=0.05,
    # xlim=(1e-5, 0.38),
    # xlim_bp=(1e-5, 0.199),
    # ylim=(0., 5.5),
    # )

#
from voronoi import Voronoi
voro = Voronoi(
    atoms_file,
    )
# voro.read_fix_phonon_log(fix_pho_log)
voro.get_A_mat()
voro.set_eigen_sets(w2, eps)
# voro.plot_LT_VDOS(
    # # nbins=200,
    # # nan_to_zero=True,
    # # plot_imaginary=True,
    # gsmear_std=0.05,
    # xlim=(1e-5, 0.38),
    # xlim_bp=(1e-5, 0.199),
    # ylim=(0., 5.5),
    # )
voro.get_atomic_VDOS_and_AM(
    freq_range = (0.375, 0.775),
    # reduced=True,
    # show_2d=True,
    )
