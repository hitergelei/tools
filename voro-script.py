#!/usr/bin/env python
import numpy as np

atoms_file = 'init.dump'
# fix_pho_log = 'dynmat.log'

from vibration_solver import VibSolver
vs = VibSolver(
    atoms_file   = atoms_file,
    calc         = 'lmp',
    displacement = 0.01,
    # cp_files     = [],
    )
# vs._produce_structures()
# vs._get_forces()
# vs.get_force_constants()
# print(vs.get_dynamical_matrix())
# vs.get_eigen_sets()
vs.plot_VDOS()

from voronoi import Voronoi
voro = Voronoi(atoms_file)
# voro.read_fix_phonon_log(fix_pho_log)
# voro.plot_VDOS()
