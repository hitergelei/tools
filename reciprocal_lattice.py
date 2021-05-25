#!/usr/bin/env python
import numpy as np

def get_reciprocal_lattice(lattice_matrix):
    """
    INPUT
    lattice_matrix ((3,3) np.array): = (a1, a2, a3)
    RETURN
    reciprocal_lattice ((3,3) np.array): = (b1, b2, b3)
    """
    return 2 *np.pi *np.linalg.inv(np.array(lattice_matrix)).T
