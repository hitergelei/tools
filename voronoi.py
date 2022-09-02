#!/usr/bin/env python
import numpy as np
from subprocess import call

class Voronoi():
    """
    ### This code is based on the Voronoi package ""
    ### Reference: 
    ### Code by Young Jae Choi @ POSTECH, Republic of Korea
    """

    def __init__(
        self,
        atoms_file,
        ):
        self.atoms_file = atoms_file
        self.dynmat = None

    def set_dynmat(self, dynmat):       
        """
        """
        self.dynmat = dynmat

    def read_fix_phonon_log(
        self,
        log_file,
        ):
        with open(log_file) as f:
            lines = f.readlines()
        raw = np.array(lines[-1][:-1].split(), dtype=float)[3:] #*1e-6 /8.61734e-05
        dynmat = raw[::2] + 1.j*raw[1::2]
        N = int(np.sqrt(dynmat.shape[0])/3.)
        # dynmat = np.transpose(dynmat.reshape(N, N, 3, 3), (0,2,1,3)).reshape(3*N, 3*N)
        dynmat = dynmat.reshape(3*N, 3*N)
        self.dynmat = dynmat

    def get_eigen_set(
        self,
        ):
        if self.dynmat is None:
            raise RuntimeError('Dynamical matrix is not set yet.')
        w2, eps = np.linalg.eig(self.dynmat)
        self.eigval = np.sqrt(w2)
        self.eigval_plot = np.real(self.eigval) - np.imag(self.eigval)
        self.eigvec = eps.T

    def plot_VDOS(
        self,
        nbins=200,
        ):
        self.get_eigen_set()
        hist, edges = np.histogram(self.eigval_plot, bins=nbins)
        mids = (edges[0:-1] + edges[1:]) /2.
        from matplotlib import pyplot as plt
        plt.plot(hist, mids, c='k')
        plt.show()
