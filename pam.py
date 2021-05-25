#!/usr/bin/env python

import numpy as np

def mode_PAM(
    eps,
    ):
    """
    Calculate mode-PAM.
    [l_\sigma (k)]_i = \hbar * [eps_\sigma (k)]_j * M_ijk * [eps_\sigma (k)]_k
       where M_ijk = I_(n*n) \kronecker_product (-i)*\levi_civita_(ijk)

    INPUT
    eps(np.array): Displacement polarization vector of size (len(k)) x (len(\sigma)) x (3n) where n is number of atoms in a unit cell.

    RETURN
    mode_l(np.array): mode resolved PAM of size (len(k)) x (len(\sigma)) x 3
    """
    
    hbar = 0.6582119569 # (eV*fs)
    n = eps.shape[-1] //3
    Mx = np.array([
        [  0,  0,  0],
        [  0,  0,-1j],
        [  0, 1j,  0],
        ])
    My = np.array([
        [  0,  0, 1j],
        [  0,  0,  0],
        [-1j,  0,  0],
        ])
    Mz = np.array([
        [  0,-1j,  0],
        [ 1j,  0,  0],
        [  0,  0,  0],
        ])
    M = np.array([
        np.kron(np.identity(n), Mx),
        np.kron(np.identity(n), My),
        np.kron(np.identity(n), Mz),
        ])
    mode_l = []
    for i in range(eps.shape[0]):
        mode_l.append([])
        for j in range(eps.shape[1]):
            l = np.tensordot(eps[i,j].conj(), M, [0, 1])
            l = np.tensordot(l, eps[i,j], [1, 0])
            mode_l[i].append(l)

    return hbar *np.array(mode_l)

# def dfdT(w, T):
    # """
    # Calculate temperature derivative of B-E distribution.

     # \round f_0     \hbar * w_\simga (k)
    # ------------ = ---------------------- * -----------------------
      # \round T           k_B * T^2           [exp(\hbar * 

    # w: Harmonic phonon frequency (=\omega_\sigma (k))

    # """
