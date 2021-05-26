#!/usr/bin/env python

import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Plot mode phonon angular momentum.
    """)
    # Positional arguments
    parser.add_argument('phonopy_pckl', type=str, help='Phonopy class object saved in pickle format.')
    # # Optional arguments
    parser.add_argument('-m', '--mesh', type=int, default=20, help='Set k-point mesh at which calculate PAM. Only (mesh, mesh, mesh) supported now. Takes integer input.')

    return parser.parse_args()

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
    # M.shape == (3, 3n, 3n)
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

    # model_l.shape == (len(q), len(sigma), 3)
    return hbar *np.array(mode_l)

# def dfdT(w, T):
    # """
    # Calculate temperature derivative of B-E distribution.

     # \round f_0     \hbar * w_\simga (k)
    # ------------ = ---------------------- * -----------------------
      # \round T           k_B * T^2           [exp(\hbar * 

    # w: Harmonic phonon frequency (=\omega_\sigma (k))

    # """

def calc_PAM():
    dlfksj

if __name__ == '__main__':
    ## Intro
    import datetime
    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ Phys. Dep. of POSTECH in Korea <<<<<'.center(120))
    print(('Code runtime : '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('Plot mode phonon angular momentum.'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()

    # fname_kappa = 'kappa-m151515.hdf5'
    # import h5py
    # f = h5py.File(fname_kappa, 'r')
    # # irre_q = np.array(f['qpoint'])
    # gamma  = np.array(f['gamma'])
    # vg     = np.array(f['group_velocity'])
    # T      = np.array(f['temperature'])
    # wei    = np.array(f['weight'])

