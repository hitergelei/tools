#!/usr/bin/env python

import numpy as np

from ase import units
hbar = 0.6582119569 # (eV*fs)
k_B = units.kB # (eV/K)

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Plot mode phonon angular momentum.
    """)
    # Positional arguments
    parser.add_argument('unitcell', type=str, help='ASE readable unitcell file.')
    parser.add_argument('phono3py_pckl', type=str, help='Phono3py class object saved in pickle format.')
    # # Optional arguments
    parser.add_argument('-t', '--tau', type=float, help='Phonon lifetime for constant lifetime approximation. In fs unit.')
    # # Optional arguments

    return parser.parse_args()

def mode_PAM(
    eps,
    ):
    """
    Calculate mode-PAM.
    [l_\sigma (k)]_i = \hbar * [eps_\sigma (k)]_j * M_ijk * [eps_\sigma (k)]_k
       where M_ijk = I_(n*n) \kronecker_product (-i)*\levi_civita_(ijk)

    INPUT
    eps(np.array): Displacement polarization vector of size (len(k), len(\sigma), 3n) where n is number of atoms in a unit cell.

    RETURN
    mode_l(np.array): mode resolved PAM of size (len(k), len(\sigma), 3)
    """
    
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

    # mode_l.shape == (len(q), len(sigma), 3)
    return hbar *np.array(mode_l)

def f_deriv(w, T):
    """
    Calculate temperature derivative of B-E distribution.

     \round f_0     \beta * \hbar * w         exp(\beta * \hbar * w)
    ------------ = ------------------- * --------------------------------
      \round T              T             [exp(\beta * \hbar * w) - 1]^2

    w (=\omega_\sigma (k)) : Harmonic phonon frequency 
    \beta = (k_B * T)^-1

    RETURN
    dfdT : shape=(len(T), len(k), len(\sigma))
    """

    # beta.shape = (len(T))
    beta = 1. /k_B /T
    # bhw.shape = (len(T), len(k), len(sigma))
    bhw = np.expand_dims(beta, axis=[1,2]) *hbar *np.expand_dims(w, axis=0) *1e-3
    # return shape = (len(T), len(k), len(sigma))
    return bhw /np.expand_dims(T, axis=[1,2]) *np.exp(bhw) /(np.exp(bhw) -1.)**2

def response(eps, w, T, V, tau, v_g):
    """
    Calculate response tensor, \alpha.
                   1    len(k)*3*len(atoms)
    \alpha_ij = - --- * (       sum       ) tau * l_i * (v_g)_j * dfdT
                   V        \k,  \sigma

    l: mode-PAM

    INPUT
    eps : Displacement polarization vector of size (len(k), len(\sigma), 3n) where n is number of atoms in a unit cell.
    w : Harmonic phonon frequencies in Thz, shape=(len(k), len(\sigma)).
    T : Temperature in kelvin, shape=(len(T))
    tau : Lifetime of each phonon mode. shape=(len(T), len(k), len(\sigma))
        or Constant lifetime approximation. (type float)
    v_g : Group velocity vector of each phonon mode. shape=(len(k), len(\sigma), 3)

    RETURN
    alpha : Response tensor. shape=(3,3)
    """

    # l.shape == (len(k), len(sigma), 3)
    l = mode_PAM(eps)
    # dfdT.shape == (len(T), len(k), len(sigma))
    dfdT = f_deriv(w, T)
    # alpha.shape == (len(k), len(sigma), 3, 3)
    alpha = np.matmul(np.expand_dims(l, axis=3), np.expand_dims(v_g, axis=2))
    alpha /= -V
    # alpha.shape == (1, len(k), len(sigma), 3, 3)
    alpha = np.expand_dims(alpha, axis=0)
    if isinstance(tau, float):
        alpha *= tau
    else:
        alpha = alpha * np.expand_dims(tau, axis=[3,4])
    alpha *= np.expand_dims(dfdT, axis=[3,4])
    # return shape = (len(T), 3, 3)
    return np.sum(np.sum(alpha ,axis=1), axis=1)

def get_w_n_eps(phonopy_obj, q):
    freq = []
    eps = []
    for i in range(len(q)):
        f, e = phonopy_obj.get_frequencies_with_eigenvectors(q[i])
        freq.append(f)
        eps.append(e.T)
    # freq.shape == (len(q), # bands)
    freq = np.array(freq)
    # eps.shape == (len(q), # bands, # bands)
    eps = np.array(eps)
    return freq, eps

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

    import pickle as pckl
    tc = pckl.load(open(args.phono3py_pckl, 'rb'))._thermal_conductivity
    mesh = tc._mesh
    # q.shape = (len(q), 3)
    q = np.array(tc._qpoints)
    # gamma.shape = (len(T), len(q), len(sigma))
    gamma = np.array(tc._gamma)[0]
    tau = 1.0 / np.where(gamma > 0, gamma, np.inf) / (2 * 2 * np.pi)
    # v_g.shape = (len(q), len(sigma), 3)
    v_g = np.array(tc._gv)
    # T.shape = (len(T))
    T = np.array(tc._temperatures)
    # w.shape = (len(k), len(sigma))
    w = np.array(tc._frequencies)
    # eps.shape = (len(k), len(sigma), len(sigma))
    eps = np.array(tc._eigenvectors)

    from ase.io import read
    V_uc = read(args.unitcell).get_volume()
    V = V_uc * mesh[0] * mesh[1] * mesh[2]

    # alpha shape=(len(T), 3, 3)
    alpha = response(eps, w, T, V, tau, v_g)
    
    for i in range(len(T)):
        print('T={}(K)'.format(T[i]))
        print(np.real(alpha[i]))
