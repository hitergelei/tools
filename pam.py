#!/usr/bin/env python

import numpy as np

from ase import units
hbar = 1e-3 *0.6582119569 # (eV*ps)
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
    parser.add_argument('phonopy_pckl', type=str, help='Phonopy class object saved in pickle format.')
    parser.add_argument('phono3py_pckl', type=str, help='Phono3py class object saved in pickle format.')
    # # Optional arguments
    parser.add_argument('-t', '--tau', type=float, default=None, help='Phonon lifetime for constant lifetime approximation. In ps unit.')
    parser.add_argument('-p', '--plot_pam', action='store_true', help='Plot mode-PAM.')
    parser.add_argument('-i', '--tau_histo', action='store_true', help='Plot lifetime histogram.')
    parser.add_argument('-v', '--vg_histo', action='store_true', help='Plot group velocity histogram.')
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
    bhw = np.expand_dims(beta, axis=[1,2]) *hbar *np.expand_dims(w, axis=0)
    # return shape = (len(T), len(k), len(sigma))
    return bhw /np.expand_dims(T, axis=[1,2]) *np.exp(bhw) /(np.exp(bhw) -1.)**2

def response(eps, w, T, V, tau, v_g, band_alpha=False):
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
    band_alpha: Return band-resolved alpha, if provided. (type boolean)

    RETURN
    alpha : Response tensor.
        shape=(            3, 3) for band_alpha=False
        shape=(len(sigma), 3, 3) for band_alpha=True
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
    alpha = alpha * np.expand_dims(dfdT, axis=[3,4])
    if band_alpha:
        # return shape = (len(T), len(sigma), 3, 3)
        return np.sum(alpha ,axis=1)
    else:
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
    tc = pckl.load(open(args.phono3py_pckl, 'rb')).thermal_conductivity
    po = pckl.load(open(args.phonopy_pckl, 'rb'))
    # mesh.shape = (3)
    mesh = tc._mesh
    # q_map.shape = (len(q))
    q_map = tc._grid_mapping_table
    # ir_q.shape = (len(ir_q), 3)
    ir_q = tc.get_qpoints()
    # q.shape = (len(q), 3)
    q = tc._grid_address[:len(q_map)] /mesh
    # T.shape = (len(T))
    T = tc.get_temperatures()

    if args.tau:
        tau = args.tau
    else:
        # ir_gamma.shape = (len(T), len(ir_q), len(sigma))
        ir_gamma = np.array(tc.get_gamma()[0])
        #
        gamma = np.zeros((len(T), len(q), ir_gamma.shape[2]))
        gamma[:,tc.get_grid_points()] = ir_gamma
        gamma = gamma[:, q_map]
        # tau.shape = (len(T), len(q), len(sigma))
        tau = 1. / np.where(gamma > 0, gamma, np.inf) / (2 * 2 * np.pi)

    # v_g = np.array(v_g)
    po.run_qpoints(q, with_eigenvectors=True, with_group_velocities=True)
    # w.shape == (len(q), len(sigma))
    w = po.qpoints.frequencies
    print('Max frequency={}Thz'.format(np.max(w)))
    # eps.shape == (len(q), len(sigma), len(sigma))
    eps = np.transpose(po.qpoints.eigenvectors, [0,2,1])
    # v_g.shape == (len(q), len(sigma), 3)
    v_g = po.qpoints.group_velocities

    from ase.io import read
    V_uc = read(args.unitcell).get_volume()
    V = V_uc * mesh[0] * mesh[1] * mesh[2]

    # band_alpha shape=(len(T), len(sigma), 3, 3)
    band_alpha = response(eps, w, T, V, tau, v_g, band_alpha=True)
    # ( eV ps / A^2 K ) to ( J s / m^2 K )
    scale = units._e * 1e-12 * 1e20 
    band_alpha *= scale

    # alpha shape=(len(T), 3, 3)
    alpha = np.sum(band_alpha, axis=1)
    
    # save
    for i in range(len(T)):
        print('T={}(K)'.format(T[i]))
        print('alpha ( J s / m^2 K ) =')
        print(np.real(np.sum(alpha[i], axis=1)))
        np.save('alpha-tau{}-qx{}{}{}-{}K.npy'.format(args.tau, *mesh, T[i]), band_alpha[i])

    # Only check purpose.
    if args.plot_pam:
        from reciprocal_lattice import get_reciprocal_lattice
        recip_latt = get_reciprocal_lattice(read(args.unitcell).get_cell())
        q_cart = np.matmul(q, recip_latt)

        mode_l = mode_PAM(eps)

        # Plot PAM
        from matplotlib import pyplot as plt
        sigma = list(range(mode_l.shape[1]))
        for s in sigma:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.quiver(
                q_cart[:, 0], q_cart[:, 1], q_cart[:, 2],
                mode_l[:, s, 0], mode_l[:, s, 1], mode_l[:, s, 2],
                length=1e3,
                )
            from ss_util import axisEqual3D
            axisEqual3D(ax)
            ax.set_title('$\sigma$={}'.format(s+1), fontsize='x-large')
            ax.set_xticks([0])
            ax.set_yticks([0])
            ax.set_xlabel(r'$k_x$', fontsize='x-large')
            ax.set_ylabel(r'$k_y$', fontsize='x-large')
            ax.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.show()

    if args.tau_histo:
        from matplotlib import pyplot as plt
        print('tau.shape={}'.format(tau.shape))
        for i in range(len(tau)):
            plt.figure()
            plt.hist(np.reshape(tau[i], -1), bins=np.logspace(np.log10(1e-5),np.log10(1e5), 50))
            plt.title('Lifetime, T={}K'.format(T[i]))
            plt.gca().set_xscale("log")
        plt.show()

    if args.vg_histo:
        from matplotlib import pyplot as plt
        print('v_g.shape={}'.format(v_g.shape))
        plt.hist(np.reshape(np.linalg.norm(v_g *100, axis=-1), -1))
        plt.title('Group velocity')
        plt.show()
