#!/usr/bin/env python

import numpy as np

from ase import units
hbar = 1e-3 *0.6582119569 # (eV*ps)
k_B = units.kB # (eV/K)
Mx_T = np.array([
    [  0,  0,  0],
    [  0,  0,  1],
    [  0,  1,  0],
    ])
My_T = np.array([
    [  0,  0,  1],
    [  0,  0,  0],
    [  1,  0,  0],
    ])
Mz_T = np.array([
    [  0,  1,  0],
    [  1,  0,  0],
    [  0,  0,  0],
    ])
M_AT = np.array([Mx_T, My_T, Mz_T])

Mx_M = np.array([
    [  0,  0,  0],
    [  0,  0,-1j],
    [  0, 1j,  0],
    ])
My_M = np.array([
    [  0,  0, 1j],
    [  0,  0,  0],
    [-1j,  0,  0],
    ])
Mz_M = np.array([
    [  0,-1j,  0],
    [ 1j,  0,  0],
    [  0,  0,  0],
    ])
M_AM = np.array([Mx_M, My_M, Mz_M])

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
    parser.add_argument('-t', '--tau', default=None, help='Phonon lifetime for constant lifetime approximation. In ps unit. "auto" will set the constant as the average value. [Default: no approximation].')
    parser.add_argument('-p', '--plot_pam', action='store_true', help='Plot mode-PAM.')
    parser.add_argument('-i', '--tau_histo', action='store_true', help='Plot lifetime histogram.')
    parser.add_argument('-v', '--vg_histo', action='store_true', help='Plot group velocity histogram.')
    parser.add_argument('--set_l_0_zero', action='store_true', help='Set PAM of 1st band as zero. Use it very carefully.')
    parser.add_argument('--temperature', type=float, nargs='+',
        help='Set temperature manually. Only works for constant lifetime approximation. Multiple temperature can be set.')
    parser.add_argument('--dim2', action='store_true', help='2-D PAM per area (not per volume) calc.')
    parser.add_argument('--plot_response_ratio_to_direc', type=str, default=None,
        help='Plot f_1 / f_0 for all phonons. T-perturbation is 1% and L=100 micrometer. Provide direction of T-gradient')

    return parser.parse_args()

def AAM(
    eps,
    AAT=False,
    ):
    r"""
    Calculate atomic angular momentum. 

    INPUT
    eps(np.array): Displacement polarization vector of size (len(q), len(\sigma), 3n) where n is number of atoms in a unit cell.

    RETURN
    l_{\sigma, i, \alpha}(q) (np.array): mode resolved AAM of size (len(q), len(\sigma), len(atoms), 3)
    """
    
    n = eps.shape[-1] //3
    #eps.shape = (len(k), len(\sigma), n, 3)
    eps = eps.reshape(eps.shape[0], eps.shape[1], n, 3)
    if AAT:
        # M.shape == (3, 3, 3)
        M = M_AT
    else:
        M = M_AM
    mode_l = []
    for i in range(3):
        l = np.expand_dims(eps.conj(), axis=3) @ np.expand_dims(M[i], axis=(0, 1, 2)) @ np.expand_dims(eps, axis=4)
        mode_l.append(l.reshape(eps.shape[0], eps.shape[1], n, 1))
    mode_l = np.concatenate(mode_l, axis=3)

    # mode_l.shape == (len(q), len(\sigma), n, 3)
    return hbar *mode_l

def mode_PAM(
    eps,
    mode_PAT=False,
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
    if mode_PAT:
        # M.shape == (3, 3n, 3n)
        M = np.array([
            np.kron(np.identity(n), Mx_T),
            np.kron(np.identity(n), My_T),
            np.kron(np.identity(n), Mz_T),
            ])
    else:
        # M.shape == (3, 3n, 3n)
        M = np.array([
            np.kron(np.identity(n), Mx_M),
            np.kron(np.identity(n), My_M),
            np.kron(np.identity(n), Mz_M),
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

def f_zero(w, T):
    """
    Calculate B-E distribution.
    f_0 = [exp(\beta * \hbar * w) - 1]**(-1)

    w (=\omega_\sigma (k)) : Harmonic phonon frequency [w.shape = (len(k), len(\sigma))]
    \beta = (k_B * T)^-1

    RETURN
    f: shape=(len(T), len(k), len(\sigma))
    """
    # beta.shape = (len(T))
    beta = 1. /k_B /T
    # bhw.shape = (len(T), len(k), len(sigma))
    bhw = np.expand_dims(beta, axis=[1,2]) *hbar *np.expand_dims(w, axis=0)
    # return shape = (len(T), len(k), len(sigma))
    return 1. /(np.exp(bhw) -1.)

def f_deriv(w, T):
    """
    Calculate temperature derivative of B-E distribution.

     \round f_0     \beta * \hbar * w         exp(\beta * \hbar * w)
    ------------ = ------------------- * --------------------------------
      \round T              T             [exp(\beta * \hbar * w) - 1]^2

    w (=\omega_\sigma (k)) : Harmonic phonon frequency [w.shape = (len(k), len(\sigma))]
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

def f_response(w, T, tau, v_g, direc, L=1e6, dT_over_T=0.01):
    """
    Calculate response of B-E distribution.
    f_response = -tau * (v_g)_i * dfdT * (dTdx_i)

    INPUT
    w : Harmonic phonon frequencies in Thz, shape=(len(k), len(\sigma)).
    T : Temperature in kelvin, shape=(len(T))
    tau : Lifetime of each phonon mode. shape=(len(T), len(k), len(\sigma))
        or Constant lifetime approximation. (type float)
    v_g : Group velocity vector of each phonon mode. shape=(len(k), len(\sigma), 3)
    direc: ('x', 'y', or 'z') (type str)
    L : Size of sample in Angstrom unit. (Default: 100 micrometer = 1e6 Angstrom, type float)
    dT_over_T: Delta_T / T (Perturbation strength, type float)

    RETURN
    f_1 : Response value of B-E distribution. f_1.shape = (len(T), len(k), len(sigma))
    """
    if direc == 'x':
        v_g_i = v_g[:,:,0]
    elif direc == 'y':
        v_g_i = v_g[:,:,1]
    elif direc == 'z':
        v_g_i = v_g[:,:,2]
    else:
        raise RuntimeError('Wrong argument for "direc" parameter.')
    # dTdx_i.shape = (len(T))
    dTdx_i = T * dT_over_T / L
    # dfdT.shape = (len(T), len(k), len(sigma))
    dfdT = f_deriv(w, T)
    # RETURN.shape = (len(T), len(k), len(sigma))
    return -1 * tau * np.expand_dims(v_g_i, axis=0) * dfdT * np.expand_dims(dTdx_i, axis=[1,2])

def response(eps, w, T, V, tau, v_g, band_alpha=False, set_l_0_zero=False):
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
    if set_l_0_zero:
        l[:,0] = 0.
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
    if args.temperature:
        T = np.array(args.temperature, dtype=float)

    if args.tau == 'auto' or args.tau == None:
        # ir_gamma.shape = (len(T), len(ir_q), len(sigma))
        ir_gamma = np.array(tc.get_gamma()[0])
        #
        gamma = np.zeros((len(T), len(q), ir_gamma.shape[2]))
        gamma[:,tc.get_grid_points()] = ir_gamma
        gamma = gamma[:, q_map]
        # tau.shape = (len(T), len(q), len(sigma))
        tau = 1. / np.where(gamma > 0, gamma, np.inf) / (2 * 2 * np.pi)
        if args.tau == 'auto':
            tau = np.mean(np.mean(tau, axis=1), axis=1)
            tau = np.expand_dims(tau, [1,2])
        const_tau = 1.
    else:
        tau = 1e12
        const_tau = float(args.tau) *1e-12

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
    if args.dim2:
        a1, a2 = read(args.unitcell).get_cell()[:2][:,:2]
        V_uc = np.linalg.norm(np.cross(a1, a2))
        V = V_uc * mesh[0] * mesh[1]
    else:
        V_uc = read(args.unitcell).get_volume()
        V = V_uc * mesh[0] * mesh[1] * mesh[2]

    if args.plot_response_ratio_to_direc:
        f_1 = f_response(w, T, tau, v_g, direc=args.plot_response_ratio_to_direc)
        f_0 = f_zero(w, T)
        from matplotlib import pyplot as plt
        max_ratio = []
        for i in range(len(T)):
            # plt.figure()
            # ratio = np.sort(np.reshape(f_1[i] / f_0[i], -1))
            # plt.plot(ratio[:len(ratio)//2:10])
            # plt.title('{}K'.format(T[i]), fontsize='x-large')
            max_ratio.append(np.amax(np.reshape(f_1[i] / f_0[i], -1)))
        plt.plot(T, max_ratio, c='k')
        plt.yscale('log')
        plt.axhline(y=1e-1, ls='--', label='0.1', c='k')
        plt.legend(fontsize='large').set_draggable(True)
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.xlabel('Temperature (K)', fontsize='x-large')
        plt.ylabel(r'max($f_1$/$f_0$)', fontsize='x-large')
        plt.subplots_adjust(left=0.20, bottom=0.20, right=0.80, top=0.80)
        plt.grid(alpha=0.2)
        plt.show()
        exit()

    # band_alpha shape=(len(T), len(sigma), 3, 3)
    band_alpha = response(eps, w, T, V, tau, v_g, band_alpha=True, set_l_0_zero=args.set_l_0_zero)
    if args.dim2:
        # ( eV ps / A K ) to ( J s / m K )
        scale = units._e * 1e-12 * 1e10 
    else:
        # ( eV ps / A^2 K ) to ( J s / m^2 K )
        scale = units._e * 1e-12 * 1e20 
    band_alpha *= scale

    # alpha shape=(len(T), 3, 3)
    alpha = np.sum(band_alpha, axis=1)
    
    # save
    for i in range(len(T)):
        print('T={}(K)'.format(T[i]))
        if args.dim2:
            print('alpha ( J s / m K ) =')
            np.save('alpha-2d-tau{}-qx{}{}{}-{}K.npy'.format(args.tau, *mesh, T[i]), band_alpha[i])
        else:
            print('alpha ( J s / m^2 K ) =')
            np.save('alpha-tau{}-qx{}{}{}-{}K.npy'.format(args.tau, *mesh, T[i]), band_alpha[i])
        print(np.real(alpha[i] * const_tau))

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
