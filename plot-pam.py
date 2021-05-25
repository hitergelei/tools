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
    parser.add_argument('-s', '--sigma', type=int, default=None, help='Phonon band index to plot. [Default: Plot all]')
    parser.add_argument('-m', '--mesh', type=int, default=40, help='Set k-point mesh at which calculate PAM. Only (mesh, mesh, 1) supported now. Takes integer input.')
    parser.add_argument('-t', '--plot_3d', action='store_true', help='If provided, plot 3d texture of PAM.')

    return parser.parse_args()

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
    fname_phonon = args.phonopy_pckl
    m = args.mesh

    # Params
    # fname_kappa = 'kappa-m151515.hdf5'
    mesh = (m, m, 1)

    # MAIN
    # # q = 
    # # Load qpoints
    # import h5py
    # f = h5py.File(fname_kappa, 'r')
    # # irre_q = np.array(f['qpoint'])
    # gamma  = np.array(f['gamma'])
    # vg     = np.array(f['group_velocity'])
    # T      = np.array(f['temperature'])
    # wei    = np.array(f['weight'])

    # Load phonon
    import pickle as pckl
    phonon = pckl.load(open(fname_phonon, 'rb'))
    phonon.run_mesh(
        mesh,
        is_time_reversal=False,
        is_mesh_symmetry=False,
        with_eigenvectors=True,
        is_gamma_center=True,
        )
    q = phonon.get_mesh_dict()['qpoints']
    from reciprocal_lattice import get_reciprocal_lattice
    recip_latt = get_reciprocal_lattice(phonon.unitcell.cell)
    q_cart = np.matmul(q, recip_latt)

    # 
    freq = []
    eps = []
    for i in range(len(q)):
        f, e = phonon.get_frequencies_with_eigenvectors(q[i])
        freq.append(f)
        eps.append(e)
    freq = np.array(freq)
    eps = np.array(eps)

    # Calc
    from pam import mode_PAM
    mode_l = mode_PAM(eps)

    # Plot PAM
    from matplotlib import pyplot as plt
    if args.sigma:
        sigma = [args.sigma-1]
    else:
        sigma = list(range(mode_l.shape[1]))
    for s in sigma:
        if args.plot_3d:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.quiver(
                q_cart[:, 0], q_cart[:, 1], q_cart[:, 2],
                mode_l[:, s, 0], mode_l[:, s, 1], mode_l[:, s, 2],
                length=0.1,
                )
            from ss_util import axisEqual3D
            axisEqual3D(ax)
        else:
            fig, ax = plt.subplots()
            ax.quiver(
                q_cart[:, 0], q_cart[:, 1],
                mode_l[:, s, 0], mode_l[:, s, 1],
                # scale=10.0,
                )
            ax.set_aspect('equal')
        ax.set_title('$\sigma$={}'.format(s), fontsize='x-large')
    plt.show()

