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
    parser.add_argument('-m', '--mesh', type=int, default=20, help='Set k-point mesh at which calculate PAM. Only (mesh, mesh, 1) supported now. Takes integer input.')
    parser.add_argument('-q', '--qpoint_npy', type=str, default=False, help='Provide q-points manually. --mesh tag is ignored if this tag is activated.')
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

    # # MAIN
    # Load phonon
    import pickle as pckl
    phonon = pckl.load(open(args.phonopy_pckl, 'rb'))
    if args.qpoint_npy:
        q = np.load(args.qpoint_npy)
    else:
        m = args.mesh
        mesh = (m, m, 1)
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
        eps.append(e.T)
    # freq.shape == (len(q), # bands)
    freq = np.array(freq)
    # eps.shape == (len(q), # bands, # bands)
    eps = np.array(eps)

    # Calc
    from pam import mode_PAM
    mode_l = mode_PAM(eps)

    #
    if not args.plot_3d:
        # Set center to be origin
        q_cart -= np.mean(q_cart, axis=0)
        #
        perp = np.zeros(3)
        for i in range(len(q_cart)):
            for j in range(len(q_cart)):
                cro = np.cross(q_cart[i], q_cart[j])
                if np.dot(cro, perp) < 0:
                    perp -= cro
                else:
                    perp += cro
        # Find new xyz-axes
        z = perp /np.linalg.norm(perp)
        x = q_cart[np.argmax(np.linalg.norm(q_cart, axis=1))]
        y = np.cross(z, x)
        y /= np.linalg.norm(y)
        x = np.cross(y, z)
        # Rotation matrix
        R = np.linalg.inv([x,y,z]).T
        q_cart_2d = np.matmul(q_cart, R.T)
        mode_l_2d = np.matmul(mode_l, R.T)

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
                q_cart_2d[:, 0], q_cart_2d[:, 1],
                mode_l_2d[:, s, 0], mode_l_2d[:, s, 1],
                # units='xy',
                # linewidths=l_size[:,s] /np.max(l_size[:,s]) *2.,
                # edgecolors='k',
                # scale=10.0,
                # minshaft=0,
                # minlength=0,
                headwidth=5,
                # headlength=1,
                pivot='mid',
                )
            ax.set_aspect('equal')
        ax.set_title('$\sigma$={}'.format(s+1), fontsize='x-large')
        ax.set_xticks([0])
        ax.set_yticks([0])
        ax.set_xlabel(r'$k_x$', fontsize='x-large')
        ax.set_ylabel(r'$k_y$', fontsize='x-large')
        ax.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.show()

