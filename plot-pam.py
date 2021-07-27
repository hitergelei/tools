#!/usr/bin/env python

import numpy as np
hbar = 1e-3 *0.6582119569 # (eV*ps)

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
    parser.add_argument('-m', '--mesh', type=int, default=20,
        help='Set k-point mesh at which calculate PAM. Only (mesh, mesh, 1) supported now. Takes integer input. [Default: 20]')
    parser.add_argument('-q', '--qpoint_npy', type=str, default=False, help='Provide q-points manually. --mesh tag is ignored if this tag is activated.')
    parser.add_argument('-x', '--set_x', type=float, nargs=3, default=None, help='Set x-axis of 2D-plot. Insert b1, b2, and b3 for 2pi/a*(b1, b2, b3) vector.')
    parser.add_argument('-t', '--plot_3d', action='store_true', help='If provided, plot 3d texture of PAM.')
    parser.add_argument('-c', '--scale', default=0.01, help='Set PAM vector scale for 2d plot. None is autoscale. [Default: 0.01]')
    parser.add_argument('-z', '--plot_l_z', action='store_true', help='If provided, plot out-of-plane PAM. Works only for 2d plot.')
    parser.add_argument('-r', '--relative', action='store_true', help='If provided, each plot will be re-scaled.')

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
    mode_t = mode_PAM(eps, mode_PAT=True)

    #
    if not args.plot_3d:
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
        if args.set_x:
            x = np.expand_dims(args.set_x, axis=0)
            x = np.matmul(x, recip_latt)[0]
        else:
            x = q_cart[np.argmax(np.linalg.norm(q_cart, axis=1))].copy()
        x /= np.linalg.norm(x)
        y = np.cross(z, x)
        # Rotation matrix
        R = np.linalg.inv([x,y,z]).T
        q_cart_2d = np.matmul(q_cart, R.T)
        mode_l_2d = np.matmul(mode_l, R.T)
        mode_t_2d = np.matmul(mode_t, R.T)

        from matplotlib import pyplot as plt
        fig, ax = plt.subplots()
        ax.scatter(q_cart_2d[:, 0], q_cart_2d[:, 1])
        ax.set_aspect('equal')

    # Plot PAM
    from matplotlib import pyplot as plt
    if args.sigma:
        sigma = [args.sigma-1]
    else:
        sigma = list(range(mode_l.shape[1]))
    for s in sigma:
        if args.plot_3d:
            # PAM
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.quiver(
                q_cart[:, 0], q_cart[:, 1], q_cart[:, 2],
                mode_l[:, s, 0], mode_l[:, s, 1], mode_l[:, s, 2],
                length=0.1,
                )
            from ss_util import axisEqual3D
            axisEqual3D(ax)
            # PAT
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.quiver(
                q_cart[:, 0], q_cart[:, 1], q_cart[:, 2],
                mode_t[:, s, 0], mode_t[:, s, 1], mode_t[:, s, 2],
                length=0.1,
                )
            from ss_util import axisEqual3D
            axisEqual3D(ax)

        elif args.plot_l_z:
            from scipy.interpolate import griddata
            intvl = (np.max(q_cart_2d[:, 0]) - np.min(q_cart_2d[:, 0])) /200.
            gridx, gridy = np.mgrid[
                np.min(q_cart_2d[:, 0]):np.max(q_cart_2d[:, 0]):intvl,
                np.min(q_cart_2d[:, 1]):np.max(q_cart_2d[:, 1]):intvl,
                ]
            gdat = griddata(
                np.transpose([q_cart_2d[:, 0], q_cart_2d[:, 1]]),
                np.real(mode_l_2d[:, s, 2]) /hbar,
                (gridx, gridy),
                method='cubic',
                )
            fig, ax = plt.subplots()
            if args.relative:
                ish = ax.imshow(gdat, cmap='seismic')
            else:
                ish = ax.imshow(gdat, cmap='seismic', vmin=-1., vmax=1.)
            # ax.contourf(
                # gridx, gridy,
                # gdat,
                # 100,
                # cmap='seismic',
                # )
            cbar = fig.colorbar(ish)
            cbar.set_label(r'$l_{\sigma,z}(k) / \hbar$', fontsize='x-large')
            cbar.ax.tick_params(labelsize='x-large')
            ax.set_aspect('equal')
            ax.set_title('PAM, $\sigma$={}, mesh={}X{}'.format(s+1, args.mesh, args.mesh), fontsize='x-large')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel(r'$k_x$', fontsize='x-large')
            ax.set_ylabel(r'$k_y$', fontsize='x-large')
            ax.tick_params(axis="both",direction="in", labelsize='x-large')
            plt.subplots_adjust(right=0.80)

        else:
            # PAM
            fig, ax = plt.subplots()
            ax.quiver(
                q_cart_2d[:, 0], q_cart_2d[:, 1],
                mode_l_2d[:, s, 0], mode_l_2d[:, s, 1],
                # units='xy',
                # linewidths=l_size[:,s] /np.max(l_size[:,s]) *2.,
                # edgecolors='k',
                scale=args.scale,
                # minshaft=0,
                # minlength=0,
                headwidth=5,
                # headlength=1,
                pivot='mid',
                )
            ax.set_aspect('equal')
            ax.set_title('PAM, $\sigma$={}'.format(s+1), fontsize='x-large')
            ax.set_xticks([0])
            ax.set_yticks([0])
            ax.set_xlabel(r'$k_x$', fontsize='x-large')
            ax.set_ylabel(r'$k_y$', fontsize='x-large')
            ax.tick_params(axis="both",direction="in", labelsize='x-large')
            # PAT
            fig, ax = plt.subplots()
            ax.quiver(
                q_cart_2d[:, 0], q_cart_2d[:, 1],
                mode_t_2d[:, s, 0], mode_t_2d[:, s, 1],
                # units='xy',
                # linewidths=l_size[:,s] /np.max(l_size[:,s]) *2.,
                # edgecolors='k',
                scale=args.scale,
                # minshaft=0,
                # minlength=0,
                headwidth=5,
                # headlength=1,
                pivot='mid',
                )
            ax.set_aspect('equal')
            ax.set_title('PAT, $\sigma$={}'.format(s+1), fontsize='x-large')
            ax.set_xticks([0])
            ax.set_yticks([0])
            ax.set_xlabel(r'$k_x$', fontsize='x-large')
            ax.set_ylabel(r'$k_y$', fontsize='x-large')
            ax.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.show()

