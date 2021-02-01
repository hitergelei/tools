#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Plot group velocity from kappa~~.hdf5 file.
    """)
    # Positional arguments
    parser.add_argument('hdf5_file', type=str, help='HDF5 file.')
    # Optional arguments
    parser.add_argument('-u', '--ymax', type=float, default=None, help='Upper limit for group velocity data.')
    parser.add_argument('-n', '--nbin', type=int, default=100, help='Number of bins for each axis')
    parser.add_argument('-c', '--cmap', type=str, default='jet', help='cmap for KDE plot.')

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
    print('Plot group velocity from kappa~~.hdf5 file.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args = argparse()
    hdf5_file = args.hdf5_file
    nbin      = args.nbin

    import h5py
    h5 = h5py.File(hdf5_file, 'r')

    #@ Main
    # Shape: (irreducible q-points)
    w = np.array(h5['weight'])
    # Shape: (irreducible q-points, bands)
    v_g = np.linalg.norm(h5['group_velocity'], axis=2)
    # Shape: (irreducible q-points, bands)
    f = np.array(h5['frequency'])

    # Flatten
    nband = v_g.shape[1]
    w   = np.repeat(w, nband)
    v_g = v_g.ravel()
    f   = f.ravel()

    #  Unit conversion
    v_g_unit = 100 # (THz*Angstrom) to (m/s)
    v_g *= v_g_unit

    #
    xmax = np.max(f) *1.05
    mask = np.array([True]*len(v_g), dtype=bool)
    if args.ymax:
        ymax = args.ymax
        mask *= v_g <= ymax
    else:
        ymax = np.max(v_g) *1.05
    #
    w   = w[mask]
    v_g = v_g[mask]
    f   = f[mask]

    #
    v_g_rep = np.repeat(v_g, w)
    f_rep = np.repeat(f, w)

    #
    from kde import run_KDE
    KDE, nbin, lims = run_KDE(
        f_rep,
        v_g_rep,
        nbin,
        xmin=0,
        xmax=xmax,
        ymin=0,
        ymax=ymax,
        )

    from matplotlib import pyplot as plt
    plt.pcolormesh(
        KDE[0][:,:nbin], KDE[1][:,:nbin],
        KDE[2][:,:nbin],
        cmap=args.cmap,
        )

    # Cbar
    cb = plt.colorbar()
    cb.ax.set_yticklabels([])
    cb.ax.get_yaxis().labelpad = 15
    cb.set_label(label='Density (arb. unit)', size='x-large', rotation=270)

    #
    plt.scatter(
        f,
        v_g,
        s=1,
        c='w',
        marker='.'
        )

    #
    plt.xlabel('Frequency (THz)', fontsize='x-large')
    plt.ylabel('Group velocity (m/s)', fontsize='x-large')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.subplots_adjust(left=0.19, bottom=0.15, right=0.96, top=0.92)
    plt.xlim(0, xmax)
    plt.ylim(0, ymax)
    plt.show()
