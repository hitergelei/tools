#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Plot cumulative LTC from kappa~~.hdf5 file.
    """)
    # Positional arguments
    parser.add_argument('hdf5_file', type=str, help='HDF5 file.')
    parser.add_argument('unitcell', type=str, help='Unit-cell structure file of Phonopy readable.')
    # # Optional arguments
    parser.add_argument('-x', '--xmax', type=float, default=None, help='Upper limit for x-axis.')
    parser.add_argument('-y', '--ymax', type=float, default=None, help='Upper limit for y-axis.')
    parser.add_argument('-n', '--nbin', type=int, default=100, help='Number of bins for frequency.')
    parser.add_argument('-p', '--partial', action='store_true', help='Plot every direction (x,y,z).')
    # parser.add_argument('-c', '--cmap', type=str, default='jet', help='cmap for KDE plot.')

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
    print('Plot cumulative LTC from kappa~~.hdf5 file.'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()
    hdf5_file = args.hdf5_file
    afile     = args.unitcell
    nbin      = args.nbin

    # Load temperature
    import h5py
    f = h5py.File(hdf5_file, 'r')
    T = np.array(f['temperature'])
    # Calc
    from subprocess import call
    call('phono3py-kaccum {} -c {} --nsp {} > kaccum.dat'.format(hdf5_file, afile, nbin), shell=True)
    data = np.loadtxt('kaccum.dat').reshape(len(T), nbin, 13)
    #
    f = data[:,:,0]
    ka = data[:,:,1:4]
    dka_df = data[:,:,7:10]

    # Plot
    from matplotlib import pyplot as plt
    if args.partial:
        labels = ['x', 'y', 'z']
        colors = ['r', 'g', 'b']
    for j in range(len(T)):
        plt.figure()
        if args.partial:
            for i in range(3):
                plt.plot(
                    f[j], ka[j,:,i],
                    label='$\kappa_{{{}}}(f)$'.format(labels[i]),
                    c=colors[i],
                    ls='--',
                    )
                plt.plot(
                    f[j], dka_df[j,:,i],
                    label='$d\kappa_{{{}}}(f)$/$df$'.format(labels[i]),
                    c=colors[i],
                    )
        else:
            plt.plot(
                f[j], np.mean(ka[j], axis=1),
                label='$\kappa(f)$',
                c='k',
                ls='--',
                )
            plt.plot(
                f[j], np.mean(dka_df[j], axis=1),
                label='$d\kappa(f)$/$df$',
                c='k',
                )

        #
        plt.xlabel('Frequency (THz)', fontsize='x-large')
        plt.ylabel('Thermal conductivity (W/mK)', fontsize='x-large')
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        # plt.subplots_adjust(left=0.19, bottom=0.15, right=0.96, top=0.92)
        plt.legend(fontsize='large')
        plt.title('{}K'.format(T[j]), fontsize='x-large')
        plt.grid(alpha=0.5)
        plt.xlim(0, args.xmax)
        plt.ylim(0, args.ymax)
    plt.show()
