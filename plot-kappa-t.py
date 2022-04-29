#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    This code draws a kappa to temperature plot.
    """)
    # Positional arguments
    parser.add_argument('hdf5_file', type=str, help='HDF5 file of phono3py rta/bte calculation.')
    # Optional arguments
    parser.add_argument('-p', '--partial', action='store_true', help='Plot partial LTC (xx, yy, zz, yz, xz, xy).')
    return parser.parse_args()

if __name__ == '__main__':
    ## Intro
    import datetime
    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ POSTECH <<<<<'.center(120))
    print(('Code runtime : '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('Exclude specified images from the ASE readable atoms-list file.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    import h5py
    f = h5py.File(args.hdf5_file, 'r')
    data_list = list(f)
    insts = []
    if 'kappa' in data_list:
        insts.append('kappa')
    if 'kappa_RTA' in data_list:
        insts.append('kappa_RTA')
    T = np.array(f['temperature'])

    data = {}
    for inst in insts:
        data[inst] = np.array(f[inst])

    from matplotlib import pyplot as plt
    for inst in insts:
        plt.figure()
        if args.partial:
            plt.plot(T, data[inst][:,0], label='$\kappa_{xx}$', c='r')
            plt.plot(T, data[inst][:,1], label='$\kappa_{yy}$', c='g')
            plt.plot(T, data[inst][:,2], label='$\kappa_{zz}$', c='b')
            plt.plot(T, data[inst][:,3], label='$\kappa_{yz}$')
            plt.plot(T, data[inst][:,4], label='$\kappa_{xz}$')
            plt.plot(T, data[inst][:,5], label='$\kappa_{xy}$')
        plt.plot(T, np.mean(data[inst][:,0:3], axis=-1), label='$\kappa$', c='k')
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.xlabel('Temperature (K)', fontsize='x-large')
        plt.ylabel('LTC (W/mK)', fontsize='x-large')
        plt.legend(fontsize='large')
        plt.title('{}, {}'.format(args.hdf5_file, inst), fontsize='x-large')
        plt.xlim(np.min(T),np.max(T))
        plt.ylim(0, None)
        plt.subplots_adjust(left=0.20, bottom=0.20, right=0.80, top=0.80)
        plt.grid(alpha=0.4)
    plt.show()
