#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    This code will give you the volume and density histogram.
    """)
    # Positional arguments
    parser.add_argument('alist_file', type=str, help='ASE readable atoms list file name.')
    # Optional arguments
    parser.add_argument('-n', '--ase_slice', type=str, default=':', help='ASE readable format of slice for atoms list file. [default: ":"]')
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
    print('This code will give you the volume and density histogram.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    from ase.io import read
    alist = read(args.alist_file, args.ase_slice)

    vol_arr = []
    n_arr = []
    for atoms in alist:
        vol_arr.append(np.linalg.det(atoms.get_cell()))
        n_arr.append(len(atoms))
    vol_arr = np.array(vol_arr)
    n_arr = np.array(n_arr)

    # Post process
    volperatom_arr = vol_arr/n_arr
    density_arr = 1./volperatom_arr
    d_average = np.mean(density_arr)
    d_std = np.std(density_arr)
    v_average = np.mean(volperatom_arr)
    v_std = np.std(volperatom_arr)
        
    ## plot
    from matplotlib import pyplot as plt
    font = {'family':'Arial'}
    plt.rc('font', **font)
    # Density
    n, bins, patches = plt.hist(density_arr, bins=100, facecolor='gray', alpha=0.70)
    max_height = np.sort(n)[-10]
    plt.title('Density Histogram (atoms/Ang^3)', fontsize='x-large')
    plt.xlabel('%d images, average = %.4e, sigma = %.4e' % (len(density_arr), d_average, d_std), fontsize='x-large')
    plt.ylabel('population', fontsize='x-large')
    plt.barh(max_height/5, d_std, height=max_height/50, left=d_average, color='black')
    plt.axvline(d_average, c='k')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.grid(alpha=0.4)
    # Volume
    plt.figure()
    n, bins, patches = plt.hist(volperatom_arr, bins=100, facecolor='gray', alpha=0.70)
    max_height = np.sort(n)[-10]
    plt.title('Volume Histogram (Ang^3/atom)', fontsize='x-large')
    plt.xlabel('%d images, average = %.4e, sigma = %.4e' % (len(volperatom_arr), v_average, v_std), fontsize='x-large')
    plt.ylabel('population', fontsize='x-large')
    plt.barh(max_height/5, v_std, height=max_height/50, left=v_average, color='black')
    plt.axvline(v_average, c='k')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.grid(alpha=0.4)
    plt.show()
