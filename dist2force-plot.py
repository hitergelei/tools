#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    This code will plot distance to force magnitude plot.
    """)
    # Positional arguments
    parser.add_argument('trajfile', type=str, help='ASE readable file.')
    # parser.add_argument('dt', type=float, help='Time step of J0Jt.dat file. Becareful about unit. (Real: fs, metal: ps)')
    # Optional arguments
    parser.add_argument('-n', '--img_slice', type=str, default=':', help='Image range following python convention. default=":" (e.g.) -n :1000:10')
    # parser.add_argument('-u', '--lammps_unit', type=str, default='metal', help='Set unit of J0Jt.dat file between metal and real. [default: metal]')
    # parser.add_argument('-l', '--dont_load', dest='load_bool', action='store_false', help="If provided, don't load the data. [default: load]")
    # parser.add_argument('-s', '--dont_save', dest='save_bool', action='store_false', help="If provided, don't save the data. [default: save]")
    # parser.add_argument('-p', '--dont_plot', dest='plot_bool', action='store_false', help="If provided, don't plot the data. [default: plot]")

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
    print('This code will plot distance to force magnitude plot.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args = argparse()

    from ase.io import read
    alist = read(args.trajfile, args.img_slice)
    if not isinstance(alist, list):
        alist = [alist]

    d = []
    f = []
    for atoms in alist:
        i = np.argsort(atoms.get_all_distances(mic=True)[0])
        d.append(atoms.get_all_distances(mic=True)[0][i])
        f.append(np.linalg.norm(atoms.get_forces(), axis=-1)[i])
    d = np.concatenate(d, axis=0)
    f = np.concatenate(f, axis=0)

    from matplotlib import pyplot as plt
    plt.yscale('log')
    plt.scatter(d, f, c='r')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.ylabel('Force, |$F_i$| (eV/$\AA$)', fontsize='x-large')
    plt.xlabel('Distance, $r_i-r_0$ ($\AA$)', fontsize='x-large')
    plt.grid(alpha=0.5)
    plt.show()
