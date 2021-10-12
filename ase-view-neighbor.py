#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    ASE GUI plot of specified atoms.
    """)
    # Positional arguments
    parser.add_argument('alist_file', type=str, help='ASE readable atoms list file.')
    parser.add_argument('atom_ind', type=int, help='Index of the center atom')
    # # Optional arguments
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='ASE understanable slice in str format. [Default: all]')
    parser.add_argument('-r', '--cutoff_radius', type=float, default=4, help='Cutoff radius for neighbor boolean. [Default: 4 Ang.')

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
    print('ASE GUI plot of specified atoms.'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()
    atom_i = args.atom_ind
    cut_r = args.cutoff_radius

    from ase.io import read
    alist = read(args.alist_file, args.image_slice)
    if not isinstance(alist, list):
        alist = [alist]

    new_alist = []
    for i in range(len(alist)):
        d = alist[i].get_distances(
            atom_i,
            np.arange(len(alist[i])).tolist(),
            mic=True,
            )
        m = d < cut_r
        new_alist.append(alist[i][m])
        print('\nimg #{} - {}'.format(i, np.arange(len(alist[i]))[m].tolist()))

    from ase.visualize import view
    view(new_alist)
