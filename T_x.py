#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Plot T(x) along the provided direction.
    """)
    # Positional arguments
    parser.add_argument('alist_ase', type=str, help='ASE readable atoms-object list.')
    parser.add_argument('direction', type=float, nargs=3, help='Direction of x variable. x, y, and z components. [Example: 1 0 0]')
    # # Optional arguments
    parser.add_argument('-b', '--nbins', type=int, default=100, help='Number of bins for the function T(x). [Default: 100]')
    parser.add_argument('-n', '--img_slice', type=str, default=':', help='ASE understandable slice [Default: ":" ]')

    return parser.parse_args()

def get_T_x(alist, direc, nbins=100,):
    """
    Return T(x) along the provided direction.
    Temperature of atoms is averaged for those in the same bin.
    x (scalar): position \dot direc

    INPUT
    alist (list): List of ASE readable atoms object.
    direc (array): Direction of x variable. shape=(3).
    nbins (int): Number of bins for the function T(x).
    """
    
    direc /= np.linalg.norm(direc)
    # cell.shape = (len(alist), 3, 3)
    cell = []
    # posi.shape = (len(alist), len(atoms), 3)
    posi  = []
    for atoms in alist:
        cell.append(atoms.get_cell())
        posi.append(atoms.get_positions())
    cell = np.array(cell)
    posi = np.array(posi)

    #
    x_max = np.max(np.tensordot(np.sum(cell, axis=1), direc, axes=(1, 0)))
    delta = x_max /nbins
    centers = np.arange(nbins) *delta +delta/2

    #
    # disp.shape = (len(alist), len(atoms))
    disp = np.tensordot(posi, direc, axes=[2, 0])
    ind = disp //delta

    # T_x.shape = (nbins, len(alist))
    T_x = []
    divider = []
    for i in range(nbins):
        T_x.append([])
        divider.append([])
        for j in range(len(alist)):
            T_x[i].append(alist[j][ind[j]==i].get_temperature())
            divider[i].append(np.sum(ind[j]==i))
            T_x[i][j] *= divider[i][j]
    T_x = np.sum(T_x, axis=1) /np.sum(divider, axis=1)
    return T_x, centers

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
    print('Plot T(x) along the provided direction.'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()

    from ase.io import read
    alist = read(args.alist_ase, args.img_slice)
    if not isinstance(alist, list):
        alist = [alist]
    T_x, centers = get_T_x(alist, args.direction, args.nbins)
    
    # Plot
    from matplotlib import pyplot as plt
    plt.plot(centers, T_x)
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.title(r'x-direction={}'.format(args.direction), fontsize='x-large')
    plt.xlabel(r'x ($\AA$)', fontsize='x-large')
    plt.ylabel(r'Temperature (K)', fontsize='x-large')
    plt.grid(alpha=0.5)
    plt.show()
