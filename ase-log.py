#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    This code will plot various quentities stored in the ASE readable file.
    """)
    # Positional arguments
    parser.add_argument('fname', type=str, help='ASE readable file.')
    # Optional arguments
    parser.add_argument('-t', '--dt', type=float, default=None, help='Time interval of provide file. (unit of ps) [default: steps unit]')
    parser.add_argument('-n', '--slice', type=str, default=':', help='ASE readable slice. [default: ":"]')
    parser.add_argument('-l', '--large_fig', action='store_true', help='Plot small figure.')

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
    print('This code will plot various quentities stored in the ASE readable file.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args  = argparse()
    fname = args.fname
    dt    = args.dt
    slc   = args.slice

    #
    from ase.io import read
    alist = read(fname, slc)
    if not isinstance(alist, list):
        alist = [alist]
    anum = len(alist[0])

    #
    anum = []
    vol = []
    ke = []
    pe = []
    T = []
    for atoms in alist:
        anum.append(len(atoms))
        vol.append(atoms.get_volume())
        ke.append(atoms.get_kinetic_energy())
        pe.append(atoms.get_potential_energy())
        T.append(atoms.get_temperature())
    anum = np.array(anum)
    den = anum /np.array(vol)
    kepa = np.array(ke) /anum
    pepa = np.array(pe) /anum
    tepa = kepa + pepa
    t = np.arange(len(den), dtype=float)
    if dt:
        t *= dt

    #
    from matplotlib import pyplot as plt
    props = [den, kepa, pepa, tepa, T]
    labels = ['Number density ($\AA^{-3}$)', 'Kinetic energy (eV/atom)', 'Potential energy (eV/atom)', 'Total energy (eV/atom)', 'Temperature (K)']
    for i in range(len(props)):
        plt.figure()
        plt.plot(t, props[i], c='k')
        if dt:
            plt.xlabel('Time (ps)', fontsize='x-large')
        else:
            plt.xlabel('Steps', fontsize='x-large')
        plt.ylabel(labels[i], fontsize='x-large')
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        if args.large_fig:
            plt.subplots_adjust(left=0.20)
        else:
            plt.subplots_adjust(left=0.25, bottom=0.25, right=0.75, top=0.75, wspace=0.2, hspace=0.2)
        plt.grid(alpha=0.4)
    plt.show()

