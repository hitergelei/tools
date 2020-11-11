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

    #
    from ase.io import read
    alist = read(fname, ':')
    if not isinstance(alist, list):
        alist = [alist]
    anum = len(alist[0])

    #
    vol = []
    ke = []
    pe = []
    T = []
    for atoms in alist:
        vol.append(atoms.get_volume())
        ke.append(atoms.get_kinetic_energy())
        pe.append(atoms.get_potential_energy())
        T.append(atoms.get_temperature())
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
    labels = ['Number density ($\AA^{-3}$)', 'Kinetic energy per atom (eV)', 'Potential energy per atom (eV)', 'Total energy per atom (eV)', 'Temperature (K)']
    for i in range(len(props)):
        plt.figure()
        plt.plot(t, props[i], c='k')
        if dt:
            plt.xlabel('Time (ps)', fontsize='x-large')
        else:
            plt.xlabel('Steps', fontsize='x-large')
        plt.ylabel(labels[i], fontsize='x-large')
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.subplots_adjust(left=0.20)
    plt.show()

