#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Get heat capacity from kappa-***.hdf5 file.
    """)
    # Positional arguments
    parser.add_argument('hdf5_file', type=str, help='HDF5 file.')
    parser.add_argument('prim_file', type=str, help='ASE readable structure of primitive unit-cell.')

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
    print('Get heat capacity from kappa-***.hdf5 file.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args = argparse()
    hdf5_file = args.hdf5_file
    prim_file = args.prim_file

    import h5py
    f = h5py.File(hdf5_file, 'r')
    from ase.io import read
    V = read(prim_file).get_volume()
    
    c_q = np.array(f['heat_capacity'])
    w_q = np.array(f['weight'])
    T = np.array(f['temperature'])

    sum_w = np.sum(w_q)
    c = []
    print('NOTE) File: {} '.format(hdf5_file))
    for i in range(len(T)):
        c.append(np.sum(c_q[i] * np.expand_dims(w_q, axis=1)) /sum_w /V)
        print('==> C({}K) = {} (eV / K Angs^3) = {} (J / K m^3)'.format(
            T[i],
            c[i],
            c[i] *1.60217662e-19 *1e30,
            ))

