#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    """)
    # Positional arguments
    parser.add_argument('hdf5_file', type=str, help='HDF5 file.')

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
    print(''.center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args = argparse()
    hdf5_file = args.hdf5_file

    import h5py
    f = h5py.File(hdf5_file, 'r')
