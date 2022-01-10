#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Shuffle the sequence of an ASE readable atoms list file.
    """)
    # Positional arguments
    parser.add_argument('alist_files', type=str, nargs='+', help='ASE reabable atoms list files. Multiple files can be provided.')
    # Optional arguments
    parser.add_argument('-f', '--file_format', type=str, default=None, help='Atoms list file format can be specified. [Default: auto detection]')

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
    print('Shuffle the sequence of an ASE readable atoms list file.'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()

    # @ Main
    from ase.io import read, write
    for alist_file in args.alist_files:
        alist = read(alist_file, ':', format=args.file_format)
        from random import shuffle
        shuffle(alist)
        write('sffld-{}'.format(alist_file), alist)
        print(' * File sffld-{} has been written.'.format(alist_file))
    print('')
