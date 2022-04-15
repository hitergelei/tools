#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    Exclude specified images from the ASE readable atoms list file.
    """)
    # Positional arguments
    parser.add_argument('inp_file', type=str, help='ASE readable atoms list file name.')
    parser.add_argument('exc_ind', type=int, nargs='+', help='IMPORTANT! Specify in python numbering (Starts from 0).')
    # Optional arguments
    parser.add_argument('-r', '--replace', type=int, help='Fill the vacancy with the provided image of index. [Default: just exclusion and no filling]')
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

    from ase.io import read, write
    alist = read(args.inp_file, ':')
    inds = list(range(len(alist)))

    # Pre-process
    exc_ind = []
    for i in args.exc_ind:
        if i < 0:
            exc_ind.append(i+len(alist))
        else:
            exc_ind.append(i)

    if args.replace < 0:
        replace = i+len(alist)
    else:
        replace = args.replace

    # Sanity check
    if replace in exc_ind:
        raise ValueError('Refilling with excluded argument is meaningless.')
    
    # MAIN
    if args.replace:
        for i in exc_ind:
            inds[i] = args.replace
    else:
        inds = list(set(inds) - set(exc_ind))

    new_alist = []
    for i in inds:
        new_alist.append(alist[i])

    write('exc-'+args.inp_file, new_alist)
