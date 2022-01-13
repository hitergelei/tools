#!/usr/bin/env python
import numpy as np

def wrap_alist(alist, log=False):
    for i in range(len(alist)):
        ## Backup results info
        if hasattr(alist[i]._calc, 'results'):
            results = alist[i]._calc.results.copy()
        else:
            results = None
        alist[i].wrap(eps=0.)
        ## Recover results info
        if results is not None:
            alist[i]._calc.results = results
            alist[i]._calc.atoms = alist[i].copy()
        #### Print every 1000 process
        if log and i % 1000 == 0: 
            print((str(i)+'-th image wrapped').center(120))
    return alist

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Wrap the ASE readable atoms list files.
    """)
    # Positional arguments
    parser.add_argument('alist_files', type=str, nargs='+', help='ASE reabable atoms list files. Multiple files can be provided.')
    # Optional arguments
    parser.add_argument('-f', '--file_format', type=str, default=None, help='Atoms list file format can be specified. [Default: same as input]')
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
    print('Wrap the ASE readable atoms list files.'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()

    # @ Main
    for alist_file in args.alist_files:
        from ase.io import read, write
        alist = read(alist_file, ':')
        if not isinstance(alist, list):
            alist = [alist]
        wrapped_alist = wrap_alist(alist, log=True)
        write('wrapped-'+alist_file, wrapped_alist, format=args.file_format)
        print(' * File wrapped-{} has been written.'.format(alist_file))
    print('')
