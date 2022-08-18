#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    View alist with specific species.
    """)
    # Positional arguments
    parser.add_argument('alist_file', type=str, help='ASE readable atoms list file name.')
    parser.add_argument('spec', type=str, nargs='+', help='Chemical symbols to view.')
    # Optional arguments
    parser.add_argument('-n', '--img_slice', type=str, default=':', help='Image range following python convention. default=":" (e.g.) -n :1000:10')
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
    print(' View alist with specific species.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    from ase.io import read
    alist = read(args.alist_file, args.img_slice)
    if not isinstance(alist, list):
        alist = [alist]

    from ss_util import screen_species
    new_alist = screen_species(alist, args.spec)

    from ase.visualize import view
    view(new_alist)
