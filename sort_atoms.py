#!/usr/bin/env python
import numpy as np

def argsort_atoms_by_chem(
    atoms,
    chem_order=None,
    ):
    chems = np.array(atoms.get_chemical_symbols())
    if not chem_order:
        chem_order = np.unique(chems)
    argsort = []
    for chem in chem_order:
        argsort.extend(np.arange(len(atoms))[chems == chem])
    return argsort

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Sort Atoms object of ASE readable structures.
    If multiple Atoms objects are provided, order is determined only by the first.
    """)
    # Positional arguments
    parser.add_argument('alist_file', type=str, help='ASE readable alist file.')
    # # Optional arguments
    parser.add_argument('-c', '--chem_order', type=str, nargs='+', default=None, help='Set order of chemical symbols. E.g. Ge Te')

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
    print('Sort Atoms object of ASE readable structures.'.center(120))
    print('If multiple Atoms objects are provided, order is determined only by the first.'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()
    alist_file = args.alist_file
    chem_order = args.chem_order

    # @ Main
    from ase.io import read, write
    alist = read(alist_file)
    if not isinstance(alist, list):
        alist = [alist]
    argsort = argsort_atoms_by_chem(alist[0], chem_order)

    new_alist = []
    for i in range(len(alist)):
        new_alist.append(alist[i][argsort])

    write('sorted-{}'.format(alist_file), new_alist)
