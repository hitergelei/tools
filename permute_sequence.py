#!/usr/bin/env python

import numpy as np
from copy import deepcopy

def permute_sequence(
    atoms,
    order,
    ):
    new_atoms = atoms[order]
    if atoms._calc is not None:
        new_atoms._calc = deepcopy(atoms._calc)
        new_atoms._calc.atoms = new_atoms.copy()
        new_atoms._calc.results['forces'] = atoms._calc.results['forces'][order]
    return new_atoms

if __name__ == '__main__':
    import sys
    import datetime

    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print("")
    print(">>>>> Code by Young Jae Choi @ POSTECH <<<<<".center(120))
    print(("code started time: "+time).center(120))
    print("")
    print("==================================================================================================".center(120))
    print("")
    print("Useage  ==> ./ase-permute.py >structure file< >order<".center(120))
    print("Example ==> ./ase-permute.py gst.traj 4 5 6 7 1 2 3 4".center(120))
    print("")
    print("==================================================================================================".center(120))
    
    ## Check input
    alist_file = sys.argv[1]
    from ase.io import read, write
    alist = read(alist_file, ':')
    if not isinstance(alist, list):
        alist = [alist]
    if len(sys.argv)-2 != len(alist[0]):
        raise ValueError('The number of arguments is wrong')
    order = np.array(sys.argv[2:], dtype=np.int)

    ## Permute
    new_alist = []
    for i in range(len(alist)):
        if i % 1000 == 999:
            print(('==== Processing {}-th image ===='.format(i+1)).center(120))
        new_atoms = permute_sequence(alist[i], order)
        new_alist.append(new_atoms)
    write('permuted_'+alist_file, new_alist)
    print()
    print(('Saved :: permuted_'+alist_file).center(120))

