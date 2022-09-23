#!/usr/bin/env python
import numpy as np

def _parse_slice(s):
    # Slice format
    if ':' in s:
        a = [int(e) if e.strip() else None for e in s.split(":")]
    # Int format
    else:
        if int(s) == -1:
            a = [-1, None]
        else:
            a = [int(s), int(s)+1]
    return slice(*a)

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    This code converts ase readable file to VASP POSCAR formats.
    """)
    # Positional arguments
    parser.add_argument('inp_file', type=str, help='ASE readable atoms list file name.')
    # Optional arguments
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='Image slice following python convention. default=":" (e.g.) -n :1000:10')
    parser.add_argument('-f', '--file_format', type=str, default=None, help='Specify file format. default=(ase)automatic')
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
    print('This code will read ASE files.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    print('  >>>> Pre-definitions <<<<')
    print('     alist   : Atoms list you gave as an object')
    print('     nlist   : len(alist)')
    print('     atoms   : alist[0]')
    print('     atom    : atoms[0]')
    print('     calc    : atoms._calc')
    print('     results : atoms._calc.results')
    print('     Traj    : ase.io.trajectory:Trajectiory == Traj')
    print('     coord   : atoms.get_positions()')
    print('     frac    : atoms.get_scaled_positions()')
    print('     forces  : atoms.get_forces()')
    print('     energy  : atoms.get_potential_energy()')
    print('     stress  : atoms.get_stress()')
    print('     chem    : atoms.get_chemical_symbols()')
    print('     view    : from ase.visualize import view')
    print('')
    print('============================================================================='.center(120))
    print('')
    args = argparse()

    ##
    from ase.io import read, write
    from ase.io.trajectory import Trajectory as Traj
    from ase.visualize import view

    alist = read(
        args.inp_file,
        index  = args.image_slice,
        format = args.file_format,
        )
    if not isinstance(alist, list):
        alist = [alist]
    try:
        nlist = len(alist)
    except:
        print('           ********* no atoms list in the file')
    try:
        atoms = alist[0]
    except:
        print('           ********* no an atoms object')
    try:
        calc = atoms._calc
    except:
        print('           ********* no a calc object')
    try:
        results = calc.results
    except:
        print('           ********* no a results object')
    try:
        forces = atoms.get_forces()
    except:
        print('           ********* no a force info')
    try:
        energy = atoms.get_potential_energy()
    except:
        print('           ********* no an energy info')
    try:
        stress = atoms.get_stress()
    except:
        print('           ********* no a stress info')
    try:
        coord = atoms.get_positions()
    except:
        print('           ********* no a position info')
    try:
        frac = atoms.get_scaled_positions()
    except:
        print('           ********* no a position info')
    try:
        atom = atoms[0]
    except:
        print('           ********* no an atom exsit')
    try:
        chem = atoms.get_chemical_symbols()
    except:
        print('           ********* no chemical symbols info')
