#!/usr/bin/env python
import numpy as np

## Functions
def to_upper_triangular_cell(cell):
    """
    cell = (3x3 array) Cell parameters. Each raw corresponds a, b, c and each column corresponds to x, y, z.
    new_axis = (3x3 array) New xyz axis in old xyz coordinate. Each raw corresponds x', y', z' and each column corresponds to x, y, z.
    """
    new_axis = np.zeros((3,3))
    new_axis[2] = cell[2]
    new_axis[0] = np.cross(cell[1], cell[2])
    new_axis[1] = np.cross(new_axis[2], new_axis[0])
    for i in range(3):
        new_axis[i] /= np.linalg.norm(new_axis[i])
    return np.matmul(cell, new_axis.T)

def to_lower_triangular_cell(cell):
    """
    cell = (3x3 array) Cell parameters. Each raw corresponds a, b, c and each column corresponds to x, y, z.
    new_axis = (3x3 array) New xyz axis in old xyz coordinate. Each raw corresponds x', y', z' and each column corresponds to x, y, z.
    """
    new_axis = np.zeros((3,3))
    new_axis[0] = cell[0]
    new_axis[1] = np.cross(cell[2], cell[0])
    new_axis[2] = np.cross(new_axis[0], new_axis[1])
    for i in range(3):
        new_axis[i] /= np.linalg.norm(new_axis[i])
    return np.matmul(cell, new_axis.T)

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    Make the cell upper(or lower) triangular cell.
    """)
    # Positional arguments
    parser.add_argument('alist_files', type=str, nargs='+', help='ASE readable files.')
    # Optional arguments
    parser.add_argument('-l', '--lower_tri', action='store_true', help='Make lower triangular cell if specified. [Default: upper triangular]')
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
    print('Make the cell upper(or lower) triangular cell.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    ## args
    alist_files = args.alist_files

    ## Main
    from ase.io import read, write
    from ase.io.formats import filetype
    for i in range(len(alist_files)):
        alist = read(alist_files[i], ':')
        new_alist = []
        if args.lower_tri:
            for atoms in alist:
                atoms.set_cell(to_lower_triangular_cell(atoms.get_cell()), scale_atoms=True)
                new_alist.append(atoms)
            write('lower-tri_{}'.format(alist_files[i]), new_alist, format=filetype(alist_files[i]))
        else:
            for atoms in alist:
                atoms.set_cell(to_upper_triangular_cell(atoms.get_cell()), scale_atoms=True)
                new_alist.append(atoms)
            write('upper-tri_{}'.format(alist_files[i]), new_alist, format=filetype(alist_files[i]))

