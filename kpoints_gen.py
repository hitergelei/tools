#!/usr/bin/env python

import numpy as np

def get_grid_num(cell_mat, precision):
    """
    precision ~ Real space cell parameter * Reciprocal lattice cell parameter
    """
    # Dtype
    precision = np.float(precision)
    if np.shape(cell_mat) == (3,3):
        cell_mat = np.linalg.norm(cell_mat, axis=-1)
    elif np.shape(cell_mat) == (3,):
        cell_mat = np.array(cell_mat)
    else:
        raise ValueError('Cell matrix shape[={}] is wrong. (3,3) or (3,) should be given.'.format(np.shape(cell_mat)))
    # Get grid numbers
    k_grids = np.around(precision / cell_mat).astype(np.int)
    return k_grids

def write_KPOINTS(k_grids):
    with open('KPOINTS', 'w') as txt:
        txt.write('KPOINTS\n0\nGamma\n{} {} {}\n0 0 0'.format(k_grids[0], k_grids[1], k_grids[2]))
        
if __name__ == '__main__':
    import sys
    import datetime

    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ POSTECH <<<<<'.center(120))
    print(('code started time: '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('')
    print('Useage  ==> ./kpoints_gen.py >One structure file< >Precision<'.center(120))
    print('Example ==> ./kpoints_gen.py gst_kooi.traj 55'.center(120))
    print('Note) Precision ~ Real space cell parameter * Reciprocal lattice cell parameter'.center(120))
    print('Note) Gamma centered mode use only'.center(120))
    print('')
    print('=================================================================================================='.center(120))

    ## Input data
    from ase.io import read
    cell_mat = read(sys.argv[1], -1).get_cell()
    prec = sys.argv[2]
    k_grids = get_grid_num(cell_mat, prec)

    ## Write KPOINTS
    write_KPOINTS(k_grids)
