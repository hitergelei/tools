#!/usr/bin/env python
import numpy as np

def unfold_positions(alist, ref_atoms=None):
    """
    This code will unfold the wrapped ASE atoms objects in PBC.
    alist: List of ase atoms object.
    ref_atoms: Unfolding refers this structure. 
    """

    # Extract position info
    r_list = []
    for i in range(len(alist)):
        r_list.append(alist[i].get_scaled_positions())
    r_list = np.array(r_list)

    # Unfold positions
    R = np.zeros(r_list.shape, int)
    R[1:] = np.add.accumulate(
        np.rint(np.array(r_list[:-1]) - np.array(r_list[1:])),
        )
    if ref_atoms:
        R = R + np.rint(ref_atoms.get_scaled_positions() - np.array(r_list[0]))
    unfolded_r_list = r_list + R

    # Update atoms object
    for i in range(len(alist)):
        alist[i].set_scaled_positions(unfolded_r_list[i])
        alist[i]._calc.atoms.set_scaled_positions(unfolded_r_list[i])

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    This code will unfold the wrapped ASE atoms objects in PBC.
    """)
    # Positional arguments
    parser.add_argument('traj_file', type=str, help='ASE readable structure file name. (Must be in current dir)')
    # Optional arguments
    parser.add_argument('-r', '--ref_file', type=str, default=None,
        help='ASE readable reference structure. Unfolding refers this structure. [Default: first atoms in traj_file]')
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='ASE understanable slice in str format. [Default: all]')
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
    print('This code will unfold the wrapped ASE atoms objects in PBC.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args = argparse()
    traj_file = args.traj_file
    image_slice = args.image_slice

    # @ Main
    from ase.io import read, write
    alist = read(traj_file, image_slice)
    if not isinstance(alist, list):
        alist = [alist]
    ref_atoms = read(args.ref_file)
    unfold_positions(alist, ref_atoms)
    write('unfold_'+traj_file, alist)
