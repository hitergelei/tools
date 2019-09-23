#!/usr/bin/env python
from ase.io.trajectory import Trajectory as Traj
from ase.io import read, write
from ase.calculators import vasp
import sys
import subprocess as sp

def _parse_slice(s):
    if ':' in s:
        a = [int(e) if e.strip() else None for e in s.split(":")]
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
    parser.add_argument('-n', '--image_slice', type=_parse_slice, default=':', help='Image slice following python convention. default=":" (e.g.) -n :1000:10')
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
    print('This code converts ase readable file to VASP POSCAR formats.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()
    inp_file = args.inp_file

    ##
    print("\nI'll extract informations from '"+inp_file+"' file.")
    print("POSCARs will be stored in folder >> 'vasp_images_"+inp_file+"'")
    alist = read(inp_file, args.image_slice)
    if not isinstance(alist, list):
        alist = [alist]
    print("\nTotally, {} POSCARs will be written\n".format(len(alist)))

    ##
    indices = range(int(1e6))[args.image_slice]
    sp.call("mkdir vasp_images_{}".format(inp_file), shell=True)
    for i in range(len(indices)):
        write("vasp_images_{}/POSCAR_{}".format(inp_file, str(indices[i]).zfill(6)), alist[i])
