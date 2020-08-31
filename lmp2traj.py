#!/usr/bin/env python

from ase import units
from ase.io.trajectory import Trajectory as Traj
from ase.io import read
import sys
import numpy as np

def read_lmp_log(
    log_file='log.lammps',
    ):
    # Load file
    with open(log_file) as f:
        lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].split()

    # Find start line
    for i in range(len(lines)):
        if len(lines[i]) > 0:
            if lines[i][0] == 'Step':
                label_line = i
                start_line = i+1

    # Find last line
    end_line = len(lines)-1 # In case, it cannot find the end line.
    for i in range(start_line, len(lines), 1):
        try:
            float(lines[i][0])
        except ValueError:
            end_line = i-1
            break
        except:
            raise RuntimeError
        else:
            pass

    # Convert as array
    arr = np.array(lines[start_line:end_line+1], float).T

    # Make dict object
    info = dict()
    for i in range(len(lines[label_line])):
        info[lines[label_line][i]] = arr[i]

    return info

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    This code will convert the out.dump file to ASE trajectory file.
    If there is atomic-energy info, potential energy can be obtained without the LAMMPS thermo log file.
    But priority for potential energy is at the thermo log file.
    """)
    # Positional arguments
        # No positional argument.
    # Optional arguments
    parser.add_argument('-i', '--dump_file', type=str, default='out.dump', help='ASE readable dump file name.')
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='ASE understanable slice in str format. [Default: all]')
    parser.add_argument('-l', '--without_log', dest='load_log', action='store_false', help='Force to do NOT read LAMMPS log file, "log.lammps". Stress and potential energy info will be obtained. If provided (even when there is atomic-energy info), it has the first priority for potential energies. [Default: log.lammps]')
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
    print('This code will convert the out.dump file to ASE trajectory file.'.center(120))
    print('If there is atomic-energy info, potential energy can be obtained without the LAMMPS thermo log file.'.center(120))
    print('But priority for potential energy is at the thermo log file.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args = argparse()
    dump_file = args.dump_file
    image_slice = args.image_slice
    load_log = args.load_log

    print('Reading dump file...'.center(120))
    dump_inp = read(dump_file, image_slice, format='lammps-dump')
    if not isinstance(dump_inp, list):
        dump_inp = [dump_inp]
    print('Successively read dump file!'.center(120))
    if load_log:
        info = read_lmp_log()
        print('Read log.lammps file.')
        print('((WARNING)) Even when there is atomic-energy info, potential energy will be read from the thermo log file.')

        if len(list(info.values())[0]) != len(dump_inp):
            print(" ***** ERROR ***** The # of images in log file and dump file do not match.")
            print("   i.e.) len(log) != len(dump_inp)")
            print("    ==> {} != {}".format(len(list(info.values())[0]), len(dump_inp)))
            raise RuntimeError()

    traj = Traj("lmp-result.traj", "w")
    for i in range(len(dump_inp)):
        if i % 1000 == 0:
            print("Writing "+str(i)+"th image")
        atoms = dump_inp[i]
        atoms._pbc = np.array([True, True, True])
        atoms._calc.atoms._pbc = np.array([True, True, True])
        if load_log:
            atoms._calc.results['energy'] = info['PotEng'][i]
            # Unit for ASE
            atoms._calc.results['stress'] = -np.array([
                info['Pxx'][i],
                info['Pyy'][i],
                info['Pzz'][i],
                info['Pyz'][i],
                info['Pxz'][i],
                info['Pxy'][i],
                ], float) * 1e-4 * units.GPa
        else:
            atoms._calc.results['energy'] = np.sum(atoms._calc.results['atomic_energies'])
        traj.write(atoms)
    traj.close()

    print("\n\n=======================================================================================")
    print("      %%%%%%%%%%% This code will covert lammps results to traj file %%%%%%%%%")
    print("         useage ==> ./lmp2traj.py 'dump file' 'energy logfile'")
    print("              e.g.) ./lmp2traj.py out.dump log.lammps")
    print("                The result file name will be lmp-result.traj")
    print("  *****NOTE***** There is some issue when ase.io.lammpsrun import dump file. *****NOTE*****")
    print("       Make sure that you revised it. (velocity-> vel /1000/units.fs, symbol issue)")
    print("=======================================================================================\n\n")

