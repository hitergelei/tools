#!/usr/bin/env python

from ase.io.trajectory import Trajectory as Traj
from ase.io import read
import sys
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    This code will convert the out.dump file to ASE trajectory file.
    If there is atomic-energy info, potential energy can be obtained without the LAMMPS thermo log file.
    """)
    # Positional arguments
        # No positional argument.
    # Optional arguments
    parser.add_argument('-i', '--dump_file', type=str, default='out.dump', help='ASE readable dump file name.')
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='ASE understanable slice in str format. [Default: all]')
    parser.add_argument('-l', '--log_file', type=str, default=None, help='Read LAMMPS log file, usually "log.lammps", to get potential energy info. If provided (even when there is atomic-energy info), it has the first priority. [Default: Do not load]')
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
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args = argparse()
    dump_file = args.dump_file
    image_slice = args.image_slice
    log_file = args.log_file

    print('Reading dump file...'.center(120))
    dump_inp = read(dump_file, image_slice, format='lammps-dump')
    if not isinstance(dump_inp, list):
        dump_inp = [dump_inp]
    print('Successively read dump file!'.center(120))
    if log_file:
        log = open(log_file, "r")
        print('Read lof file "{}".'.format(log_file))
        print('((WARNING)) Even when there is atomic-energy info, potential energy will be read from the thermo log file.')
        epot = []

        ## Get first PotEng line number
        while True:
            line = log.readline()
            field = line.split()
            if 'PotEng' in field:
                poteng_nword = field.index('PotEng') +2
                epot.append(field[poteng_nword])
                break

        ## Get delta(PotEng line number)
        poteng_dline =1
        while True:
            line = log.readline()
            if not line:
                poteng_dline = False # there is only 1 image
                break
            field = line.split()
            if 'PotEng' in field:
                epot.append(field[poteng_nword])
                break
            else:
                poteng_dline +=1

        ## Get rest of PotEng from now on
        if poteng_dline:
            while True:
                for i in range(poteng_dline):
                    line = log.readline()
                if not line:
                    break
                field = line.split()
                if field[poteng_nword -2] == 'PotEng':
                    epot.append(line.split()[poteng_nword])
                else:
                    break

        if len(epot) != len(dump_inp):
            print(" ***** ERROR ***** The # of images in log file and dump file do not match\n")
            print("   i.e.) len(epot) != len(dump_inp)\n")
            print("   "+str(len(epot))+" != "+str(len(dump_inp))+"\n")
            sys.exit(1)

    traj = Traj("lmp-result.traj", "w")
    for i in range(len(dump_inp)):
        if i % 1000 == 0:
            print("Writing "+str(i)+"th image")
        atoms = dump_inp[i]
        atoms._pbc = np.array([True, True, True])
        atoms._calc.atoms._pbc = np.array([True, True, True])
        if log_file:
            atoms._calc.results['energy'] = epot[i]
        else:
            atoms._calc.results['energy'] = np.sum(atoms._calc.results['atomic_energies'])
        traj.write(atoms)
    traj.close()


    print("\n\n#######################################################################################")
    print("      %%%%%%%%%%% This code will covert lammps results to traj file %%%%%%%%%")
    print("         useage ==> ./lmp2traj.py 'dump file' 'energy logfile'")
    print("              e.g.) ./lmp2traj.py out.dump log.lammps")
    print("                The result file name will be lmp-result.traj")
    print("  *****NOTE***** There is some issue when ase.io.lammpsrun import dump file. *****NOTE*****")
    print("       Make sure that you revised it. (velocity-> vel /1000/units.fs, symbol issue)")
    print("#######################################################################################\n\n")

