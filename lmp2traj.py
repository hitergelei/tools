#!/usr/bin/env python

from ase import units
from ase.io.trajectory import Trajectory as Traj
from ase.io import read
import sys
import numpy as np

def read_lmp_log(
    log_file='log.lammps',
    image_slice=':',
    ):
    # Load file
    with open(log_file, errors='ignore') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].split()

    # Find start line
    start_lines = []
    for i in range(len(lines)):
        if len(lines[i]) > 0:
            if lines[i][0] == 'Step':
                start_lines.append(i+1)
    print("Found start lines: {}".format(start_lines))

    # Find last line
    end_lines = []
    for j in range(len(start_lines)):
        end_line = len(lines)-1 # In case, it cannot find the end line.
        for i in range(start_lines[j], len(lines), 1):
            try:
                float(lines[i][0])
            except (IndexError, ValueError):
                end_line = i-1
                break
            except:
                raise RuntimeError
            else:
                pass
        end_lines.append(end_line)

    print("Found end line: {}".format(end_lines))

    # Convert as array
    arr = []
    for i in range(len(start_lines)):
        arr.append(np.array(lines[start_lines[i]:end_lines[i]+1], float))

    # Image slicing
    num_img = []
    for i in range(len(arr)):
        num_img.append(len(arr[i]))

    img_bool = np.array([False] *np.sum(num_img))
    from ss_util import parse_slice
    img_bool[parse_slice(image_slice)] = True

    # Partitioning img_bool
    s = 0
    i_b = []
    for i in range(len(num_img)):
        i_b.append(img_bool[s:s+num_img[i]])
        s += num_img[i]
    img_bool = i_b

    # Screening arr
    a = []
    for i in range(len(num_img)):
        a.append(arr[i][img_bool[i]].T)
    arr = a

    # Make dict object
    info = []
    for j in range(len(start_lines)):
        info.append(dict())
        for i in range(len(lines[start_lines[j]])):
            info[j][lines[start_lines[j]-1][i]] = arr[j][i]

    return info

# def read_lmp_dump(
    # dump_file='out.dump',
    # ):

    # # Load file
    # f = open(dump_file)
    # lines = f.readlines()

    # # Info
    # len_atoms = int(lines[3].split()[0])
    # len_img = lines//(len_atoms+9)
    # aid = []
    # element = []
    # for i in range(10,10+len_atoms,1):
        # aid.append(int(lines[i][0]))
        # element.append((lines[i][1]))
    # # re-order
    # order = np.argsort(aid)
    # element = np.array(element)[order]
    # tilt_items = lines[4].split()[3:]
    # if len(tilt_items) == 6:
        # tilt_bool = True
    # elif len(tilt_items) == 3:
        # tilt_bool = False
    # else:
        # raise NotImplementedError

    # # Iteration loop
    # from ase.atoms import Atoms
    # for img_i in range(len_img):
        # lo = []
        # hi = []
        # tilt = []
        # id = []
        # types = []
        # positions = []
        # element = [] ## ssrokyz
        # pe = [] ## ssrokyz
        # scaled_positions = []
        # velocities = []
        # forces = []
        # quaternions = []
        # # @ BOX
        # # save labels behind "ITEM: BOX BOUNDS" in
        # # triclinic case (>=lammps-7Jul09)
        # for l in range(img_i+5,img_i+8,1):
            # fields = lines[l].split()
            # lo.append(float(fields[0]))
            # hi.append(float(fields[1]))
            # if (len(fields) >= 3):
                # tilt.append(float(fields[2]))

        # # determine cell tilt (triclinic case!)
        # if (len(tilt) >= 3):
            # # for >=lammps-7Jul09 use labels behind
            # # "ITEM: BOX BOUNDS" to assign tilt (vector) elements ...
            # if (len(tilt_items) >= 3):
                # xy = tilt[tilt_items.index('xy')]
                # xz = tilt[tilt_items.index('xz')]
                # yz = tilt[tilt_items.index('yz')]
            # # ... otherwise assume default order in 3rd column
            # # (if the latter was present)
            # else:
                # xy = tilt[0]
                # xz = tilt[1]
                # yz = tilt[2]
        # else:
            # xy = xz = yz = 0
        # xhilo = (hi[0] - lo[0]) - (xy**2)**0.5 - (xz**2)**0.5
        # yhilo = (hi[1] - lo[1]) - (yz**2)**0.5
        # zhilo = (hi[2] - lo[2])
        # if xy < 0:
            # if xz < 0:
                # celldispx = lo[0] - xy - xz
            # else:
                # celldispx = lo[0] - xy
        # else:
            # celldispx = lo[0]
        # celldispy = lo[1]
        # celldispz = lo[2]

        # cell = [[xhilo, 0, 0], [xy, yhilo, 0], [xz, yz, zhilo]]
        # celldisp = [[celldispx, celldispy, celldispz]]
                
        # line = f.readline()

        # alist.append(Atoms(
            # symbols=element,
            # positions=positions,
            # celldisp=celldisp,
            # cell=cell,
            # ))

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
    dump_inp = read(dump_file, image_slice, format='lammps-dump-text')
    if not isinstance(dump_inp, list):
        dump_inp = [dump_inp]
    print('Successively read dump file!'.center(120))
    if load_log:
        info = read_lmp_log(image_slice=image_slice)
        print('Read log.lammps file.')
        print('((WARNING)) Even when there is atomic-energy info, potential energy will be read from the thermo log file.')

        num_img = []
        for i in range(len(info)):
            num_img.append(len(info[i][list(info[i].keys())[0]]))
        if np.sum(num_img) != len(dump_inp):
            print(" ***** ERROR ***** The # of images in log file and dump file do not match.")
            print("   i.e.) len(log) != len(dump_inp)")
            print("    ==> {} != {}".format(len(list(info.values())[0]), len(dump_inp)))
            raise RuntimeError()
    else:
        num_img=[len(dump_inp)]

    # Partitioning dump_inp
    s = 0
    d_i = []
    for i in range(len(num_img)):
        d_i.append(dump_inp[s:s+num_img[i]])
        s += num_img[i]
    dump_inp = d_i

    for d in range(len(num_img)):
        traj = Traj("lmp-results-{}-{}.traj".format(image_slice, d), "w")
        for i in range(len(dump_inp[d])):
            if i % 1000 == 0:
                print("Writing "+str(i)+"th image")
            atoms = dump_inp[d][i]
            atoms._pbc = np.array([True, True, True])
            atoms._calc.atoms._pbc = np.array([True, True, True])
            if load_log:
                atoms._calc.results['energy'] = info[d]['PotEng'][i]
                # Unit for ASE
                if 'Pxx' in info[d].keys():
                    atoms._calc.results['stress'] = -np.array([
                        info[d]['Pxx'][i],
                        info[d]['Pyy'][i],
                        info[d]['Pzz'][i],
                        info[d]['Pyz'][i],
                        info[d]['Pxz'][i],
                        info[d]['Pxy'][i],
                        ], float) * 1e-4 * units.GPa
            else:
                atoms._calc.results['energy'] = np.sum(atoms._calc.results['energies'])
            traj.write(atoms)
        traj.close()

    print("\n\n=======================================================================================")
    print("      %%%%%%%%%%% This code will covert lammps results to traj file %%%%%%%%%")
    print("         useage ==> ./lmp2traj.py 'dump file' 'energy logfile'")
    print("              e.g.) ./lmp2traj.py out.dump log.lammps")
    print("                The result file name will be lmp-results.traj")
    print("  *****NOTE***** There is some issue when ase.io.lammpsrun import dump file. *****NOTE*****")
    print("       Make sure that you revised it. (velocity-> vel /1000/units.fs, symbol issue)")
    print("=======================================================================================\n\n")

