#!/usr/bin/env python

import numpy as np
import sys

def alist2numpy(atoms_list):
    try:
        atoms_list[0].get_potential_energy()
    except:
        e_info = False
    else:
        e_info = True
    try:
        atoms_list[0].get_forces()
    except:
        f_info = False
    else:
        f_info = True
    try:
        atoms_list[0].get_stress()
    except:
        s_info = False
    else:
        s_info = True
    box = []
    coord = []
    energy = []
    force = []
    stress = []
    if e_info and f_info and s_info:
        for atoms in atoms_list:
            box.append(np.array(atoms.get_cell(), dtype='float64'))
            coord.append(np.array(atoms.get_positions(), dtype='float64'))
            energy.append(atoms.get_potential_energy())
            force.append(np.array(atoms.get_forces(), dtype='float64'))
            stress.append(np.array(atoms.get_stress(), dtype='float64'))
    elif e_info and f_info and not s_info:
        for atoms in atoms_list:
            box.append(np.array(atoms.get_cell(), dtype='float64'))
            coord.append(np.array(atoms.get_positions(), dtype='float64'))
            energy.append(atoms.get_potential_energy())
            force.append(np.array(atoms.get_forces(), dtype='float64'))
    elif e_info and not f_info and s_info:
        for atoms in atoms_list:
            box.append(np.array(atoms.get_cell(), dtype='float64'))
            coord.append(np.array(atoms.get_positions(), dtype='float64'))
            energy.append(atoms.get_potential_energy())
            stress.append(np.array(atoms.get_stress(), dtype='float64'))
    elif e_info and not f_info and not s_info:
        for atoms in atoms_list:
            box.append(np.array(atoms.get_cell(), dtype='float64'))
            coord.append(np.array(atoms.get_positions(), dtype='float64'))
            energy.append(atoms.get_potential_energy())
    elif not e_info and f_info and s_info:
        for atoms in atoms_list:
            box.append(np.array(atoms.get_cell(), dtype='float64'))
            coord.append(np.array(atoms.get_positions(), dtype='float64'))
            force.append(np.array(atoms.get_forces(), dtype='float64'))
            stress.append(np.array(atoms.get_stress(), dtype='float64'))
    elif not e_info and f_info and not s_info:
        for atoms in atoms_list:
            box.append(np.array(atoms.get_cell(), dtype='float64'))
            coord.append(np.array(atoms.get_positions(), dtype='float64'))
            force.append(np.array(atoms.get_forces(), dtype='float64'))
    elif not e_info and not f_info and s_info:
        for atoms in atoms_list:
            box.append(np.array(atoms.get_cell(), dtype='float64'))
            coord.append(np.array(atoms.get_positions(), dtype='float64'))
            stress.append(np.array(atoms.get_stress(), dtype='float64'))
    elif not e_info and not f_info and not s_info:
        for atoms in atoms_list:
            box.append(np.array(atoms.get_cell(), dtype='float64'))
            coord.append(np.array(atoms.get_positions(), dtype='float64'))
    box = np.array(box, dtype='float64')
    coord = np.array(coord, dtype='float64')
    if e_info:
        energy = np.array(energy, dtype='float64')
    if f_info:
        force = np.array(force, dtype='float64')
    if s_info:
        stress = np.array(stress, dtype='float64')
    return box, coord, energy, force, stress

if __name__ == '__main__':
    print("\n\n")
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(100))
    print("            ___________________________           ".center(100))
    print(" __________|  C o d e  b y  Y.J. Choi  |_________ ".center(100))
    print("|______________ ssrokyz@gmail.com _______________|".center(100))
    print("")
    print("*******   This code will generate npy files from trajectory file   *******".center(100))
    print("useage ==> ./traj2npy.py 'file'".center(100))
    print("EXAMPLE) ./traj2npy.py GST_ran.traj".center(100))
    print("")
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(100))
    print("")
    if len(sys.argv) is 2:
        print(("The Number of arguments(= %d) is correct." %(len(sys.argv)-1)).center(100))
        print("\n")
    else:
        print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(100))
        print("\n")
        sys.exit(1)

    traj_file = sys.argv[1]
    # t2v_ratio = float(sys.argv[2])
    # if sys.argv[3] == 'o':
        # shuffle = True
    # elif sys.argv[3] == 'x':
    shuffle = False
    # else:
        # raise ValueError('Shuffle argument you gave is somehow wrong. It should be o or x. Please check.')

    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(100))
    print('')
    print(('file name: '+traj_file).center(100))
    # print(('training/validation set ratio: '+str(t2v_ratio)).center(100))
    # print(('shuffle: '+str(shuffle)).center(100))
    print('')
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(100))
    print('')

    from time import time
    time_i = time()

    import subprocess as sp
    sp.call(['rm -rf old_npys'], shell=True)
    sp.call(['mv npys old_npys'], shell=True)
    sp.call(['mkdir -p npys'], shell=True)

    from ase.io import read
    atoms_list = read(traj_file, index=':', format='traj')
    image_num = len(atoms_list)
    # valid_num = int(image_num/(t2v_ratio+1))
    # train_num = image_num - valid_num

    log_f = open('npys/log.txt', 'w')
    log_f.write('Made by file, ../'+traj_file+'\n')
    log_f.write('Number of total images: '+str(image_num)+'\n')
    # log_f.write('Ratio of training/validation sets: ('+str(t2v_ratio)+' : 1)\n')
    # log_f.write('Shuffle: '+str(shuffle)+'\n')
    # log_f.write('Number of training sets:   '+str(train_num)+'\n')
    # log_f.write('Number of validation sets: '+str(valid_num)+'\n')
    # if shuffle:
        # from random import shuffle as sffl
        # sffl(atoms_list)
    # else:
        # log_f.write('################################## Caution ####################################\n')
        # log_f.write("     You didn't have order shuffled. Please be aware of what you're doing!     \n")
        # log_f.write('################################## Caution ####################################\n\n')

    box, coord, energy, force, stress = alist2numpy(atoms_list)

    np.save('npys/box.npy', box[:])
    np.save('npys/coord.npy', coord[:])
    if len(energy) != 0:
        np.save('npys/energy.npy', energy[:])
    else:
        log_f.write(' *** energy information not exist *** \n')
        print(' *** energy information not exist *** '.center(100))
    if len(force) != 0:
        np.save('npys/force.npy', force[:])
    else:
        log_f.write(' *** forces information not exist *** \n')
        print(' *** forces information not exist *** '.center(100))
    if len(stress) != 0:
        np.save('npys/stress.npy', stress[:])
    else:
        log_f.write(' *** stress information not exist *** \n')
        print(' *** stress information not exist *** '.center(100))
        print('')
    log_f.close()

    type_txt = open('npys/type.txt', 'w')
    from ss_util import list2numlist as l2nl
    symbols = atoms_list[0].get_chemical_symbols()
    if atoms_list[-1].get_chemical_symbols() != symbols:
        raise ValueError("Chemical symbols seem to be not consistent btw images. Please check")
    symbols_num = l2nl(list(symbols))
    for nums in symbols_num:
        type_txt.write(str(nums)+" ")
    type_txt.write('\n')
    for symbol in symbols:
        type_txt.write(str(symbol)+" ")
    type_txt.close()
    # sp.call(['cp npys/type.txt npys/training/'], shell=True)
    # sp.call(['cp npys/type.txt npys/validation/'], shell=True)

    # from ase.io.trajectory import Trajectory as Traj   
    # train_traj = Traj('npys/training/training_set.traj', 'w')
    # valid_traj = Traj('npys/validation/validation_set.traj', 'w')
    # for i in range(train_num):
        # train_traj.write(atoms_list[i])
    # for i in range(train_num,train_num+valid_num):
        # valid_traj.write(atoms_list[i])

    time_f = time()
    time_d = time_f - time_i
    print(('Total time used: '+str(time_d)+' sec ').center(100))
    print('\n\n')
    log_f = open('npys/log.txt', 'a')
    log_f.write('Total time used: '+str(time_d)+' sec ')
    log_f.close()
