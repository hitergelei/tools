#!/usr/bin/env python

import numpy as np

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
    print("Useage  ==> ./ase-force-rectify.py >structure file< >cutoff force(eV/Ang)<".center(120))
    print("Example ==> ./ase-force-rectify.py gst.traj 10.".center(120))
    print("")
    print("Return ------------> blacklist_(original name).npy, rectified_(original name).traj".center(120))
    print("")
    print("==================================================================================================".center(120))
    

    ## Preprocess
    cutoff_force = np.array(sys.argv[2], dtype=np.float)
    from ase.io import read, write
    alist = read(sys.argv[1], ':')

    ## Gather max forces of each system
    blacklist = []
    for i in range(len(alist)):
        forces = np.linalg.norm(alist[i].get_forces(), axis=-1)
        atom_ind = np.arange(len(alist[i]))[forces > cutoff_force]
        if len(atom_ind) != 0:
            blacklist.extend(np.transpose([[i]*len(atom_ind), atom_ind, forces[atom_ind]]))
            # print(blacklist[-1])
    blacklist = np.array(blacklist)

    ## Save
    # Numpy save
    np.save('blacklist_cut{}_{}'.format(int(cutoff_force), sys.argv[1]), blacklist)
    # Alist save
    img2remove = np.unique(blacklist[:,0]).astype(np.int)
    select_bool = np.ones(len(alist), dtype=np.bool)
    select_bool[img2remove] = False
    new_alist = []
    for i in range(len(alist)):
        if select_bool[i]:
            new_alist.append(alist[i])
    write('rectified_'+sys.argv[1], new_alist, 'traj')

    ## Print
    print("")
    print(" >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RESULTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ".center(120))
    print("{} file contains {} atoms object images.".format(sys.argv[1], len(alist)).center(120))
    print("{} images contain {} atomic-forces higher than cutoff {} (eV/Ang).".format(len(img2remove), len(blacklist), sys.argv[2]).center(120))
    print("{} images will be excluded. ( {} % of original )".format(len(img2remove), len(img2remove)/len(alist)*100.).center(120))
    print("")

