#!/usr/bin/env python

from ase.io import read
from ase.io.trajectory import Trajectory as Traj
import sys

if __name__ == '__main__':
    print("\n\n")
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(100))
    print("            ___________________________           ".center(100))
    print(" __________|  C o d e  b y  Y.J. Choi  |_________ ".center(100))
    print("|______________ ssrokyz@gmail.com _______________|".center(100))
    print("")
    print("*******   This code will give you the shuffled traj file   *******".center(100))
    print("useage ==> ./ase-shuffle.py 'file' ('format')".center(100))
    print("EXAMPLE) ./ase-shuffle.py vasprun.xml (vasp-xml)".center(100))
    print("")
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(100))
    print("")
    if len(sys.argv) == 2 or len(sys.argv) == 3:
        print(("The Number of arguments(= %d) is correct." %(len(sys.argv)-1)).center(100))
        print("\n")
    else:
        print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(100))
        print("\n")
        sys.exit(1)
    read_file = sys.argv[1]
    alist = read(
        read_file,
        index  = ':',
        format = None if len(sys.argv) == 2 else sys.argv[2],
        )
    from random import shuffle
    shuffle(alist)
    in_format = sys.argv[1].split('.')[-1]
    traj_sfile = 'sffld_'+sys.argv[1] if in_format == 'traj' else 'sffld_'+sys.argv[1]+'.traj'
    traj_save = Traj(traj_sfile, 'w')
    for atoms in alist:
        traj_save.write(atoms)
