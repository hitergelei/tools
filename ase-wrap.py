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
    print('')
    print("==================================================================================================".center(120))
    print('')
    print("Useage  ==> ./ase-wrap.py 'file'".center(120))
    print("Example ==> ./ase-wrap.py gst.traj".center(120))
    print('')
    print("==================================================================================================".center(120))
    
    if len(sys.argv) == 2:
        obj_file = sys.argv[1]
    else:
        raise ValueError("*****ERROR***** The number of arguments is not correct *****ERROR*****")

    from ase.io import read, write
    obj = read(obj_file, ':')
    if not isinstance(obj, list):
        obj = [obj]

    # from ase.io.trajectory import Trajectory as Traj
    # if obj_file.split('.')[-1] == 'traj':
        # traj = Traj('wrapped-'+obj_file, 'w')
    # else:
        # traj = Traj('wrapped-'+obj_file+'.traj', 'w')
    alist = []

    #### Check the available info
    atoms = obj[0]
    try:
        atoms.get_potential_energy()
    except:
        print('No energy info'.center(120))
        E = False
    else:
        E = True
    try:
        atoms.get_forces()
    except:
        print('No force info'.center(120))
        F = False
    else:
        F = True
    try:
        atoms.get_stress()
    except:
        print('No stress info'.center(120))
        S = False
    else:
        S = True

    i=0
    for atoms in obj:
        ## Backup results info
        if atoms._calc is not None:
            result = atoms._calc.results.copy()
        atoms.wrap(eps=0.)
        ## Recover results info
        if atoms._calc is not None:
            atoms._calc.result = result.copy()
            atoms._calc.atoms = atoms.copy()
        # traj.write(atoms)
        alist.append(atoms)
        #### Print every 1000 process
        i+=1
        if i % 1000 == 0: 
            print((str(i)+'-th image written').center(120))

    if obj_file.split('.')[-1] == 'traj':
        write('wrapped-'+obj_file, alist)
    else:
        write('wrapped-'+obj_file+'.traj', alist)
