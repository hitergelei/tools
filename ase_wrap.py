#!/usr/bin/env python

import numpy as np

def wrap_alist(alist, log=False):
    for i in range(len(alist)):
        ## Backup results info
        if hasattr(alist[i]._calc, 'results'):
            results = alist[i]._calc.results.copy()
        else:
            results = None
        alist[i].wrap(eps=0.)
        ## Recover results info
        if results is not None:
            alist[i]._calc.results = results
            alist[i]._calc.atoms = alist[i].copy()
        #### Print every 1000 process
        if log and i % 1000 == 0: 
            print((str(i)+'-th image wrapped').center(120))
    return alist

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
        alist_file = sys.argv[1]
    else:
        raise ValueError("*****ERROR***** The number of arguments is not correct *****ERROR*****")

    from ase.io import read, write
    alist = read(alist_file, ':')
    if not isinstance(alist, list):
        alist = [alist]

    wrapped_alist = wrap_alist(alist, log=True)
    write('wrapped-'+alist_file, wrapped_alist, format='traj')
