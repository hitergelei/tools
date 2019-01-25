#!/usr/bin/env python

import numpy as np

def D4_2_D2(D4_arr):
    shape = list(D4_arr.shape)
    if len(shape) != 4:
        raise ValueError('Array dimension is not 4')
    D2_arr = np.zeros((shape[0]*shape[2], shape[1]*shape[3]))
    for a in range(shape[0]):
        for b in range(shape[1]):
            for c in range(shape[2]):
                for d in range(shape[3]):
                    D2_arr[3*a+c, 3*b+d] = D4_arr[a, b, c, d]
    return D2_arr

if __name__ == '__main__':
    import sys
    import datetime

    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print("")
    print(">>>>> Code by Youngjae Choi @ POSTECH <<<<<".center(120))
    print(("code started time: "+time).center(120))
    print("")
    print("===================================================================================".center(120))
    print("")
    print("Useage ==> ./fc_mapping.py 'phonopy phonon object pckl file' 'vmin' 'vmax'".center(120))
    print("    Or ==> ./fc_mapping.py 'phonopy phonon object pckl file'              ".center(120))
    print("Note) set 'vmin' or 'vmax' as None to automatic setting".center(120))
    print("*** Note) Phonon object MUST HAVE self._force_constants array already ***".center(120))
    print("")
    print("===================================================================================".center(120))
    if len(sys.argv) is 2:
        print("")
        vrange = False
        vmin = None
        vmax = None
    elif len(sys.argv) is 4:
        vmin = sys.argv[2]
        vmax = sys.argv[3]
        if vmin == 'None':
            vmin = None
        if vmax == 'None':
            vmax = None
        vrange = True
    else:
        print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(120))
        print("")
        print("")
        sys.exit(1)

    ph_pckl_file = sys.argv[1]

    #### Load object
    import pickle as pckl
    with open(ph_pckl_file, 'rb') as f:
        phonon_obj = pckl.load(f)
    ## Dimension conversion
    D2_fc = D4_2_D2(phonon_obj._force_constants)

    #### Plot
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    ## Color scheme
    cm = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', ['k', 'r', 'r', 'r', 'y', 'y', 'y', 'w', 'w', 'w'], 256)
    ## Plot figure
    img = plt.imshow(np.absolute(D2_fc), cmap=cm, vmin=vmin, vmax=vmax)
    plt.colorbar(img, cmap=cm)
    ## Draw grid
    # Grid interval
    anum = len(phonon_obj._unitcell.symbols)
    for y in range(-1, len(D2_fc), len(D2_fc) // anum):
        plt.axhline(y+0.5, color='w', linewidth=0.5)
    for x in range(-1, len(D2_fc[0]), len(D2_fc[0]) // anum):
        plt.axvline(x+0.5, color='w', linewidth=0.5)
    plt.xlim(0, None)
    plt.ylim(None, 0)
    plt.show()

