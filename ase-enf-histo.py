#!/usr/bin/env python

from matplotlib import pyplot as plt
import sys
import numpy as np

if __name__ == '__main__':
    print("\n\n")
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(120))
    print("            ___________________________           ".center(120))
    print(" __________|  C o d e  b y  Y.J. Choi  |_________ ".center(120))
    print("|______________ ssrokyz@gmail.com _______________|".center(120))
    print("")
    print("*******   This code will show you a histograms of energies and forces from atoms list   *******".center(120))
    print("useage ==> ./ase-histogram-traj.py 'atoms list file' 'data type(e/f)'                    ".center(120))
    print("    or ==> ./ase-histogram-traj.py 'atoms list file' 'data type(e/f)' 'bin min' 'bin max'".center(120))
    print("EXAMPLE) ./ase-histogram-traj.py vasprun.xml e          ".center(120))
    print("      or ./ase-histogram-traj.py vasprun.xml e -250 -200".center(120))
    print("")
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(120))
    print("")
    if len(sys.argv) is 3 or len(sys.argv) is 5 :
        print(("The Number of arguments(= %d) is correct." %(len(sys.argv)-1)).center(120))
        print("\n")
    else:
        print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(120))
        print("\n")
        sys.exit(1)
    alist_file = sys.argv[1]
    from ase.io import read
    alist = read(alist_file, ':')
    d_type = sys.argv[2]
    if len(sys.argv) == 5:
        bin_range = True
        bin_min = float(sys.argv[3])
        bin_max = float(sys.argv[4])
    else:
        bin_range = False
        bin_min = None
        bin_max = None

    if d_type == 'e':
        from ss_util import E_fromAlist
        data = E_fromAlist(alist)
    elif d_type == 'f':
        from ss_util import F_fromAlist
        data = np.linalg.norm(F_fromAlist(alist).reshape((-1,3)), axis=-1)
    else:
        raise ValueError('data type is wrong!!! Must be "e" or "f"')

    average = np.average(data)
    std = np.std(data)
    if bin_range:
        n, bins, patches = plt.hist(data, bins=200, range=[bin_min, bin_max], facecolor='purple', alpha=0.70)
    else:
        n, bins, patches = plt.hist(data, bins=200, facecolor='purple', alpha=0.70)
    max_height = np.amax(n)
    
    #### plot
    plt.title('Histogram')
    plt.xlabel('%d data, average = %.3f, sigma = %.3f' % (len(data), average, std))
    plt.ylabel('population')
    plt.barh(max_height/5, std, height=max_height/100, left=average, color='black')
    plt.grid(True)
    plt.show()

