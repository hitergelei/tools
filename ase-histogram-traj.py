#!/usr/bin/env python

from matplotlib import pyplot as plt
import sys
import numpy as np

if __name__ == '__main__':
    print("\n\n")
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(100))
    print("            ___________________________           ".center(100))
    print(" __________|  C o d e  b y  Y.J. Choi  |_________ ".center(100))
    print("|______________ ssrokyz@gmail.com _______________|".center(100))
    print("")
    print("*******   This code will show you a histogram from traj file   *******".center(100))
    print("useage ==> ./ase-histogram-traj.py 'traj file' 'data type(e/f)'                    ".center(100))
    print("    or ==> ./ase-histogram-traj.py 'traj file' 'data type(e/f)' 'bin min' 'bin max'".center(100))
    print("EXAMPLE) ./ase-histogram-traj.py gst.traj e          ".center(100))
    print("      or ./ase-histogram-traj.py gst.traj e -200 -250".center(100))
    print("")
    print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$".center(100))
    print("")
    if len(sys.argv) is 3 or len(sys.argv) is 5 :
        print(("The Number of arguments(= %d) is correct." %(len(sys.argv)-1)).center(100))
        print("\n")
    else:
        print("*****ERROR***** The number of arguments is not correct *****ERROR*****".center(100))
        print("\n")
        sys.exit(1)
    traj_file = sys.argv[1]
    from ase.io.trajectory import Trajectory as Traj
    traj = Traj(traj_file)
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
        from ss_util import E_FromTraj
        data = E_FromTraj(traj)
    elif d_type == 'f':
        from ss_util import F_FromTraj
        data = np.reshape(F_FromTraj(traj),-1)
    else:
        raise ValueError('data type is wrong!!! Must be "e" or "f"')

    average = np.average(data)
    std = np.std(data)

    if bin_range:
        n, bins, patches = plt.hist(data, bins=200, range=[bin_min, bin_max], facecolor='purple', alpha=0.70)
    else:
        n, bins, patches = plt.hist(data, bins=200, facecolor='purple', alpha=0.70)
    
    #### plot
    plt.title('histogram')
    plt.xlabel('%d data, average = %.3f, sigma = %.3f' % (len(data), average, std))
    plt.ylabel('population')
    plt.grid(True)
    plt.show()

