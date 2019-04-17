#!/usr/bin/env python

import numpy as np
import sys

print("\n")
print("#######################################################################################".center(120))
print("")
print("useage ==> ./ase-volume-histo.py 'atoms list file'".center(120))
print("           EXAMPLE) ./ase-volume-histo.py vasprun.xml".center(120))
print("")
print("#######################################################################################".center(120))
print("")
if len(sys.argv) is 2:
    print("The Number of arguments is correct.".center(120))
    print("\n")
else:
    print(">>>>> ERROR <<<<<     The number of arguments is not correct     >>>>> ERROR <<<<<".center(120))
    print("\n")
    sys.exit(1)

from ase.io import read
alist = read(sys.argv[1], ':')

vol_list = []
for atoms in alist:
    vol_list.append(np.linalg.det(atoms.get_cell()))
# Post process
anum = len(alist[0])
vol_list = np.array(vol_list)/anum
average = np.mean(vol_list)
std = np.std(vol_list)
    
## plot
from matplotlib import pyplot as plt
n, bins, patches = plt.hist(vol_list, bins=200, facecolor='purple', alpha=0.70)
max_height = np.sort(n)[-10]
plt.title('Volume Histogram (Ang^3/atom)')
plt.xlabel('%d images, average = %.3f, sigma = %.3f' % (len(vol_list), average, std))
plt.ylabel('population')
plt.barh(max_height/5, std, height=max_height/50, left=average, color='black')
plt.grid(True)
plt.show()
