import sys
import numpy as np

print('\n\n')
print('================================================================================================================='.center(120))
print('            ___________________________           '.center(120))
print(' __________|  C o d e  b y  Y.J. Choi  |_________ '.center(120))
print('|______________ ssrokyz@gmail.com _______________|'.center(120))
print('')
print('*******   This code will read numpy files (npy, npz)  *******'.center(120))
print('useage ==> python -i np-read.py >file<'.center(120))
print('EXAMPLE) python -i np-read.py data.npz'.center(120))
print('')
print('  >>>> Pre-definition <<<<')
print('      np     : Numpy class')
print('     npyz    : Loaded npy or npz object')
print('     dic     : npyz.__dict__')
print('     keys    : npyz.__dict__.keys()')
print('')
print('================================================================================================================='.center(120))
print('')
if len(sys.argv) is 2:
    print(('The Number of arguments(= %d) is correct.' %(len(sys.argv)-1)).center(120))
    print('\n')
else:
    print('*****ERROR***** The number of arguments is not correct *****ERROR*****'.center(120))
    print('\n')
    sys.exit(1)

f_name = sys.argv[1]
npyz   = np.load(f_name)
dic    = npyz.__dict__
keys   = dic.keys()

print(dic)
