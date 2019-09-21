#!/usr/bin/env python

import numpy as np

### Global params
n = 5
freq_range = (0.,2.5)
slices = [[0.6, 1.0], [1.0,1.45], [1.45,2.2]]
file_range = [0,333]
delta_T = 20. # (ps)

### Main part
if   delta_T == 35.:
    image_range = '1500:'
elif delta_T == 20.:
    image_range = '3000:'
elif delta_T == 5.:
    image_range = '4500:'
else:
    raise ValueError
DOS = []
for i in range(*file_range):
    npz = np.load('lj-Ar-prim-50K-10fs-{}x{}x{}-NVE-many-{}.traj_dt0.01_img{}.npz'.format(n,n,n,i,image_range))
    ados = npz['ADOS']
    DOS.append(np.sum(ados.reshape(-1,ados.shape[2]), axis=0))
f = npz['f']
## SLICE
# Get slice bool
slice_bool = []
for i in range(len(slices)):
    slice_bool.append(((f>=slices[i][0]).astype(int) * (f<slices[i][1]).astype(int)).astype(bool))
# Slice DOS
DOS_slice = []
for i in range(len(DOS)):
    bins = np.zeros(len(slice_bool))
    for j in range(len(slice_bool)):
        bins[j] = np.sum(DOS[i][slice_bool[j]])
    DOS_slice.append(bins)
# Print integrated variance
df = 1./delta_T
print('Means: {}'.format(np.mean(DOS_slice, axis=0)*df))
print('Variances: {}'.format(np.var(DOS_slice, axis=0)*df))
print('Mean of variances: {}'.format(np.mean(np.var(DOS_slice, axis=0)*df)))


## Plot continueous variance
variance = np.var(DOS, axis=0)
with open('tmp.txt', 'w') as txt:
    for i in range(len(f)):
        txt.write('{:10.5f}'.format(f[i]))
        txt.write('   ')
        txt.write('{:10.5f}'.format(variance[i]))
        txt.write('\n')

from matplotlib import pyplot as plt
freq_bool = ((f>=freq_range[0]).astype(int) * (f<=freq_range[1]).astype(int)).astype(bool)
f = f[freq_bool]
plt.plot(variance[freq_bool], f, label='Statistic')
anal_sol = 2.*delta_T/3./n**3 *np.mean(DOS, axis=0) / 10
plt.plot(anal_sol[freq_bool], f, label='Analytic')
# ratio = variance / anal_sol
# plt.plot(ratio[freq_bool], f)
plt.legend()
plt.ylim(0.,2.5)
plt.xlim(0,None)
plt.subplots_adjust(left=0.35, bottom=0.15, right=0.60, top=0.95, wspace=0.2, hspace=0.2)
plt.tick_params(axis="both",direction="in", labelsize='xx-large')
plt.show()



