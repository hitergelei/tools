#!/usr/bin/env python

# @ Params ======================================================================================
load_range = range(30)
dt = 0.01
corr_len = 100.
avg_intvl = 10
img_slice = ':'
two_term = False
# sub_bind = False
sub_bind = True

# @ Main ========================================================================================
import numpy as np
from ase import units
dt_ase = dt *1e3 *units.fs

from ase.io import read
hfacf = []
temp = []
volu = []
for i in load_range:
    fname = 'job-{}/kappa-whole/lmp-results.traj_dt{}_c{}_n{}_a{}_t{}_b{}.npy'.format(i, dt, corr_len, img_slice, avg_intvl, two_term, sub_bind)
    hfacf.append(np.load(fname))
    temp.append(np.load('{}-temp.npy'.format(fname)))
    volu.append(np.load('{}-volu.npy'.format(fname)))
hfacf = np.mean(hfacf, axis=0)
mean_temp = np.mean(temp)
mean_volu = np.mean(volu)

#
norm_hfacf = hfacf / hfacf[0]
kappa = np.add.accumulate(hfacf) *dt_ase /units.kB /mean_temp**2 /mean_volu *units._e *units.second *1e10

avg_norm_hfacf = []
for i in range(len(hfacf)):
    tmp = []
    for j in range(3):
        tmp.append(norm_hfacf[i,j,j])
    avg_norm_hfacf.append(np.mean(tmp))

avg_kappa = []
for i in range(len(kappa)):
    tmp = []
    for j in range(3):
        tmp.append(kappa[i,j,j])
    avg_kappa.append(np.mean(tmp))

from matplotlib import pyplot as plt
start = dt/2.
t = np.arange(len(kappa), dtype=float) *dt +start
fig, ax1 = plt.subplots(3,3)
for i in range(3):
    for j in range(3):
        ax2 = ax1[i,j].twinx()
        ax1[i,j].plot(t, norm_hfacf[:,i,j], c='b')
        ax2.plot(t, kappa[:,i,j], c='r')
        #
        ax1[i,j].set_xlabel('Time (ps)', fontsize='x-large')
        ax1[i,j].set_ylabel('Scaled HFACF (Arb. Unit)', fontsize='x-large', color='b')
        ax2.set_ylabel('$\kappa_{}$$_{}$ (W/mK)'.format(i+1, j+1), fontsize='x-large', color='r')
        ax1[i,j].tick_params(axis="x",direction="in", labelsize='x-large')
        ax1[i,j].tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='b')
        ax2.tick_params(axis="y",direction="in", labelsize='x-large',colors='r',  labelcolor='r')
        ax1[i,j].grid(alpha=0.5)
        ax2.grid(alpha=0.5)
        plt.title('IS={}'.format(img_slice), fontsize='x-large')
plt.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.95, wspace=0.80, hspace=0.40)

# Average
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(t, avg_norm_hfacf[:], c='b')
ax2.plot(t, avg_kappa[:], c='r')
#
ax1.set_xlabel('Time (ps)', fontsize='x-large')
ax1.set_ylabel('Scaled HFACF (Arb. Unit)', fontsize='x-large', color='b')
ax2.set_ylabel('$\kappa$ (W/mK)', fontsize='x-large', color='r')
ax1.tick_params(axis="x",direction="in", labelsize='x-large')
ax1.tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='b')
ax2.tick_params(axis="y",direction="in", labelsize='x-large',colors='r',  labelcolor='r')
ax1.grid(alpha=0.5)
ax2.grid(alpha=0.5)
plt.title('M) IS={}, AI={}, dt={}, T={:.2f}'.format(img_slice, avg_intvl, dt, mean_temp))
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.85, top=0.90)

plt.show()

