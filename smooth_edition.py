#!/usr/bin/env python

import numpy as np

# Params
r_cs = 0.5
r_c = 1.0
nbin = 500
x_max = 1.5

# Main
dx = x_max /nbin
X = np.arange(nbin+1, dtype=float) /nbin *x_max

def DPMD(X, r_c):
    Y = []
    for x in X:
        if x < r_c:
            Y.append(1.)
        else:
            Y.append(0.)
    return np.array(Y)

def DPMD_SE(X, r_c, r_cs):
    Y = []
    for x in X:
        if x < r_cs:
            Y.append(1.)
        elif x <= r_c:
            Y.append(0.5*np.cos(np.pi*(x-r_cs)/(r_c-r_cs)) +0.5)
        else:
            Y.append(0.)
    return np.array(Y)

from matplotlib import pyplot as plt
# DPMD
dpmd = DPMD(X, r_c)
fig, ax1 = plt.subplots()
ax1.plot(X, dpmd, c='r', label='s(x)', lw=3)
# ax1.set_ylim(-0.2, 1.2)
ax2 = ax1.twinx()
ax2.plot(X[:-1] +dx/2., (dpmd[1:] - dpmd[:-1]) /dx, c='b', label='ds(r)/dr', ls=(0, (5, 5)), lw=3)
plt.title('DPMD', fontsize='x-large')
ax1.tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='r')
ax2.tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='b')
ax1.tick_params(axis="x",direction="in", labelsize='x-large')
ax1.set_ylabel('s(r)', fontsize='x-large', color='r')
ax2.set_ylabel('ds(r)/dr', fontsize='x-large', color='b')
plt.xticks((0, r_c), (0, '$r_{c}$'))
plt.subplots_adjust(right=0.85)
ax1.grid(alpha=0.2)

# DPMD-SE
dpmd_se = DPMD_SE(X, r_c, r_cs)
fig, ax1 = plt.subplots()
ax1.plot(X, dpmd_se, c='r', label='s(x)', lw=3)
# ax1.set_ylim(-0.2, 1.2)
ax2 = ax1.twinx()
ax2.plot(X[:-1] +dx/2., (dpmd_se[1:] - dpmd_se[:-1]) /dx, c='b', label='ds(r)/dr', ls=(0, (5, 5)), lw=3)
plt.title('DPMD-SE', fontsize='x-large')
ax1.tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='r')
ax2.tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='b')
ax1.tick_params(axis="x",direction="in", labelsize='x-large')
ax1.set_ylabel('s(r)', fontsize='x-large', color='r')
ax2.set_ylabel('ds(r)/dr', fontsize='x-large', color='b')
plt.xticks((0, r_cs, r_c), (0, '$r_{cs}$', '$r_{c}$'))
plt.subplots_adjust(right=0.85)
ax1.grid(alpha=0.2)

#
plt.show()
