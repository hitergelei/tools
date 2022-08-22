#!/usr/bin/env python
import numpy as np

#
t = [0.01, 0.1, 1, 10, 100]
b_iso = [
    [1.78301662, 1.74180694, 1.72566072, 1.71975117, 1.69731809],
    [1.39936596, 1.38654058, 1.35694338, 1.35475313, 1.33475433],
    [1.58471990, 1.48531105, 1.41507355, 1.35524102, 1.29332159],
    ]
b_iso_std = [
    [0.04073385, 0.03435311, 0.03383517, 0.03565933, 0.03410886],
    [0.03522830, 0.03212352, 0.02969538, 0.03298201, 0.02851918],
    [0.03756802, 0.02545169, 0.02558674, 0.02891549, 0.02343095],
    ]
spec = ['Ge', 'Sb', 'Te']
color_list = ['b', 'r', 'g']

#
b_iso = np.array(b_iso)
b_iso_std = np.array(b_iso_std)

from matplotlib import pyplot as plt
for i in range(b_iso.shape[0]):
    plt.errorbar(t, b_iso[i], yerr=b_iso_std[i], fmt='s', c=color_list[i], mfc='none', ecolor=color_list[i], capsize=3)
    plt.plot(t, b_iso[i], c=color_list[i], label=spec[i])
plt.xscale('log')
plt.xlabel('Time (ns)', fontsize='x-large')
plt.ylabel(r'$B_{iso}\ (\AA^2)$', fontsize='x-large')
plt.tick_params(axis="both",direction="in", labelsize='x-large')
plt.subplots_adjust(left=0.11, bottom=0.25, right=0.98, top=0.75, wspace=0.2, hspace=0.2)
plt.legend(fontsize='large').set_draggable(True)
plt.grid(alpha=0.5)
plt.show()
