#!/usr/bin/env python
import numpy as np

# Hyper-parameters
dt = 0.01 # ns
alist_file = 'quench-to-cryst.traj'
img_slice = ':'
bin_edges = np.load('peaks/peak_Te.npz')['pos']

#
from ase.io import read
alist = read(alist_file, img_slice)

# pos_z.shape = (len(alist), len(atoms))
pos_z = np.zeros((len(alist), len(alist[0])))
for i in range(len(alist)):
    pos_z[i] = alist[i].get_positions()[:,2]

# 
chem = np.array(alist[0].get_chemical_symbols())
mask_ge = chem == 'Ge'
mask_sb = chem == 'Sb'

#
hist_ge = np.zeros((len(bin_edges)-1, len(pos_z)))
hist_sb = np.zeros((len(bin_edges)-1, len(pos_z)))
for i in range(len(pos_z)):
    hist_ge[:, i], _ = np.histogram(pos_z[i][mask_ge], bin_edges)
    hist_sb[:, i], _ = np.histogram(pos_z[i][mask_sb], bin_edges)
hist_ge = (hist_ge - np.expand_dims([128, 0, 0, 0, 128] *4, axis=1)) /512.
hist_sb = (hist_sb - np.expand_dims([0, 128, 0, 128, 0] *4, axis=1)) /512.

def plot_figs(
    plt,
    data,
    dt,
    labels,
    colors,
    title,
    ):
    plt.figure()
    for i in range(len(data)):
        plt.plot(np.arange(len(data[i])) *dt, data[i], label=labels[i], c=colors[i])
    plt.xlabel('Time (ns)', fontsize='x-large')
    plt.ylabel('Scaled population', fontsize='x-large')
    plt.legend(fontsize='large').set_draggable(True)
    plt.title(title, fontsize='x-large')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.subplots_adjust(left=0.25, right=0.75, bottom=0.25, top=0.75)
    plt.xlim((0, (len(data[0])-1) *dt))
    plt.ylim((0, None))
    plt.grid(alpha=0.4)

#
from matplotlib import pyplot as plt
# Layer-wise plot
plot_figs(
    plt,
    hist_ge,
    dt,
    ['Floor {}'.format(i) for i in range(len(bin_edges)-1)],
    [None]*(len(bin_edges)-1),
    'Ge',
    )
plot_figs(
    plt,
    hist_sb,
    dt,
    ['Floor {}'.format(i) for i in range(len(bin_edges)-1)],
    [None]*(len(bin_edges)-1),
    'Sb',
    )

# Class by site
class_addr = np.array([3, 4, 0, 1, 2,] *4, dtype=int)
unique_class = np.unique(class_addr)
hist_ge_site = np.zeros((len(unique_class), len(pos_z)))
hist_sb_site = np.zeros((len(unique_class), len(pos_z)))
for i in range(len(unique_class)):
    hist_ge_site[i] = np.sum(hist_ge[class_addr==unique_class[i]], axis=0)
    hist_sb_site[i] = np.sum(hist_sb[class_addr==unique_class[i]], axis=0)
plot_figs(
    plt,
    hist_ge_site,
    dt,
    ['Class {}'.format(i) for i in unique_class],
    [None]*(len(unique_class)),
    'Ge classified by site',
    )
plot_figs(
    plt,
    hist_sb_site,
    dt,
    ['Class {}'.format(i) for i in unique_class],
    [None]*(len(unique_class)),
    'Sb classified by site',
    )

# Class by symmetry
class_addr = np.array([2, 1, 0, 1, 2,] *4, dtype=int)
unique_class = np.unique(class_addr)
hist_ge_sym = np.zeros((len(unique_class), len(pos_z)))
hist_sb_sym = np.zeros((len(unique_class), len(pos_z)))
for i in range(len(unique_class)):
    hist_ge_sym[i] = np.sum(hist_ge[class_addr==unique_class[i]], axis=0)
    hist_sb_sym[i] = np.sum(hist_sb[class_addr==unique_class[i]], axis=0)
plot_figs(
    plt,
    hist_ge_sym,
    dt,
    ['Vacancy site', 'Sb site', 'Ge site'],
    ['k', 'r', 'b'],
    'Ge classified by symmetry',
    )
plot_figs(
    plt,
    hist_sb_sym,
    dt,
    ['Vacancy site', 'Sb site', 'Ge site'],
    ['k', 'r', 'b'],
    'Sb classified by symmetry',
    )

plt.show()
