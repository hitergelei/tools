#!/usr/bin/env python
import numpy as np

# Task
# calculate = True
calculate = False
# rotate    = True
rotate    = False

# params
uc_file    = '../GaN-prim.vasp'
ph2_pckl   = 'vasp-x442_d0.010_symTrue-NACTrue-fc2.bin'
ph3_pckl   = '3pho_vasp_fc2x442_fc3x332-nacTrue-qx{}{}{}-rmltc-bte.pckl'
q_range    = range(15,16)
q          = (15,15,15)
T_list     = np.arange(50, 121, 5, dtype=float) # (K)
T          = 30.
ij         = (0,1)
tau        = None
color    = ['r', 'g', 'b', 'c', 'm', 'y']
# color      = ['r', 'g', 'b']
band_group = (range(0,3), range(3,12))
# band_group = None

def rot_alpha(alpha, new_z):
    # Rotate alpha's z-axis to (111)
    new_x = np.array([1.,0.,0.])
    new_y = np.cross(new_z, new_x)
    new_y /= np.linalg.norm(new_y)
    new_x = np.cross(new_y, new_z)
    # R.shape = (1, 1, 3, 3)
    R = np.expand_dims(np.array([new_x, new_y, new_z]).T, axes=[0, 1])
    alpha = np.matmul(np.transpose(R, [0,1,3,2]), alpha)
    alpha = np.matmul(alpha, R)
    return alpha

if calculate:
    from subprocess import call
    for i in q_range:
        cmd = 'pam.py {} {} {}'.format(uc_file, ph2_pckl, ph3_pckl).format(i, i, i)
        if tau is not None:
            cmd += ' -t {}'.format(tau)
        call(cmd, shell=True)

# Load
# alpha_q.shape = (len(q_range), len(sigma), 3, 3)
alpha_q = []
for i in q_range:
    alpha_q.append(np.load('alpha-tau{}-qx{}{}{}-{}K.npy'.format(tau,i,i,i,T)))
alpha_q = np.array(alpha_q)

# alpha_T.shape = (len(T_list), len(sigma), 3, 3)
alpha_T = []
for i in T_list:
    alpha_T.append(np.load('alpha-tau{}-qx{}{}{}-{}K.npy'.format(tau,*q,i)))
alpha_T = np.array(alpha_T)

if rotate:
    from ase.io import read
    cell = read(uc_file).get_cell()
    new_z = np.sum(cell, axis=0)
    new_z /= np.linalg.norm(new_z)
    alpha_q = rot_alpha(alpha_q, new_z)
    alpha_T = rot_alpha(alpha_T, new_z)

from matplotlib import pyplot as plt
plt.plot(q_range, np.sum(alpha_q, axis=1)[:, ij[0], ij[1]])
plt.tick_params(axis="both",direction="in", labelsize='x-large')
plt.xlabel('q-mesh ($q^3$)', fontsize='x-large')
plt.ylabel(r'$\alpha_{{{}{}}}$ $( J S / m^2 K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
plt.legend(fontsize='large')
plt.title(r'At {} K, $\tau$={} ps'.format(T, tau), fontsize='x-large')
plt.xlim(np.min(q_range),np.max(q_range))
# plt.ylim(0,130)
plt.subplots_adjust(left=0.20, bottom=0.20, right=0.80, top=0.80)
plt.grid(alpha=0.5)

print('q={}x{}x{}'.format(*q), '\n', alpha_T)
from matplotlib import pyplot as plt
len_sigma = alpha_T.shape[1]
for s in range(len_sigma //len(color)):
    plt.figure()
    plt.plot(
        T_list,
        np.sum(alpha_T, axis=1)[:, ij[0], ij[1]],
        c='k',
        label='Total',
        )
    for i in range(len(color)):
        plt.plot(
            T_list,
            alpha_T[:, len(color)*s+i, ij[0], ij[1]],
            label=r'$\sigma=${}'.format(len(color)*s+i+1),
            c=color[i],
            )
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xlabel('Temperature (K)', fontsize='x-large')
    plt.ylabel(r'$\alpha_{{{}{}}}$ $( J S / m^2 K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
    plt.legend(fontsize='large').set_draggable(True)
    plt.title(r'q-mesh={}X{}X{}, $\tau$={} ps'.format(*q, tau), fontsize='x-large', pad=20)
    plt.xlim(np.min(T_list),np.max(T_list))
    # plt.ylim(0,130)
    plt.subplots_adjust(left=0.20, bottom=0.20, right=0.80, top=0.80)
    # plt.yscale('log')
    plt.grid(alpha=0.5)
plt.show()

if band_group is not None:
    from matplotlib import pyplot as plt
    plt.plot(
        T_list,
        np.sum(alpha_T, axis=1)[:, ij[0], ij[1]],
        c='k',
        label='Total',
        )
    for i in range(len(band_group)):
        plt.plot(
            T_list,
            np.sum(alpha_T[:, band_group[i], ij[0], ij[1]], axis=1),
            label='Group {}\n{}'.format(i+1, band_group[i]),
            c=color[i],
            )
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xlabel('Temperature (K)', fontsize='x-large')
    plt.ylabel(r'$\alpha_{{{}{}}}$ $( J S / m^2 K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
    plt.legend(fontsize='large').set_draggable(True)
    plt.title(r'q-mesh={}X{}X{}, $\tau$={} ps'.format(*q, tau), fontsize='x-large', pad=20)
    plt.xlim(np.min(T_list),np.max(T_list))
    # plt.ylim(0,130)
    plt.subplots_adjust(left=0.18, bottom=0.20, right=0.88, top=0.80)
    # plt.yscale('log')
    plt.grid(alpha=0.5)
    plt.show()
