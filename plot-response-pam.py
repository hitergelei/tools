#!/usr/bin/env python
import numpy as np
from ase.io import read

# PARAMS
# calculate = True
calculate = False
#
# rotate    = True
rotate    = False
#
uc_file    = '../BN-mono-prim.vasp'
ph2_pckl   = '../vasp-x441_d0.010_symTrue-NACTrue-fc2.bin'
ph3_pckl   = '3pho_vasp_fc2x441_fc3x331-nacTrue-qx{}{}{}-rmltc-rta.pckl'
#
dim2  = True
new_z = None
new_x = np.array([np.cos(1./3.*np.pi), np.sin(1./3.*np.pi), 0.])
#
# dim2  = False
# new_x = None
# cell  = read(uc_file).get_cell()
# new_z = np.sum(cell, axis = 0)

# params
q_range    = range(60,65,5)
q          = (60,60,1)
T_list     = list(np.arange(30, 63, 3, dtype=float)) + [100., 120., 200., 300.] # (K)
# T_list     = list(np.arange(5, 35, 5, dtype=float)) # (K)
T          = 30.
ij         = (2,0)
const_tau  = False
# const_tau  = True
# color    = ['r', 'g', 'b', 'c', 'm', 'y']
color      = ['r', 'g', 'b']
band_group = (range(0,3), range(3,6))
# band_group = None

def rot_alpha(alpha, new_z):
    # Rotate alpha's z-axis to (new_z)
    new_z /= np.linalg.norm(new_z)
    new_x = np.array([1.,0.,0.])
    new_y = np.cross(new_z, new_x)
    new_y /= np.linalg.norm(new_y)
    new_x = np.cross(new_y, new_z)
    # R.shape = (1, 1, 3, 3)
    R = np.expand_dims(np.array([new_x, new_y, new_z]).T, axis=[0, 1])
    alpha = np.matmul(np.transpose(R, [0,1,3,2]), alpha)
    alpha = np.matmul(alpha, R)
    return alpha

def rot_alpha_xy(alpha, new_x):
    # Rotate alpha's x-axis to (new_x)
    new_x /= np.linalg.norm(new_x)
    new_z = np.array([0.,0.,1.])
    new_y = np.cross(new_z, new_x)
    new_y /= np.linalg.norm(new_y)
    new_z = np.cross(new_x, new_y)
    # R.shape = (1, 1, 3, 3)
    R = np.expand_dims(np.array([new_x, new_y, new_z]).T, axis=[0, 1])
    alpha = np.matmul(np.transpose(R, [0,1,3,2]), alpha)
    alpha = np.matmul(alpha, R)
    # np.save('R.npy', R[0,0])
    return alpha

if const_tau:
    tau = 1e12
else:
    tau = None

if calculate:
    from subprocess import call
    for i in q_range:
        if dim2:
            cmd = 'pam.py {} {} {} --dim2'.format(uc_file, ph2_pckl, ph3_pckl).format(i, i, 1)
        else:
            cmd = 'pam.py {} {} {}'.format(uc_file, ph2_pckl, ph3_pckl).format(i, i, i)
        if tau is not None:
            cmd += ' -t {}'.format(tau)
        call(cmd, shell=True)
    if dim2:
        cmd = 'pam.py {} {} {} --dim2'.format(uc_file, ph2_pckl, ph3_pckl).format(*q)
    else:
        cmd = 'pam.py {} {} {}'.format(uc_file, ph2_pckl, ph3_pckl).format(*q)
    if tau is not None:
        cmd += ' -t {}'.format(tau)
    call(cmd, shell=True)

# Load
# alpha_q.shape = (len(q_range), len(sigma), 3, 3)
alpha_q = []
for i in q_range:
    if dim2:
        alpha_q.append(np.load('alpha-2d-tau{}-qx{}{}{}-{}K.npy'.format(tau,i,i,1,T)))
    else:
        alpha_q.append(np.load('alpha-tau{}-qx{}{}{}-{}K.npy'.format(tau,i,i,i,T)))
alpha_q = np.array(alpha_q)

# alpha_T.shape = (len(T_list), len(sigma), 3, 3)
alpha_T = []
for i in T_list:
    if dim2:
        alpha_T.append(np.load('alpha-2d-tau{}-qx{}{}{}-{}K.npy'.format(tau,*q,i)))
    else:
        alpha_T.append(np.load('alpha-tau{}-qx{}{}{}-{}K.npy'.format(tau,*q,i)))
alpha_T = np.array(alpha_T)

if const_tau:
    alpha_T[0] = 0.
    alpha_q[0] = 0.

if rotate:
    if new_x is not None or False:
        alpha_q = rot_alpha_xy(alpha_q, new_x)
        alpha_T = rot_alpha_xy(alpha_T, new_x)
        if new_z is not None or False:
            raise ValueError('Both new_x and new_z are provided.')
    elif new_z is not None or False:
        alpha_q = rot_alpha(alpha_q, new_z)
        alpha_T = rot_alpha(alpha_T, new_z)
    else:
        raise ValueError('None of new_x or new_z are provided.')

if len(alpha_q) > 0:
    from matplotlib import pyplot as plt
    plt.plot(q_range, np.sum(alpha_q, axis=1)[:, ij[0], ij[1]])
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xlabel('q-mesh ($q^3$)', fontsize='x-large')
    if const_tau:
        if dim2:
            plt.ylabel(r'$\alpha_{{{}{}}}$/$\tau$ $\rm ( J / m K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
        else:
            plt.ylabel(r'$\alpha_{{{}{}}}$/$\tau$ $\rm ( J / m^2 K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
    else:
        if dim2:
            plt.ylabel(r'$\alpha_{{{}{}}}$ $\rm ( J s / m K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
        else:
            plt.ylabel(r'$\alpha_{{{}{}}}$ $\rm ( J s / m^2 K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
    plt.title(r'At {} K'.format(T), fontsize='x-large')
    plt.legend(fontsize='large')
    plt.xlim(np.min(q_range),np.max(q_range))
    # plt.ylim(0,130)
    plt.subplots_adjust(left=0.20, bottom=0.20, right=0.80, top=0.80)
    plt.grid(alpha=0.5)

print('q={}x{}x{}'.format(*q), '\n', np.real(np.sum(alpha_T, axis=1)))
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
    if const_tau:
        if dim2:
            plt.ylabel(r'$\alpha_{{{}{}}}$/$\tau$ $\rm ( J / m K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
        else:
            plt.ylabel(r'$\alpha_{{{}{}}}$/$\tau$ $\rm ( J / m^2 K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
    else:
        if dim2:
            plt.ylabel(r'$\alpha_{{{}{}}}$ $\rm ( J s / m K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
        else:
            plt.ylabel(r'$\alpha_{{{}{}}}$ $\rm ( J s / m^2 K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
    plt.title(r'q-mesh={}X{}X{}'.format(*q), fontsize='x-large', pad=20)
    plt.legend(fontsize='large').set_draggable(True)
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
            label='Band {}-{}'.format(band_group[i][0]+1, band_group[i][-1]+1),
            c=color[i],
            )
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xlabel('Temperature (K)', fontsize='x-large')
    if const_tau:
        if dim2:
            plt.ylabel(r'$\alpha_{{{}{}}}$/$\tau$ $\rm ( J / m K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
        else:
            plt.ylabel(r'$\alpha_{{{}{}}}$/$\tau$ $\rm ( J / m^2 K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
    else:
        if dim2:
            plt.ylabel(r'$\alpha_{{{}{}}}$ $\rm ( J s / m K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
        else:
            plt.ylabel(r'$\alpha_{{{}{}}}$ $\rm ( J s / m^2 K )$'.format(ij[0]+1, ij[1]+1), fontsize='x-large')
    plt.title(r'q-mesh={}X{}X{}'.format(*q), fontsize='x-large', pad=20)
    plt.legend(fontsize='large').set_draggable(True)
    plt.xlim(np.min(T_list),np.max(T_list))
    # plt.ylim(0,130)
    plt.subplots_adjust(left=0.18, bottom=0.20, right=0.88, top=0.80)
    # plt.yscale('log')
    plt.grid(alpha=0.5)
    plt.show()
