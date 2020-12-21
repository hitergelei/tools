#!/usr/bin/env python
import numpy as np

## MLP part
# Params
disp_min = -0.25
disp_max = 0.25
disp_nbins = 500
direc = [1, 0, 0]

# Main
dx = (disp_max-disp_min) /disp_nbins
disp = np.array(range(disp_nbins+1)) /(disp_nbins) *(disp_max-disp_min) +disp_min
direc /= np.linalg.norm(direc)
from ase.io import read, write
e = []
f = []
from subprocess import call
from lmp2traj import read_lmp_log
for i in range(len(disp)):
    folder = 'calc-mlp/{}_{}'.format(i, disp[i])
    e.append(read_lmp_log('{}/log.lammps'.format(folder))['PotEng'][0])
    atoms = read('{}/out.dump'.format(folder))
    # e.append(np.sum(atoms.get_potential_energies()))
    f.append(atoms.get_forces()[0])
# Shape = (disp_nbins, 3)
f = np.squeeze(np.dot(np.expand_dims(f, 1), np.expand_dims(direc, [0,2])))
# Shape = (disp_nbins, 3)
e = np.array(e)
dedx = (e[1:] - e[:-1]) /dx

## DFT part
# Params
dft_disp_min = -0.25
dft_disp_max = 0.25
dft_disp_nbins = 8

# Main
dft_dx = (dft_disp_max-dft_disp_min) /dft_disp_nbins
dft_disp = np.array(range(dft_disp_nbins+1)) /(dft_disp_nbins) *(dft_disp_max-dft_disp_min) +dft_disp_min
direc /= np.linalg.norm(direc)
from ase.io import read, write
dft_e = []
dft_f = []
from subprocess import call
for i in range(len(dft_disp)):
    dft_folder = 'calc-dft/{}_{}'.format(i, dft_disp[i])
    dft_atoms = read('{}/vasprun.xml'.format(dft_folder))
    dft_e.append(dft_atoms.get_potential_energy())
    dft_f.append(dft_atoms.get_forces()[0])
# Shape = (disp_nbins, 3)
dft_f = np.squeeze(np.dot(np.expand_dims(dft_f, 1), np.expand_dims(direc, [0,2])))
# Shape = (disp_nbins, 3)
dft_e = np.array(dft_e)
    
#
from matplotlib import pyplot as plt
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(disp, e-np.min(e), label='E(MLP)', c='k')
ax1.scatter(
    dft_disp,
    dft_e-np.min(dft_e),
    label='E(DFT)',
    s=50,
    facecolors='none',
    edgecolors='k',
    )
ax2.plot(disp, np.abs(f), label='F(MLP)', c='r')
ax2.scatter(
    dft_disp,
    np.abs(dft_f),
    label='F(DFT)',
    s=50,
    facecolors='none',
    edgecolors='r',
    )
# ax2.plot(disp[:-1] +dx/2., np.abs(dedx), label='|dE/dx|', c='b', ls=':')

#
ax1.set_xlabel('dx ($\AA$)', fontsize='x-large')
ax1.set_ylabel('Energy (eV)', fontsize='x-large', color='k')
ax2.set_ylabel('Force (eV/$\AA$)', fontsize='x-large', color='r')
ax1.tick_params(axis="x",direction="in", labelsize='x-large')
ax1.tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='k')
ax2.tick_params(axis="y",direction="in", labelsize='x-large',colors='r',  labelcolor='r')
ax1.grid(alpha=0.4)
ax1.legend(bbox_to_anchor=(0.2,1), loc='upper left', fontsize='large')
ax2.legend(bbox_to_anchor=(0.6,1), loc='upper left', fontsize='large')
plt.title('Potential Energy Surface (PES)', fontsize='x-large')
plt.subplots_adjust(left=0.14, bottom=0.13, right=0.86, top=0.93)
plt.show()



