#!/usr/bin/env python
import numpy as np

#### Global params
radius = float(8) # (Ang)
dt = float(1) # (ps)
temp_ini = float(500) # (K)
temp_fin = float(1100)

#### Main
from ase.io import read
alist = read('lmp-result.traj', ':')
atoms = alist[0]
cell = atoms.get_cell()
center = np.sum(cell, axis=0) /2.

## Count atoms inside cutoff
counts = np.zeros(len(alist))
for i in range(len(alist)):
    counts[i] = np.sum(np.linalg.norm(alist[i].get_positions() - center, axis=-1) <= radius)

# Potential energy
epot = []
for atoms in alist:
    epot.append(atoms.get_potential_energy())
epot = np.array(epot) /len(atoms)

## post-process
volume = 4./3 * np.pi * radius **3
num_density = counts / volume
eff_vol = 1./ num_density
norm_eff_vol = eff_vol / eff_vol[0]
rel_norm_eff_vol = norm_eff_vol - 1.

#### Plot
temp_arr = np.arange(temp_ini, temp_fin, (temp_fin-temp_ini) / len(alist))
avg_rel_norm_eff_vol = rel_norm_eff_vol.copy()
# Averaging
avg_number = 40
for i in range(1,1+avg_number):
    for j in range(avg_number, len(temp_arr)-avg_number):
        avg_rel_norm_eff_vol[j] += rel_norm_eff_vol[j+i] + rel_norm_eff_vol[j-i]
avg_rel_norm_eff_vol /= float(1+avg_number*2)
# # Gaussian smearing
# from scipy.ndimage.filters import gaussian_filter1d
# gsmear = 20
# avg_rel_norm_eff_vol = gaussian_filter1d(avg_rel_norm_eff_vol, gsmear)

from matplotlib import pyplot as plt
font = {'family':'Arial'}
plt.rc('font', **font)
fig, ax1 = plt.subplots()
ax1.plot(temp_arr, rel_norm_eff_vol *100., c='0.7')
ax1.plot(temp_arr[avg_number:-avg_number], avg_rel_norm_eff_vol[avg_number:-avg_number] *100.,c='k', label='$\pm${} average'.format(avg_number))
# ax1.plot(temp_arr, avg_rel_norm_eff_vol *100., c='r', label='Gaussian smearing: {}'.format(gsmear))
ax1.set_xlabel('Temperature (K)', fontsize='x-large')
ax1.set_ylabel('Volume expansion $\Delta V/ V_0$ (%)', fontsize='x-large')
ax1.tick_params(axis="both",direction="in", labelsize='x-large')
ax2 = ax1.twinx()
ax2.plot(temp_arr, epot, c='r')
ax2.set_ylabel('Energy per atom (eV)', fontsize='x-large', c='r')
ax2.tick_params(axis="y",direction="in", labelsize='x-large',colors='r',  labelcolor='r')
plt.grid(alpha=0.2)
plt.axvline(850., ls='--', linewidth=1, c='k')
plt.axvline(910., ls='--', linewidth=1, c='k')
plt.xlim(np.amin(temp_arr), np.amax(temp_arr))
plt.subplots_adjust(right=0.85)
plt.legend(fontsize='large')
plt.legend().remove()
plt.show()


