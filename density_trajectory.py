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

## post-process
volume = 4./3 * np.pi * radius **3
num_density = counts / volume
eff_vol = 1./ num_density
norm_eff_vol = eff_vol / eff_vol[0]
rel_norm_eff_vol = norm_eff_vol - 1.

#### Plot
temp_arr = np.arange(temp_ini, temp_fin, (temp_fin-temp_ini) / len(alist))
avg_rel_norm_eff_vol = rel_norm_eff_vol.copy()
avg_number = 40
for i in range(1,1+avg_number):
    for j in range(avg_number, len(temp_arr)-avg_number):
        avg_rel_norm_eff_vol[j] += rel_norm_eff_vol[j+i] + rel_norm_eff_vol[j-i]
avg_rel_norm_eff_vol /= float(1+avg_number*2)

from matplotlib import pyplot as plt
plt.plot(temp_arr, rel_norm_eff_vol *100., c='0.7')
plt.plot(temp_arr[avg_number:-avg_number], avg_rel_norm_eff_vol[avg_number:-avg_number] *100., c='r', label='$\pm${} average'.format(avg_number))
plt.xlabel('Temperature (K)', fontsize='x-large')
plt.ylabel('$\Delta V/ V_0$ (%)', fontsize='x-large')
plt.tick_params(axis="both",direction="in", labelsize='x-large')
plt.grid(alpha=0.2)
plt.axvline(850., ls='--', linewidth=1, c='k')
plt.axvline(900., ls='--', linewidth=1, c='k')
plt.legend(fontsize='large')
plt.show()


