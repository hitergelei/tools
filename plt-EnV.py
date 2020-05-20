#!/usr/bin/env python
import numpy as np

#### Global params
dt = float(1) # (ps)

#### Main
from ase.io import read
alist = read('lmp-result.traj', '::100')
atoms = alist[0]

# Get volume and energy
epot = []
vol  = []
for atoms in alist:
    epot.append(atoms.get_potential_energy())
    vol.append(np.linalg.det(atoms.get_cell()))
epot = np.array(epot) /len(atoms)
vol = np.array(vol) /len(atoms)

#### Plot
t = np.arange(len(alist), dtype=float) *dt
from matplotlib import pyplot as plt
font = {'family':'Arial'}
plt.rc('font', **font)
fig, ax1 = plt.subplots()
ax1.plot(t, 1./vol, c='k')
ax1.set_xlabel('Time (ps)', fontsize='x-large')
ax1.set_ylabel('Number density (atoms/$\AA^3$)', fontsize='x-large')
# ax1.set_ylim(None,28.6)
ax1.tick_params(axis="both",direction="in", labelsize='x-large')
ax2 = ax1.twinx()
ax2.plot(t, epot, c='r')
ax2.set_ylabel('Energy per atom (eV)', fontsize='x-large', c='r')
# ax2.set_ylim(None,-4.0)
ax2.tick_params(axis="y",direction="in", labelsize='x-large',colors='r',  labelcolor='r')
plt.axvline(100, ls='--', linewidth=1, c='k')
plt.grid(alpha=0.2)
plt.xlim(np.amin(t), np.amax(t))
plt.subplots_adjust(right=0.85)
plt.legend(fontsize='large')
plt.legend().remove()
plt.show()


