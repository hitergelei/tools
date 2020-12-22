import numpy as np

from ase.io import read
atoms = read('vasprun.xml', -1)

d = np.sort(atoms.get_all_distances(mic=True)[0])
f = np.sort(np.linalg.norm(atoms.get_forces(), axis=-1))[::-1]

from matplotlib import pyplot as plt
plt.plot(d, f, c='r')
plt.yscale('log')
plt.tick_params(axis="both",direction="in", labelsize='x-large')
plt.ylabel('Force, |$F_i$| (eV/$\AA$)', fontsize='x-large')
plt.xlabel('Distance, $r_i-r_0$ ($\AA$)', fontsize='x-large')
plt.grid(alpha=0.5)
plt.show()
