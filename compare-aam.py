#!/usr/bin/env python
import numpy as np

tps = ['10ps', '100ps', '1ns', '10ns', '100ns']
paths = ['{}/phonon/rvdos-from0.375to0.775THz-{}.vasp.traj'.format(t, t) for t in tps]
colors = ['r', 'g', 'b', 'm', 'c']

from ase.io import read

hist = []
mids = []
aam_norm = []
for i in range(len(paths)):
    atoms = read(paths[i], 0)
    aam_norm_i = np.linalg.norm(atoms.calc.results['aam'], axis=-1)
    aam_norm.append(aam_norm_i)
    aam_norm_log_i = np.log(aam_norm_i)
    hist_log, bins_log = np.histogram(aam_norm_log_i, 50)
    mids_log = (bins_log[1:] + bins_log[:-1]) /2.
    hist.append(hist_log)
    mids.append(mids_log)

from matplotlib import pyplot as plt
for i in range(len(hist)):
    plt.plot(np.exp(mids[i]), hist[i], colors[i], label=tps[i])
    plt.axvline(np.mean(aam_norm[i]), c=colors[i], lw=1)
    plt.xscale('log')
    plt.xlabel('Atomic angular momentum (eV ps)', fontsize='x-large')
    plt.ylabel('Population', fontsize='x-large')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.subplots_adjust(left=0.20, bottom=0.20, right=0.80, top=0.80, wspace=0.2, hspace=0.2)
    plt.legend(fontsize='large').set_draggable(True)
    plt.grid(alpha=0.5)
plt.show()
