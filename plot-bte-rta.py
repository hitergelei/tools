#!/usr/bin/env python
import numpy as np

methods = ['rta', 'bte']
plot_range = range(1,16)
xlim = None
# xlim = (0,20)
ylim = None
# ylim = (0,130)

data = {}
for method in methods:
    data[method] = []
    for i in plot_range:
        q_mesh = (i,i,i)
        import h5py
        with h5py.File('{}/kappa-m{}{}{}.hdf5'.format(method, *q_mesh), 'r') as f:
            data[method].append(np.mean(f['kappa'][0][0:3]))

from matplotlib import pyplot as plt
i = list(plot_range)
plt.plot(i, data['bte'], label='BTE', c='k')
plt.plot(i, data['rta'], label='RTA', c='r')
plt.tick_params(axis="both",direction="in", labelsize='x-large')
plt.xlabel('q-mesh ($x^3$)', fontsize='x-large')
plt.ylabel('LTC (W/mK)', fontsize='x-large')
plt.legend(fontsize='large')
plt.xlim(xlim)
plt.ylim(ylim)
plt.subplots_adjust(left=0.20, bottom=0.20, right=0.80, top=0.80)
plt.grid(alpha=0.5)
plt.show()
