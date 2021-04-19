#!/usr/bin/env python
import numpy as np

insts = ['kappa', 'kappa_RTA']
q_range = range(1,15)

data = {}
for inst in insts:
    data[inst] = []
    for i in q_range:
        temp = []
        q_mesh = (i,i,i)
        import h5py
        with h5py.File('kappa-m{}{}{}.hdf5'.format(*q_mesh), 'r') as f:
            data[inst].append(np.array(f[inst]))
            temp = np.array(f['temperature'])
    data[inst] = np.array(data[inst])


from matplotlib import pyplot as plt
for i in range(len(temp)):
    #
    q = list(q_range)
    plt.figure()
    plt.plot(q, np.mean(data['kappa'][:,i,:3], axis=1), label='BTE', c='k')
    plt.plot(q, np.mean(data['kappa_RTA'][:,i,:3], axis=1), label='RTA', c='r')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xlabel('q-mesh ($x^3$)', fontsize='x-large')
    plt.ylabel('LTC (W/mK)', fontsize='x-large')
    plt.legend(fontsize='large')
    plt.title('At {} K'.format(temp[i]), fontsize='x-large')
    plt.xlim(np.min(q),np.max(q))
    plt.ylim(0,130)
    plt.subplots_adjust(left=0.20, bottom=0.20, right=0.80, top=0.80)
    plt.grid(alpha=0.5)
plt.show()
