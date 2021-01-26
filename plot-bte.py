#!/usr/bin/env python
import numpy as np

method = 'bte'
q_range = range(1,16)
ymax    = 3.1

data = []
for i in q_range:
    temp = []
    q_mesh = (i,i,i)
    import h5py
    with h5py.File('{}/kappa-m{}{}{}.hdf5'.format(method, *q_mesh), 'r') as f:
        data.append(np.array(f['kappa']))
        temp = np.array(f['temperature'])
data = np.array(data)


from matplotlib import pyplot as plt
for i in range(len(temp)):
    q = list(q_range)
    plt.figure()
    plt.plot(q, data[:,i,0], label='$\kappa_x$={:.2f}(W/mK)'.format(data[-1,i,0]), c='r')
    plt.plot(q, data[:,i,1], label='$\kappa_y$={:.2f}(W/mK)'.format(data[-1,i,1]), c='g')
    plt.plot(q, data[:,i,2], label='$\kappa_z$={:.2f}(W/mK)'.format(data[-1,i,2]), c='b')
    plt.plot(q, np.mean(data[:,i,:3], axis=1), label='$\kappa_{{avg}}$={:.2f}(W/mK)'.format(np.mean(data[-1,i,:3])), c='k')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xlabel('q-mesh ($x^3$)', fontsize='x-large')
    plt.ylabel('LTC (W/mK)', fontsize='x-large')
    plt.legend(fontsize='large')
    plt.title('{} at {} K'.format(method, temp[i]), fontsize='x-large')
    plt.xlim(np.min(q),np.max(q))
    plt.ylim(top=ymax)
    plt.subplots_adjust(left=0.20, bottom=0.20, right=0.80, top=0.80)
    plt.grid(alpha=0.5)
plt.show()
