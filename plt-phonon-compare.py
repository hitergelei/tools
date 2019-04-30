#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def load(file, mul=1):
    line = file.readline()
    line = file.readline()
    line = file.readline()
    x = []
    y = []
    band_x = np.asarray([])
    band_y = np.asarray([])
    a=0
    tick_list = set([0])
    while True:
        while True:
            line = file.readline()
            if not line: break
            lists = line.split()
            if len(lists) == 0:
                tick_list.add(float(a))
                continue
            a, b = lists
            ########## new band line
            if float(a) == 0:
                x.append(band_x)
                y.append(band_y)
                band_x = np.asarray([])
                band_y = np.asarray([])
                band_x = np.append(band_x, float(a))
                band_y = np.append(band_y, float(b) * mul)
            else:
                band_x = np.append(band_x, float(a))
                band_y = np.append(band_y, float(b) * mul)
        if not line: break
    return x, y, tick_list

######### read inputs
band_1_f = open('dft-band.in', 'r')
band_2_f = open('dpmd-band.in', 'r')
# band_1_x, band_1_y, tick_list = load(band_1_f, 4.13567)
band_1_x, band_1_y, tick_list = load(band_1_f)
band_2_x, band_2_y, tick_list = load(band_2_f)
fig, ax = plt.subplots()

########### plot
for i in range(len(band_1_x)):
    if len(band_1_x[i]) is not 0:
        ax.plot(
            band_1_x[i],
            band_1_y[i],
            'C0',
            # color = '#112482',
            )
        # ax.scatter(band_2_x[i][::10], band_2_y[i][::10], s=40, facecolors='none', edgecolors='r')
    ######### making error bar
    error = band_1_y[i] * 0.05
    if len(band_2_x[i]) is not 0:
        ax.errorbar(
            band_2_x[i][::15],
            band_2_y[i][::15],
            yerr              = error[::15],
            fmt               = 'ro',
            # color           = '#eb2d01',
            capsize           = 3,
            elinewidth        = 1,
            markeredgewidth   = 1,
            # markerfacecolor = 'none',
            markersize        = 5,
            )

######### legend option
ax.plot(
    band_1_x[0],
    band_1_y[0],
    'C0',
    # color = '#112482',
    label = 'DFT',
    )
ax.errorbar(
    band_2_x[i][::15],
    band_2_y[i][::15],
    yerr              = error[::15],
    fmt               = 'ro',
    capsize           = 3,
    elinewidth        = 1,
    markeredgewidth   = 1,
    # markerfacecolor = 'none',
    markersize        = 5,
    label             = 'NNP('+r'$\pm$'+'5%)',
    )
ax.legend()

########### plot
ax.grid(color='0.8')
plt.xticks(list(sorted(tick_list)), (r'$\Gamma$', 'X', 'K', r'$\Gamma$', 'L'))
plt.show()
