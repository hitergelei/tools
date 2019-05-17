#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Global variables
band_1_f = open('dft-band.in', 'r')
band_2_f = open('dpmd-band.in', 'r')
unit = 'THz'
scatter_num = 40
# errorbar_frac = 0.00
# errorbar_num = 20
legend_bool = False
scatter_size = 10
ylim_low = None
ylim_up = None

def load_data(file):
    # Dump first three lines
    line = file.readline()
    line = file.readline()
    line = file.readline()
    # Initialize variables
    x_list = []
    y_list = []
    piece_x = []
    piece_y = []
    x=0.
    tick_set = set([0.])
    # Readline iteration
    while True:
        line = file.readline()
        if not line:
            x_list.append(piece_x)
            y_list.append(piece_y)
            break
        lists = line.split()
        # Add tick
        if len(lists) == 0:
            tick_set.add(float(x))
        else:
            x, y = lists
            # New band line
            if float(x) == 0 and len(piece_x) != 0:
                x_list.append(piece_x)
                y_list.append(piece_y)
                # Reinitialize
                piece_x = []
                piece_y = []
            piece_x.append(np.float(x))
            piece_y.append(np.float(y))
    return np.array(x_list), np.array(y_list), tick_set

def get_scatter_tick(x_arr, scatter_num):
    interval = x_arr[-1] / scatter_num
    ideal_p_list = np.arange(scatter_num) * interval + interval/2
    real_p_list = []
    for ideal_p_i in ideal_p_list:
        real_p_list.append(x_arr[np.abs(x_arr - ideal_p_i).argmin()])
    return real_p_list

def get_scatter_bool_arr(x_arr, scatter_tick_list):
    bool_i_list = []
    for scatter_tick_i in scatter_tick_list:
        bool_i_list.append(x_arr == scatter_tick_i)
    return np.sum(bool_i_list, axis=0).astype(np.bool)


if __name__ == '__main__':
    ## read inputs
    # band_1_x, band_1_E, tick_set = load_data(band_1_f, 4.13567)
    band_1_x, band_1_E, tick_set = load_data(band_1_f)
    band_2_x, band_2_E, tick_set = load_data(band_2_f)
    ## Scale energy
    if unit == 'meV':
        band_1_E = band_1_E * 4.13567
        band_2_E = band_2_E * 4.13567
    elif unit != 'THz':
        raise ValueError('Provided unit is wrong')

    fig, ax = plt.subplots()
    scatter_tick_list = get_scatter_tick(band_1_x[0], scatter_num)
    scatter_bool_arr = get_scatter_bool_arr(band_1_x[0], scatter_tick_list)

    ########### plot
    for i in range(len(band_1_x)):
        ax.plot(
            band_1_x[i],
            band_1_E[i],
            'C0',
            # color = '#112482',
            )
    for i in range(len(band_1_x)):
        ax.scatter(
            band_2_x[i][scatter_bool_arr],
            band_2_E[i][scatter_bool_arr],
            s          = scatter_size,
            # facecolors = 'none',
            facecolors = 'r',
            edgecolors = 'r',
            # alpha      = 0.5,
            zorder     = 10,
            )
        # making error bar
        # error = band_1_E[i] * errorbar_frac
        # ax.plot(
            # band_2_x[i],
            # band_2_E[i],
            # 'C1',
            # # color = '#112482',
            # )
        # ax.errorbar(
            # band_2_x[i][ebar_bool_arr],
            # band_2_E[i][ebar_bool_arr],
            # yerr              = error[ebar_bool_arr],
            # fmt               = 'ro',
            # # color           = '#eb2d01',
            # capsize           = 0,
            # elinewidth        = 1,
            # markeredgewidth   = 1,
            # # markerfacecolor = 'none',
            # markersize        = markersize,
            # )

    ######### legend option
    if legend_bool:
        ax.plot(
            band_1_x[0],
            band_1_E[0],
            'C0',
            # color = '#112482',
            label = 'DFT',
            )
        ax.scatter(
            band_2_x[0][scatter_bool_arr][0],
            band_2_E[0][scatter_bool_arr][0],
            s          = scatter_size,
            facecolors = 'r',
            edgecolors = 'r',
            label      = 'ML-Pot',
            zorder     = 10,
            )
        # ax.errorbar(
            # band_2_x[i][ebar_bool_arr],
            # band_2_E[i][ebar_bool_arr],
            # yerr              = 0,
            # fmt               = 'ro',
            # capsize           = 0,
            # elinewidth        = 0,
            # markeredgewidth   = 0,
            # # markerfacecolor = 'none',
            # markersize        = markersize,
            # label             = 'ML-Pot(5%)',
            # )
        ax.legend(loc='upper right', fontsize='large')

    ## Normalized Distance
    err_std = np.std(band_1_E - band_2_E)
    print('Root Mean Square Error = {:.3f} {}'.format(err_std, unit))
    

    ########### plot
    ax.grid(color='0.8')
    plt.xticks(list(sorted(tick_set)), ('L', '$\Gamma$', 'X', 'U|K', '$\Gamma$'), fontsize='xx-large')
    # plt.xticks(list(sorted(tick_set)), ('K', '$\Gamma$', 'M'), fontsize='xx-large')
    plt.yticks(fontsize='xx-large')
    plt.xlabel('Frequency RMSE = {:.3f} ({})'.format(err_std, unit), fontsize='large')
    plt.xlim(0., np.amax(band_1_x))
    plt.ylim(ylim_low, ylim_up)
    plt.ylabel('Frequency({})'.format(unit), fontsize='xx-large')
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.85, top=0.95, wspace=0.2, hspace=0.2)
    plt.show()

