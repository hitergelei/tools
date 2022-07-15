#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    DESCRIPTION
    """)
    # Positional arguments
    parser.add_argument('band_1_f', help='First band to compare (phonopy format).')
    parser.add_argument('band_2_f', help='Second band to compare (phonopy format).')
    # Optional arguments
    parser.add_argument('-t', '--tick_labels', type=str, nargs='+', default=None, help='Tick labels for band plot')
    parser.add_argument('-f', '--unit', type=str, choices=['THz', 'meV'], default='THz', help='Frequency unit to plot (THz or meV) [default: THz].')
    parser.add_argument('-b', '--labels', type=str, nargs=2, default=['Band 1', 'Band 2'], help='Labels for band 1 and 2 [default: Band 1, Band 2].')
    parser.add_argument('-u', '--upper_lim', type=float, default=None, help='Upper limit for frequency axis [default: automatic]')
    parser.add_argument('-l', '--lower_lim', type=float, default=None, help='Lower limit for frequency axis [default: automatic]')
    return parser.parse_args()

def load_data(file):
    # Dump first line
    lines = file.readlines()
    # 
    tick_list = [float(t) for t in lines[1][1:-1].split()]
    print('tick list: {}'.format(tick_list))

    # Extract data
    x = []
    y = []
    # Band iteration iteration
    line_i = 1
    while True:
        x_i = []
        y_i = []
        # tick interation
        for j in range(len(tick_list)-1):
            x_ij = []
            y_ij = []
            while True:
                line_i += 1
                line = lines[line_i][:-1].split()
                if len(line) == 0:
                    break
                else:
                    x_ij.append(float(line[0]))
                    y_ij.append(float(line[1]))
            x_i.append(x_ij)
            y_i.append(y_ij)
        line_i += 1
        x.append(x_i)
        y.append(y_i)
        if line_i >= len(lines)-1:
            break
    return np.array(tick_list, dtype=float), np.array(x, dtype=float), np.array(y, dtype=float)

if __name__ == '__main__':
    ## Intro
    import datetime
    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ POSTECH <<<<<'.center(120))
    print(('Code runtime : '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('DESCRIPTION'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    # Global variables
    band_1_f = open(args.band_1_f, 'r')
    band_2_f = open(args.band_2_f, 'r')
    unit = args.unit
    ylim_low = args.lower_lim
    ylim_up = args.upper_lim
    tick_labels = args.tick_labels

    ## read inputs
    tick_list_1, band_1_x, band_1_E = load_data(band_1_f)
    tick_list_2, band_2_x, band_2_E = load_data(band_2_f)
    tick_list = tick_list_1
    if not (tick_list_1 == tick_list_2).all():
        print('WARNING) Ticks for band 1 and 2 are different')

    ## Scale energy
    if unit == 'meV':
        band_1_E = band_1_E * 4.13567
        band_2_E = band_2_E * 4.13567

    import matplotlib.pyplot as plt
    font = {'family':'Arial'}
    plt.rc('font', **font)

    ########### plot
    for i in range(len(band_1_x)):
        for j in range(len(tick_list)-1):
            if i == 0 and j == 0:
                labels=args.labels
            else:
                labels=[None, None]

            plt.plot(
                band_1_x[i, j],
                band_1_E[i, j],
                c='k',
                label=labels[0],
                )

            plt.plot(
                band_2_x[i, j],
                band_2_E[i, j],
                ls = (0, (3,3)),
                c = 'r',
                label=labels[1],
                )

    ######### legend option
    plt.legend(fontsize='large').set_draggable(True)

    ## Normalized Distance
    err_std = np.std(band_1_E - band_2_E)
    print('Standard deviation of frequency difference = {:.3e} {}'.format(err_std, unit))
    

    ########### plot
    plt.grid()
    if tick_labels:
        plt.xticks(list(tick_list), tick_labels)
    else:
        plt.xticks(list(tick_list), ['' for i in range(len(tick_list))])
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xlabel('Frequency RMSE = {:.3e} ({})'.format(err_std, unit), fontsize='large')
    plt.xlim(tick_list[0], tick_list[-1])
    plt.ylim(ylim_low, ylim_up)
    plt.ylabel('Frequency ({})'.format(unit), fontsize='x-large')
    plt.subplots_adjust(left=0.25, bottom=0.25, right=0.75, top=0.75, wspace=0.2, hspace=0.2)
    plt.show()

