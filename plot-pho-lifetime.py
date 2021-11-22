#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Plot phonon lifetime. (Scatter & average line)
    """)
    # Positional arguments
    parser.add_argument('phono3py_pckl', type=str, help='Phono3py class object saved in pickle format.')
    # # Optional arguments
    parser.add_argument('-b', '--nbins', type=int, default=100, help='Number of bins for average plot. [Default: 100]')
    parser.add_argument('-t', '--temperature', type=float, nargs='+',
        help='Set temperature manually. The pckl file must include data at the temperatures. Multiple temperature can be set.')
    parser.add_argument('-c', '--color_list', type=str, nargs='+', help='Color list for the plot. [pyplot convention]')
    parser.add_argument('-x', '--xlim', type=float, nargs=2, help='Provide xmin and xmax [e.g. -x 1 1e4]')
    parser.add_argument('-y', '--ylim', type=float, nargs=2, help='Provide ymin and ymax [e.g. -x 0 50]')

    return parser.parse_args()

if __name__ == '__main__':
    ## Intro
    import datetime
    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ Phys. Dep. of POSTECH in Korea <<<<<'.center(120))
    print(('Code runtime : '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('Plot phonon lifetime. (Scatter & average line)'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()
    nbins = args.nbins

    # @ Main
    import pickle as pckl
    tc = pckl.load(open(args.phono3py_pckl, 'rb')).thermal_conductivity

    # T.shape = (len(T))
    T = tc.get_temperatures()
    if args.temperature:
        T_mask = []
        for j in range(len(T)):
            T_mask.append(T[j] in args.temperature)
        T = T[T_mask]
    # w.shape = (len(ir_q), len(sigma))
    w = tc.get_frequencies()
    w_flat = w.flatten()
    nbands = w.shape[1]
    # ir_gamma.shape = (len(T), len(ir_q), len(sigma))
    ir_gamma = np.array(tc.get_gamma()[0])
    if args.temperature:
        ir_gamma = ir_gamma[T_mask]
    ir_gamma_flat = ir_gamma.reshape(len(T), -1)
    # c.shape = (len(ir_q))
    c = tc.get_grid_weights()
    c_rep = np.repeat(c, w.shape[1])
    # tau.shape = (len(T), len(ir_q), len(sigma))
    tau = 1. / np.where(ir_gamma > 0, ir_gamma, np.inf) / (2 * 2 * np.pi)
    tau_flat = tau.reshape(len(T), -1)

    # @ Averaging
    w_max = np.max(w)
    bin_width = w_max /nbins
    # bins.shape = (nbins+1)
    bins = np.arange(0, w_max+bin_width/2., bin_width, dtype=float)
    # bin_centers.shape = (nbins)
    bin_centers = bins[:-1] + bin_width/2.

    #
    cls = np.digitize(w_flat, bins) -1

    tau_avg = []
    for j in range(len(T)):
        tau_avg.append([])
        for i in range(nbins):
            mask = cls == i
            tau_avg[j].append(np.sum(tau_flat[j, mask] * c_rep[mask]) / np.sum(c_rep[mask]))
    tau_avg = np.array(tau_avg)

    # Table
    w_sigma_mean = []
    tau_sigma_mean = []
    from ss_util import bcolors
    for i in range(nbands):
        w_sigma_mean.append(np.sum(w[:, i] * c) / np.sum(c))
        tau_sigma_mean.append([])
        print('\n{}Band # = {}, mean[w]= {:10.4f} (THz){}'.format(bcolors.okblue, i+1, w_sigma_mean[i], bcolors.endc))
        for j in range(len(T)):
            tau_sigma_mean[i].append(np.sum(tau[j, :, i] * c) / np.sum(c))
            print('mean[tau({}K)]= {:10.4f} (ps)'.format(T[j], tau_sigma_mean[i][j]))
        print('{} * tau({}K) / tau({}K) = {:6.3f}{}'.format(bcolors.okgreen, T[0], T[-1], tau_sigma_mean[i][0] / tau_sigma_mean[i][-1], bcolors.endc))
    tau_sigma_mean = np.array(tau_sigma_mean)
    print()
    for j in range(len(T)):
        print('(All bands) mean[tau({}K)]= {:10.4f} (ps)'.format(T[j], np.mean(tau_sigma_mean[:, j])))


    # Plot
    from matplotlib import pyplot as plt
    for j in range(len(T)):
        if args.color_list is not None:
            plt.plot(bin_centers, tau_avg[j], label='{}K'.format(T[j]), c=args.color_list[j])
        else:
            plt.plot(bin_centers, tau_avg[j], label='{}K'.format(T[j]))
    if args.xlim is not None:
        plt.xlim(args.xlim)
    if args.ylim is not None:
        plt.ylim(args.ylim)
    plt.legend(fontsize='large').set_draggable(True)
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xlabel('Frequency (THz)', fontsize='x-large')
    plt.ylabel('Lifetime (ps)', fontsize='x-large')
    plt.title('{}, n(bins)={}'.format(args.phono3py_pckl, nbins), fontsize='small')
    #
    tau_flat = np.reshape(tau_avg, [-1])
    min_tau = np.min(tau_flat[tau_flat > 0])
    max_tau = np.max(tau_flat[np.logical_not(np.isnan(tau_flat))])
    plt.yscale('log')
    plt.xticks(range(0,25,5))
    plt.yticks([1, 10, 100, 1000, 10000, 100000])
    plt.subplots_adjust(
        0.30, # left
        0.25, # bottom
        0.70, # right
        0.75, # top
        0.10, # wspace
        0.20, # hspace
        )
    plt.grid(alpha=0.5)
    plt.show()

