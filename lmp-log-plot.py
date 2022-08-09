#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Lammps logfile plot.
    """)
    # Optional arguments
    parser.add_argument('-f', '--log_file', type=str, default='log.lammps', help='Specify the log file to plot. Default: log.lammps')
    parser.add_argument('-a', '--init_pos', type=str, default=None, help='Provide initial Atoms structure file to get number of atoms info. [default: None]')
    parser.add_argument('-t', '--ns_unit', action='store_true', help='Plot with time unit of ns. [default: ps]')
    parser.add_argument('-x', '--xticks', type=float, nargs='+', default=None, help='Set x-ticks manually.')
    parser.add_argument('-y', '--yticks', type=float, nargs='+', default=None, help='Set y-ticks manually.')
    parser.add_argument('-l', '--xlog', action='store_true', help='Set x-axis as log-scale')
    parser.add_argument('-m', '--ylog', action='store_true', help='Set y-axis as log-scale')
    parser.add_argument('-s', '--sigma', type=int, help='The half number of steps for standard deviation plot')
    return parser.parse_args()

if __name__ == '__main__':
    # > Intro
    import datetime
    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ POSTECH <<<<<'.center(120))
    print(('Code runtime : '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('This code will plot the lammps log file for you.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    # > Read input params
    # params
    log_file = args.log_file

    # > Main
    if args.init_pos:
        from ase.io import read
        anum = len(read(args.init_pos, 0))
    else:
        anum = None

    #
    from lmp2traj import read_lmp_log
    info = read_lmp_log(log_file)
    with open(log_file) as f:

        # Read time interval
        while 1:
            words = f.readline().split()
            if words == []:
                pass
            elif words[0] == 'timestep':
                if words[1][0] == '$':
                    words = f.readline().split()
                dt = float(words[1])
                break
            elif words[0] == 'Time':
                dt = float(words[3])
                break

    for j in range(len(info)):
        # > Post process
        t = info[j]['Step'] *dt
        if args.ns_unit:
            t /= 1000.
        if args.xlog:
            t += t[1] - t[0]
        # t = np.arange(len(info[j][list(info[j].keys())[0]]), dtype=float) * dt

        # Remove dummy info.
        avail_info = []
        for (key, value) in info[j].items():
            if np.std(np.array(value, dtype=float)) > 0.:
                avail_info.append(key)

        #
        extensive_keys = ['TotEng', 'PotEng', 'KinEng', 'Volume']
        if anum is not None:
            for key in extensive_keys:
                if key in info[j].keys():
                    avail_info.append(key+' per atom')
                    info[j][key+' per atom'] = info[j][key] /anum
            if 'Volume' in info[j].keys():
                avail_info.append('Density')
                info[j]['Density'] = anum /info[j]['Volume']



        # > Plot
        from matplotlib import pyplot as plt
        for i in range(len(avail_info)):
            fig, ax1 = plt.subplots()
            plt.plot(t, info[j][avail_info[i]], c='k')
            # plt.xlim((t[0], t[-1]))
            if args.xlog:
                plt.xscale('log')
            if args.ylog:
                plt.yscale('log')
            plt.title(avail_info[i], fontsize='x-large')
            plt.ylabel(avail_info[i], fontsize='x-large')
            if args.ns_unit:
                plt.xlabel('Time (ns)', fontsize='x-large')
            else:
                plt.xlabel('Time (ps)', fontsize='x-large')
            if args.xticks is not None:
                plt.xticks(args.xticks)
            if args.yticks is not None:
                plt.yticks(args.yticks)
            plt.tick_params(axis="both",direction="in", labelsize='x-large')
            if args.sigma is not None:
                if len(t) > args.sigma*2 +1:
                    ax2 = ax1.twinx()
                    std = np.zeros(len(t), dtype=float)
                    for k in range(args.sigma, len(t)-args.sigma):
                        std[k] = np.std(info[j][avail_info[i]][k-args.sigma:k+args.sigma+1])
                    ax2.plot(t, std, c='r')
                    ax2.set_ylabel('Standard deviation', fontsize='x-large', c='r')
                    ax2.tick_params(axis="y",direction="in", labelsize='x-large', colors='r', labelcolor='r')
            plt.subplots_adjust(left=0.30, bottom=0.25, right=0.70, top=0.75, wspace=0.2, hspace=0.2)
            plt.grid(alpha=0.5)
    plt.show()



