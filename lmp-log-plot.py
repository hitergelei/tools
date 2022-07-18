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
            plt.figure()
            plt.plot(t, info[j][avail_info[i]], c='k')
            # plt.xlim((t[0], t[-1]))
            plt.title(avail_info[i], fontsize='x-large')
            plt.ylabel(avail_info[i], fontsize='x-large')
            if args.ns_unit:
                plt.xlabel('Time (ns)', fontsize='x-large')
            else:
                plt.xlabel('Time (ps)', fontsize='x-large')
            plt.tick_params(axis="both",direction="in", labelsize='x-large')
            plt.subplots_adjust(left=0.30, bottom=0.25, right=0.70, top=0.75, wspace=0.2, hspace=0.2)
            plt.grid(alpha=0.5)
    plt.show()



