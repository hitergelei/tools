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
        t = np.arange(len(info[j][list(info[j].keys())[0]]), dtype=float) * dt

        # Remove dummy info.
        avail_info = []
        for (key, value) in info[j].items():
            if np.std(np.array(value, dtype=float)) > 1e-10:
                avail_info.append(key)

        # > Plot
        from matplotlib import pyplot as plt
        for i in range(len(avail_info)):
            plt.figure()
            plt.plot(t, info[j][avail_info[i]], c='k')
            plt.title(avail_info[i], fontsize='x-large')
            plt.tick_params(axis="both",direction="in", labelsize='x-large')
            plt.grid(alpha=0.2)
    plt.show()



