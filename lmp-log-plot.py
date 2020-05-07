#!/usr/bin/env python

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
    import numpy as np
    with open(log_file) as f:

        # Read time interval
        while 1:
            words = f.readline().split()
            if words == []:
                pass
            elif words[0] == 'timestep':
                dt = float(words[1])
                break
            elif words[0] == 'Time':
                dt = float(words[3])
                break

        # Check the kinds of info we have.
        info_name = []
        while 1:
            line = f.readline()
            if not line:
                raise RuntimeError('No an information exists.')
            words = line.split()
            if words[1] == 'Step':
                while 1:
                    words = f.readline().split()
                    if words[1] == 'Step':
                        break
                    for i in words[::3]:
                        info_name.append(i)
                break

        # Start again
        f.seek(0)
        while 1:
            words = f.readline().split()
            if len(words) < 2:
                pass
            elif words[1] == 'Step':
                t = [int(words[2])]
                break

        # Gather info
        info = []
        new_info = []
        while 1:
            line = f.readline()
            if not line:
                if new_info:
                    info.append(new_info)
                break
            words = line.split()
            if len(words) < 2:
                break
            elif words[1] == 'Step':
                t.append(int(words[2]))
                info.append(new_info)
                new_info = []
            else:
                new_info.extend(words[2::3])

    # Check if the info is complete.
    if len(info[-1]) != len(info_name):
        del(info[-1])
        del(t[-1])
    if len(info) != len(t):
        del(t[-1])
    if len(info) != len(t):
        raise ValueError("What's wrong? T-T")

    # > Post process
    t = np.array(t) * dt
    info = np.array(info, dtype=np.float64).T

    # Remove dummy info.
    mask_bool = np.std(info, axis=1) != 0.
    info = info[mask_bool]
    info_name = np.array(info_name)[mask_bool]

    # > Plot
    from matplotlib import pyplot as plt
    for i in range(len(info_name)):
        plt.figure()
        plt.plot(t, info[i], c='k')
        plt.title(info_name[i], fontsize='x-large')
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.grid(alpha=0.2)
    plt.show()



