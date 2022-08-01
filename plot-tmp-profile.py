#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    Plot tmp.profile file.
    """)
    # Positional arguments
    parser.add_argument('thermo_file', type=str, help='LAMMPS thermo file including heat flux data labeled as f_hf')
    parser.add_argument('tmp_file', type=str, help='LAMMPS format tmp_file.')
    parser.add_argument('dt', type=float, help='Variable dt for LAMMPS simulation (ps unit)')
    parser.add_argument('a', type=float, help='Lattice parameter along the direction of interest in Angst. unit')
    # Optional arguments
    parser.add_argument('-t', '--dump_frame', type=int, help='Dump the specified number of frames from head of the data.')
    return parser.parse_args()

def first_line(line):
    t, nbins, natoms = line[:-1].split()
    return float(t), int(nbins), int(natoms)

def data_line(line):
    chunk_i, coord, ncount, temp = line[:-1].split()
    return int(chunk_i), float(coord), float(ncount), float(temp)

def convergence_test(plt, t, T):
    """
    t (array of shape=(ndat))
        - Time array in ps unit
    T (array of shape=(ndat, nbins))
        - Temperature data
    """
    plt.figure()
    plt.plot(t, T[:, 0], c='b', label='0-th')
    halfway = len(T[0])//2
    plt.plot(t, T[:, halfway], c='r', label='{}-th'.format(halfway))
    plt.xlabel('Time (ps)', fontsize='x-large')
    plt.ylabel('Temperature (K)', fontsize='x-large')
    plt.title('Convergence test', fontsize='x-large', pad=10.)
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.subplots_adjust(left=0.12, bottom=0.25, right=0.99, top=0.75, wspace=0.2, hspace=0.2)
    plt.legend(fontsize='large').set_draggable(True)
    plt.grid(alpha=0.5)

def T_plot(plt, x, T):
    plt.figure()
    plt.plot(x, np.mean(T, axis=0), c='k')
    plt.xlabel(r'Coordinate ($\AA$)', fontsize='x-large')
    plt.ylabel('Temperature (K)', fontsize='x-large')
    plt.title('Temperature plot', fontsize='x-large', pad=10.)
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.subplots_adjust(left=0.12, bottom=0.25, right=0.99, top=0.75, wspace=0.2, hspace=0.2)
    plt.grid(alpha=0.5)
    
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
    print('Plot tmp.profile file.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    from lmp2traj import read_lmp_log
    thermo_info = read_lmp_log(args.thermo_file, )

    with open(args.tmp_file, 'r') as f:
        lines = f.readlines()

    _, nbins, natoms = first_line(lines[3])

    ndat = (len(lines)-3) // (nbins+1)
    t = np.zeros(ndat, dtype=float)
    T = np.zeros((ndat, nbins), dtype=float)
    coord = np.zeros(nbins, dtype=float)
    for i in range(ndat):
        t[i], _1, _2 = first_line(lines[3 + i*(nbins+1)])
        for j in range(nbins):
            _1, coord[j], _3, T[i][j] = data_line(lines[3 + i*(nbins+1) + 1 + j])
    # convert to ps unit
    t *= args.dt 
    # convert to Angst.
    x = coord *args.a

    from matplotlib import pyplot as plt
    # convergence_test(plt, t, T)
    if args.init_frame is not None:
        T = T[args.dump_frame:]
    T_plot(plt, x, T)
    plt.show()
