#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    Plot tmp.profile file.
    """)
    # Positional arguments
    parser.add_argument('thermo_file', type=str, help='LAMMPS thermo file including heat transfer data labeled as f_hf')
    parser.add_argument('tmp_file', type=str, help='LAMMPS format tmp_file.')
    parser.add_argument('dt', type=float, help='Variable dt for LAMMPS simulation (ps unit)')
    parser.add_argument('a', type=float, help='Lattice parameter along the direction of interest in Angst. unit')
    parser.add_argument('area', type=float, help='Area of intersection perpendicular to direction of the interest in Angst^2 unit.')
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
    from scipy.stats import linregress as lr
    halfway = len(x) //2
    fit_l = lr(x[3:halfway-2], np.mean(T[:, 3:halfway-2], axis=0))
    fit_r = lr(x[halfway+3:-2], np.mean(T[:, halfway+3:-2], axis=0))
    # gradT in unit of ( K / Angst )
    gradT = (fit_l.slope - fit_r.slope) /2.

    plt.figure()
    plt.errorbar(x, np.mean(T, axis=0), yerr=np.std(T, axis=0), fmt='s', c='k', mfc='none', ecolor='k', capsize=3)
    plt.plot(x[3:halfway-2], x[3:halfway-2] *fit_l.slope + fit_l.intercept, c='r', label=r'$\nabla T$={:.2f}'.format(fit_l.slope))
    plt.plot(x[halfway+3:-2], x[halfway+3:-2] *fit_r.slope + fit_r.intercept, c='b', label=r'$\nabla T$={:.2f}'.format(fit_r.slope))
    plt.xlabel(r'Coordinate ($\rm \AA$)', fontsize='x-large')
    plt.ylabel('Temperature (K)', fontsize='x-large')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.subplots_adjust(left=0.25, bottom=0.25, right=0.75, top=0.75, wspace=0.2, hspace=0.2)
    plt.legend(fontsize='large').set_draggable(True)
    plt.grid(alpha=0.5)
    return gradT
    
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
    heat_trans = thermo_info[0]['f_hf'][1:]

    with open(args.tmp_file, 'r') as f:
        lines = f.readlines()

    _, nbins, natoms = first_line(lines[3])

    ndat = (len(lines)-3) // (nbins+1)
    if ndat != len(heat_trans):
        raise RuntimeError('{} file and {} file, the number of time point mismatch'.format(args.thermo_file, args.tmp_file))
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

    #
    if args.dump_frame is not None:
        t          = t[args.dump_frame:]
        T          = T[args.dump_frame:]
        heat_trans = heat_trans[args.dump_frame:]

    # heat_flux in unit of ( eV / ps Angst^2 )
    J = (heat_trans[-1] - heat_trans[0]) /(t[-1]-t[0]) /args.area /2.

    from matplotlib import pyplot as plt
    gradT = T_plot(plt, x, T)

    # kappa in unit of ( eV / ps Angst K )
    kappa = J / gradT
    from ase import units
    # kappa_si in unit of ( W / m K )
    kappa_si = kappa *units._e *1e12 *1e10

    plt.title(r'$\kappa$={:.4f} (W/mK)'.format(kappa_si), fontsize='x-large', pad=10.)
    plt.show()
