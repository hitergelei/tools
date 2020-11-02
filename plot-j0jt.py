#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    This code will plot HFACF and lattice thermal conductivity from lammps J0Jt.dat file.
    Run this code in the folder containing 'J0Jt.dat' and 'out.dump' file.
    """)
    # Positional arguments
    parser.add_argument('temp', type=float, help='Temperature of equilibrium MD simulation. Unit in Kelvin.')
    parser.add_argument('dt', type=float, help='Time step of J0Jt.dat file. Becareful about unit. (Real: fs, metal: ps)')
    # Optional arguments
    parser.add_argument('-u', '--lammps_unit', type=str, default='metal', help='Set unit of J0Jt.dat file between metal and real. [default: metal]')
    parser.add_argument('-a', '--avg_intvl', type=int, default=1, help='"tau" sample interval. [default: 1]')
    parser.add_argument('-b', '--bin_intvl', type=int, default=1, help='"t" sample interval. [default: 1]')
    parser.add_argument('-c', '--corr_len', type=float, default=None, help='Set max correlation time length in ps unit. i.e. x-axis length. [default: maximal]')
    parser.add_argument('-l', '--dont_load', dest='load_bool', action='store_false', help="If provided, don't load the data. [default: load]")
    parser.add_argument('-s', '--dont_save', dest='save_bool', action='store_false', help="If provided, don't save the data. [default: save]")
    parser.add_argument('-p', '--dont_plot', dest='plot_bool', action='store_false', help="If provided, don't plot the data. [default: plot]")

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
    print('This code will plot HFACF and lattice thermal conductivity from lammps J0Jt.dat file.'.center(120))
    print('Run this code in the folder containing "J0Jt.dat" and "out.dump" file.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args = argparse()
    temp = args.temp
    dt = args.dt
    lammps_unit = args.lammps_unit
    avg_intvl = args.avg_intvl
    bin_intvl = args.bin_intvl
    corr_len = args.corr_len
    load_bool = args.load_bool
    save_bool = args.save_bool
    plot_bool = args.plot_bool

    # Load
    from os import stat
    mod_date = '{}-{}'.format(stat('out.dump')[8], stat('J0Jt.dat')[8])
    fname = 'j0jt-plot/T{}_dt{}_avg-intvl{}_bin-intvl{}_corr-len{}.npy'.format(temp, dt, avg_intvl, bin_intvl, corr_len)
    if load_bool:
        try:
            hfacf = np.load(fname)
            len_t = np.load('j0jt-plot/len_t.npy'.format(temp, dt))
            len_tau = np.load('j0jt-plot/len_tau.npy'.format(temp, dt))
        except:
            load_bool = False
            print('Failed to load "{}" file.'.format(fname))
        else:
            if corr_len is None:
                c_l = len_t
            else:
                c_l = int(corr_len/ dt)
            print('Successfully load "{}" file.'.format(fname))

    # Main
    if not load_bool:
        with open('J0Jt.dat') as f:
            lines = f.readlines()
        len_t = int(lines[3].split()[1])
        len_tau = (len(lines)-3)//(len_t+1)

        if corr_len is None:
            c_l = len_t
        else:
            c_l = int(corr_len/ dt)

        hfacf = np.zeros((c_l//bin_intvl,3))
        for tau in range(c_l,len_tau,avg_intvl):
            # print('Processing with tau={}'.format(tau))
            hfacf += np.loadtxt(lines[tau*(len_t+1)+4:tau*(len_t+1)+4+c_l:bin_intvl], usecols=(3,4,5))
        hfacf = np.mean(hfacf, axis=1) /((len_tau-c_l) //avg_intvl)

        # Save
        if save_bool:
            from subprocess import call
            call('mkdir j0jt-plot', shell=True)
            np.save(fname, hfacf)
            np.save('j0jt-plot/len_t.npy', len_t)
            np.save('j0jt-plot/len_tau.npy', len_tau)

    #
    scale_hfacf = hfacf / hfacf[0]
    # Scale
    from ase.io import read
    atoms = read('out.dump', 0)
    vol = np.linalg.det(atoms.get_cell())
    if lammps_unit == 'metal':
        scale =  1.60218e-19 *1e22 /8.61733034e-5 / temp**2 /vol *dt *bin_intvl
    elif lammps_unit == 'real':
        scale = (4186/6.02214e23)**2 *1e25 /1.3806504e-23/ temp**2 /vol *dt *bin_intvl
    #
    kappa = np.add.accumulate(hfacf) *scale


    # Plot
    if plot_bool:
        #
        t = np.arange(len(kappa), dtype=float) *dt *bin_intvl
        # set time axis unit to ps.
        if lammps_unit == 'real':
            t /= 1e3

        #
        from matplotlib import pyplot as plt
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.plot(t, scale_hfacf, c='b')
        ax2.plot(t, kappa, c='r')
        #
        ax1.set_xlabel('Time (ps)', fontsize='x-large')
        ax1.set_ylabel('Scaled HFACF (Arb. Unit)', fontsize='x-large', color='b')
        ax2.set_ylabel('$\kappa$ (W/mK)', fontsize='x-large', color='r')
        ax1.tick_params(axis="x",direction="in", labelsize='x-large')
        ax1.tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='b')
        ax2.tick_params(axis="y",direction="in", labelsize='x-large',colors='r',  labelcolor='r')
        plt.title('T={}K, c_l={}, l_a={}, dt={:.3f}, ai={}, bi={}'.format(temp, c_l, len_tau-c_l, len_tau, dt, avg_intvl, bin_intvl))
        plt.subplots_adjust(left=0.15, bottom=0.15, right=0.85, top=0.90)

        plt.show()
