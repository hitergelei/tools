#!/usr/bin/env python

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    dE/dt plot.
    """)
    # Optional arguments
    parser.add_argument('traj_file', type=str, help='Trajectory file to analyze.')
    parser.add_argument('dt', type=float, help='Time interval of traj_file in ps.')
    parser.add_argument('-s', '--gauss_sigma', type=int, default=0, help='Specify the sigma of smearing in unit of steps (int), if you wanna do Gaussian smearing to dE/dt plot. Default: No smearing')
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
    print('This code will plot dE/dt.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    # > Read input params
    dt          = args.dt
    traj_file   = args.traj_file
    gauss_sigma = args.gauss_sigma

    # > Main
    import numpy as np
    from ase.io import read
    alist = read(traj_file, ':')
    len_atoms = len(alist[0])
    epot = []
    for i in range(len(alist)):
        epot.append(alist[i].get_potential_energy()/ len_atoms)
    epot_0 = epot[:-2]
    epot_1 = epot[1:-1]
    epot_2 = epot[2:]
    epot_l = np.mean([epot_0, epot_1], axis=0)
    epot_r = np.mean([epot_1, epot_2], axis=0)
    dEdt = (epot_r - epot_l) /dt
    if gauss_sigma != 0:
        from scipy.ndimage import gaussian_filter
        dEdt = gaussian_filter(dEdt, sigma=gauss_sigma)
    t = np.arange(len(alist)) *dt

    # > Plot
    from matplotlib import pyplot as plt
    # E
    fig, ax1 = plt.subplots()
    ax1.plot(t, epot, c='k')
    ax1.set_xlabel('Time (ps)', fontsize='x-large')
    ax1.set_ylabel('Potential Energy (eV/atom)', fontsize='x-large')
    ax1.tick_params(axis="both",direction="in", labelsize='x-large')
    # dEdt
    ax2 = ax1.twinx()
    ax2.plot(t[1:-1], dEdt*1e3, c='r', label='$\Delta t$={}(ps)\nGaussian smearing $\sigma$={} (ps)'.format(dt, gauss_sigma *dt))
    ax2.set_ylabel('$\Delta E$/$\Delta t$ (meV/ps)', fontsize='x-large', c='r')
    ax2.tick_params(axis="y",direction="in", labelsize='x-large', colors='r', labelcolor='r')
    ax1.xaxis.grid(alpha=0.2)
    ax2.grid(alpha=0.2)
    plt.legend(loc=(0.00, 1.05), fontsize='x-large')
    plt.subplots_adjust(left=0.19, right=0.85, top=0.80)
    plt.show()
