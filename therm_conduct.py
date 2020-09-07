#!/usr/bin/env python
import numpy as np

def thermal_current(
    positions,
    atomic_energies,
    dt,
    ):
    """

    Get thermal current J.

    INPUT
    - positions (array): Atomic positions of cartesian coordinate in unit of "Angstrom". shape=(# of structures, len(atoms), space dimension(=3))
    - atomic_energies (array): Atomic energies(=kinetic+potential) in unit of "eV". shape=(# of structure, len(atoms))
    - dt (float): Time interval between two successive structures in unit of "ps".

    OUTPUT
    - A thermal current vector of shape=(# of structures-1, space dimension(=3))
    Note) The output function (thermal current vector, J) has time domain right-shifted by dt/2 from the original domain.

    """

    # @ Get thermal charge vector
    # -> shape = (len(alist), len(atoms), 1)
    atomic_energies = np.expand_dims(atomic_energies, 2)
    # -> shape = (len(alist), 3)
    thm_chg_vec = np.sum(positions * atomic_energies, 1)
    # -> shape = (len(alist)-1, 3)
    J = (thm_chg_vec[1:] - thm_chg_vec[:-1]) /dt
    return J

def heat_flux_ACF(
    thermal_current,
    num_avg_steps,
    ):
    """

    INPUT
    - thermal_current (array): 
    - num_avg_steps (int): Must be a positive integer.

    OUPUT
    - heat flux autocorrelation fuction (HFACF) matrix (array)
    Add formula

    """

    # J_tau
    # -> shape = (num_avg_steps, 3, 1)
    j_tau = np.expand_dims(thermal_current[:num_avg_steps], 2)

    # J_t+tau
    j_t_tau = []
    for i in range(num_avg_steps):
        l = len(thermal_current)-num_avg_steps+i
        j_t_tau.append(thermal_current[i:l])
    # -> shape = (num_avg_steps, len(thermal_current)-num_avg_steps, 1, 3)
    j_t_tau = np.expand_dims(j_t_tau, 2)

    # @ HFACF
    hfacf = []
    for i in range(len(j_t_tau)):
        hfacf.append(np.matmul(j_tau[i], j_t_tau[i]))
    # -> shape = (len(thermal_current)-num_avg_steps, 3, 3)
    hfacf = np.mean(hfacf, 0)

    return hfacf

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    This code will calculate thermal conductivity and heat-flux autocorrelation function (HFACF).
    Simulation must be about isothermal-isobaric MD.
    """)
    # Positional arguments
    parser.add_argument('traj_file', type=str, help='ASE readable atoms list file name.')
    parser.add_argument('dt', type=float, help='Time interval between steps selected in unit of picosec.')
    # Optional arguments
    parser.add_argument('-t', '--num_avg_steps', type=int, default=10, help='Number of steps for the time average in HFACF function. [Default: 10]')
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='ASE readable slice. default=":" (e.g.) -n :1000:10')
    parser.add_argument('-g', '--gauss_sigma', type=int, default=100, help='Number of steps for sigma of the Gaussian-smearing plot. [Default: 100]')
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
    print('This code will calculate thermal conductivity and heat-flux autocorrelation function (HFACF).'.center(120))
    print('Simulation must be about isothermal-isobaric MD.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args = argparse()
    traj_file     = args.traj_file
    dt            = args.dt
    num_avg_steps = args.num_avg_steps
    img_slice     = args.image_slice
    gauss_sigma   = args.gauss_sigma
    from ss_util import parse_slice
    real_slice = parse_slice(img_slice)

    from ase.io import read
    alist = read(traj_file, img_slice)
    posi = []
    a_pe = []
    a_ke = []
    temp = []
    volu = []
    for atoms in alist:
        posi.append(atoms.get_positions())
        a_pe.append(atoms.get_atomic_energies())
        a_ke.append(np.linalg.norm(atoms.get_momenta(), axis=-1)**2 /2 /atoms.get_masses())
        temp.append(atoms.get_temperature())
        volu.append(atoms.get_volume())
    posi = np.array(posi)
    a_ke = np.array(a_ke)
    a_pe = np.array(a_pe)
    # -> shape = (len(alist), len(atoms))
    a_te = a_pe + a_ke

    # from matplotlib import pyplot as plt
    # # for i in range(0,len(a_te[0]),1000):
    # for i in range(1):
        # plt.plot(a_te[:,i], c='k')
        # plt.plot(np.array(a_pe)[:,i])
        # print(np.std(np.array(a_pe)[:,i]))
        # # plt.plot(np.array(a_ke)[:,i])
    # # fig = plt.figure()
    # # for i in range(0,len(a_te[0]),1000):
        # # plt.plot(np.array(a_pe)[:,i])

    # plt.figure()
    # plt.plot(np.sum(a_te, 1), c='k')
    # plt.plot(np.sum(a_pe, 1))
    # print(np.std(np.sum(a_pe,1)))
    # # plt.plot(a_te[0])
    # # fig = plt.figure()
    # # plt.plot(a_pe[0])
    # # fig = plt.figure()
    # # plt.plot(a_ke[0])
    # plt.show()
    # exit(0)

    #
    j = thermal_current(posi, a_te, dt)
    hfacf = heat_flux_ACF(j, num_avg_steps)
    norm_hfacf = hfacf / hfacf[0]

    # from matplotlib import pyplot as plt
    # plt.plot(j[:,0])
    # plt.figure()
    # plt.plot(np.sum(a_te,1))
    # plt.figure()
    # plt.plot(np.sum(a_pe,1))
    # plt.figure()
    # plt.plot(np.sum(a_ke,1))
    # plt.figure()
    # plt.plot(np.sum(a_pe+a_ke,1))
    # plt.figure()
    # plt.plot(a_te[:,0])
    # plt.figure()
    # plt.plot(test[:,0])
    # plt.figure()
    # plt.plot(hfacf[:,0,0])
    # plt.show()
    # exit(0)

    # @ Thermal conductivity
    from ase import units
    # Unit for kappa: Joul *sec /(meter *Kelvin)
    kappa = np.add.accumulate(hfacf) /units.kB /np.mean(temp)**2 /np.mean(volu)**2 *units._e /units.second *1e10

    # det_hfacf = []
    # for i in range(len(hfacf)):
        # det_hfacf.append(np.linalg.det(hfacf[i]))

    # det_kappa = []
    # for i in range(len(kappa)):
        # det_kappa.append(np.linalg.det(kappa[i]))

    from scipy.ndimage import gaussian_filter1d as gf
    gf_norm_hfacf = gf(norm_hfacf, gauss_sigma, 0)
    gf_kappa = gf(kappa, gauss_sigma, 0)
    # gf_det_hfacf = gf(det_hfacf, gauss_sigma, 0)
    # gf_det_kappa = gf(det_kappa, gauss_sigma, 0)

    from matplotlib import pyplot as plt
    if real_slice.start:
        start = real_slice.start *dt + dt/2.
    else:
        start = dt/2.
    t = np.arange(len(kappa), dtype=float) *dt +start
    fig, ax1 = plt.subplots(3,3)
    for i in range(3):
        for j in range(3):
            ax2 = ax1[i,j].twinx()
            ax1[i,j].plot(t, norm_hfacf[:,i,j], alpha=0.5, c='b')
            ax2.plot(t, kappa[:,i,j], alpha=0.5, c='r')
            #
            ax1[i,j].plot(t, gf_norm_hfacf[:,i,j], c='b')
            ymin = np.amin(gf_norm_hfacf[:,i,j])
            ymax = np.amax(gf_norm_hfacf[:,i,j])
            ax1[i,j].set_ylim(1.1*ymin-0.1*ymax, 1.1*ymax-0.1*ymin)
            # 
            ax2.plot(t, gf_kappa[:,i,j], c='r')
            ymin = np.amin(gf_kappa[:,i,j])
            ymax = np.amax(gf_kappa[:,i,j])
            ax2.set_ylim(1.1*ymin-0.1*ymax, 1.1*ymax-0.1*ymin)
            #
            ax1[i,j].set_xlabel('Time (ps)', fontsize='x-large')
            ax1[i,j].set_ylabel('Scaled HFACF (Arb. Unit)', fontsize='x-large', color='b')
            ax2.set_ylabel('$\kappa_{}$$_{}$ (W/mK)'.format(i+1, j+1), fontsize='x-large', color='r')
            ax1[i,j].tick_params(axis="x",direction="in", labelsize='x-large')
            ax1[i,j].tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='b')
            ax2.tick_params(axis="y",direction="in", labelsize='x-large',colors='r',  labelcolor='r')
            plt.title('NAS={}, IS={}, $\sigma$={}'.format(num_avg_steps, img_slice, gauss_sigma), fontsize='x-large')
    # # Determinants
    # fig, ax1 = plt.subplots()
    # ax2 = ax1.twinx()
    # ax1.plot(t, det_hfacf[:], alpha=0.5, c='b')
    # ax2.plot(t, det_kappa[:], alpha=0.5, c='r')
    # #
    # ax1.plot(t, gf_det_hfacf[:], c='b')
    # ymin = np.amin(gf_det_hfacf[:])
    # ymax = np.amax(gf_det_hfacf[:])
    # ax1.set_ylim(1.1*ymin-0.1*ymax, 1.1*ymax-0.1*ymin)
    # #
    # ax2.plot(t, gf_det_kappa[:], c='r')
    # ymin = np.amin(gf_det_kappa[:])
    # ymax = np.amax(gf_det_kappa[:])
    # ax2.set_ylim(1.1*ymin-0.1*ymax, 1.1*ymax-0.1*ymin)
    # #
    # ax1.set_xlabel('Time (ps)', fontsize='x-large')
    # ax1.set_ylabel('Number density (atoms/$\AA^3$)', fontsize='x-large', color='b')
    # ax2.set_ylabel('Energy per atom (eV)', fontsize='x-large', color='r')
    # ax1.tick_params(axis="x",direction="in", labelsize='x-large')
    # ax1.tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='b')
    # ax2.tick_params(axis="y",direction="in", labelsize='x-large',colors='r',  labelcolor='r')
    # plt.title('$\kappa$, NAS={}, IS={}, $\sigma$={}'.format(num_avg_steps, img_slice, gauss_sigma), fontsize='x-large')
    plt.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.95, wspace=0.80, hspace=0.40)
    plt.show()

