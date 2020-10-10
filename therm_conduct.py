#!/usr/bin/env python
import numpy as np

def classify_atoms(
    posi,
    lattice_params,
    ):
    """
    Devide box into 8 pieces to avoid problems from periodic boundary condition.
    Each atom will be classfied to the set for the nearest corner.
    The classification will be done independently for each ensemble case.

    INPUT
    - posi (array): Cartesian coordinate of the initial structure for one ensemble of len_t long (in unit of "Angstrom").
                         Shape = (# of ensemble average, len(atoms), space dimension(=3))
    - lattice_params (array): 
                       Shape = (# of ensemble average, 3, 3)

    OUTPUT
    - Class of each atom. Kinds of class = range(8)
      Shape = (# of ensemble average, len(atoms))
    """

    # @ Get 8 corners
    corners = []
    for n in range(len(lattice_params)):
        tmp = []
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    tmp.append(np.sum(lattice_params[n]* np.array([[i],[j],[k]]), axis=0))
        corners.append(tmp)
    # Shape = (len(posi), 1, 8, 3)
    corners = np.expand_dims(corners, 1)
    # Shape = (len(posi), len(atoms), 1, 3)
    posi = np.expand_dims(posi, 2)
    
    # @ Classify
    # Shape = (len(posi), len(atoms), 8)
    dist = np.linalg.norm(corners - posi, axis=3)

    return np.argmin(dist, axis=2)

def energy_moment(
    positions,
    atomic_energies,
    ):
    """

    Get thermal flux J.

    INPUT
    - positions (array): Atomic positions of cartesian coordinate in ASE unit (Angstrom). shape=(# of structures, len(atoms), space dimension(=3))
    - atomic_energies (array): Atomic energies(=kinetic+potential) in unit of "eV". shape=(# of structure, len(atoms))

    OUTPUT
    - A thermal charge vector of shape=(# of structures, space dimension(=3))

    """

    # @ Get thermal charge vector
    # -> shape = (len_t, len(atoms), 1)
    atomic_energies = np.expand_dims(atomic_energies, 2)
    # -> shape = (len_t, 3)
    energy_moment = np.sum(positions * atomic_energies, 1)
    return energy_moment

def einstein_eq(
    energy_moment,
    dt,
    ):
    """

    INPUT
    - energy_moment (array): A thermal charge vector of shape=(# of structures, space dimension(=3))

    OUTPUT
    - 
    add equation 

    """

    # shape of (len_t, 3)
    e_mom = energy_moment - np.expand_dims(energy_moment[0], axis=0)
    # shape of (len_t, 3, 3)
    mom_mat = np.matmul(
        np.expand_dims(e_mom, axis=2),
        np.expand_dims(e_mom, axis=1),
        )
    t = np.expand_dims((np.arange(len(energy_moment), dtype=float)) *dt, [1,2])
    t[0,0,0] = 1e3

    kappa = mom_mat /2. /(dt*len(t)) /units.kB /np.mean(temp)**2 /np.mean(volu) *units._e *units.second *1e10
    return kappa

def thermal_flux(
    positions,
    velocities,
    atomic_energies,
    dt,
    ):
    """

    Get thermal flux J.

    INPUT
    - positions (array): Atomic positions of cartesian coordinate in ASE unit (Angstrom). shape=(# of structures, len(atoms), space dimension(=3))
    - velocities (array): Atomic velocities in ASE unit. shape=(# of structures, len(atoms), space dimension(=3))
    - atomic_energies (array): Atomic energies(=kinetic+potential) in unit of "eV". shape=(# of structure, len(atoms))
    - dt (float): Time interval between two successive structures in unit of ASE.

    OUTPUT
    - A thermal flux vector of shape=(# of structures-1, space dimension(=3))
    Note) The output function (thermal flux vector, J) has time domain right-shifted by dt/2 from the original domain.
    >>add equation

    """

    # @ first term
    first = np.sum(velocities * np.expand_dims(atomic_energies, axis=2), axis=1)
    # -> shape = (# of structures-1, 3)
    first = (first[1:] + first[:-1]) /2.

    # @ second term
    # -> shape = (# of structures-1, len(atoms))
    dEdt = (atomic_energies[1:] - atomic_energies[:-1]) /dt
    # -> shape = (# of structures-1, 3)
    second = np.sum(positions[1:] * np.expand_dims(dEdt, axis=2), axis=1)
    
    # -> shape = (# of structures-1, 3)
    J = first + second
    return J

def heat_flux_ACF(
    thermal_flux,
    # num_avg_steps,
    # avg_intvl,
    ):
    """

    INPUT
    - thermal_flux (array): 
    # - num_avg_steps (int): Must be a positive integer.
    # - avg_intvl (int): Interval between HFACF average in unit of steps

    OUPUT
    - heat flux autocorrelation fuction (HFACF) matrix (array)
    Add formula

    """

    # J_tau
    # -> shape = (3, 1)
    j_tau = np.expand_dims(thermal_flux[0], 1)

    # J_t+tau
    # -> shape = (len_t-1, 1, 3)
    j_t_tau = np.expand_dims(thermal_flux, 1)

    # @ HFACF
    # -> shape = (len_t-1, 3, 3)
    hfacf = np.matmul(j_tau, j_t_tau)

    # @ Take symmetric part only.
    hfacf = (hfacf + np.transpose(hfacf, [0,2,1]))/2.

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
    parser.add_argument('-a', '--avg_intvl', type=int, default=1, help='Interval between HFACF average in unit of steps')
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='ASE readable slice. default=":" (e.g.) -n :1000:10')
    parser.add_argument('-g', '--gauss_sigma', type=int, default=100, help='Number of steps for sigma of the Gaussian-smearing plot. [Default: 100]')
    # parser.add_argument('-d', '--dont_devide', dest=devide_box, action='store_false', help='[Default: devide box into 8 pieces.]')
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
    avg_intvl     = args.avg_intvl
    # devide_box    = args.devide_box

    from ss_util import parse_slice
    real_slice = parse_slice(img_slice)

    from ase import units
    dt *= 1e3 *units.fs
    from ase.io import read
    alist = read(traj_file, img_slice)
    posi = []
    velo = []
    a_pe = []
    a_ke = []
    temp = []
    volu = []
    latt = []
    for atoms in alist:
        posi.append(atoms.get_positions())
        velo.append(atoms.get_velocities())
        a_pe.append(atoms.get_potential_energies())
        a_ke.append(np.linalg.norm(atoms.get_momenta(), axis=-1)**2 /2 /atoms.get_masses())
        temp.append(atoms.get_temperature())
        volu.append(atoms.get_volume())
        latt.append(atoms.get_cell())
    # -> shape = (len(alist), len(atoms), 3)
    posi = np.array(posi)
    velo = np.array(velo)
    a_ke = np.array(a_ke)
    a_pe = np.array(a_pe)
    # -> shape = (len(alist), 3, 3)
    latt = np.array(latt)
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


    # # @ Devide box
    # # Shape = (len(alist), len(atoms))
    # masker = classify_atoms(posi, latt)
        
    # # @ Cut ensembles
    # for i in range(num_avg_steps):
        # l = len(thermal_flux)-(num_avg_steps-i-1)*avg_intvl
        # j_t_tau.append(thermal_flux[i*avg_intvl:l])

    # @ Get HFACF
    len_t = len(posi)-(num_avg_steps-1)*avg_intvl
    hfacf = []
    for i in range(num_avg_steps):
        J = thermal_flux(
            posi[i*avg_intvl:i*avg_intvl+len_t],
            velo[i*avg_intvl:i*avg_intvl+len_t],
            a_te[i*avg_intvl:i*avg_intvl+len_t],
            dt,
            )
        hfacf.append(heat_flux_ACF(J))
    hfacf = np.mean(hfacf, axis=0)
    norm_hfacf = hfacf / hfacf[0]
    kappa = np.add.accumulate(hfacf) *dt /units.kB /np.mean(temp)**2 /np.mean(volu) *units._e *units.second *1e10

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

    # # @ Thermal conductivity
    # # Unit for kappa: Joul /(sec *meter *Kelvin)
    # len_t = len(posi)-(num_avg_steps-1)*avg_intvl
    # kappa = []
    # for i in range(num_avg_steps):
        # e_mom = energy_moment(posi[i*avg_intvl:i*avg_intvl+len_t], a_te[i*avg_intvl:i*avg_intvl+len_t])
        # kappa.append(einstein_eq(e_mom, dt))
    # kappa = np.mean(kappa, axis=0)
    # hfacf = kappa[1:] - kappa[:-1]
    # norm_hfacf = hfacf / hfacf[0]

    # det_hfacf = []
    # for i in range(len(hfacf)):
        # det_hfacf.append(np.linalg.det(hfacf[i]))

    # det_kappa = []
    # for i in range(len(kappa)):
        # det_kappa.append(np.linalg.det(kappa[i]))

    avg_norm_hfacf = []
    for i in range(len(hfacf)):
        tmp = []
        for j in range(3):
            tmp.append(norm_hfacf[i,j,j])
        avg_norm_hfacf.append(np.mean(tmp))

    avg_kappa = []
    for i in range(len(kappa)):
        tmp = []
        for j in range(3):
            tmp.append(kappa[i,j,j])
        avg_kappa.append(np.mean(tmp))

    from scipy.ndimage import gaussian_filter1d as gf
    gf_norm_hfacf = gf(norm_hfacf, gauss_sigma, 0)
    gf_kappa = gf(kappa, gauss_sigma, 0)
    # gf_det_hfacf = gf(det_hfacf, gauss_sigma, 0)
    # gf_det_kappa = gf(det_kappa, gauss_sigma, 0)
    gf_avg_norm_hfacf = gf(avg_norm_hfacf, gauss_sigma, 0)
    gf_avg_kappa = gf(avg_kappa, gauss_sigma, 0)

    from matplotlib import pyplot as plt
    if real_slice.start:
        start = real_slice.start * args.dt + args.dt/2.
    else:
        start = args.dt/2.
    t = np.arange(len(kappa), dtype=float) *args.dt +start
    fig, ax1 = plt.subplots(3,3)
    for i in range(3):
        for j in range(3):
            ax2 = ax1[i,j].twinx()
            ax1[i,j].plot(t, norm_hfacf[:,i,j], alpha=0.5, c='b')
            ax2.plot(t, kappa[:,i,j], alpha=0.5, c='r')
            #
            ax1[i,j].plot(t, gf_norm_hfacf[:,i,j], c='b')
            # ymin = np.amin(gf_norm_hfacf[:,i,j])
            # ymax = np.amax(gf_norm_hfacf[:,i,j])
            ymin = np.amin(norm_hfacf[:,i,j])
            ymax = np.amax(norm_hfacf[:,i,j])
            ax1[i,j].set_ylim(1.1*ymin-0.1*ymax, 1.1*ymax-0.1*ymin)
            # 
            ax2.plot(t, gf_kappa[:,i,j], c='r')
            # ymin = np.amin(gf_kappa[:,i,j])
            # ymax = np.amax(gf_kappa[:,i,j])
            ymin = np.amin(kappa[:,i,j])
            ymax = np.amax(kappa[:,i,j])
            ax2.set_ylim(1.1*ymin-0.1*ymax, 1.1*ymax-0.1*ymin)
            #
            ax1[i,j].set_xlabel('Time (ps)', fontsize='x-large')
            ax1[i,j].set_ylabel('Scaled HFACF (Arb. Unit)', fontsize='x-large', color='b')
            ax2.set_ylabel('$\kappa_{}$$_{}$ (W/mK)'.format(i+1, j+1), fontsize='x-large', color='r')
            ax1[i,j].tick_params(axis="x",direction="in", labelsize='x-large')
            ax1[i,j].tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='b')
            ax2.tick_params(axis="y",direction="in", labelsize='x-large',colors='r',  labelcolor='r')
            plt.title('NAS={}, IS={}, $\sigma$={}'.format(num_avg_steps, img_slice, gauss_sigma), fontsize='x-large')
    plt.subplots_adjust(left=0.10, bottom=0.05, right=0.90, top=0.95, wspace=0.80, hspace=0.40)

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

    # Average
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(t, avg_norm_hfacf[:], alpha=0.5, c='b')
    ax2.plot(t, avg_kappa[:], alpha=0.5, c='r')
    #
    ax1.plot(t, gf_avg_norm_hfacf[:], c='b')
    ymin = np.amin(avg_norm_hfacf[:])
    ymax = np.amax(avg_norm_hfacf[:])
    ax1.set_ylim(1.1*ymin-0.1*ymax, 1.1*ymax-0.1*ymin)
    #
    ax2.plot(t, gf_avg_kappa[:], c='r')
    ymin = np.amin(avg_kappa[:])
    ymax = np.amax(avg_kappa[:])
    ax2.set_ylim(1.1*ymin-0.1*ymax, 1.1*ymax-0.1*ymin)
    #
    ax1.set_xlabel('Time (ps)', fontsize='x-large')
    ax1.set_ylabel('Scaled HFACF (Arb. Unit)', fontsize='x-large', color='b')
    ax2.set_ylabel('$\kappa$ (W/mK)', fontsize='x-large', color='r')
    ax1.tick_params(axis="x",direction="in", labelsize='x-large')
    ax1.tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='b')
    ax2.tick_params(axis="y",direction="in", labelsize='x-large',colors='r',  labelcolor='r')
    plt.title('NM) IS={}, NAS={}, AI={}, $\sigma$={}, dt={}, T={:.2f}'.format(img_slice, num_avg_steps, avg_intvl, gauss_sigma, args.dt, np.mean(temp)))
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.85, top=0.90)

    plt.show()

