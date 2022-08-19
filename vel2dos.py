#!/usr/bin/env python

import numpy as np

def velocity2phononPDOS(atomic_mass_arr, average_temp, velocity_arr, dt):
    """
    atomic_mass_arr (array) : Array of atomic masses of all atoms in unit cell in amu (atomic mass unit).
    average_temp (float)    : Average temperature of provided ensemble trajectory in kelvin.
    velocity_arr (array)    : All of velocities in trajectory of all atoms. Array shape == ( len(traj), len(atoms), 3 )
    dt (float)             : Time interval of images provided in (ps) unit (ps^(-1) --> THz)

    ================================================
                        Theory                      
    ================================================
    This fomula gives phonon density of states from 
    velocities of atoms in trajectory.

                      1         N             _             
                     --- * (summation) m_n * |v_n(f)|^2     
      1               2        n=1                          
     ---- * g(f) = --------------------------------------   
     M*df                        3                          
                     M*df * M * --- * NkT                   
                                 2                          
    _
    v_n(f) : velocity vector in frequency space
    m_n    : Atomic mass of atom n
    N      : Number of atoms in unit cell
    df     : Magnitude of frequency bin
    M      : Number of time-step (= 2 times of number of positive frequency bins)
    n      : Atom index
    k      : Boltzmann constant
    T      : Temperature
    g(f)   : Phonon density of state

    """
    # Get params
    image_num = len(velocity_arr)
    natoms = len(atomic_mass_arr)
    t_init = 0.
    t_fin = (image_num-1) * dt
    # Get discrete time domain
    t = np.arange(t_init, t_fin+dt, dt)
    # Get frequency domain
    f = np.fft.rfftfreq(image_num) / dt
    d_f = f[1] - f[0]
    # Reshape
    velocity_arr = np.reshape(velocity_arr, (image_num, natoms*3)).T
    # Fourier transform for each DOF
    v_f_arr = []
    for DOF_i in range(len(velocity_arr)):
        # Fourier transform and normalize
        v_f = np.fft.rfft(velocity_arr[DOF_i], norm=None)
        v_f_arr.append(v_f)
    v_f_arr = np.array(v_f_arr)
    # calculate the formula (First factor 2 is from throwing away negative frequencies) (Or 0~inf integral of g(f) equals to be 1/2)
    from ase import units
    ADOS = 2. * (np.repeat(atomic_mass_arr, 3) / 3. / len(t)**2 / natoms / units.kB / average_temp / d_f * np.square(np.abs(v_f_arr)).T).T

    return f, np.reshape(ADOS, (natoms, 3, -1))
        

def plot_total_DOS(
    plt,
    f,
    ADOS,
    unit='THz',
    freqlim_low=None,
    freqlim_up=None,
    doslim_low=None,
    doslim_up=None,
    legend_bool=True,
    boson_peak=False,
    ):
    """
    f (arr)
    ADOS (arr)
    unit (str) : 'THz' and 'meV' are implimented.
    freqlim_low (float or None) : 
    freqlim_up  (float or None) : 
    """
    ## Preprocess
    DOS = np.sum(np.sum(ADOS, axis=0), axis=0)

    ## Plot
    fig, ax = plt.subplots()
    ax.plot(DOS, f, c='k')
    ax.set_ylim((freqlim_low, freqlim_up))
    if boson_peak:
        ax.set_xlabel('$g(\omega)/\omega^2$', fontsize='x-large')
    else:
        ax.set_xlabel('$g(\omega)$', fontsize='x-large')
    # ax.set_xticklabels([])
    ax.set_ylabel('Frequency ({})'.format(unit), fontsize='x-large')
    ax.set_xlim((doslim_low,doslim_up))
    plt.subplots_adjust(left=0.40, bottom=0.25, right=0.60, top=0.752, wspace=0.2, hspace=0.2)
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    ax.xaxis.set_major_locator(plt.MaxNLocator(1))
    plt.grid(alpha=0.4)

def plot_direc_DOS(
    plt,
    f,
    ADOS,
    unit='THz',
    freqlim_low=None,
    freqlim_up=None,
    doslim_low=None,
    doslim_up=None,
    legend_bool=True,
    boson_peak=False,
    ):
    """
    f (arr)
    ADOS (arr)
    unit (str) : 'THz' and 'meV' are implimented.
    freqlim_low (float or None) : 
    freqlim_up  (float or None) : 
    """
    ## Preprocess
    DOS = np.sum(ADOS, axis=0)

    ## Plot
    fig, ax = plt.subplots()
    ax.plot(DOS[0], f, c='g', label='x')
    ax.plot(DOS[1], f, c='r', label='y')
    ax.plot(DOS[2], f, c='b', label='z')
    ax.plot(np.sum(DOS, axis=0), f, c='k')
    ax.set_ylim((freqlim_low, freqlim_up))
    if boson_peak:
        ax.set_xlabel('$g(\omega)/\omega^2$', fontsize='x-large')
    else:
        ax.set_xlabel('$g(\omega)$', fontsize='x-large')
    # ax.set_xticklabels([])
    ax.set_ylabel('Frequency ({})'.format(unit), fontsize='x-large')
    ax.set_xlim((doslim_low,doslim_up))
    plt.subplots_adjust(left=0.40, bottom=0.25, right=0.60, top=0.752, wspace=0.2, hspace=0.2)
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    ax.xaxis.set_major_locator(plt.MaxNLocator(1))
    plt.grid(alpha=0.4)
    if legend_bool:
        plt.legend(fontsize='large', handlelength=1).set_draggable(True)
    else:
        plt.legend().set_visible(False)

def plot_chem_DOS(
    plt,
    f,
    ADOS,
    species_arr,
    unit='THz',
    freqlim_low=None,
    freqlim_up=None,
    doslim_low=None,
    doslim_up=None,
    lc_list=None,
    legend_bool=True,
    boson_peak=False,
    ):
    """
    f (arr)
    ADOS (arr)
    species_arr (arr) :
    unit (str) : 'THz' and 'meV' are implimented.
    freqlim_low (float or None) : 
    freqlim_up  (float or None) : 
    """
    ## Preprocess
    species_arr = np.array(species_arr)
    # Species
    unique_spec = np.unique(species_arr)
    print(unique_spec)
    atomic_index = np.arange(len(species_arr))
    pdos_indices = []
    for spec in unique_spec:
        pdos_indices.append(atomic_index[species_arr==spec])
    # Gather along species PDOS_list.shape=(len(unique_spec), 3, len(f))
    PDOS_list = []
    for spec_i in range(len(unique_spec)):
        PDOS_list.append(np.sum(ADOS[pdos_indices[spec_i]], axis=0))
    if not lc_list:
        lc_list = [None]*len(unique_spec)
    elif len(lc_list) != len(unique_spec):
        raise ValueError('len of    lc_list=={}   and    unique_spec=={}    must be same.'.format(lc_list, unique_spec))
    PDOS_list = np.array(PDOS_list)


    ## Chem plot
    fig, ax = plt.subplots()
    #
    for spec_i in range(len(unique_spec)):
        ax.plot(np.sum(PDOS_list[spec_i], axis=0), f, label=unique_spec[spec_i], c=lc_list[spec_i])
    ax.plot(np.sum(np.sum(PDOS_list, axis=0), axis=0), f, color='k')
    ax.set_ylim((freqlim_low, freqlim_up))
    if boson_peak:
        ax.set_xlabel('$g(\omega)/\omega^2$', fontsize='x-large')
    else:
        ax.set_xlabel('$g(\omega)$', fontsize='x-large')
    # ax.set_xticklabels([])
    ax.set_ylabel('Frequency ({})'.format(unit), fontsize='x-large')
    # ax.fill_between(np.sum(PDOS_list,axis=0), f, color='k', alpha=0.3)
    ax.set_xlim((doslim_low, doslim_up))
    plt.subplots_adjust(left=0.40, bottom=0.25, right=0.60, top=0.752, wspace=0.2, hspace=0.2)
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    ax.xaxis.set_major_locator(plt.MaxNLocator(1))
    plt.grid(alpha=0.4)
    if legend_bool:
        plt.legend(fontsize='large', handlelength=1).set_draggable(True)
    else:
        plt.legend().set_visible(False)

    ## Chemical-directional plot
    for spec_i in range(len(unique_spec)):
        fig, ax = plt.subplots()
        ax.plot(PDOS_list[spec_i, 0], f, label='x', c='g')
        ax.plot(PDOS_list[spec_i, 1], f, label='y', c='r')
        ax.plot(PDOS_list[spec_i, 2], f, label='z', c='b')
        ax.plot(np.sum(PDOS_list[spec_i], axis=0), f, c='k')
        ax.set_ylim((freqlim_low, freqlim_up))
        if boson_peak:
            ax.set_xlabel('$g(\omega)/\omega^2$', fontsize='x-large')
        else:
            ax.set_xlabel('$g(\omega)$', fontsize='x-large')
        # ax.set_xticklabels([])
        ax.set_ylabel('Frequency ({})'.format(unit), fontsize='x-large')
        # ax.fill_between(np.sum(PDOS_list,axis=0), f, color='k', alpha=0.3)
        ax.set_xlim((doslim_low, doslim_up))
        plt.subplots_adjust(left=0.40, bottom=0.25, right=0.60, top=0.752, wspace=0.2, hspace=0.2)
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        ax.xaxis.set_major_locator(plt.MaxNLocator(1))
        plt.grid(alpha=0.4)
        if legend_bool:
            plt.legend(fontsize='large', handlelength=1).set_draggable(True)
        else:
            plt.legend().set_visible(False)
        plt.title(unique_spec[spec_i])

    ## Directional-chemical plot
    direcs = ['x', 'y', 'z']
    for direc_i in range(len(direcs)):
        fig, ax = plt.subplots()
        #
        for spec_i in range(len(unique_spec)):
            ax.plot(PDOS_list[spec_i, direc_i], f, label=unique_spec[spec_i], c=lc_list[spec_i])
        ax.plot(np.sum(PDOS_list[:, direc_i], axis=0), f, c='k')
        ax.set_ylim((freqlim_low, freqlim_up))
        if boson_peak:
            ax.set_xlabel('$g(\omega)/\omega^2$', fontsize='x-large')
        else:
            ax.set_xlabel('$g(\omega)$', fontsize='x-large')
        # ax.set_xticklabels([])
        ax.set_ylabel('Frequency ({})'.format(unit), fontsize='x-large')
        # ax.fill_between(np.sum(PDOS_list,axis=0), f, color='k', alpha=0.3)
        ax.set_xlim((doslim_low, doslim_up))
        plt.subplots_adjust(left=0.40, bottom=0.25, right=0.60, top=0.752, wspace=0.2, hspace=0.2)
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        ax.xaxis.set_major_locator(plt.MaxNLocator(1))
        plt.grid(alpha=0.4)
        if legend_bool:
            plt.legend(fontsize='large', handlelength=1).set_draggable(True)
        else:
            plt.legend().set_visible(False)
        plt.title(direcs[direc_i])

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    This code will give you the phonon-(partial/total)DOS from MD trajectory.
    """)
    # Positional arguments
    parser.add_argument('inp_file_list', type=str, nargs='+', help='ASE readable atoms list file name. When multiple input files are provided, DOS will averaged.')
    parser.add_argument('dt', type=float, help='Time interval between images selected in unit of picosec.')
    # Optional arguments
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='Image range following python convention. default=":" (e.g.) -n :1000:10')
    parser.add_argument('-p', '--chemical_DOS', action='store_true', help='Plot species-wise partial DOS.')
    parser.add_argument('-g', '--gsmear', type=float, default=0., help='Width(simga, STD) of Gaussian smearing in frequency unit. Zero means no smearing. [default: 0]')
    parser.add_argument('-l', '--freqlim_low', type=float, default=0., help='Set frequency lower limit for plot. [default: 0]')
    parser.add_argument('-u', '--freqlim_up', type=float, default=None, help='Set frequency upper limit for plot. Auto detect as default.')
    parser.add_argument('-m', '--doslim_low', type=float, default=1e-5, help='Set DOS lower limit for plot. [default: 1e-5]')
    parser.add_argument('-v', '--doslim_up', type=float, default=None, help='Set DOS upper limit for plot. Auto detect as default.')
    parser.add_argument('-s', '--dont_save', dest='save_bool', action='store_false', help='If provided, ADOS arrays will not be saved. Default: Save array')
    parser.add_argument('-o', '--dont_load', dest='load_bool', action='store_false', help='If provided, ADOS arrays will not be loaded. Default: Load if possible')
    parser.add_argument('-t', '--dont_plot', dest='plot_bool', action='store_false', help='Do not plot, if provided. [default: Plot].')
    parser.add_argument('-f', '--DOS_factor', type=float, default=1., help='DOS multiply factor. As default, integral of total DOS is 1. (cf. In case of phonopy, 3N, where N is number of atoms in a primitive cell.)')
    parser.add_argument('-b', '--no_legend', dest='legend_bool', action='store_false', help='No legend plot. [default: True for partial DOS].')
    parser.add_argument('-c', '--lc_list', type=str, nargs='+', default=None, help='Line color list. For partial DOS, len(lc_list) == len(np.unique(chem)) [default: automatic].')
    parser.add_argument('-k', '--boson_peak', action='store_true', help='Plot g(f)/f**2 to seek boson peak. [default: False].')
    parser.add_argument('-e', '--meV_unit', action='store_true', help='Use meV unit [default: THz]')
    return parser.parse_args()

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
    print('This code will give you the phonon-(partial/total)DOS from MD trajectory.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    ## Read input params
    dt          = args.dt
    cdos_bool   = args.chemical_DOS
    gsmear_std  = args.gsmear
    freqlim_low = args.freqlim_low
    freqlim_up  = args.freqlim_up
    doslim_low     = args.doslim_low
    doslim_up      = args.doslim_up
    plot_bool   = args.plot_bool
    if args.meV_unit:
        unit = 'meV'
    else:
        unit = 'THz'
    from ss_util import str_slice_to_list
    slice_list = str_slice_to_list(args.image_slice)

    ## Main loop
    from ase.io import read
    ADOS_list = []
    for i in range(len(args.inp_file_list)):
        f_name = args.inp_file_list[i].split('/')[-1]
        f_path = args.inp_file_list[i][:-(len(f_name)+1)]
        npz_name = './{}/vel2dos-saved/{}_dt{}_img{}-{}-{}.npz'.format(f_path, f_name, dt, slice_list[0], slice_list[1], slice_list[2])
        try:
            args.load_bool = True
            npz = np.load(npz_name)
            f    = npz['f']
            ADOS = npz['ADOS']
        except:
            # log
            if i % 10 == 0:
                print("Getting {}-th file's VACF".format(i))
            # Read file
            alist = read(args.inp_file_list[i], args.image_slice)
            # Get atomic masses
            atomic_mass_arr = alist[0].get_masses()
            # Gather velocities and temperatures
            v_arr = []
            temp_arr = []
            for atoms in alist:
                v_arr.append(atoms.get_velocities())
                temp_arr.append(atoms.get_temperature())
            # Get average temperature
            average_temp = np.mean(temp_arr)
            # Get VACF
            f, ADOS = velocity2phononPDOS(atomic_mass_arr, average_temp, v_arr, dt)
            if args.save_bool:
                from subprocess import call
                call('mkdir -p ./{}/vel2dos-saved'.format(f_path), shell=True)
                np.savez(npz_name, f=f, ADOS=ADOS)
        else:
            if i % 100 == 0:
                print("Loading {}-th file's PDOS".format(i))
        ADOS_list.append(ADOS)
        # # For debugging
        # if np.shape(np.array(ADOS)) != (8, 3, 1750):
            # print(np.shape(np.array(ADOS)))
            # print(args.inp_file_list[i])
    # ## write txt file
    # with open('tmp.txt', 'w') as txt:
        # sum = np.sum(np.sum(ADOS, axis=0), axis=0)
        # for j in range(ADOS.shape[-1]):
            # txt.write('{:10.5f}'.format(f[j]))
            # txt.write('        ')
            # txt.write('{:10.5e}'.format(sum[j]))
            # txt.write('\n')
    # # Print area (Normalization test)
    # d_f = 1. / (dt * (len(f)-1)*2)
    # print(np.sum(ADOS) * d_f)
    #
    # Frequency scaling
    if unit is 'THz':
        pass
    elif unit is 'meV':
        from phonopy import units
        f *= units.THztoEv * 1e3
    d_f = f[1] - f[0]
    ADOS_list = np.array(ADOS_list)

    if not gsmear_std == 0:
        print(' Gaussian smearing...')
        from scipy.ndimage.filters import gaussian_filter1d
        ADOS_list = gaussian_filter1d(ADOS_list, gsmear_std /d_f)

    if plot_bool:
        ## Averaging ADOS shape=(len(atoms), 3, len(f))
        ADOS = np.mean(ADOS_list, axis=0) * args.DOS_factor
        ## Boson peak
        if args.boson_peak:
            ADOS /= f**2
        ## Plot
        from matplotlib import pyplot as plt
        font = {'family':'Arial'}
        plt.rc('font', **font)
        plot_total_DOS(
            plt,
            f,
            ADOS,
            unit='THz',
            freqlim_low=freqlim_low,
            freqlim_up=freqlim_up,
            doslim_low=doslim_low,
            doslim_up=doslim_up,
            legend_bool=args.legend_bool,
            boson_peak=args.boson_peak,
            )
        plot_direc_DOS(
            plt,
            f,
            ADOS,
            unit='THz',
            freqlim_low=freqlim_low,
            freqlim_up=freqlim_up,
            doslim_low=doslim_low,
            doslim_up=doslim_up,
            legend_bool=args.legend_bool,
            boson_peak=args.boson_peak,
            )
        if cdos_bool:
            atoms = read(args.inp_file_list[i], 0)
            plot_chem_DOS(
                plt,
                f,
                ADOS,
                atoms.get_chemical_symbols(),
                unit='THz',
                freqlim_low=freqlim_low,
                freqlim_up=freqlim_up,
                doslim_low=doslim_low,
                doslim_up=doslim_up,
                lc_list=args.lc_list,
                legend_bool=args.legend_bool,
                boson_peak=args.boson_peak,
                )
        plt.show()
