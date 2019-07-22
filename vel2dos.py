#!/usr/bin/env python

import numpy as np

def velocity2phononPDOS(atomic_mass_arr, average_temp, velocity_arr, d_t):
    """
    atomic_mass_arr (array) : Array of atomic masses of all atoms in unit cell in amu (atomic mass unit).
    average_temp (float)    : Average temperature of provided ensemble trajectory in kelvin.
    velocity_arr (array)    : All of velocities in trajectory of all atoms. Array shape == ( len(traj), len(atoms), 3 )
    d_t (float)             : Time interval of images provided in (ps) unit (ps^(-1) --> THz)

    ================================================
                        Theory                      
    ================================================
    This fomula gives phonon density of states from 
    velocities of atoms in trajectory.

                        N             _             
             1/2 * (summation) m_n * |v_n(f)|^2     
                       n=1                          
    g(f) = --------------------------------------   
                        3/2 * NkT                   
                                                    
    _
    v_n(f) : velocity vector in frequency space
    m_n    : Atomic mass of atom n
    N      : Number of atoms in unit cell
    n      : Atom index
    k      : Boltzmann constant
    T      : Temperature
    g(f)   : Phonon density of state

    """
    # Get params
    image_num = len(velocity_arr)
    natoms = len(atomic_mass_arr)
    t_init = 0.
    t_fin = (image_num-1) * d_t
    # Get discrete time domain
    t = np.arange(t_init, t_fin+d_t, d_t)
    # Get frequency domain (Use Arfken representation)
    f = np.fft.fftfreq(image_num) / d_t
    # Get positive part only
    f = f[:int(round(len(t)/2.))]
    # Reshape
    velocity_arr = np.reshape(velocity_arr, (image_num, natoms*3)).T
    # Fourier transform for each DOF
    v_f_arr = []
    for DOF_i in range(len(velocity_arr)):
        v_f = np.fft.fft(velocity_arr[DOF_i])
        # Get positive part only
        v_f = v_f[:int(round(len(t)/2.))]
        v_f_arr.append(v_f)
    v_f_arr = np.array(v_f_arr)
    # calculate the formula
    from ase import units
    # ADOS =  (np.repeat(atomic_mass_arr, 3) / 3. / natoms / units.kB / average_temp * np.square(np.abs(v_f_arr)).T).T
    ADOS = np.square(np.abs(v_f_arr)) / 3. / natoms / units.kB / average_temp

    return f, ADOS
        

def plot_total_DOS(
    f,
    DOS,
    unit='THz',
    freqlim_low=None,
    freqlim_up=None,
    flip_xy=True,
    ):
    """
    f (arr)
    DOS (arr)
    unit (str) : 'THz' and 'meV' are implimented.
    freqlim_low (float or None) : 
    freqlim_up  (float or None) : 
    flip_xy (bool)              :
    """
    ## Preprocess
    # Frequency scaling
    if unit is 'THz':
        pass
    elif unit is 'meV':
        from phonopy import units
        f *= units.THztoEv * 1e3

    ## Plot
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots()
    if flip_xy:
        ax.plot(DOS, f)
        ax.set_ylim((freqlim_low, freqlim_up))
        ax.set_xlabel('DOS (a.u.)', fontsize='x-large')
        ax.set_ylabel('Frequency ({})'.format(unit), fontsize='x-large')
        ax.set_xlim((0,None))
    else:
        ax.plot(f, DOS)
        ax.set_xlim((freqlim_low, freqlim_up))
        ax.set_xlabel('Frequency ({})'.format(unit), fontsize='x-large')
        ax.set_ylabel('DOS (a.u.)', fontsize='x-large')
        ax.set_ylim((0,None))
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.grid(alpha=0.2)
    plt.subplots_adjust(left=0.35, bottom=0.15, right=0.60, top=0.95, wspace=0.2, hspace=0.2)
    plt.show()

def plot_partial_DOS(
    f,
    ADOS,
    species_arr,
    unit='THz',
    freqlim_low=None,
    freqlim_up=None,
    flip_xy=True,
    ):
    """
    f (arr)
    ADOS (arr)
    species_arr (arr) :
    unit (str) : 'THz' and 'meV' are implimented.
    freqlim_low (float or None) : 
    freqlim_up  (float or None) : 
    flip_xy (bool)              :
    """
    ## Preprocess
    species_arr = np.array(species_arr)
    # Frequency scaling
    if unit is 'THz':
        pass
    elif unit is 'meV':
        from phonopy import units
        f *= units.THztoEv * 1e3
    # Species
    unique_spec = np.unique(species_arr)
    atomic_index = np.arange(len(species_arr))
    pdos_indices = []
    for spec in unique_spec:
        pdos_indices.append(atomic_index[species_arr==spec])
    # Gather along species
    PDOS_list = []
    for spec_i in range(len(unique_spec)):
        PDOS_list.append(np.sum(ADOS[pdos_indices[spec_i]],axis=0))

    ## Plot
    from matplotlib import pyplot as plt
    fig, ax = plt.subplots()
    if flip_xy:
        for spec_i in range(len(unique_spec)):
            ax.plot(PDOS_list[spec_i], f, label=unique_spec[spec_i])
        ax.set_ylim((freqlim_low, freqlim_up))
        ax.set_xlabel('DOS (a.u.)', fontsize='x-large')
        ax.set_ylabel('Frequency ({})'.format(unit), fontsize='x-large')
        ax.fill_between(np.sum(PDOS_list,axis=0), f, color='k', alpha=0.3)
        ax.set_xlim((0,None))
    else:
        for spec_i in range(len(unique_spec)):
            ax.plot(f, PDOS_list[spec_i], label=unique_spec[spec_i])
        ax.set_xlim((freqlim_low, freqlim_up))
        ax.set_xlabel('Frequency ({})'.format(unit), fontsize='x-large')
        ax.set_ylabel('DOS (a.u.)', fontsize='x-large')
        ax.fill_betweenx(np.sum(PDOS_list,axis=0), f, color='k', alpha=0.3)
        ax.set_ylim((0,None))
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.grid(alpha=0.2)
    plt.legend()
    plt.subplots_adjust(left=0.35, bottom=0.15, right=0.60, top=0.95, wspace=0.2, hspace=0.2)
    plt.show()

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    This code will give you the phonon-(partial/total)DOS from MD trajectory.
    """)
    # Positional arguments
    parser.add_argument('inp_file', type=str, help='ASE readable atoms list file name.')
    parser.add_argument('d_t', type=float, help='Time interval between images selected in unit of picosec.')
    # Optional arguments
    parser.add_argument('-n', '--image_range', type=str, default=':', help='Image range following python convention. default=":" (e.g.) -n :1000:10')
    parser.add_argument('-p', '--partial_DOS', action='store_true', help='If activated, return partial DOS. (If not, total DOS as default)')
    parser.add_argument('-l', '--freqlim_low', type=float, default=0., help='Set frequency lower limit for plot. Zero as default.')
    parser.add_argument('-u', '--freqlim_up', type=float, default=None, help='Set frequency upper limit for plot. Auto detect as default.')
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
    d_t = args.d_t
    pdos_bool = args.partial_DOS
    freqlim_low = args.freqlim_low
    freqlim_up = args.freqlim_up
    # inp_file
    from ase.io import read
    alist = read(args.inp_file, args.image_range)

    ## Preprocess
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
    # Main
    f, ADOS = velocity2phononPDOS(atomic_mass_arr, average_temp, v_arr, d_t)
    if pdos_bool:
        plot_partial_DOS(
            f,
            ADOS,
            atoms.get_chemical_symbols(),
            unit='THz',
            freqlim_low=freqlim_low,
            freqlim_up=freqlim_up,
            flip_xy=True,
            )
    else:
        plot_total_DOS(
            f,
            np.sum(ADOS, axis=0),
            unit='THz',
            freqlim_low=freqlim_low,
            freqlim_up=freqlim_up,
            flip_xy=True,
            )
