        #####    Code by YJ Choi of CNPL, Dep. of Phys. POSTECH, Korea    #####
        #####                    ssrokyz@postech.ac.kr                    #####
        ##### Some codes are copied from phonopy. Code for Phonopy v2.1.2 #####

import numpy as np
from subprocess import call

#def get_fc(directions, delta):
#    """ Get forces, FCs """
#    import subprocess as sp
#
#    image_num = len(directions)
#    pwd = (sp.check_output("pwd"))[:-1]
#    for i in range(image_num):
#        ndir = str(delta) + "-" + str(directions[i][0]) + "-" + \
#               str(directions[i][1]) + str(directions[i][2]) + str(directions[i][3])

def bu_and_mkdir(calc_dir, ndir):
    call(['rm -rf '+calc_dir+'/bu-'+ndir], shell=True)
    call(['mv '+calc_dir+'/'+ndir+' '+calc_dir+'/bu-'+ndir], shell=True)
    call(['mkdir -p '+calc_dir+'/'+ndir], shell=True)

def get_subdir_name(order, disp): 
    sign = []
    for i in range(3):
        num = np.sign(disp[order][i+1], dtype=np.float)
        if num == 1.:
            sign.append('+')
        elif num == 0.:
            sign.append('0')
        elif num == -1.:
            sign.append('-')
        else:
            raise ValueError
    return "pos"+str(order).zfill(3) +"_atom"+str(disp[order][0]).zfill(3) +"_direc"+sign[0]+sign[1]+sign[2]

def calc_vasp(phonon, disp, calc_dir, F_0_correction, amp_calc):
    """ Calculate Force Constant with Vasp """
    forces = []
    ndir = ''
    # import io
    from phonopy.interface.vasp import VasprunxmlExpat
    # Make KPOINTS file
    from kpoints_gen import get_grid_num
    kgrid = get_grid_num(phonon.get_supercell().cell, precision=55.)
    from ase.io import read
    for i in range(len(disp)):
        print(' >>> Starting {:d}-th image calculation <<< '.format(i).center(120))
        ndir_prev = ndir
        ndir = get_subdir_name(i, disp)
        bu_and_mkdir(calc_dir, ndir)
        call(['cp INCAR POTCAR '+calc_dir+'/poscars/POSCAR-'+str(i).zfill(3)+' '+calc_dir+'/'+ndir], shell=True)
        with open(calc_dir+'/'+ndir+'/KPOINTS', 'w') as txt:
            txt.write('KPOINTS\n0\nGamma\n{} {} {}\n0 0 0'.format(kgrid[0], kgrid[1], kgrid[2]))
        call(['cp POSCAR-'+str(i).zfill(3)+' POSCAR'], cwd=calc_dir+'/'+ndir, shell=True)
        call(['cp WAVECAR CHGCAR ../'+ndir], cwd=calc_dir+'/'+ndir_prev, shell=True)
        call(['mpiexec.hydra -np $NSLOTS vasp_std > out'], cwd = calc_dir+'/'+ndir, shell=True)
        result = read(calc_dir+'/'+ndir+'/vasprun.xml', 0)
        # with io.open(calc_dir+'/'+ndir+'/vasprun.xml', 'rb') as xmlfile:
            # vasprun = VasprunxmlExpat(xmlfile)
            # if vasprun.parse():
                # force_now = (vasprun.get_forces()).tolist()
#                   epsilon = vasprun.get_epsilon()
#                   borns = vasprun.get_born()
#                   lattice = vasprun.get_lattice()[-1]
#                   points = vasprun.get_points()[-1]
#                   symbols = vasprun.get_symbols()
#                   unit_cell = PhonopyAtoms(symbols=symbols,
#                                            scaled_positions=points,
#                                            cell=lattice)
                # forces.extend(force_now)
        if i != 0:
            forces.extend(result.get_forces())
    phonon.set_forces(np.array(forces))

def calc_dpmd(phonon, disp, calc_dir, F_0_correction, amp_calc):
    """ Calculate Force Constant with DPMD """
    forces = []
    from ase.io.lammpsrun import read_lammps_dump as read_dump
    from ase.io import write
    for i in range(len(disp)):
        print(' >>> Starting {:d}-th image calculation <<< '.format(i).center(120))
        ndir = get_subdir_name(i, disp)
        bu_and_mkdir(calc_dir, ndir)
        call(['cp frozen_model.pb input.in '+calc_dir+'/poscars/POSCAR-'+str(i).zfill(3)+' '+calc_dir+'/'+ndir], shell=True)
        call(['lmp-pos2lmp.awk POSCAR-'+str(i).zfill(3)+' > structure.in'], cwd = calc_dir+'/'+ndir, shell = True)
        call(['mpiexec.hydra -np $NSLOTS lmp_mpi -in input.in > out'], cwd = calc_dir+'/'+ndir, shell = True)
        atoms = read_dump(calc_dir+'/'+ndir+'/out.dump', index=0, order=True)
        if i == 0:
            F_0 = atoms.get_forces(apply_constraint=False)
        else:
            if F_0_correction:
                atoms._calc.results['forces'] -= F_0
            forces.append(atoms.get_forces(apply_constraint=False).tolist())
        write(calc_dir+'/'+ndir+'/result.traj', atoms, 'traj')
    phonon.set_forces(np.array(forces))

def calc_amp(phonon, disp, calc_dir, F_0_correction, amp_calc):
    """ Calculate Force Constant with AMP """
    numeric_F_dx=0.001
    parallel=True
    forces = []
    from ase.io import read, write
    for i in range(len(disp)):
        print(' >>> Starting {:d}-th image calculation <<< '.format(i).center(120))
        ndir = get_subdir_name(i, disp)
        bu_and_mkdir(calc_dir, ndir)
        call(['cp '+calc_dir+'/poscars/POSCAR-'+str(i+1).zfill(3)+' '+calc_dir+'/'+ndir], shell=True)

        ########### calculate forces & atomic energies with amp ############
        atoms = read(calc_dir+'/'+ndir+'/POSCAR-'+str(i+1).zfill(3), format = 'vasp')
        atoms.set_pbc(True)
        atoms.set_calculator(amp_calc)

        #********** numerical force must precede ***********
        # force_now = [calc.calculate_numerical_forces(
            # atoms,
            # d = numeric_F_dx,
            # parallel = parallel,
            # )] # alternative
        if i!=0:
            forces.append(atoms.get_forces(apply_constraint=False).tolist()) # alternative
        write(calc_dir+'/'+ndir+'/result.traj', atoms, 'traj')
    phonon.set_forces(np.array(forces))

def calc_amp_tf(phonon, disp, calc_dir, F_0_correction, amp_calc):
    """ Calculate Force Constant with AMP with tensorflow """
    numeric_F_dx=0.001
    parallel=True
    forces = []
    from ase.io import read, write
    for i in range(len(disp)):
        print(' >>> Starting {:d}-th image calculation <<< '.format(i).center(120))
        ndir = get_subdir_name(i, disp)
        bu_and_mkdir(calc_dir, ndir)
        call(['cp '+calc_dir+'/poscars/POSCAR-'+str(i+1).zfill(3)+' '+calc_dir+'/'+ndir], shell=True)

        ########### calculate forces & atomic energies with amp ############
        atoms = read(calc_dir+'/'+ndir+'/POSCAR-'+str(i+1).zfill(3), format = 'vasp')
        atoms.set_pbc(True)
        atoms.set_calculator(amp_calc)

        #********** numerical force must precede ***********
        if i!=0:
            forces.append(calc.calculate_numerical_forces(atoms, d = numeric_F_dx, parallel = parallel))
        atoms._calc.results['forces'] = forces[-1]
        write(calc_dir+'/'+ndir+'/result.traj', atoms, 'traj')
    phonon.set_forces(np.array(forces))

def calc_amp_tf_bunch(phonon, disp, calc_dir, F_0_correction, amp_calc):
    """ Calculate Force Constant with AMP with tensorflow (fast version) """
    numeric_F_dx=0.001
    parallel=True
    forces = []
    from ase.io import read, write
    for i in range(len(disp)):
        print(' >>> Starting {:d}-th image calculation <<< '.format(i).center(120))
        ndir = get_subdir_name(i, disp)
        bu_and_mkdir(calc_dir, ndir)
        call(['cp '+calc_dir+'/poscars/POSCAR-'+str(i+1).zfill(3)+' '+calc_dir+'/'+ndir], shell=True)

        ########### calculate forces & atomic energies with amp ############
        atoms = read(calc_dir+'/'+ndir+'/POSCAR-'+str(i+1).zfill(3), format = 'vasp')
        atoms.set_pbc(True)
        atoms.set_calculator(amp_calc)

        #force_now = [atoms.get_forces(apply_constraint=False).tolist()] # alternative
        #stress_now = calc.calculate_numerical_stress(atoms) # just in case
        #********** energy calculation ***********
        energy_now = atoms.get_potential_energy()
        energies_now = atoms.get_potential_energies()
        #********** force information restore ***********
        atoms._calc.results['forces'] = np.asarray(force_now[0])
        if verbose:
            print(energies_now)
            print(atoms._calc.results['forces'])
        if i!=0:
            forces.append(calc.calculate_numerical_forces(atoms, d = numeric_F_dx, parallel = parallel))
        atoms._calc.results['forces'] = forces[-1]
        write(calc_dir+'/'+ndir+'/result.traj', atoms, 'traj')
    phonon.set_forces(np.array(forces))

def calc_phonon(calculator, phonon, acoustic_sum_rule=True, F_0_correction=False, verbose=False, amp_calc=None):
    """
    calc -- Specify calculator. One of these. [vasp, dpmd, amp, amp_tf, amp_tf_bunch]
    phonon -- Phonopy phonon object.
    """
    import sys
    import pickle as pckl
    if verbose:
        np.set_printoptions(threshold=sys.maxsize)
    delta    = np.linalg.norm(phonon.get_displacements()[0][1:4])
    disp     = [[0,0,0,0]] + phonon.get_displacements()
    job_name = 'x{}{}{}_d{:5.3f}_sym{}'.format(
        phonon.get_supercell_matrix()[0][0],
        phonon.get_supercell_matrix()[1][1],
        phonon.get_supercell_matrix()[2][2],
        delta,
        phonon._is_symmetry,
        )
    # Define names
    pckl_name = job_name+'.bin'
    calc_dir = './calcs/'+job_name
    # Environment preset
    call(['rm -rf '+calc_dir+'/poscars'], shell=True)
    call(['mkdir -p '+calc_dir+'/poscars'], shell=True)
    call(['cp SPOSCAR POSCAR-000'], shell=True)
    call(['mv POSCAR-* SPOSCAR '+calc_dir+'/poscars/'], shell=True)

    # Load saved pickle file
    print('')
    try:
        phonon = pckl.load(open(pckl_name, 'rb'))
        if phonon.get_force_constants() is None:
            raise ValueError(' Error:: Phonon object in pickle file is not containing force constants info')
    except:
        print('<<  CAUTION  >>  Fail to load pickle file.                <<  CAUTION  >>'.center(120))
        print(('<<  CAUTION  >>'+':: Expected file name ::'.center(43)+'<<  CAUTION  >>').center(120))
        print(('<<  CAUTION  >>'+pckl_name.center(43)+'<<  CAUTION  >>').center(120))
        do_calc = True
    else:
        print('>>>>>>> Pickle file "{}" has been loaded. <<<<<<<<'.format(pckl_name).center(120))
        do_calc = False
    print('')
    if do_calc:
        if calculator == 'vasp': 
            calc = calc_vasp
        elif calculator == 'dpmd':
            calc = calc_dpmd
        elif calculator == 'amp':
            calc = calc_amp
        elif calculator == 'amp_tf':
            calc = calc_amp_tf
        elif calculator == 'amp_tf_bunch':
            calc = calc_amp_tf_bunch
        else:
            raise ValueError('Unknown calculator ({}) has been provided.'.format(calculator))
        print(">>>>>>> {} calculation will be carried out. <<<<<<<<".format(calculator).center(120))
        calc(phonon, disp, calc_dir, F_0_correction, amp_calc)
            
        if verbose:
            print('\n\n'+'forces'+'\n'+str(phonon.forces))
        # Produce fc
        phonon.produce_force_constants()
        if acoustic_sum_rule:
            phonon.symmetrize_force_constants()

        if verbose:
            print('\n\ndisplacement_dataset =>\n\n')
            print(phonon.get_displacement_dataset())
            print('\n\nforce_constants =>\n\n')
            print(phonon.get_force_constants())
        pckl.dump(phonon, open(pckl_name, 'wb'), protocol=2)
    return phonon

def make_band(path, N_q):
    bands = []
    for i in range(len(path)):
        q_start  = np.array(path[i][0])
        q_end    = np.array(path[i][1])
        band = []
        for j in range(N_q+1):
            band.append(q_start + (q_end - q_start) / N_q * j)
        bands.append(band)
    return bands

# def plot_band(phonon, labels=None):
    # import matplotlib.pyplot as plt
    # if labels:
        # from matplotlib import rc
        # rc('font',**{'family':'serif','sans-serif':['Times']})
        # rc('text', usetex=False)

    # fig, ax = plt.subplots()
    # ax.xaxis.set_ticks_position('both')
    # ax.yaxis.set_ticks_position('both')
    # ax.xaxis.set_tick_params(which='both', direction='in')
    # ax.yaxis.set_tick_params(which='both', direction='in')

    # phonon._band_structure.plot(plt, labels=labels)
    # return plt

def set_projection(phonon, proj_eigvec):
    self = phonon._band_structure
    self._projections = {}
    self._proj_freq = {}

    for key in proj_eigvec.keys():
        proj_tmp = []
        freq_tmp = []
        for _path in self._paths:
            proj_tmp2 = []
            freq_tmp2 = []
            for _q in _path:
                freq, eigvec = phonon.get_frequencies_with_eigenvectors(_q)
                eigvec = eigvec.T   # numpy.linalg.eigh returns transposed eigen vector
                #### Gather
                proj_tmp2.append(eigvec)
                freq_tmp2.append(freq)
            #### Gather
            proj_tmp.append(proj_tmp2)
            freq_tmp.append(freq_tmp2)
        #### Transform to array
        proj_tmp = np.array(proj_tmp)
        freq_tmp = np.array(freq_tmp)
        #### calculate similarity
        proj_tmp = np.square(np.absolute(np.squeeze(np.matmul(np.conj(proj_tmp), np.expand_dims(proj_eigvec[key], axis=-1)))))
        self._projections[key] = proj_tmp
        self._proj_freq[key]   = freq_tmp

def bs_plot(self, plt, ax, proj_size_factor, proj_colors, proj_alpha, reverse_seq, plot_legend, labels=None):
    if self._projections is not None:
        #### Define key list
        key_list = list(self._projections.keys())
        if reverse_seq:
            key_list.reverse()
        #### Pick colors
        proj_colors = proj_colors[len(proj_colors)-len(key_list):]
        if reverse_seq:
            proj_colors.reverse()
        #### 
        legend = []
        #### Iter for projector eigenvectors
        for key in key_list:
            #### Iter for q_path fragments
            # for distances, frequencies, projections, proj_freq in zip(self._distances,
            for distances, projections, freq in zip(self._distances,
                                                         self._projections[key],
                                                         self._proj_freq[key]):
                #### Iter for band lines
                for i in range(len(freq.T)):
                    plt.plot(distances, freq.T[i], 'k-')
                    legend_tmp = plt.scatter(
                        distances,
                        freq.T[i],
                        proj_size_factor * projections.T[i],
                        proj_colors[-1],
                        alpha=proj_alpha,
                        edgecolors='none',
                        label=key,
                        )
            #### Gather just the one sample legend
            legend.append(legend_tmp)
            #### Throw away used color
            proj_colors.pop()
        #### Legend plot
        if reverse_seq:
            legend.reverse(); key_list.reverse()
        if plot_legend:
            plt.legend(legend, key_list, scatterpoints = 1, fontsize='xx-large', loc='upper right')
    else:
        for distances, frequencies in zip(self._distances,
                                          self._frequencies):
            for freqs in frequencies.T:
                plt.plot(distances, freqs, 'k-')

    plt.ylabel('Frequency')
    plt.xlabel('Wave vector')

    if labels and len(labels) == len(self._special_points):
        plt.xticks(self._special_points, labels) # ssrokyz
    else:
        plt.xticks(self._special_points, [''] * len(self._special_points))
    plt.xlim(0, self._distance)
    plt.axhline(y=0, linestyle=':', linewidth=0.5, color='k')

def plot_band_and_dos(
    phonon,
    pdos_indices     = None,
    labels           = None,
    unit             = 'THz',
    proj_eigvec      = None,
    proj_size_factor = 400.,
    proj_colors      = ['r', 'g', 'b', 'c', 'm', 'y'],
    proj_alpha       = 0.5,
    plot_legend      = False,
    ylim_lower       = None,
    ylim_upper       = None,
    reverse_seq      = False,
    ):
    """
    proj_eigvec = (dict) = Eigenvectors that will be used for projection. Keys of dict will be used as label for pyplot.
    proj_colors = (list) = Use colors sequently. Duplication is possible.
    reverse_seq = (bool) = Change sequence of plotting scatters. (Decides which one goes up. Only esthetic purpose.)
    """

    #### Variable setting
    self = phonon
    proj_colors.reverse()

    import matplotlib.pyplot as plt
    if labels:
        from matplotlib import rc
        rc('font',**{'family':'serif','sans-serif':['Times']})
        rc('text', usetex=True)

    import matplotlib.gridspec as gridspec
    plt.figure(figsize=(20, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[7,4,1])
    ax2 = plt.subplot(gs[0, 2])
    if pdos_indices is None:
        self._total_dos.plot(ax2,
                             ylabel="",
                             draw_grid=True,
                             flip_xy=True)
    else:
        self._pdos.plot(ax2,
                        indices=pdos_indices,
                        ylabel="",
                        draw_grid=True,
                        flip_xy=True)
    ax2.set_xlim((0, None))
    # ax2.set_title('DOS', fontsize=35)
    ax2.set_xlabel('')
    ax2.set_xticklabels([])
    plt.setp(ax2.get_yticklabels(), visible=False)

    ax1 = plt.subplot(gs[0, 1], sharey=ax2)

    #### projection
    self._band_structure._projections = None
    if proj_eigvec is not None:
        # dtype managing
        for key in proj_eigvec.keys():
            proj_eigvec[key] = np.array(proj_eigvec[key], dtype=np.complex128)
        set_projection(self, proj_eigvec)
    bs_plot(self._band_structure, plt, ax1, proj_size_factor, proj_colors, proj_alpha, reverse_seq, plot_legend, labels=labels)
    if unit == 'meV':
        plt.ylabel('Frequency(meV)', fontsize=35)
    elif unit == 'THz':
        plt.ylabel('Frequency(THz)', fontsize=35)
    # ax1.set_title('Phonon dispersion', fontsize=35)
    ax1.set_xlabel('')

    plt.yticks(fontsize=35)
    plt.xticks(fontsize=35)
    plt.grid(True)
    plt.subplots_adjust(wspace=0.0)
    plt.ylim(ylim_lower, ylim_upper)
    plt.tight_layout()

    return plt

    ################# old version backup (phonopy 1.13.2)
    # import matplotlib.pyplot as plt
    # import matplotlib.gridspec as gridspec
    # if labels:
        # from matplotlib import rc
        # rc('font',**{'family':'serif','sans-serif':['Times']})
        # rc('text', usetex=False)

    # #### Variable setting
    # proj_colors.reverse()

    # plt.figure(figsize=(10, 6))
    # gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1])
    # ax1 = plt.subplot(gs[0, 0])
    # ax1.xaxis.set_ticks_position('both')
    # ax1.yaxis.set_ticks_position('both')
    # ax1.xaxis.set_tick_params(which='both', direction='in')
    # ax1.yaxis.set_tick_params(which='both', direction='in')

    # #### projection
    # phonon._band_structure._projections = None
    # if proj_eigvec is not None:
        # # dtype managing
        # for key in proj_eigvec.keys():
            # proj_eigvec[key] = np.array(proj_eigvec[key], dtype=np.complex128)
        # set_projection(phonon, proj_eigvec)

    # bs_plot(phonon._band_structure, plt, ax1, proj_size_factor, proj_colors, proj_alpha, reverse_seq, labels=labels)
    # if unit == 'meV':
        # plt.ylabel('Frequency(meV)', fontsize=22)
    # elif unit == 'THz':
        # plt.ylabel('Frequency(THz)', fontsize=22)
    # plt.xlabel('')
    # plt.grid(True)
    # plt.title('Phonon dispersion', fontsize=24)
    # plt.yticks(fontsize=20)

    # ax2 = plt.subplot(gs[0, 1], sharey=ax1)
    # ax2.xaxis.set_ticks_position('both')
    # ax2.yaxis.set_ticks_position('both')
    # ax2.xaxis.set_tick_params(which='both', direction='in')
    # ax2.yaxis.set_tick_params(which='both', direction='in')
    # plt.subplots_adjust(wspace=0.08)
    # plt.setp(ax2.get_yticklabels(), visible=False)

    # if pdos_indices is None:
        # phonon._total_dos.plot(plt,
                               # ylabel="",
                               # draw_grid=True,
                               # flip_xy=True)
    # else:
        # phonon._pdos.plot(plt,
                          # indices=pdos_indices,
                          # ylabel="",
                          # draw_grid=True,
                          # flip_xy=True)

    # ax2.set_xlim((0, None))
    # plt.title('DOS', fontsize=24)
    # plt.xlabel('')
    # plt.xticks([])
    # plt.ylim(ylim_lower, ylim_upper)

    # return plt

def two_dos_plot(
    phonon_dict,
    color_dict = None,
    ylim_upper = None,
    ylim_lower = 0,
    ):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_tick_params(which='both', direction='in')
    ax.yaxis.set_tick_params(which='both', direction='in')

    for key in phonon_dict.keys():
        if color_dict is not None:
            phonon_dict[key]._total_dos.plot(plt, flip_xy=True, draw_grid=False, color=color_dict[key], legend=key)
        else:
            phonon_dict[key]._total_dos.plot(plt, flip_xy=True, draw_grid=False, color=None, legend=key)

    ax.set_ylim((ylim_lower, ylim_upper))
    plt.legend()

    return plt
