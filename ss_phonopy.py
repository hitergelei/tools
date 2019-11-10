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

def get_subdir_name(order, disp, prec=1e-5): 
    sign = []
    for i in range(3):
        num = disp[order][i+1]
        if num > prec:
            sign.append('+')
        elif num < -1 * prec:
            sign.append('-')
        else:
            sign.append('0')
    return "pos"+str(order).zfill(3) +"_atom"+str(disp[order][0]).zfill(3) +"_direc"+sign[0]+sign[1]+sign[2]

def calc_vasp(phonon, disp, calc_dir, F_0_correction, ase_calc):
    """ Calculate Force Constant with Vasp """
    forces = []
    ndir = ''
    # import io
    from phonopy.interface.vasp import VasprunxmlExpat
    # Make KPOINTS file
    from kpoints_gen import get_grid_num
    kgrid = get_grid_num(phonon.get_supercell().cell, precision=55.)
    from ase.io import read
    for i in range(1,len(disp)):
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
        # if i != 0:
        forces.append(result.get_forces())
    phonon.set_forces(np.array(forces))

def calc_dpmd(phonon, disp, calc_dir, F_0_correction, ase_calc):
    """ Calculate Force Constant with DPMD """
    forces = []
    from ase.io import read, write
    for i in range(len(disp)):
        print(' >>> Starting {:d}-th image calculation <<< '.format(i).center(120))
        ndir = get_subdir_name(i, disp)
        bu_and_mkdir(calc_dir, ndir)
        call(['cp frozen_model.pb input.in '+calc_dir+'/poscars/POSCAR-'+str(i).zfill(3)+' '+calc_dir+'/'+ndir], shell=True)
        call(['lmp-pos2lmp.awk POSCAR-'+str(i).zfill(3)+' > structure.in'], cwd = calc_dir+'/'+ndir, shell = True)
        call(['mpiexec.hydra -np $NSLOTS lmp_mpi -in input.in > out'], cwd = calc_dir+'/'+ndir, shell = True)
        atoms = read(calc_dir+'/'+ndir+'/out.dump', index=0, format='lammps-dump', order=True)
        if i == 0:
            F_0 = atoms.get_forces(apply_constraint=False)
        else:
            if F_0_correction:
                atoms._calc.results['forces'] -= F_0
            forces.append(atoms.get_forces(apply_constraint=False).tolist())
        write(calc_dir+'/'+ndir+'/result.traj', atoms, 'traj')
    phonon.set_forces(np.array(forces))

def calc_ase_calc(phonon, disp, calc_dir, F_0_correction, ase_calc):
    """ Calculate Force Constant with any ase implemented calculator """
    # numeric_F_dx=0.001
    # parallel=True
    forces = []
    from ase.io import read, write
    for i in range(len(disp)):
        print(' >>> Starting {:d}-th image calculation <<< '.format(i).center(120))
        ndir = get_subdir_name(i, disp)
        bu_and_mkdir(calc_dir, ndir)
        call(['cp '+calc_dir+'/poscars/POSCAR-'+str(i).zfill(3)+' '+calc_dir+'/'+ndir], shell=True)

        ########### calculate forces & atomic energies with amp ############
        atoms = read(calc_dir+'/'+ndir+'/POSCAR-'+str(i).zfill(3), format = 'vasp')
        atoms.set_pbc(True)
        atoms.set_calculator(ase_calc)

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

def calc_amp_tf(phonon, disp, calc_dir, F_0_correction, ase_calc):
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
        atoms.set_calculator(ase_calc)

        #********** numerical force must precede ***********
        if i!=0:
            forces.append(calc.calculate_numerical_forces(atoms, d = numeric_F_dx, parallel = parallel))
        atoms._calc.results['forces'] = forces[-1]
        write(calc_dir+'/'+ndir+'/result.traj', atoms, 'traj')
    phonon.set_forces(np.array(forces))

def calc_amp_tf_bunch(phonon, disp, calc_dir, F_0_correction, ase_calc):
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
        atoms.set_calculator(ase_calc)

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

def calc_phonon(calculator, phonon, acoustic_sum_rule=True, F_0_correction=False, verbose=False, ase_calc=None):
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
        elif calculator == 'ase_calc':
            calc = calc_ase_calc
        elif calculator == 'amp_tf':
            calc = calc_amp_tf
        elif calculator == 'amp_tf_bunch':
            calc = calc_amp_tf_bunch
        else:
            raise ValueError('Unknown calculator ({}) has been provided.'.format(calculator))
        print(">>>>>>> {} calculation will be carried out. <<<<<<<<".format(calculator).center(120))
        calc(phonon, disp, calc_dir, F_0_correction, ase_calc)
            
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

def bs_plot(self, plt, ax, scatter_max_size, proj_facecolors, proj_edgecolors, proj_alpha, reverse_seq, legend_bool, labels=None, scatter_interval=1):
    if self._projections is not None:
        #### Define key list
        key_list = list(self._projections.keys())
        if reverse_seq:
            key_list.reverse()
        #### Pick colors
        proj_facecolors = proj_facecolors[len(proj_facecolors)-len(key_list):]
        proj_edgecolors = proj_edgecolors[len(proj_edgecolors)-len(key_list):]
        if reverse_seq:
            proj_facecolors.reverse()
            proj_edgecolors.reverse()
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
                        distances[::scatter_interval],
                        freq.T[i][::scatter_interval],
                        (scatter_max_size * projections.T[i])[::scatter_interval],
                        alpha=proj_alpha,
                        facecolors=proj_facecolors[-1],
                        edgecolors=proj_edgecolors[-1],
                        label=key,
                        )
            proj_facecolors.pop()
            proj_edgecolors.pop()
            #### Gather just the one sample legend
            legend.append(legend_tmp)
        #### Legend plot
        if reverse_seq:
            legend.reverse(); key_list.reverse()
        if legend_bool:
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
    scatter_max_size = 400.,
    proj_facecolors  = ['r', 'g', 'b', 'c', 'm', 'y'],
    proj_edgecolors  = ['r', 'g', 'b', 'c', 'm', 'y'],
    proj_alpha       = 0.5,
    legend_bool      = False,
    ylim_lower       = None,
    ylim_upper       = None,
    reverse_seq      = False,
    scatter_interval = 1,
    ):
    """
    proj_eigvec = (dict) = Eigenvectors that will be used for projection. Keys of dict will be used as label for pyplot.
    proj_facecolors = (list) = Use colors sequently. Duplication is possible.
    reverse_seq = (bool) = Change sequence of plotting scatters. (Decides which one goes up. Only esthetic purpose.)
    """

    #### Variable setting
    self = phonon
    proj_facecolors.reverse()
    proj_edgecolors.reverse()

    from matplotlib import pyplot as plt
    font = {'family':'Arial'}
    plt.rc('font', **font)
    # if labels:
        # from matplotlib import rc
        # rc('font',**{'family':'serif','sans-serif':['Times']})
        # rc('text', usetex=True)

    import matplotlib.gridspec as gridspec
    font = {'family':'Arial'}
    plt.rc('font', **font)
    # plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4,1])
    ax2 = plt.subplot(gs[0, 1])
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

    ax1 = plt.subplot(gs[0, 0], sharey=ax2)

    #### projection
    self._band_structure._projections = None
    if proj_eigvec is not None:
        # dtype managing
        for key in proj_eigvec.keys():
            proj_eigvec[key] = np.array(proj_eigvec[key], dtype=np.complex128)
        set_projection(self, proj_eigvec)
    bs_plot(self._band_structure, plt, ax1, scatter_max_size, proj_facecolors, proj_edgecolors, proj_alpha, reverse_seq, legend_bool, labels, scatter_interval)
    plt.ylabel('Frequency ({})'.format(unit), fontsize=24)
    # ax1.set_title('Phonon dispersion', fontsize=35)
    ax1.set_xlabel('')

    plt.tick_params(axis="both",direction="in", labelsize=24)
    plt.grid(alpha=0.2)
    plt.subplots_adjust(left=0.20, bottom=0.08, right=0.80, top=0.95, wspace=0.10, hspace=0.20)
    plt.ylim(ylim_lower, ylim_upper)
    # plt.tight_layout()

    return plt

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

def calc_dos(
    phonon,
    mode_projection,
    pdos_precision,
    dos_tetrahedron,
    ):
    ## Get name
    delta    = np.linalg.norm(phonon.get_displacements()[0][1:4])
    job_name = 'x{}{}{}_d{:5.3f}_sym{}'.format(
        phonon.get_supercell_matrix()[0][0],
        phonon.get_supercell_matrix()[1][1],
        phonon.get_supercell_matrix()[2][2],
        delta,
        phonon._is_symmetry,
        )
    dir_name = job_name + '.pdos/'
    call('mkdir {}'.format(dir_name), shell=True)
    from kpoints_gen import get_grid_num
    k_grids = get_grid_num(phonon.get_supercell().cell, precision=pdos_precision)
    pckl_name = dir_name + 'p-{:d}_k-{:d}-{:d}-{:d}_tetra-{}.bin'.format(pdos_precision, k_grids[0],k_grids[1],k_grids[2], dos_tetrahedron)
    ## 
    import pickle as pckl
    try:
        phonon._pdos = pckl.load(open(pckl_name, 'rb'))
    except:
        print('Failed to load PDOS pickle file. Expected file name: {}'.format(pckl_name).center(120))
        print('PDOS calculation will be carried out'.center(120))
        phonon.run_mesh(
            [k_grids[0], k_grids[1], k_grids[2]],
            is_mesh_symmetry=False,
            with_eigenvectors=True,
            is_gamma_center=True,
            )
        phonon.run_projected_dos(
            use_tetrahedron_method=dos_tetrahedron,
            )
        pckl.dump(phonon._pdos, open(pckl_name, 'wb'), protocol=2)
    else:
        print('Successfully load {}'.format(pckl_name).center(120))

    if mode_projection:
        mode_dir_name = pckl_name[:-3]+'mode_proj/'
        call('mkdir {}'.format(mode_dir_name), shell=True)
        phonon._pdos_mode = dict()
        for mode in list(mode_projection.keys()):
            mode_pckl_name = mode_dir_name + '{}.bin'.format(mode)
            try:
                phonon._pdos_mode[mode] = pckl.load(open(mode_pckl_name, 'rb'))
            except:
                print('Failed to load mode projected PDOS pickle file. Expected file name: {}'.format(mode_pckl_name).center(120))
                print('Mode projected PDOS calculation will be carried out ({} mode)'.format(mode).center(120))
                from copy import deepcopy
                phonon._pdos_mode[mode] = deepcopy(phonon._pdos)
                sim = []
                for i in range(len(phonon._pdos_mode[mode]._eigenvectors)):
                    sim.append(np.square(np.abs(np.dot(np.transpose(phonon._pdos_mode[mode]._eigenvectors[i]), np.reshape(mode_projection[mode], (-1))))))
                phonon._pdos_mode[mode]._weights = np.array(sim)
                phonon._pdos_mode[mode]._run_tetrahedron_method()
                pckl.dump(phonon._pdos_mode[mode], open(mode_pckl_name, 'wb'), protocol=2)
            else:
                print('Successfully load {}'.format(mode_pckl_name).center(120))

def plot_pdos(
    phonon,
    chemical_pdos,
    pdos_colors,
    total_dos_bool,
    mode_projection,
    proj_multiple_coef,
    proj_colors,
    unit,
    freqlim_low,
    freqlim_up,
    doslim_low,
    doslim_up,
    flip_pdos_xy,
    legend_bool,
    save_svg,
    ):
    if chemical_pdos:
        atoms = phonon._unitcell
        from ss_util import get_chem_ind_arr
        unique_chem, pdos_indices = get_chem_ind_arr(atoms.get_chemical_symbols())
    else:
        (unique_chem, pdos_indices) = (None, None)
    plt, ax = phonon.plot_projected_dos(
        pdos_indices,
        pdos_colors,
        unique_chem,
        total_dos_bool,
        flip_xy=flip_pdos_xy,
        )
    if mode_projection:
        if proj_colors:
            proj_colors.reverse()
        else:
            proj_colors = [None]*len(mode_projections.keys())

        from phonopy.phonon.dos import get_pdos
        mode_pdos_list = []
        for mode in list(mode_projection.keys()):
            mode_pdos_list.append(get_pdos(
                ax,
                phonon._pdos_mode[mode]._partial_dos,
                pdos_indices,
                ))
        mode_pdos = proj_multiple_coef *np.sum(mode_pdos_list, axis=1)
        for pdos, mode in zip(mode_pdos, list(mode_projection.keys())):
            if flip_pdos_xy:
                ax.plot(-pdos, phonon._pdos._frequency_points, linewidth=1, label=mode, c=proj_colors.pop())
            else:
                ax.plot(phonon._pdos._frequency_points, -pdos, linewidth=1, label=mode, c=proj_colors.pop())
        mode_pdos_sum = np.sum(mode_pdos, axis=0)
        if flip_pdos_xy:
            ax.fill_between(-mode_pdos_sum, phonon._pdos._frequency_points, color='k', alpha=0.15)
        else:
            ax.fill_betweenx(-mode_pdos_sum, phonon._pdos._frequency_points, color='k', alpha=0.15)

    if flip_pdos_xy:
        plt.xlabel('PDOS (arb. units)',fontsize=24)
        plt.ylabel('Frequency ({})'.format(unit), fontsize=24)
        plt.subplots_adjust(left=0.30, bottom=0.15, right=0.75, top=0.97, wspace=0.2, hspace=0.2)
        plt.ylim(freqlim_low, freqlim_up)
        plt.xlim(doslim_low, doslim_up)
    else:
        plt.xlabel('Frequency ({})'.format(unit), fontsize=24)
        plt.ylabel('PDOS (arb. units)',fontsize=24)
        plt.subplots_adjust(left=0.12, bottom=0.30, right=0.99, top=0.70, wspace=0.2, hspace=0.2)
        plt.xlim(freqlim_low, freqlim_up)
        plt.ylim(doslim_low, doslim_up)
    from math import ceil, floor
    ax.set_yticks(range(int(ceil(freqlim_low)),int(floor(freqlim_up)+1),1))
    ax.tick_params(axis="both",direction="in", labelsize=24)
    plt.legend(fontsize='xx-large')
    if not legend_bool:
        plt.legend().remove()
    plt.grid(alpha=0.2)
    if save_svg:
        plt.savefig('figure.svg', format='svg')
    plt.show()
