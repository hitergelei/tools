#!/usr/bin/env python

from subprocess import call
from ase.io import read
import numpy as np
import pickle as pckl

def _read_atoms(
    wdir,
    calc,
    ):
    if calc == 'vasp':
        atoms = read('{}/vasprun.xml'.format(wdir))
    elif calc == 'lmp':
        atoms = read('{}/out.dump'.format(wdir))
    else:
        raise ValueError('Unknown calculator ({}) has been provided.'.format(calculator))
    return atoms

def calc_vasp(
    wdir,
    ):
    call('mpirun -np $NSLOTS vasp_std > out', shell=True, cwd=wdir)
    return _read_atoms(wdir, 'vasp')

def calc_lmp(
    wdir,
    ):
    call('lmp-pos2lmp.awk POSCAR > structure.in', shell=True, cwd=wdir)
    call('lmp_mpi -in input-phonon.in -sf intel > out', shell=True, cwd=wdir)
    return _read_atoms(wdir, 'lmp')

def _calc_forces(
    wdir,
    calc,
    cp_files=None,
    ):
    if calc == 'vasp':
        calculator = calc_vasp
    elif calc == 'lmp':
        calculator = calc_lmp
    else:
        raise ValueError('Unknown calculator ({}) has been provided.'.format(calculator))

    #
    if cp_files:
        cp_str = ''
        for f in cp_files:
            cp_str += f+' '
        call('cp {} {}'.format(cp_str, wdir), shell=True)

    #
    atoms = calculator(wdir)
    return atoms.get_forces()

def calc_fcs(
    phono3py,
    calc,
    unitcell_f='Unknown',
    cp_files=None,
    sym_fc=True,
    fc_calc=None,
    r_cut=None,
    ):
    """
    fc_calc: None or 'alm'
    r_cut: None or float. Use hiPhive if provided.
    """

    # Check if structure is lower triangular cell
    for c in ((0,1), (0,2), (1,2)):
        if phono3py.primitive.get_cell()[c[0],c[1]] != 0. and calc == 'lmp':
            raise ValueError('Please provide lower triangular cell.')

    #
    phono3py.generate_displacements()
    fc2_snd = phono3py.get_phonon_supercells_with_displacements()
    fc3_snd = phono3py.get_supercells_with_displacements()

    fc2_job_name = '{}-x{}{}{}_d0.010_sym{}-fc2'.format(
        calc,
        *np.diag(phono3py.get_phonon_supercell_matrix()),
        phono3py._is_symmetry,
        )
    fc3_job_name = '{}-x{}{}{}_d0.030_sym{}-fc3'.format(
        calc,
        *np.diag(phono3py.get_supercell_matrix()),
        phono3py._is_symmetry,
        )
    fc2_npy_name = '{}-forces.npy'.format(fc2_job_name)
    fc3_npy_name = '{}-forces.npy'.format(fc3_job_name)

    calc_dir = 'calcs'
    fc2_path = '{}/{}'.format(calc_dir, fc2_job_name)
    fc3_path = '{}/{}'.format(calc_dir, fc3_job_name)

    #
    from phonopy.interface.vasp import write_vasp
    call('mkdir -p {}/poscars'.format(fc2_path), shell=True)
    write_vasp('{}/poscars/SPOSCAR'.format(fc2_path), phono3py.get_phonon_supercell())
    for i in range(len(fc2_snd)):
        write_vasp('{}/poscars/POSCAR-{:03d}'.format(fc2_path, i+1), fc2_snd[i])
    #
    call('mkdir -p {}/poscars'.format(fc3_path), shell=True)
    write_vasp('{}/poscars/SPOSCAR'.format(fc3_path), phono3py.get_supercell())
    for i in range(len(fc3_snd)):
        write_vasp('{}/poscars/POSCAR-{:05d}'.format(fc3_path, i+1), fc3_snd[i])

    try:
        fc2_forces = np.load(fc2_npy_name)
    except:
        calc_fc2 = True
        print('\n*** NOTE) Failed to load {} file. ***\n'.format(fc2_npy_name))
    else:
        calc_fc2 = False
        print('\n=== NOTE) Loaded: {} ===\n'.format(fc2_npy_name))

    if calc_fc2:
        print('=== NOTE) Starting fc2 calculations! ===\n')
        fc2_forces = []
        for i in range(len(fc2_snd)):
            folder = 'disp-{:03d}'.format(i+1)
            wdir = '{}/{}'.format(fc2_path, folder)
            budir = '{}/bu-{}'.format(fc2_path, folder)
            call('rm -rf {}'.format(budir), shell=True)
            call('mv {} {}'.format(wdir, budir), shell=True)
            call('mkdir -p {}'.format(wdir), shell=True)
            write_vasp('{}/POSCAR'.format(wdir), fc2_snd[i])
            fc2_forces.append(_calc_forces(wdir, calc, cp_files))
            print(' == Progress: {}/{}'.format(i+1, len(fc2_snd)))
        np.save(fc2_npy_name, fc2_forces)

    try:
        fc3_forces = np.load(fc3_npy_name)
    except:
        calc_fc3 = True
        print('\n*** NOTE) Failed to load {} file. ***\n'.format(fc3_npy_name))
    else:
        calc_fc3 = False
        print('\n=== NOTE) Loaded: {} ===\n'.format(fc3_npy_name))
            
    if calc_fc3:
        print('=== NOTE) Starting fc3 calculations! ===\n')
        fc3_forces = []
        for i in range(len(fc3_snd)):
            folder = 'disp-{:05d}'.format(i+1)
            wdir = '{}/{}'.format(fc3_path, folder)
            budir = '{}/bu-{}'.format(fc3_path, folder)
            call('rm -rf {}'.format(budir), shell=True)
            call('mv {} {}'.format(wdir, budir), shell=True)
            call('mkdir -p {}'.format(wdir), shell=True)
            write_vasp('{}/POSCAR'.format(wdir), fc3_snd[i])
            fc3_forces.append(_calc_forces(wdir, calc, cp_files))
            print(' == Progress: {}/{}'.format(i+1, len(fc3_snd)))
        np.save(fc3_npy_name, fc3_forces)

    #
    phono3py.phonon_forces = np.array(fc2_forces)
    phono3py.forces = np.array(fc3_forces)

    #
    phono3py.produce_fc3(
        symmetrize_fc3r=sym_fc,
        fc_calculator=fc_calc,
        )
    phono3py.produce_fc2(
        symmetrize_fc2=sym_fc,
        fc_calculator=fc_calc,
        )
    if r_cut:
        #
        from ase.io import read
        from ase.calculators.singlepoint import SinglePointDFTCalculator
        fc2_disp_ds = phono3py.get_phonon_displacement_dataset()['first_atoms']
        fc2_sposcar = read('{}/poscars/SPOSCAR'.format(fc2_path))
        fc2_supers = []
        for i in range(len(fc2_disp_ds)):
            displacements = np.zeros(fc2_sposcar.get_positions().shape, dtype=float)
            displacements[fc2_disp_ds[i]['number']] += fc2_disp_ds[i]['displacement']
            fc2_supers.append(fc2_sposcar.copy())
            fc2_supers[-1].new_array(
                'displacements',
                displacements,
                )
            fc2_supers[-1]._calc = SinglePointDFTCalculator(fc2_supers[-1])
            fc2_supers[-1]._calc.results['forces'] = fc2_disp_ds[i]['forces']
        # rotational sum rule
        from hiphive import ForceConstants, ClusterSpace, StructureContainer
        cs = ClusterSpace(fc2_sposcar, [r_cut])
        sc = StructureContainer(cs)
        for s in fc2_supers:
            sc.add_structure(s)
        from hiphive.fitting import Optimizer
        opt = Optimizer(sc.get_fit_data(), train_size=1.0)
        opt.train()
        print(opt)
        parameters = opt.parameters
        from hiphive import ForceConstantPotential, enforce_rotational_sum_rules
        parameters_rot = enforce_rotational_sum_rules(cs, parameters, ['Huang', 'Born-Huang'])
        fcp_rot = ForceConstantPotential(cs, parameters_rot)
        fcs = fcp_rot.get_force_constants(fc2_sposcar).get_fc_array(order=2)
        phono3py.set_fc2(fcs)

        # # Not implemented yet in hiPhive
        # fc3_disp_ds = phono3py.get_displacement_dataset()['first_atoms']
        # fc3_sposcar = read('{}/poscars/SPOSCAR'.format(fc3_path))
        # fc3_supers = []
        # for i in range(len(fc3_disp_ds)):
            # displacements_1st = np.zeros(fc3_sposcar.get_positions().shape, dtype=float)
            # displacements_1st[fc3_disp_ds[i]['number']] += fc3_disp_ds[i]['displacement']
            # for j in range(len(fc3_disp_ds[i]['second_atoms'])):
                # displacements_2nd = np.zeros(fc3_sposcar.get_positions().shape, dtype=float)
                # displacements_2nd[fc3_disp_ds[i]['second_atoms'][j]['number']] += fc3_disp_ds[i]['second_atoms'][j]['displacement']
                # displacements = displacements_1st + displacements_2nd
                # fc3_supers.append(fc3_sposcar.copy())
                # fc3_supers[-1].new_array(
                    # 'displacements',
                    # displacements,
                    # )
                # fc3_supers[-1]._calc = SinglePointDFTCalculator(fc3_supers[-1])
                # fc3_supers[-1]._calc.results['forces'] = fc3_disp_ds[i]['second_atoms'][j]['forces']

    return phono3py

def plot_fc3_gruneisen_band(data, labels, g_max, g_min, f_max, f_min):
    # Imported from Togo's Phono3py code.
    import matplotlib.pyplot as plt

    d = []
    g = []
    f = []
    distance = 0.0

    ticks = [0.0]
    for path in data['path']:
        for q in path['phonon']:
            d.append(q['distance'] + distance)
            g.append([band['gruneisen'] for band in q['band']])
            f.append([band['frequency'] for band in q['band']])
        distance += path['phonon'][-1]['distance']
        ticks.append(distance)

    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(d, g, '-')
    ax2.plot(d, f, '-')
    ax1.set_xticks(ticks)
    ax2.set_xticks(ticks)
    if labels:
        ax1.set_xticklabels(labels)
        ax2.set_xticklabels(labels)

    if g_max is not None:
        ax1.set_ylim(ymax=g_max)
    if g_min is not None:
        ax1.set_ylim(ymin=g_min)
    if f_max is not None:
        ax2.set_ylim(ymax=f_max)
    if f_min is not None:
        ax2.set_ylim(ymin=f_min)

    ax1.set_xlim(ticks[0], ticks[-1])
    ax2.set_xlim(ticks[0], ticks[-1])
    ax1.tick_params(axis="both",direction="in", labelsize='x-large')
    ax2.tick_params(axis="both",direction="in", labelsize='x-large')
    ax1.grid(alpha=0.5)
    ax2.grid(alpha=0.5)

    return plt

def plot_fc3_gruneisen_yaml(
    labels=None,
    g_max=None,
    g_min=None,
    f_max=None,
    f_min=None,
    ):
    import yaml
    try:
        from yaml import CLoader as Loader
    except ImportError:
        from yaml import Loader
    with open('gruneisen.yaml') as f:
        data = yaml.load(f.read(), Loader=Loader)
    plt = plot_fc3_gruneisen_band(
        data,
        labels,
        g_max,
        g_min,
        f_max,
        f_min,
        )
    plt.show()
    return plt

