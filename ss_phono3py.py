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
    call('mpiexec.hydra -machinefile $TMPDIR/machines -np $NSLOTS vasp_std > out', shell=True, cwd=wdir)
    return _read_atoms(wdir, 'vasp')

def calc_lmp(
    wdir,
    ):
    call('lmp-pos2lmp.awk POSCAR > structure.in', shell=True, cwd=wdir)
    call('lmp_mpi -in input-phonon.in > out', shell=True, cwd=wdir)
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

def calc_forces(
    phono3py,
    calc,
    unitcell_f='Unknown',
    cp_files=None,
    save=True,
    load=True,
    ):

    # Check if structure is lower triangular cell
    for c in ((0,1), (0,2), (1,2)):
        if phono3py.primitive.get_cell()[c[0],c[1]] != 0. and calc == 'lmp':
            raise ValueError('Please provide lower triangular cell.')

    #
    job_name = '3pho_{}_{}_sc2-{}-{}-{}_sc3-{}-{}-{}'.format(
        calc,
        unitcell_f,
        *np.diag(phono3py.phonon_supercell_matrix),
        *np.diag(phono3py.supercell_matrix),
        )
    fc2_name = '{}_{}_sc2-{}-{}-{}'.format(
        calc,
        unitcell_f,
        *np.diag(phono3py.phonon_supercell_matrix),
        )
    fc3_name = '{}_{}_sc3-{}-{}-{}'.format(
        calc,
        unitcell_f,
        *np.diag(phono3py.supercell_matrix),
        )

    if load:
        try:
            phono3py = pckl.load(open('saved-pckls/{}-forces.pckl'.format(job_name), 'rb'))
        except:
            load=False
            load_fc2_calc=True
            load_fc3_calc=True
            print('\n\n\n*** NOTE) Failed to load saved-pckls/{}-forces.pckl file. ***\n\n\n'.format(job_name))
        else:
            print('\n\n\n=== NOTE) Loaded: saved-pckls/{}-forces.pckl ===\n\n\n'.format(job_name))
            return phono3py

    if not load:
        calc_dir = 'calcs'
        fc2_path = '{}/{}'.format(calc_dir, fc2_name)
        fc2_bu_path = '{}/bu-{}'.format(calc_dir, fc2_name)
        fc3_path = '{}/{}'.format(calc_dir, fc3_name)
        fc3_bu_path = '{}/bu-{}'.format(calc_dir, fc3_name)

        #
        phono3py.generate_displacements()
        sc2 = phono3py.get_phonon_supercells_with_displacements()
        sc3 = phono3py.get_supercells_with_displacements()

        #
        if load_fc2_calc:
            try:
                sc2_forces = []
                for i in range(len(sc2)):
                    wdir = '{}/disp-{:05d}'.format(fc2_path, i+1)
                    sc2_forces.append(_read_atoms(wdir, calc).get_forces())
            except:
                print('*** NOTE) Failed to load previous fc2 calc results. ***\n')
                load_fc2_calc=False
                #
                call('rm -rf {}'.format(fc2_bu_path), shell=True)
                call('mv {} {}'.format(fc2_path, fc2_bu_path), shell=True)
                call('mkdir -p {}'.format(fc2_path), shell=True)
            else:
                print('*** NOTE) Loaded the previous fc2 calc results. Be aware! ***\n')
        if not load_fc2_calc:
            print('=== NOTE) Start fc2 calculations! ===\n')
            from phonopy.interface.vasp import write_vasp
            sc2_forces = []
            for i in range(len(sc2)):
                wdir = '{}/disp-{:05d}'.format(fc2_path, i+1)
                call('mkdir -p {}'.format(wdir), shell=True)
                write_vasp('{}/POSCAR'.format(wdir), sc2[i])
                sc2_forces.append(_calc_forces(wdir, calc, cp_files))
                print(' == Progress: {}/{}'.format(i+1, len(sc2)))

        if load_fc3_calc:
            try:
                sc3_forces = []
                for i in range(len(sc3)):
                    wdir = '{}/disp-{:05d}'.format(fc3_path, i+1)
                    sc3_forces.append(_read_atoms(wdir, calc).get_forces())
            except:
                print('*** NOTE) Failed to load previous fc3 calc results. ***\n')
                load_fc3_calc=False
                #
                call('rm -rf {}'.format(fc3_bu_path), shell=True)
                call('mv {} {}'.format(fc3_path, fc3_bu_path), shell=True)
                call('mkdir -p {}'.format(fc3_path), shell=True)
            else:
                print('*** NOTE) Loaded the previous fc3 calc results. Be aware! ***\n')
        if not load_fc3_calc:
            print('=== NOTE) Start fc3 calculations! ===\n')
            from phonopy.interface.vasp import write_vasp
            sc3_forces = []
            for i in range(len(sc3)):
                wdir = '{}/disp-{:05d}'.format(fc3_path, i+1)
                call('mkdir -p {}'.format(wdir), shell=True)
                write_vasp('{}/POSCAR'.format(wdir), sc3[i])
                sc3_forces.append(_calc_forces(wdir, calc, cp_files))
                print(' == Progress: {}/{}'.format(i+1, len(sc3)))

        #
        phono3py.phonon_forces = np.array(sc2_forces)
        phono3py.forces = np.array(sc3_forces)

        #
        if save:
            call('mkdir saved-pckls', shell=True)
            pckl.dump(
                phono3py,
                open('saved-pckls/{}-forces.pckl'.format(job_name), 'wb'),
                protocol=3,
                )
            print('=== {}-forces.pckl file has been saved. ==='.format(job_name))

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

