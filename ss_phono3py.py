#!/usr/bin/env python

from subprocess import call
from ase.io import read
import numpy as np
import pickle as pckl

def calc_vasp(
    wdir,
    calc,
    ):
    call('mpiexec.hydra -machinefile $TMPDIR/machines -np $NSLOTS vasp_std > out', shell=True, cwd=wdir)
    return read('{}/vasprun.xml'.format(wdir))

def calc_lmp(
    wdir,
    calc,
    ):
    call('lmp-pos2lmp.awk POSCAR > structure.in', shell=True, cwd=wdir)
    call('lmp_mpi -in input-phonon.in > out', shell=True, cwd=wdir)
    return read('{}/out.dump'.format(wdir))

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
    atoms = calculator(wdir, calc)
    return atoms.get_forces()

def calc_forces(
    phono3py,
    calc,
    unitcell_f='Unknown',
    cp_files=None,
    save=True,
    load=True,
    ):

    #
    job_name = '3pho_{}_{}_sc2-{}-{}-{}_sc3-{}-{}-{}'.format(
        calc,
        unitcell_f,
        *np.diag(phono3py.phonon_supercell_matrix),
        *np.diag(phono3py.phonon_supercell_matrix),
        )

    if load:
        try:
            phono3py = pckl.load(open('saved-pckls/{}-forces.pckl'.format(job_name), 'rb'))
        except:
            load=False
            print('\n\n\n*** NOTE) Failed to load saved-pckls/{}-forces.pckl file. ***\n\n\n'.format(job_name))
        else:
            print('\n\n\n=== NOTE) Loaded: saved-pckls/{}-forces.pckl ===\n\n\n'.format(job_name))
            return phono3py

    if not load:
        print('=== Starting force calculations. ==='.format(job_name))
        phono3py._job_name = job_name
        calc_dir = 'calcs'
        path = '{}/{}'.format(calc_dir, job_name)
        bu_path = '{}/bu-{}'.format(calc_dir, job_name)

        #
        call('rm -rf {}'.format(bu_path), shell=True)
        call('mv {} {}'.format(path, bu_path), shell=True)
        call('mkdir -p {}'.format(path), shell=True)

        #
        phono3py.generate_displacements()
        sc2 = phono3py.get_phonon_supercells_with_displacements()
        sc3 = phono3py.get_supercells_with_displacements()

        #
        sc2_forces = []
        from phonopy.interface.vasp import write_vasp
        for i in range(len(sc2)):
            wdir = '{}/fc2/disp-{:05d}'.format(path, i+1)
            call('mkdir -p {}'.format(wdir), shell=True)
            write_vasp('{}/POSCAR'.format(wdir), sc2[i])
            sc2_forces.append(_calc_forces(wdir, calc, cp_files))
        phono3py.phonon_forces = np.array(sc2_forces)

        sc3_forces = []
        for i in range(len(sc3)):
            wdir = '{}/fc3/disp-{:05d}'.format(path, i+1)
            call('mkdir -p {}'.format(wdir), shell=True)
            write_vasp('{}/POSCAR'.format(wdir), sc3[i])
            sc3_forces.append(_calc_forces(wdir, calc, cp_files))
        phono3py.forces = np.array(sc3_forces)

        if save:
            call('mkdir saved-pckls', shell=True)
            pckl.dump(
                phono3py,
                open('saved-pckls/{}-forces.pckl'.format(job_name), 'wb'),
                protocol=3,
                )
            print('=== {}-forces.pckl file has been saved. ==='.format(job_name))

        return phono3py

