#!/usr/bin/env python
import numpy as np

# 
prim_cell  = 'Si-diamond-prim.vasp'
N          = 5
NNN        = ((N,0,0),(0,N,0),(0,0,N))
T          = 300 #K
seed_range = range(0,30,1)
calc       = 'lmp'
cp_files   = ('frozen_model.pb',)

def get_eigen_set(
    force_constants,
    masses,
    ):
    """
    force_constants (array): Force constants of supercell in phonopy convention.
    masses (list): Mass list in ASE convension (atomic unit)
    """
    masses = np.array(masses)
    fc = np.reshape(
        np.transpose(
            force_constants,
            axes=(0,2,1,3),
            ),
        (3*len(masses), 3*len(masses)),
        )
    # Build dynamical matrix
    rminv = (masses ** -0.5).repeat(3)
    dynamical_matrix = fc * rminv[:, None] * rminv[None, :]
    # Solve eigenvalue problem to compute phonon spectrum and eigenvectors
    eigen_set = np.linalg.eigh(dynamical_matrix)
    return eigen_set

def get_phonon_distrib_alist(
    prim_cell,
    NNN,
    T,
    seed_range,
    calc,
    cp_files=None,
    plus_minus=True,
    delta=0.01,
    ):
    """
    prim_cell (str)
        ASE readable primitive unitcell structure file.
    NNN (list of int, shape: (3,3))
        Supercell matrix used to calculate force constants.
    T (float)
        Temperature in Kelvin unit.
    seed_range (range or list of int)
        Seeds for random number generation.
    calc (str, choices: ('lmp', 'vasp'))
        Calculator for force calculation.
    cp_files (list of str)
        List of input files for force calculation.
    plus_minus (bool)
        Option handed to ASE function.
    delta (float)
        Displacement in Angstrom for finite displacement method.
    """

    # Main
    from ase.io import read, write
    atoms = read(prim_cell)
    from phonopy import Phonopy
    phonon = Phonopy(
        atoms,
        NNN,
        )
    phonon.generate_displacements(
        delta,
        )
    pho_super = phonon.get_supercell()
    pho_disp  = phonon.get_supercells_with_displacements()
    from phonopy.interface import vasp
    vasp.write_supercells_with_displacements(
        pho_super,
        pho_disp,
        )
    import ss_phonopy as ssp
    phonon = ssp.calc_phonon(
        calc,
        phonon,
        cp_files=cp_files,
        )
    masses = pho_super.masses
    job_name = 'eig_x{}{}{}_d{:5.3f}_sym{}'.format(NNN[0][0], NNN[1][1], NNN[2][2], delta, phonon._is_symmetry)
    try:
        e_s = np.load('{}.npz'.format(job_name))
    except:
        print('Failed to load {}.npz file. Start to solve eigen problem.'.format(job_name))
        eigen_set = get_eigen_set(
            phonon.get_force_constants(),
            masses,
            )
        np.savez('{}.npz'.format(job_name), w2=eigen_set[0], D=eigen_set[1])
    else:
        print('Successfully loaded {}.npz file!'.format(job_name))
        eigen_set = (e_s['w2'], e_s['D'])

    from ase.md.velocitydistribution import phonon_harmonics
    from ase import Atoms, units
    alist = []
    for i in seed_range:
        d_ac, v_ac = phonon_harmonics(
            None,
            masses,
            T *units.kB,
            eigen_set=eigen_set,
            # rng=np.random.rand,
            seed=i,
            quantum=True,
            plus_minus=plus_minus,
            # return_eigensolution=False,
            # failfast=True,
            )
        new_super = Atoms(
            cell       = pho_super.get_cell(),
            symbols    = pho_super.get_chemical_symbols(),
            positions  = pho_super.get_positions() + d_ac,
            velocities = v_ac,
            pbc        = True,
            )
        alist.append(new_super)
    return(alist)
