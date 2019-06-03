#!/usr/bin/env python

import numpy as np

    ## Global params
# calc = 'dpmd'
calc = 'vasp'
from phonopy.interface import vasp
atoms = vasp.read_vasp('POSCAR_tmp')
N                  = 1
NNN                = [[N,0,0],[0,N,0],[0,0,1]]
delta              = 0.050
# primitive_matrix   = [[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]]
primitive_matrix   = [[1,0,0],[0,1,0],[0,0,1]]
symmetry           = True
# symmetry           = '+-'
# phonon_or_pdos     = 'phonon'
phonon_or_pdos     = 'pdos'
freqlim_up         = None
freqlim_low        = None
unit               = 'THz'
# unit               = 'meV'
legend_bool        = True
    ## PDOS arguments
pdos_precision     = 100
chemical_pdos      = True
flip_pdos_xy       = True
dos_tetrahedron    = True

#
if symmetry is True:
    is_symmetry = True
    is_plusminus = 'auto'
elif symmetry == '+-':
    is_symmetry = True
    is_plusminus = True
elif symmetry is False:
    is_symmetry = False
    is_plusminus = 'auto'
else:
    raise ValueError('symmetry parameter "{}" is unknown.'.format(symmetry))

#
from phonopy import Phonopy, units
if unit == 'THz':
    factor = units.VaspToTHz,
elif unit == 'meV':
    factor = units.VaspToEv * 1000,
else:
    raise ValueError('Unit parameter, "{}" is unknown'.format(unit))
phonon = Phonopy(
    atoms,
    # [[N,0,0],[0,N,0],[0,0,N]],
    NNN,
    # primitive_matrix = [[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]],
    primitive_matrix = primitive_matrix,
    factor           = factor,
    is_symmetry      = is_symmetry,
    )
phonon.generate_displacements(
    distance = delta,
    is_plusminus = is_plusminus,
    )

pho_disp = phonon.get_supercells_with_displacements()

vasp.write_supercells_with_displacements(
    phonon.get_supercell(),
    pho_disp,
    )


######### get new phonon object ##########
# amp_calc = Amp.load('es_class-checkpoint.amp', label='es_class')
import ss_phonopy as ssp
phonon = ssp.calc_phonon(
    calc,
    phonon,
    # acoustic_sum_rule=False,
    # F_0_correction=True,
    # verbose=True,
    # amp_calc=amp_calc,
    )

######### Band structure ##########
from ase.dft.kpoints import ibz_points, bandpath
points = ibz_points['hexagonal']
GM = points['Gamma']
M = points['M']
K = points['K']
A = points['A']
L = points['L']
H = points['H']

#point_names = ['$\Gamma$', 'X', 'K', '$\Gamma$', 'L']
path = [[K, GM], [GM, M]]
# path = [[GM, X], [X, K], [K, GM], [GM, L]]
N_q = 100

bands = ssp.make_band(path, N_q)

phonon.set_band_structure(
    bands,
    is_eigenvectors=True,
    )


######### eigenvectors #########
# freq, eigvec = phonon.get_frequencies_with_eigenvectors([0.00000333,0.00000333,0.])
# freq, eigvec = phonon.get_frequencies_with_eigenvectors([0.,0.000001,0.])
freq, eigvec = phonon.get_frequencies_with_eigenvectors(M)
eigvec = eigvec.T
np.savez('freqNeigvec', freq=freq, eigvec=eigvec)

########## Plot ################
#ssp.plot_band(phonon, labels = ['GM', 'X', 'U|K', 'GM', 'L']).show()
yaml_name = 'gnuplot_phonon.yaml'
phonon.write_yaml_band_structure(filename=yaml_name)
from subprocess import call
call(['phonopy-bandplot --gnuplot '+yaml_name+' > band-'+calc+'.in'], shell=True)
# eu1 = np.load('eu1.npy')
# eu2 = np.load('eu2.npy')

#### Band plot
if phonon_or_pdos == 'phonon':
    from kpoints_gen import get_grid_num
    k_grids = get_grid_num(phonon.get_supercell().cell, precision=55.)
    phonon.run_mesh(
        [k_grids[0], k_grids[1], k_grids[2]],
        # is_mesh_symmetry=False,
        with_eigenvectors=True,
        )
    phonon.run_total_dos()

    # eu1 = eigvec[3]
    # eu2 = eigvec[4]
    ssp.plot_band_and_dos(
        phonon,
        labels           = ['K', '$\Gamma$', 'M'],
        unit             = unit,
        # proj_eigvec      = {'Eu-A':eu1, 'Eu-B':eu2},
        # proj_size_factor = 400.,
        # proj_colors      = ['g', 'r'],
        # proj_alpha       = 0.1,
        # proj_legend      = True,
        ylim_lower       = freqlim_low,
        ylim_upper       = freqlim_up,
        # reverse_seq      = True,
        ).show()

#### Partial DOS plot
if phonon_or_pdos == 'pdos':
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
            )
        phonon.run_projected_dos(
            use_tetrahedron_method=dos_tetrahedron,
            )
        pckl.dump(phonon._pdos, open(pckl_name, 'wb'), protocol=2)
    else:
        print('Successfully load {}'.format(pckl_name).center(120))

    phonon.write_projected_dos('pdos-'+calc+'.in')
    if chemical_pdos:
        chem_arr = np.array(atoms.get_chemical_symbols())
        unique_chem = np.unique(chem_arr)
        ind = np.arange(len(chem_arr))
        pdos_indices = []
        for chem in unique_chem:
            pdos_indices.append(ind[chem_arr==chem])
        if not legend_bool:
            unique_chem = None

    else:
        (unique_chem, pdos_indices) = (None, None)
    plt = phonon.plot_projected_dos(
        pdos_indices,
        unique_chem,
        flip_xy=flip_pdos_xy,
        )
    if flip_pdos_xy:
        plt.xlabel('PDOS (a.u.)',fontsize='xx-large')
        plt.ylabel('Frequency ({})'.format(unit), fontsize='xx-large')
        plt.subplots_adjust(left=0.35, bottom=0.15, right=0.65, top=0.97, wspace=0.2, hspace=0.2)
        plt.ylim(freqlim_low, freqlim_up)
    else:
        plt.xlabel('Frequency ({})'.format(unit), fontsize='xx-large')
        plt.ylabel('PDOS (a.u.)',fontsize='xx-large')
        plt.subplots_adjust(left=0.12, bottom=0.30, right=0.99, top=0.70, wspace=0.2, hspace=0.2)
        plt.xlim(freqlim_low, freqlim_up)
    plt.xticks(fontsize='xx-large')
    plt.yticks(fontsize='xx-large')
    plt.grid(alpha=0.2)
    plt.show()

