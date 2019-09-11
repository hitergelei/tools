#!/usr/bin/env python

import numpy as np
from ase import units as ase_units

    ## Global params
calc = 'dpmd'
# calc = 'vasp'
# calc = 'ase_calc'
## ASE calc
# ase_calc = Amp.load('es_class-checkpoint.amp', label='es_class')
# from ase.calculators.lj import LennardJones as LJ
# ase_calc = LJ(epsilon=120 *ase_units.kB, sigma=0.34 *ase_units.nm)
## Params
from phonopy.interface import vasp
atoms = vasp.read_vasp('POSCAR_rlx')
N                  = 4
NNN                = [[N,0,0],[0,N,0],[0,0,1]]
delta              = 0.050
# primitive_matrix   = [[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]]
# primitive_matrix   = [[0.25,0.25,0],[0,0.25,0.25],[0.25,0,0.25]]
primitive_matrix   = [[1,0,0],[0,1,0],[0,0,1]]
symmetry           = True
# symmetry           = '+-'
# phonon_or_pdos     = 'phonon'
phonon_or_pdos     = 'pdos'
freqlim_up         = 6.0
freqlim_low        = -0.5
unit               = 'THz'
# unit               = 'meV'
legend_bool        = True
plot_bool          = True
# mode_projection    = {'g1':np.load('g1.npy'), 'g2':np.load('g2.npy'), 'g3':np.load('g3.npy')}
# mode_projection    = None
mode_projection    = {'Eu-A':np.load('eu1.npy'), 'Eu-B':np.load('eu2.npy')}
save_svg           = True
    ## PDOS arguments
pdos_precision     = 250
chemical_pdos      = True
flip_pdos_xy       = True
dos_tetrahedron    = True
total_dos_bool     = True
doslim_up          = None
doslim_low         = None
    ## Phonon arguments
# reverse_seq        = True
reverse_seq        = False

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
from phonopy import Phonopy, units as phonopy_units
if unit == 'THz':
    factor = phonopy_units.VaspToTHz,
elif unit == 'meV':
    factor = phonopy_units.VaspToEv * 1000,
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
#
import ss_phonopy as ssp
phonon = ssp.calc_phonon(
    calc,
    phonon,
    # acoustic_sum_rule=False,
    # F_0_correction=True,
    # verbose=True,
    # ase_calc=ase_calc,
    )

######### Band structure ##########
from ase.dft.kpoints import ibz_points, bandpath
points = ibz_points['hexagonal']
G = points['Gamma']
M = points['M']
K = points['K']
A = points['A']
L = points['L']
H = points['H']
# points = ibz_points['fcc']
# G = points['Gamma']
# X = points['X']
# W = points['W']
# K = points['K']
# U = points['U']
# L = points['L']

#point_names = ['$\Gamma$', 'X', 'K', '$\Gamma$', 'L']
path = [[K, G], [G, M]]
# path = [[G, X], [X, U], [K, G], [G, L]]
N_q = 100

bands = ssp.make_band(path, N_q)

phonon.set_band_structure(
    bands,
    is_eigenvectors=True,
    )


######### eigenvectors #########
# freq, eigvec = phonon.get_frequencies_with_eigenvectors([0.00000333,0.00000333,0.])
# freq, eigvec = phonon.get_frequencies_with_eigenvectors([0.,0.000001,0.])
freq, eigvec = phonon.get_frequencies_with_eigenvectors(G)
eigvec = eigvec.T
np.savez('freqNeigvec', freq=freq, eigvec=eigvec)

########## Plot ################
from subprocess import call
#### Band plot
if phonon_or_pdos == 'phonon':
    from kpoints_gen import get_grid_num
    k_grids = get_grid_num(phonon.get_supercell().cell, precision=100.)
    phonon.run_mesh(
        [k_grids[0], k_grids[1], k_grids[2]],
        # is_mesh_symmetry=False,
        with_eigenvectors=True,
        is_gamma_center=True,
        )
    phonon.run_total_dos()

    # eu1 = eigvec[3]
    # eu2 = eigvec[4]
    yaml_name = 'gnuplot_phonon.yaml'
    phonon.write_yaml_band_structure(filename=yaml_name)
    call(['phonopy-bandplot --gnuplot '+yaml_name+' > band-'+calc+'.in'], shell=True)
    ssp.plot_band_and_dos(
        phonon,
        labels           = ['K', '$\Gamma$', 'M'],
        # labels           = ['$\Gamma$', 'X', 'U|K', '$\Gamma$', 'L'],
        unit             = unit,
        proj_eigvec      = mode_projection,
        proj_size_factor = 400.,
        proj_colors      = ['g', 'r', 'b'],
        proj_alpha       = 0.5,
        legend_bool      = legend_bool,
        ylim_lower       = freqlim_low,
        ylim_upper       = freqlim_up,
        reverse_seq      = reverse_seq,
        ).show()
    # Only band plot
    #ssp.plot_band(phonon, labels = ['G', 'X', 'U|K', 'G', 'L']).show()

#### Partial DOS plot
if phonon_or_pdos == 'pdos':
    ssp.calc_dos(
        phonon,
        mode_projection,
        pdos_precision,
        dos_tetrahedron,
        )

    phonon.write_projected_dos('pdos-'+calc+'.in')
    if plot_bool:
        ssp.plot_pdos(
            phonon,
            chemical_pdos,
            total_dos_bool,
            mode_projection,
            unit,
            freqlim_low,
            freqlim_up,
            doslim_low,
            doslim_up,
            flip_pdos_xy,
            legend_bool,
            save_svg,
            )


