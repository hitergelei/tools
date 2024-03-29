#!/usr/bin/env python

from os import environ
environ['CUDA_VISIBLE_DEVICES']=''
import numpy as np
from ase import units as ase_units

## Global params
calc = 'lmp'
# calc = 'vasp'
# calc = 'ase_calc'
## ASE calc
# ase_calc = Amp.load('es_class-checkpoint.amp', label='es_class')
# from ase.calculators.lj import LennardJones as LJ
# ase_calc = LJ(epsilon=120 *ase_units.kB, sigma=0.34 *ase_units.nm)
# cp_files = None
# cp_files = ['Si.tersoff',]
cp_files = ['frozen_model.pb',]
acou_sum_rule = True
# acou_sum_rule = False
# rot_sum_rule = True
rot_sum_rule = False
# r_cut = 4.5
r_cut = None
# masses = [1., 1.2]
masses = None

## Params
from os import environ
environ['CUDA_VISIBLE_DEVICES'] = ''
from phonopy.interface import vasp
atoms = vasp.read_vasp('POSCAR_4_FerroGeTe-wo-vdw')
N                  = 4
NNN                = [[N,0,0],[0,N,0],[0,0,1]]
delta              = 0.010
# primitive_matrix   = [[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]]
# primitive_matrix   = [[0.25,0.25,0],[0,0.25,0.25],[0.25,0,0.25]]
primitive_matrix   = [[1,0,0],[0,1,0],[0,0,1]]
symmetry           = True
# symmetry           = '+-'
# symmetry           = False
# nac                = True
nac                = False
run_mode           = 'phonon'
# run_mode           = 'pdos'
# run_mode           = 'jdos'
# num_freq_points    = 500
# temp               = np.concatenate([
    # np.arange( 0,11,1, dtype=float),
    # np.arange(12,60,3, dtype=float),
    # np.arange(60,400,10, dtype=float),
    # ]).tolist() # (K)
freqlim_up         = 5.5
# freqlim_up         = None
freqlim_low        = 0
# freqlim_low        = None
unit               = 'THz'
# unit               = 'meV'
legend_bool        = False
plot_bool          = True
# mode_projection    = {'g1':np.load('g1.npy'), 'g2':np.load('g2.npy'), 'g3':np.load('g3.npy')}
# mode_projection    = None
mode_projection    = {'Eu-A':np.load('eu2.npy'), 'Eu-B':np.load('eu1.npy')}
scatter_max_size   = 100
proj_alpha         = 0.2
save_svg           = True
scatter_interval   = 1
proj_facecolors    = ['g', 'r', 'b']
proj_edgecolors    = ['g', 'r', 'b']
    ## PDOS arguments
pdos_mesh          = [60,60,60]
chemical_pdos      = True
proj_multiple_coef = 8.
pdos_colors        = ['teal','firebrick','olive']
flip_pdos_xy       = True
dos_tetrahedron    = True
total_dos_bool     = True
doslim_up          = None
doslim_low         = None
    ## Phonon arguments
# reverse_seq        = True
reverse_seq        = False
# plot_adjust        = None
plot_adjust        = (
    0.25, # left
    0.40, # bottom
    0.75, # right
    0.85, # top
    0.10, # wspace
    0.20, # hspace
    )

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

if nac:
    from phonopy.interface.vasp import get_born_vasprunxml
    born_chg, eps, _ = get_born_vasprunxml(
        is_symmetry=False,
        symmetrize_tensors=True,
        )
    print(born_chg, eps)
    from phonopy.interface.calculator import get_default_physical_units
    nac_factor = get_default_physical_units('vasp')['nac_factor']
    nac_params = {
        'born': born_chg,
        'dielectric':eps,
        'factor':nac_factor,
        # 'method':'wang',
        }
else:
    nac_params = None

#
from phonopy import Phonopy, units as phonopy_units
if unit == 'THz':
    factor = phonopy_units.VaspToTHz
elif unit == 'meV':
    factor = phonopy_units.VaspToEv * 1000
else:
    raise ValueError('Unit parameter, "{}" is unknown'.format(unit))
phonon = Phonopy(
    atoms,
    # [[N,0,0],[0,N,0],[0,0,N]],
    NNN,
    # primitive_matrix = [[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]],
    primitive_matrix = primitive_matrix,
    nac_params       = nac_params,
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
    acoustic_sum_rule=acou_sum_rule,
    rot_sum_rule=rot_sum_rule,
    r_cut=r_cut,
    # F_0_correction=True,
    # verbose=True,
    # ase_calc=ase_calc,
    cp_files=cp_files,
    )
if masses is not None:
    phonon.masses = masses
if nac:
    phonon.dynamical_matrix.show_nac_message()

######### Band structure ##########
from ase.dft.kpoints import ibz_points
# points = ibz_points['orthorhombic']
# G = points['Gamma']
# R = points['R']
# S = points['S']
# T = points['T']
# U = points['U']
# X = points['X']
# Y = points['Y']
# Z = points['Z']
# path = [[G, X], [X, S], [S, Y], [Y, G]]
# labels = ['$\Gamma$', 'X', 'S', 'Y', '$\Gamma$']

# points = ibz_points['hexagonal2']
points = ibz_points['hexagonal']
G = points['Gamma']
M = points['M']
K = points['K']
A = points['A']
L = points['L']
H = points['H']
path = [[K, G], [G, M], [M, K]] #, [G, A], [A, H], [H, L], [L, A],]
labels = ['K', '$\Gamma$', 'M', 'K']#, 'A', 'H', 'L', 'A']

# points = ibz_points['fcc']
# G = points['Gamma']
# X = points['X']
# W = points['W']
# K = points['K']
# U = points['U']
# L = points['L']
# path = [[G, X], [X, U], [K, G], [G, L]]
# labels = ['$\Gamma$', 'X', 'U|K', '$\Gamma$', 'L']

# points = {
    # 'Gamma': [0.,0.,0.],
    # 'X':[1/2., 1/2., 0.],
    # 'U':[0.6301369863, 0.6301369863, 0.2397260274],
    # 'K':[0.7602739726, 0.3698630137, 0.3698630137],
    # 'L':[1/2., 1/2., 1/2.],
    # }
# G = points['Gamma']
# X = points['X']
# U = points['U']
# K = points['K']
# L = points['L']
# path = [[G, X], [X, U], [K, G], [G, L]]
# labels = ['$\Gamma$', 'X', 'U|K', '$\Gamma$', 'L']

# DOS grid
# from kpoints_gen import get_grid_num
# k_grids = get_grid_num(phonon.get_primitive().cell, precision=pdos_precision)
phonon.run_mesh(
    pdos_mesh,
    # is_mesh_symmetry=False,
    with_eigenvectors=True,
    is_gamma_center=True,
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
if run_mode == 'phonon':
    #
    N_q = 100
    bands = ssp.make_band(path, N_q)

    phonon.run_band_structure(
        bands,
        with_eigenvectors=True,
        )

    phonon.run_total_dos()

    # eu1 = eigvec[3]
    # eu2 = eigvec[4]
    yaml_name = 'gnuplot_phonon.yaml'
    phonon.write_yaml_band_structure(filename=yaml_name)
    call(['phonopy-bandplot --gnuplot '+yaml_name+' > band-'+calc+'.in'], shell=True)
    ssp.plot_band_and_dos(
        phonon,
        labels           = labels,
        unit             = unit,
        proj_eigvec      = mode_projection,
        scatter_max_size = scatter_max_size,
        proj_facecolors  = proj_facecolors,
        proj_edgecolors  = proj_edgecolors,
        proj_alpha       = proj_alpha,
        legend_bool      = legend_bool,
        ylim_lower       = freqlim_low,
        ylim_upper       = freqlim_up,
        reverse_seq      = reverse_seq,
        scatter_interval = scatter_interval,
        plot_adjust      = plot_adjust,
        ).show()
    # Only band plot
    #ssp.plot_band(phonon, labels = labels).show()

#### Partial DOS plot
if run_mode == 'pdos':
    ssp.calc_dos(
        phonon,
        mode_projection,
        250,
        dos_tetrahedron,
        )

    phonon.write_projected_dos('pdos-'+calc+'.in')
    if plot_bool:
        ssp.plot_pdos(
            phonon,
            chemical_pdos,
            pdos_colors,
            total_dos_bool,
            mode_projection,
            proj_multiple_coef,
            proj_edgecolors,
            unit,
            freqlim_low,
            freqlim_up,
            doslim_low,
            doslim_up,
            flip_pdos_xy,
            legend_bool,
            save_svg,
            )

if run_mode == 'jdos':
    np.save('ir_grid_points-{}{}{}.npy'.format(*pdos_mesh), phonon.mesh.get_ir_grid_points())
    np.save('grid_addresses-{}{}{}.npy'.format(*pdos_mesh), phonon.mesh.get_grid_address())
    np.save('ir_grid_weights-{}{}{}.npy'.format(*pdos_mesh), phonon.mesh.get_weights())
    from phono3py import Phono3pyJointDos as JDOS
    jdos = JDOS(
        phonon.get_supercell(),
        phonon.get_primitive(),
        pdos_mesh,
        phonon.get_force_constants(),
        nac_params=nac_params,
        temperatures=temp,
        num_frequency_points=num_freq_points,
        )
    jdos.run(
        grid_points=phonon.mesh.get_ir_grid_points(),
        write_jdos=True,
        )
