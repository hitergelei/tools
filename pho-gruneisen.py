#!/usr/bin/env python

from os import environ
environ['CUDA_VISIBLE_DEVICES']=''
import numpy as np
from ase import units as ase_units

    ## Global params
# calc = 'lmp'
calc = 'vasp'
# calc = 'ase_calc'
## ASE calc
# ase_calc = Amp.load('es_class-checkpoint.amp', label='es_class')
# from ase.calculators.lj import LennardJones as LJ
# ase_calc = LJ(epsilon=120 *ase_units.kB, sigma=0.34 *ase_units.nm)
# cp_files = ['frozen_model.pb',]
cp_files = None

## Params
from phonopy.interface import vasp
# atoms = vasp.read_vasp('POSCAR_GeTe_conv')
atoms = vasp.read_vasp('gete-alpha-prim.vasp')
# n_snapshots        = 100 # Enable --> alamode FC fit function
n_snapshots        = None
N                  = 4
NNN                = [[N,0,0],[0,N,0],[0,0,N]]
delta              = 0.010
# primitive_matrix   = [[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]]
# primitive_matrix   = [[0.25,0.25,0],[0,0.25,0.25],[0.25,0,0.25]]
primitive_matrix   = [[1,0,0],[0,1,0],[0,0,1]]
symmetry           = True
# symmetry           = '+-'
# nac                = True
nac                = False
objective          = 'phonon'
unit               = 'THz'
# unit               = 'meV'
legend_bool        = False
plot_bool          = True
#
# g_max = 2.0
# g_min = -3.0
g_max = None
g_min = None
# f_max = 15.5
# f_min = -0.25
f_max = None
f_min = None

#
if n_snapshots:
    fc_calc = 'alm'
else:
    fc_calc = None

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
from phonopy import Phonopy, PhonopyGruneisen, units as phonopy_units
if unit == 'THz':
    factor = phonopy_units.VaspToTHz,
elif unit == 'meV':
    factor = phonopy_units.VaspToEv * 1000,
else:
    raise ValueError('Unit parameter, "{}" is unknown'.format(unit))

#
d_vol = [0.99, 1.00, 1.01]
phonons = []
for i in range(3):
    new_atoms = atoms.copy()
    new_atoms.set_cell(atoms.get_cell()*d_vol[i]**(1/3.))
    phonons.append(Phonopy(
        new_atoms,
        # [[N,0,0],[0,N,0],[0,0,N]],
        NNN,
        # primitive_matrix = [[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]],
        primitive_matrix = primitive_matrix,
        nac_params       = nac_params,
        factor           = factor,
        is_symmetry      = is_symmetry,
        ))
    phonons[i].generate_displacements(
        distance = delta,
        is_plusminus = is_plusminus,
        number_of_snapshots = n_snapshots,
        )

    pho_disp = phonons[i].get_supercells_with_displacements()

    vasp.write_supercells_with_displacements(
        phonons[i].get_supercell(),
        pho_disp,
        )


    ######### get new phonon object ##########
    #
    import ss_phonopy as ssp
    phonons[i] = ssp.calc_phonon(
        calc,
        phonons[i],
        # acoustic_sum_rule=False,
        # F_0_correction=True,
        # verbose=True,
        # ase_calc=ase_calc,
        cp_files=cp_files,
        subscript=i,
        fc_calc=fc_calc,
        )
    if nac:
        phonons[i].dynamical_matrix.show_nac_message()

gru_pho = PhonopyGruneisen(
    phonons[1],
    phonons[2],
    phonons[0],
    )

######### Band structure ##########

# points = ibz_points['hexagonal2']
# # points = ibz_points['hexagonal']
# G = points['Gamma']
# M = points['M']
# K = points['K']
# A = points['A']
# L = points['L']
# H = points['H']
# path = [[K, G], [G, M]]
# labels = ['K', '$\Gamma$', 'M',]
# path = [[G, K], [K, M], [M, G], [G, A], [A, H], [H, L], [L, A],]
# labels = ['$\Gamma$', 'K', 'M', '$\Gamma$', 'A', 'H', 'L', 'A']

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

X = [0.3700000000,    0.0         ,   -0.3700000000]
G = [0.0         ,    0.0         ,    0.0         ]
Z = [0.5         ,    0.5         ,    0.5         ]
F = [0.5         ,    0.5         ,    0.0         ]
L = [0.5         ,    0.0         ,    0.0         ]
path = [[X, G], [G, Z], [Z, F], [F, L], [L, G]]
labels = ['X', '$\Gamma$', 'T', 'F', 'L', '$\Gamma$']

N_q = 100
bands = ssp.make_band(path, N_q)
gru_pho.set_band_structure(bands)
gru_pho.get_band_structure()
gru_pho.write_yaml_band_structure()

########## Plot ################

# Only band plot
from ss_phonopy import plot_gruneisen_band
plt = plot_gruneisen_band(
    gru_pho,
    g_max = g_max,
    g_min = g_min,
    f_max = f_max,
    f_min = f_min,
    labels = labels,
    )
plt.subplots_adjust(left=0.15, right=0.90)
plt.show()

