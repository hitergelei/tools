#!/usr/bin/env python

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
cp_files = ['frozen_model.pb',]

## Params
from phonopy.interface import vasp
atoms = vasp.read_vasp('Si-diamond-prim.vasp')
N                  = 4
NNN                = [[N,0,0],[0,N,0],[0,0,N]]
delta              = 0.050
# primitive_matrix   = [[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]]
# primitive_matrix   = [[0.25,0.25,0],[0,0.25,0.25],[0.25,0,0.25]]
primitive_matrix   = [[1,0,0],[0,1,0],[0,0,1]]
symmetry           = True
# symmetry           = '+-'
objective          = 'phonon'
# freqlim_up         = 6.0
freqlim_up         = None
# freqlim_low        = -0.5
freqlim_low        = None
unit               = 'THz'
# unit               = 'meV'
legend_bool        = False
plot_bool          = True

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
        factor           = factor,
        is_symmetry      = is_symmetry,
        ))
    phonons[i].generate_displacements(
        distance = delta,
        is_plusminus = is_plusminus,
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
        subscript=i
        )

gru_pho = PhonopyGruneisen(
    phonons[1],
    phonons[2],
    phonons[0],
    )

######### Band structure ##########
from ase.dft.kpoints import ibz_points, bandpath
# points = ibz_points['hexagonal']
# G = points['Gamma']
# M = points['M']
# K = points['K']
# A = points['A']
# L = points['L']
# H = points['H']
points = ibz_points['fcc']
G = points['Gamma']
X = points['X']
W = points['W']
K = points['K']
U = points['U']
L = points['L']

# path = [[K, G], [G, M]]
path = [[G, X], [X, U], [K, G], [G, L]]
N_q = 100
bands = ssp.make_band(path, N_q)
gru_pho.set_band_structure(bands)
gru_pho.get_band_structure()

########## Plot ################

# Only band plot
plt = gru_pho.plot_band_structure()
# ssp.plot_band(
    # phonon,
    # labels = ['$\Gamma$', 'X', 'U|K', '$\Gamma$', 'L'],
    # ).show()
plt.show()

