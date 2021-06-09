#!/usr/bin/env python
import numpy as np

# Params
q_range    = range( 1,16)
NNN2       = [4, 4, 4]
NNN3       = [3, 3, 3]
prim_mat   = [[1,0,0],[0,1,0],[0,0,1]]
unitcell_f = 'gete-alpha-prim-dpmd.vasp'
calc       = 'lmp'
cp_files   = ['frozen_model.pb', 'input-phonon.in']
# calc       = 'vasp'
# cp_files   = ['INCAR', 'POTCAR',]
# run_mode   = 'onle f'
# run_mode   = 'ltc-rta'
run_mode   = 'ltc-bte'
# run_mode   = 'self-e'
# run_mode   = 'gruneisen'
temp       = (10, 50, 100, 200, 300, 500, 700,) # (K)
#
sym_fc     = True
# sym_fc     = False
nac        = True
# nac        = False
#
# fc_calc    = 'alm'
fc_calc    = None
#
# g_max = 5.5
# g_min = -12.4
g_max = None
g_min = None
# f_max = 5.5
# f_min = -0.25
f_max = None
f_min = None

if nac:
    from phonopy.interface.vasp import get_born_vasprunxml
    born_chg, eps, _ = get_born_vasprunxml(
        is_symmetry=False,
        symmetrize_tensors=True,
        )
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

from os import environ
environ['CUDA_VISIBLE_DEVICES'] = ''
from phonopy.interface import vasp
atoms = vasp.read_vasp(unitcell_f)

from phono3py import Phono3py
pho = Phono3py(
    unitcell                = atoms,
    supercell_matrix        = NNN3,
    primitive_matrix        = prim_mat,
    phonon_supercell_matrix = NNN2,
    # masses                  = None,
    # mesh                    = None,
    # band_indices            = None,
    # sigmas                  = None,
    # sigma_cutoff            = None,
    # cutoff_frequency        = 1e-4,
    # frequency_factor_to_THz = VaspToTHz,
    # is_symmetry             = True,
    # is_mesh_symmetry        = True,
    symmetrize_fc3q         = sym_fc,
    # symprec                 = 1e-5,
    # calculator              = None,
    # log_level               = 0,
    # lapack_zheev_uplo       = 'L',
    )

if nac:
    pho.set_nac_params(nac_params)

from ase.dft.kpoints import ibz_points
# points = ibz_points['hexagonal']
# G = points['Gamma']
# M = points['M']
# K = points['K']
# A = points['A']
# L = points['L']
# H = points['H']
# path = [[K, G], [G, M]]
# labels = ['K', '$\Gamma$', 'M',]

# points = ibz_points['fcc']
# G = points['Gamma']
# X = points['X']
# W = points['W']
# K = points['K']
# U = points['U']
# L = points['L']
# path = [[G, X], [X, U], [K, G], [G, L]]
# labels = ['$\Gamma$', 'X', 'U|K', '$\Gamma$', 'L']

points = {
    'Gamma': [0.,0.,0.],
    'X':[1/2., 1/2., 0.],
    'U':[0.6301369863, 0.6301369863, 0.2397260274],
    'K':[0.7602739726, 0.3698630137, 0.3698630137],
    'L':[1/2., 1/2., 1/2.],
    }
G = points['Gamma']
X = points['X']
U = points['U']
K = points['K']
L = points['L']
path = [[G, X], [X, U], [K, G], [G, L]]
labels = ['$\Gamma$', 'X', 'U|K', '$\Gamma$', 'L']

N_q = 100
import ss_phonopy as ssp
bands = ssp.make_band(path, N_q)

# # Useless part
# # print(pho.dataset.keys(), pho.phonon_dataset.keys())
# len_disp = len(pho.dataset['first_atoms']) *len(pho.dataset['first_atoms'][0]['second_atoms']) + len(pho.phonon_dataset['first_atoms'])
# # print(pho.dataset, pho.phonon_dataset)
# sc3 = pho.get_supercells_with_displacements()
# sc2 = pho.get_phonon_supercells_with_displacements()
# # print(sc2)
# # len_disp == len(sc) # sc is including harmonic supercells.
# # print(len_disp)


from ss_phono3py import calc_forces
pho = calc_forces(
    pho,
    calc,
    unitcell_f,
    cp_files,
    )
pho.produce_fc3(
    symmetrize_fc3r=sym_fc,
    fc_calculator=fc_calc,
    )
pho.produce_fc2(
    symmetrize_fc2=sym_fc,
    fc_calculator=fc_calc,
    )

for i in q_range:
    q_mesh     = [i,i,i]
    pho.mesh_numbers = q_mesh
    pho.init_phph_interaction()

    grid_points=list(range(len(pho._interaction._grid_address)))
    if run_mode == 'only f':
        continue

    elif run_mode == 'ltc-rta':
        pho.run_thermal_conductivity(
            is_LBTE=False,
            temperatures=temp,
            # is_isotope=False,
            # mass_variances=None,
            grid_points=grid_points,
            # boundary_mfp=None,  # in micrometre
            # solve_collective_phonon=False,
            # use_ave_pp=False,
            # gamma_unit_conversion=None,
            # mesh_divisors=None,
            # coarse_mesh_shifts=None,
            # is_reducible_collision_matrix=False,
            # is_kappa_star=True,
            # gv_delta_q=None,  # for group velocity
            # is_full_pp=False,
            # pinv_cutoff=1.0e-8,  # for pseudo-inversion of collision matrix
            # pinv_solver=0,  # solver of pseudo-inversion of collision matrix
            # write_gamma=False,
            # read_gamma=False,
            # is_N_U=False,
            write_kappa=True,
            # write_gamma_detail=False,
            # write_collision=False,
            # read_collision=False,
            # write_pp=False,
            # read_pp=False,
            # write_LBTE_solution=False,
            # compression="gzip",
            # input_filename=None,
            # output_filename=None,
            )
        # import h5py
        # with h5py.File('kappa-m{}{}{}.hdf5'.format(*q_mesh), 'r') as f:
            # kappa = np.array(f['kappa'])
        # for j in range(len(temp)):
            # print(' >> q-mesh: {}x{}x{} <<'.format(*q_mesh))
            # print(' ==> kappa(T={}K) = {}(W/mK)'.format(temp[j], kappa[j]))

    elif run_mode == 'ltc-bte':
        pho.run_thermal_conductivity(
            is_LBTE=True,
            temperatures=temp,
            # is_isotope=False,
            # mass_variances=None,
            grid_points=grid_points,
            # boundary_mfp=None,  # in micrometre
            # solve_collective_phonon=False,
            # use_ave_pp=False,
            # gamma_unit_conversion=None,
            # mesh_divisors=None,
            # coarse_mesh_shifts=None,
            # is_reducible_collision_matrix=False,
            # is_kappa_star=True,
            # gv_delta_q=None,  # for group velocity
            # is_full_pp=False,
            # pinv_cutoff=1.0e-8,  # for pseudo-inversion of collision matrix
            # pinv_solver=0,  # solver of pseudo-inversion of collision matrix
            # write_gamma=False,
            # read_gamma=False,
            # is_N_U=False,
            write_kappa=True,
            # write_gamma_detail=False,
            # write_collision=False,
            # read_collision=False,
            # write_pp=False,
            # read_pp=False,
            # write_LBTE_solution=False,
            # compression="gzip",
            # input_filename=None,
            # output_filename=None,
            )
        # import h5py
        # with h5py.File('kappa-m{}{}{}.hdf5'.format(*q_mesh), 'r') as f:
            # kappa = np.array(f['kappa'])
        # for j in range(len(temp)):
            # print(' >> q-mesh: {}x{}x{} <<'.format(*q_mesh))
            # print(' ==> kappa(T={}K) = {}(W/mK)'.format(temp[j], kappa[j]))

    elif run_mode == 'self-e':
        pts, delta = pho.run_real_self_energy(
            grid_points=grid_points,
            temperatures=temp,
            # run_on_bands=False,
            # frequency_points=None,
            # frequency_step=None,
            # num_frequency_points=None,
            # epsilons=None,
            # write_txt=False,
            # write_hdf5=False,
            # output_filename=None,
            )
        pts, gamma = pho.run_imag_self_energy(
            grid_points=grid_points,
            temperatures=temp,
            # frequency_points=None,
            # frequency_step=None,
            # num_frequency_points=None,
            # scattering_event_class=None,
            # write_txt=False,
            # write_gamma_detail=False,
            # keep_gamma_detail=False,
            # output_filename=None,
            )

    elif run_mode == 'gruneisen':
        if NNN2 != NNN3:
            raise ValueError('Supercell for fc2 and fc3 should be same.')
        from phono3py.phonon3.gruneisen import run_gruneisen_parameters
        from phonopy.units import VaspToTHz
        run_gruneisen_parameters(
            pho.fc2,
            pho.fc3,
            pho.supercell,
            pho.primitive,
            band_paths=np.array(bands),
            mesh=None,
            rotations=None,
            qpoints=None,
            nac_params=nac_params,
            # nac_q_direction=None,
            # ion_clamped=True,
            factor=VaspToTHz,
            # symprec=1e-5,
            # output_filename=None,
            # log_level=1,
            )
        from ss_phono3py import plot_fc3_gruneisen_yaml
        plot_fc3_gruneisen_yaml(
            labels=labels,
            g_max=g_max,
            g_min=g_min,
            f_max=f_max,
            f_min=f_min,
            )
    pckl_name = '3pho_{}_fc2x{}{}{}_fc3x{}{}{}-qx{}{}{}.pckl'.format(calc, *NNN2, *NNN3, *q_mesh)
    import pickle as pckl
    pckl.dump(pho, open(pckl_name, 'wb'))

