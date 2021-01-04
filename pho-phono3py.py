#!/usr/bin/env python

# Params
NNN2       = [3, 3, 3]
NNN3       = [3, 3, 3]
prim_mat   = [[1,0,0],[0,1,0],[0,0,1]]
unitcell_f = 'Si-diamond-prim.vasp'
# calc       = 'lmp'
# cp_files   = ['frozen_model.pb', 'input-phonon.in']
calc       = 'vasp'
cp_files   = ['INCAR', 'POTCAR', 'WAVECAR', 'CHGCAR']
# run_mode   = 'ltc-rta'
run_mode   = 'ltc-bte'
# run_mode   = 'self-e'
temp       = 300 # (K)
save       = True
load       = True

for i in range(1,21):
    q_mesh     = [i,i,i]
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
        # symmetrize_fc3q         = False,
        # symprec                 = 1e-5,
        # calculator              = None,
        # log_level               = 0,
        # lapack_zheev_uplo       = 'L',
        )

    # # Band path
    # from ase.dft.kpoints import ibz_points, bandpath
    # # points = ibz_points['hexagonal']
    # # G = points['Gamma']
    # # M = points['M']
    # # K = points['K']
    # # A = points['A']
    # # L = points['L']
    # # H = points['H']
    # points = ibz_points['fcc']
    # G = points['Gamma']
    # X = points['X']
    # W = points['W']
    # K = points['K']
    # U = points['U']
    # L = points['L']

    # # path = [[K, G], [G, M]]
    # path = [[G, X], [X, U], [K, G], [G, L]]
    # N_q = 100

    # from ss_phonopy import make_band
    # band_path = make_band(path, N_q)

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
        save,
        load,
        )
    pho.produce_fc3()
    pho.produce_fc2()
    pho.mesh_numbers = q_mesh
    pho.init_phph_interaction()

    if run_mode == 'ltc-rta':
        pho.run_thermal_conductivity(
            is_LBTE=False,
            temperatures=(temp,),
            # is_isotope=False,
            # mass_variances=None,
            # grid_points=None,
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
        import h5py
        with h5py.File('kappa-m{}{}{}.hdf5'.format(*q_mesh), 'r') as f:
            print(f['kappa'][0])

    elif run_mode == 'ltc-bte':
        pho.run_thermal_conductivity(
            is_LBTE=True,
            temperatures=(temp,),
            # is_isotope=False,
            # mass_variances=None,
            # grid_points=None,
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
        import h5py
        with h5py.File('kappa-m{}{}{}.hdf5'.format(*q_mesh), 'r') as f:
            print(f['kappa'][0])

    elif run_mode == 'self-e':
        grid_points = list(range(10))
        pts, delta = pho.run_real_self_energy(
            grid_points=grid_points,
            temperatures=(temp,),
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
            temperatures=(temp,),
            # frequency_points=None,
            # frequency_step=None,
            # num_frequency_points=None,
            # scattering_event_class=None,
            # write_txt=False,
            # write_gamma_detail=False,
            # keep_gamma_detail=False,
            # output_filename=None,
            )
