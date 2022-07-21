#!/usr/bin/env python

import numpy as np

# Global parameters
# angle_cutoff = 145.
# bond_cutoff = 3.7
angle_cutoff = 135.
bond_cutoff = 4.0
bonding_rules = [['Ge','Sb'], ['Te']]
# bonding_rules = [['Ge','Sb','X'], ['Te']]
alist_file = 'out.dump'
# alist_file = 'pre-out.dump'
# alist_file = 'vacancy/out.dump/4.0-GeSb-Te/135.0/dv3.007-ru-2.2-rv-2.0-nt-Te/vac-8000:8001:1.traj'
# alist_file = 'gst-cubic-md-last50.traj'
# alist_file = 'test.traj'
# alist_slice = ':'
# alist_slice = '0:200'
# alist_slice = '::100'
# alist_slice = '5000:15000:10'
# alist_slice = '120000:175000:100'
# alist_slice = '350'
from sys import argv
alist_slice = argv[1]

from structure_analyses import Structure_analyses
sa = Structure_analyses(
    alist_file,
    alist_slice,
    # dt=0.01,
    dt=10.,
    # dt=100,
    )

# sa.plot_avg_coord_nums(
    # bond_cutoff,
    # bonding_rules,
    # # load_bool=False,
    # # save_bool=False,
    # )

# sa.plot_chain_length_stat(
    # bond_cutoff,
    # angle_cutoff,
    # bonding_rules,
    # # inf_as_zero=True,
    # # deriv_sigma=5,
    # # load_bool=False,
    # # save_bool=False,
    # )

# sa.plot_chain_length_histo(
    # bond_cutoff,
    # angle_cutoff,
    # bonding_rules,
    # # load_bool=False,
    # # save_bool=False,
    # )

# sa.view_chains(
    # bond_cutoff,
    # angle_cutoff,
    # bonding_rules,
    # chain_length=0,
    # )

# sa.plot_positional_deviation(
    # in_num_avg=1,
    # out_num_avg=1,
    # out_avg_dn=1,
    # # return_intvl=10,
    # )

# sa.plot_chain_pieces_stat(
    # bond_cutoff,
    # angle_cutoff,
    # bonding_rules,
    # num_bins=100,
    # # load_bool=False,
    # # save_bool=False,
    # )

# sa.plot_terminal_stat(
    # bond_cutoff,
    # angle_cutoff,
    # 2,
    # bonding_rules,
    # ['b','r','y'],
    # # load_bool=False,
    # # save_bool=False,
    # )

# sa.plot_atomic_energy_histo(
    # bond_cutoff,
    # angle_cutoff,
    # bonding_rules,
    # unique_seq = {'Te': [
        # ['Ge', 'Te', 'Ge'],
        # ['Ge', 'Te', 'Sb'],
        # ['Sb', 'Te', 'Sb'],
        # ['Ge', 'Te', 'V'],
        # ['Sb', 'Te', 'V'],
        # # ['Ge', 'Te', 'X'],
        # # ['Sb', 'Te', 'X'],
        # ]},
    # # include_terminal=False,
    # # num_bins=200,
    # # load_bool=False,
    # # save_bool=False,
    # )

# sa.view_terminal_atoms(
    # bond_cutoff,
    # angle_cutoff,
    # nth_term    = 1,
    # plot_chem   = ['Te'],
    # bond_rules  = bonding_rules,
    # # load_bool = False,
    # # save_bool = False,
    # )

sa.replace_vacancies_w_X(
    bond_cutoff,
    angle_cutoff,
    bond_rules   = bonding_rules,
    vac_dist     = 3.007,
    # vac_dist     = 1.500,
    unite_cutoff = 2.2,
    vac_radius   = 2.0,
    next_to      = ['Te'],
    # wrap         = False,
    # view_X       = True,
    # load_bool    = False,
    # save_bool    = False,
    # file_name    = None,
    )
