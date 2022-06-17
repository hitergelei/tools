#!/usr/bin/env python

import numpy as np

# Global parameters
angle_cutoff = 145.
bond_cutoff = 3.7
bonding_rules = [['Ge','Sb'], ['Te']]
alist_file = 'out.dump'
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
    # angle_cutoff,
    # bond_cutoff,
    # bonding_rules,
    # # inf_as_zero=True,
    # # deriv_sigma=15,
    # # load_bool=False,
    # # save_bool=False,
    # )

# sa.plot_chain_length_histo(
    # angle_cutoff,
    # bond_cutoff,
    # bonding_rules,
    # # load_bool=False,
    # # save_bool=False,
    # )

# sa.view_chains(
    # angle_cutoff,
    # bond_cutoff,
    # bonding_rules,
    # chain_length=8,
    # )

# sa.plot_positional_deviation(
    # in_num_avg=1,
    # out_num_avg=1,
    # out_avg_dn=1,
    # # return_intvl=10,
    # )

# sa.plot_chain_pieces_stat(
    # angle_cutoff,
    # bond_cutoff,
    # bonding_rules,
    # num_bins=100,
    # # load_bool=False,
    # # save_bool=False,
    # )

sa.plot_terminal_hist(
    angle_cutoff,
    bond_cutoff,
    2,
    bonding_rules,
    ['b','r','y'],
    # load_bool=False,
    # save_bool=False,
    )
