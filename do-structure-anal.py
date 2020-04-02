#!/usr/bin/env python

import numpy as np

# Global parameters
angle_cutoff = 145.
bond_cutoff = 3.7
bonding_rules = [['Ge','Sb'], ['Te']]
# alist_file = 'lmp-result.traj'
# alist_file = 'gst-cubic-md-last50.traj'
alist_file = 'test.traj'
alist_slice = ':'
# alist_slice = '::1000'
# alist_slice = ':2000:300'
# alist_slice = ':1000'
# alist_slice = '180000'

from structure_analyses import Structure_analyses
sa = Structure_analyses(
    alist_file,
    alist_slice,
    dt=0.01,
    )

sa.plot_avg_coord_nums(
    bond_cutoff,
    bonding_rules,
    # load_bool=False,
    # save_bool=False,
    )

# sa.plot_chain_length_stat(
    # angle_cutoff,
    # bond_cutoff,
    # bonding_rules,
    # # inf_as_zero=True,
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
    # # chain_length=10,
    # )

# sa.plot_positional_deviation(
    # in_num_avg=1,
    # out_num_avg=10,
    # out_avg_dn=1,
    # # return_intvl=10,
    # )
