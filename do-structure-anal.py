#!/usr/bin/env python

import numpy as np

# Global parameters
angle_cutoff = 145.
bond_cutoff = 3.7
bonding_rules = [['Ge','Sb'], ['Te']]
step_intvl = 1
alist_ind_list = list(range(0,1000,step_intvl))
# alist_ind_list = list(range(0,2000,step_intvl))
# alist_ind_list = list(range(10))
# alist_ind_list = [180000]
time_intvl = 10 * step_intvl / 1000 # (ps)

from structure_analyses import Structure_analyses
sa = Structure_analyses(
    # 'lmp-result.traj',
    'GST-melting.traj'
    # 'gst-cubic-md-last50.traj',
    # 'test.traj',
    )

# sa.plot_avg_coord_nums(
    # bond_cutoff,
    # bonding_rules,
    # alist_ind_list,
    # time_intvl,
    # # load_bool=False,
    # # save_bool=False,
    # )

# sa.plot_chain_length_stat(
    # angle_cutoff,
    # bond_cutoff,
    # bonding_rules,
    # alist_ind_list,
    # time_intvl,
    # # inf_as_zero=True,
    # # load_bool=False,
    # # save_bool=False,
    # )

sa.plot_chain_length_histo(
    angle_cutoff,
    bond_cutoff,
    bonding_rules,
    alist_ind_list,
    # [15000]
    # load_bool=False,
    # save_bool=False,
    )

# sa.view_chains(
    # angle_cutoff,
    # bond_cutoff,
    # bonding_rules,
    # # [15000],
    # alist_ind_list,
    # chain_length=2,
    # )