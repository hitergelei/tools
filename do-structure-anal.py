#!/usr/bin/env python

import numpy as np

# Global parameters
angle_cutoff = 135.
bond_cutoff = 3.7
bonding_rules = [['Ge','Sb'], ['Te']]
step_intvl = 100
alist_ind_list = list(range(0,180001,step_intvl))
# alist_ind_list = list(range(0,2000,step_intvl))
# alist_ind_list = list(range(10))
# alist_ind_list = [180000]
time_intvl = 10 * step_intvl / 1000 # (ps)

from structure_analyses import Structure_analyses
sa = Structure_analyses(
    'lmp-result.traj',
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
    # # load_bool=False,
    # # save_bool=False,
    # )

# sa.plot_chain_length_histo(
    # 15000,
    # angle_cutoff,
    # bond_cutoff,
    # bonding_rules,
    # # load_bool=False,
    # # save_bool=False,
    # )

sa.view_chains(
    15000,
    angle_cutoff,
    bond_cutoff,
    bonding_rules,
    )
