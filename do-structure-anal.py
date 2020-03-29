#!/usr/bin/env python

import numpy as np

# Global parameters
bond_cutoff = 3.7
bonding_rules = [['Ge','Sb'], ['Te']]
step_intvl = 100000
image_ind_list = list(range(3,180000,step_intvl))
time_intvl = 10 * step_intvl / 1000 # (ps)

from structure_analyses import Structure_analyses
sa = Structure_analyses(
    'lmp-result.traj',
    image_ind_list,
    )
sa.plot_avg_coord_nums(
    sa.get_avg_coord_nums(
        bond_cutoff,
        bonding_rules,
        # load_bool=False,
        # save_bool=False,
        ),
    time_intvl,
    )
