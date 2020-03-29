#!/usr/bin/env python

import numpy as np

def get_neighbors(
    alist,
    bond_cutoff,
    bonding_rules=None,
    dtype='float64',
    log_fname=None,
    ):
    """
    alist (List of ase atoms objects.)
        - If memory issue occurs, use this function iteratively with smaller size of atoms list.
    bond_cutoff (Float)
        - Cutoff distance for neighbor.
    bonding_rules (None or list including two groups.)
        - Only the elements in the other group will be considered as neighbor.
            e.g.) [(['Ge','Sb'],['Te'])]
        - None --> Ignore chemical kinds.
    """
    len_alist = len(alist)
    len_atoms = len(alist[0])
    types = np.array(alist[0].get_chemical_symbols())
    if bonding_rules is not None:
        group1_bool = np.sum([types == t for t in bonding_rules[0]], axis=0).astype(bool)
        group2_bool = np.sum([types == t for t in bonding_rules[1]], axis=0).astype(bool)
    box   = []
    coord = []
    for i in range(len_alist):
        box.append(alist[i].get_cell())
        coord.append(alist[i].get_scaled_positions())
    #--> shape of (len_alist, 3, 3)
    box   = np.array(box, dtype=dtype)
    #--> shape of (len_alist, len_atoms, 3)
    coord = np.array(coord, dtype=dtype)

    # Get relative coords
    #--> shape of (len_alist, len_atoms, len_atoms, 3)
    rel_coord = ((np.tile(np.expand_dims(coord, axis=1), [1, len_atoms, 1, 1]) - np.expand_dims(coord, axis=2) + \
        np.array([0.5, 0.5, 0.5], dtype=dtype)) % 1.0 - np.array([0.5, 0.5, 0.5], dtype=dtype)) @ np.reshape(box, [len_alist, 1, 3, 3])
    #--> shape of (len_alist, len_atoms, len_atoms, 5)
    concat = np.concatenate(
        (
            np.expand_dims(np.tile(np.arange(len_atoms), [len_alist, len_atoms, 1]), axis=3),
            np.linalg.norm(rel_coord, axis=3, keepdims=True),
            rel_coord,
            ),
        axis=3,
        )

    #
    indices = []
    lengths = []
    directions = []
    if log_fname:
        log = open(log_fname, 'w')
    for i in range(len_alist):
        if log_fname:
            log.write('Step {}\n'.format(i))
        indices_i = []
        lengths_i = []
        directions_i = []
        for j in range(len_atoms):
            # Cut out wrong bonds and atoms far away.
            type_mask = np.tile([True], len_atoms)
            if bonding_rules is not None:
                if group1_bool[j]:
                    type_mask = group2_bool.copy()
                else:
                    type_mask = group1_bool.copy()
            # Cut out atoms far away.
            self_false = np.tile([True], len_atoms)
            self_false[j] = False
            tmp = concat[i,j][(concat[i,j][:,1] < bond_cutoff) * self_false * type_mask]
            # Sort.
            tmp = tmp[np.argsort(tmp[:,1])]
            indices_i.append(tmp[:,0].astype(int))
            lengths_i.append(tmp[:,1])
            directions_i.append(tmp[:,2:5] / np.expand_dims(tmp[:,1], axis=1))
        indices.append(indices_i)
        lengths.append(lengths_i)
        directions.append(directions_i)
    return indices, lengths, directions

class Structure_analyses(object):
    """
    """
    def __init__(
        self,
        alist_file,
        alist_ind_list=None,
        ):
        """
        alist_file (str)
            - ASE readable alist file name including path.
            - All atoms objects must have same sequencies of chemcial symbols along the atomic indices.
        alist_ind_list (list or None)
            - List of indices of atoms object in alist to load.
        """
        #
        self.alist_file = alist_file
        if alist_ind_list is None:
            from ase.io import read
            alist_ind_list = np.range(len(read(alist_file, ':')))
        self.alist_ind_list = alist_ind_list
        #
        from ase.io import read
        atoms = read(alist_file, 0)
        self.types = np.array(atoms.get_chemical_symbols())
        self.types_unique = np.unique(self.types)
        self.types_dict = {}
        for ty in self.types_unique:
            self.types_dict[ty] = self.types == ty
        self.len_alist = len(alist_ind_list)
        self.len_atoms = len(self.types)

    def produce_neighbor_info(
        self,
        bond_cutoff,
        bonding_rules=None,
        returns=['indices', 'lengths', 'directions'],
        load_bool=True,
        save_bool=True,
        ):
        """
        """
        #
        self.bond_cutoff = bond_cutoff
        if bonding_rules is not None:
            self.bonding_rules_str = ['', '']
            for i in range(2):
                for j in range(len(bonding_rules[i])):
                    self.bonding_rules_str[i] += bonding_rules[i][j]
        else:
            self.bonding_rules_str = ['all', 'all']
        #
        info = {}
        for ret in ['indices', 'lengths', 'directions']:
            info[ret] = []
        #
        for alist_ind in self.alist_ind_list:
            # File
            from ase.io import read
            atoms = read(self.alist_file, alist_ind)

            # Calc
            import pickle as pckl
            #
            path = 'neigh_saved/{}.d/{}-{}-{}/{}'.format(self.alist_file, bond_cutoff, self.bonding_rules_str[0], self.bonding_rules_str[1], alist_ind)
            try:
                assert load_bool == True
                print('Trying to load pckl files')
                indices_i   = pckl.load(open('{}/indices.pckl'   .format(path), 'rb'))
                lengths_i   = pckl.load(open('{}/lengths.pckl'   .format(path), 'rb'))
                directions_i = pckl.load(open('{}/directions.pckl'.format(path), 'rb'))
            except:
                print('Failed to load pckl files')
                indices_i, lengths_i, directions_i = get_neighbors([atoms], bond_cutoff, bonding_rules, dtype='float32')
                indices_i    = indices_i[0]
                lengths_i    = lengths_i[0]
                directions_i = directions_i[0]
                #save
                if save_bool:
                    from subprocess import call
                    call('mkdir -p {}'.format(path), shell=True)
                    pckl.dump(indices_i,    open('{}/indices.pckl'   .format(path), 'wb'))
                    pckl.dump(lengths_i,    open('{}/lengths.pckl'   .format(path), 'wb'))
                    pckl.dump(directions_i, open('{}/directions.pckl'.format(path), 'wb'))
                    print('Saved pckl files')
            else:
                print('Loaded pckl files')
            # Gather
            info['indices']   .append(indices_i)
            info['lengths']   .append(lengths_i)
            info['directions'].append(directions_i)
        return [info[ret] for ret in returns]

    def get_avg_coord_nums(
        self,
        bond_cutoff,
        bonding_rules=None,
        load_bool=True,
        save_bool=True,
        ):

        indices = self.produce_neighbor_info(
            bond_cutoff,
            bonding_rules,
            returns=['indices'],
            load_bool=load_bool,
            save_bool=save_bool,
            )[0]

        # Count num_bonds
        num_bonds = {}
        for ty in self.types_unique:
            num_bonds[ty] = []
        for i in range(self.len_alist):
            for ty in self.types_unique:
                num_bonds[ty].append(0)
            for j in range(self.len_atoms):
                num_bonds[self.types[j]][-1] += len(indices[i][j])

        return num_bonds

    def plot_avg_coord_nums(
        self,
        num_bonds,
        time_intvl=1.,
        ):
        """
        num_bonds (dict)
            - Can be obtained from function, self.get_avg_coord_nums.
        time_intvl (float)
            - This is used only in plots.
            - Time interval in picoseconds between two atoms objects included in alist_ind_list.
            - Default is 1. (x-axis will be just bin index.)
        """

        from matplotlib import pyplot as plt
        font = {'family':'Arial'}
        plt.rc('font', **font)
        colors = ['r','b','k','g','m','y']
        colors=colors[::-1]
        for ty in self.types_unique:
            plt.plot(np.arange(self.len_alist) *time_intvl, np.array(num_bonds[ty])/np.sum(self.types_dict[ty]), c=colors.pop(), label=ty)
        plt.legend()
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.title('cut={} $\AA$, bond={}-{}, t-intvl={}'.format(self.bond_cutoff, self.bonding_rules_str[0], self.bonding_rules_str[1], time_intvl), fontsize='x-large')
        plt.xlabel('Time (ps)', fontsize='x-large')
        plt.ylabel('Average coordination number', fontsize='x-large')
        plt.grid(alpha=0.5)
        plt.show()

