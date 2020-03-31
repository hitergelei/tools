#!/usr/bin/env python

"""
    Young-Jae Choi
    Computational Nano-Physics Laboratory
    Physics Dept. POSTECH, Pohang, South Korea.
"""

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
        ):
        """
        alist_file (str)
            - ASE readable alist file name including path.
            - All atoms objects must have same sequencies of chemcial symbols along the atomic indices.
        """
        #
        self.alist_file = alist_file
        #
        from ase.io import read
        atoms = read(alist_file, 0)
        self.types = np.array(atoms.get_chemical_symbols())
        self.types_unique = np.unique(self.types)
        self.types_dict = {}
        for ty in self.types_unique:
            self.types_dict[ty] = self.types == ty
        self.len_atoms = len(self.types)

    def produce_neighbor_info(
        self,
        bond_cutoff,
        bonding_rules=None,
        alist_ind_list=None,
        returns=['indices', 'lengths', 'directions'],
        load_bool=True,
        save_bool=True,
        ):
        """
        alist_ind_list (list or None)
            - List of indices of atoms object in alist to load.
        """
        #
        if alist_ind_list is None:
            from ase.io import read
            alist_ind_list = np.arange(len(read(self.alist_file, ':')))
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
        from ase.io import read
        import pickle as pckl
        from subprocess import call
        for alist_ind in alist_ind_list:
            # File
            atoms = read(self.alist_file, alist_ind)

            # Calc
            #
            path = 'neigh_saved/{}.d/{}-{}-{}/{}'.format(
                self.alist_file,
                bond_cutoff,
                self.bonding_rules_str[0],
                self.bonding_rules_str[1],
                alist_ind,
                )
            try:
                assert load_bool == True
                print('Trying to load pckl files.')
                indices_i    = pckl.load(open('{}/indices.pckl'   .format(path), 'rb'))
                lengths_i    = pckl.load(open('{}/lengths.pckl'   .format(path), 'rb'))
                directions_i = pckl.load(open('{}/directions.pckl'.format(path), 'rb'))
            except:
                print('Failed to load pckl files from {}.'.format(path))
                indices_i, lengths_i, directions_i = get_neighbors([atoms], bond_cutoff, bonding_rules, dtype='float32')
                indices_i    = indices_i[0]
                lengths_i    = lengths_i[0]
                directions_i = directions_i[0]
                #save
                if save_bool:
                    call('mkdir -p {}'.format(path), shell=True)
                    pckl.dump(indices_i,    open('{}/indices.pckl'   .format(path), 'wb'))
                    pckl.dump(lengths_i,    open('{}/lengths.pckl'   .format(path), 'wb'))
                    pckl.dump(directions_i, open('{}/directions.pckl'.format(path), 'wb'))
                    print('Saved pckl files at {}.'.format(path))
            else:
                print('Loaded pckl files.')
            # Gather
            info['indices']   .append(indices_i)
            info['lengths']   .append(lengths_i)
            info['directions'].append(directions_i)
        return [info[ret] for ret in returns]

    def get_chain_set(
        self,
        angle_cutoff,
        bond_cutoff,
        bonding_rules=None,
        alist_ind_list=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        angle_cutoff (float)
            - In degree.
        alist_ind_list (list or None)
            - List of indices of atoms object in alist to load.
        """
        #
        if alist_ind_list is None:
            from ase.io import read
            alist_ind_list = np.arange(len(read(self.alist_file, ':')))
        #
        indices, directions = self.produce_neighbor_info(
            bond_cutoff,
            bonding_rules,
            alist_ind_list,
            returns=['indices', 'directions'],
            load_bool=load_bool,
            save_bool=save_bool,
            )

        import pickle as pckl
        from subprocess import call
        chain_set = []
        for i in range(len(alist_ind_list)):
            alist_ind = alist_ind_list[i]
            print('Getting chain set of image #{}'.format(alist_ind))
            path = 'neigh_saved/{}.d/{}-{}-{}/{}/{}'.format(
                self.alist_file,
                bond_cutoff,
                self.bonding_rules_str[0],
                self.bonding_rules_str[1],
                alist_ind,
                angle_cutoff,
                )
            try:
                assert load_bool == True
                print('Trying to load chain_set.pckl file.')
                chain_set_i = pckl.load(open('{}/chain_set.pckl'.format(path), 'rb'))
            except:
                print('Failed to load {}/chain_set.pckl file.'.format(path))
            else:
                print('Loaded {}/chain_set.pckl file.'.format(path))
                chain_set.append(chain_set_i)
                continue

            # Gather Chain pieces.
            chain_pieces = []
            for j in range(self.len_atoms):
                # Get upper triangle due to the duplications
                #--> shape of [number of bonds for atom j]*2
                cosines = np.triu(directions[i][j] @ directions[i][j].T, 1)
                chain_sites = np.transpose(np.where(cosines < np.cos(angle_cutoff / 180 * np.pi)))
                for inds in chain_sites:
                    #--> shape of (3,)
                    chain_pieces.append(list(np.concatenate((indices[i][j][inds], [j]))[[0,2,1]]))
            # Classify chain pieces.
            chain_set_i = []
            while chain_pieces:
                new_chain = chain_pieces.pop()
                # Extend head or tail.
                while 1:
                    changed = False
                    head = new_chain[:2]
                    tail = new_chain[-2:]
                    for p in range(len(chain_pieces)):
                        if chain_pieces[p][:2] == tail:
                            new_chain.append(chain_pieces[p][2])
                            del(chain_pieces[p])
                            changed = True
                            break
                        elif chain_pieces[p][:2] == head[::-1]:
                            new_chain = new_chain[::-1]
                            new_chain.append(chain_pieces[p][2])
                            del(chain_pieces[p])
                            changed = True
                            break
                        elif chain_pieces[p][-2:] == head:
                            new_chain = new_chain[::-1]
                            new_chain.append(chain_pieces[p][0])
                            del(chain_pieces[p])
                            changed = True
                            break
                        elif chain_pieces[p][-2:] == tail[::-1]:
                            new_chain.append(chain_pieces[p][0])
                            del(chain_pieces[p])
                            changed = True
                            break
                    if not changed:
                        chain_set_i.append(new_chain)
                        break
            chain_set.append(chain_set_i)
            #save
            if save_bool:
                call('mkdir {}'.format(path), shell=True)
                pckl.dump(chain_set_i, open('{}/chain_set.pckl'.format(path), 'wb'))
                print('Saved {}/chain_set.pckl file.'.format(path))
        return chain_set

    def get_chain_lengths(
        self,
        angle_cutoff,
        bond_cutoff,
        bonding_rules=None,
        alist_ind_list=None,
        inf_as_zero=False,
        load_bool=True,
        save_bool=True,
        ):
        """
        alist_ind_list (list or None)
            - List of indices of atoms object in alist to load.
        inf_as_zero (bool)
            - If true, return the length of infinite chain as zero.
        return (list)
            - List of chain lengths.
            - Note) Chain length of 0 means infinite chain (due to periodicity). 
                    They have a tail same as its head.
        """
        chain_set = self.get_chain_set(
            angle_cutoff,
            bond_cutoff,
            bonding_rules,
            alist_ind_list,
            load_bool,
            save_bool,
            )

        # Get chain lenghts.
        lengths = []
        for i in range(len(chain_set)):
            lengths_i = []
            for j in range(len(chain_set[i])):
                if chain_set[i][j][:2] == chain_set[i][j][-2:] and inf_as_zero:
                    lengths_i.append(0)
                else:
                    lengths_i.append(len(chain_set[i][j])-1)
            lengths.append(lengths_i)
        return lengths

    def plot_chain_length_stat(
        self,
        angle_cutoff,
        bond_cutoff,
        bonding_rules=None,
        alist_ind_list=None,
        time_intvl=1.,
        inf_as_zero=False,
        load_bool=True,
        save_bool=True,
        ):
        """
        alist_ind_list (list or None)
            - List of indices of atoms object in alist to load.
        time_intvl (float)
            - This is used only in plots.
            - Time interval in picoseconds between two atoms objects included in alist_ind_list.
            - Default is 1. (x-axis will be just bin index.)
        inf_as_zero (bool)
            - If true, let the length of infinite chain as zero.
        """

        lengths = self.get_chain_lengths(
            angle_cutoff,
            bond_cutoff,
            bonding_rules,
            alist_ind_list,
            inf_as_zero,
            load_bool,
            save_bool,
            )

        # Get info
        num_chain = []
        sum_chain = []
        # max_chain = []
        for i in range(len(lengths)):
            num_chain.append(len(lengths[i]))
            sum_chain.append(np.sum(lengths[i]))
            # max_chain.append(np.max(lengths[i]))
        num_chain = np.array(num_chain, dtype='int')
        sum_chain = np.array(sum_chain, dtype='int')
        # max_chain = np.array(max_chain, dtype='int')

        # Plot
        time_arr = np.arange(len(lengths)) * time_intvl
        from matplotlib import pyplot as plt
        font = {'family':'Arial'}
        plt.rc('font', **font)
        fig, ax1 = plt.subplots()
        ax1.plot(
            time_arr,
            sum_chain,
            label='Sum of lengths',
            c='k',
            )
        # ax1.plot(
            # time_arr,
            # max_chain,
            # label='Max. of lengths',
            # c='k',
            # )
        ax2 = ax1.twinx()
        # ax2.plot(
            # time_arr,
            # num_chain,
            # label='Num. of chains',
            # c='b',
            # )
        ax2.plot(
            time_arr,
            sum_chain / num_chain,
            label='Mean of lengths',
            c='r',
            )
        ax1.tick_params(axis="both",direction="in", labelsize='x-large')
        ax2.tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='r')
        plt.title('cut={} $\AA$ & {} deg, bond={}-{}, t-intvl={}'.format(
            bond_cutoff,
            angle_cutoff,
            self.bonding_rules_str[0],
            self.bonding_rules_str[1],
            time_intvl,
            ), fontsize='x-large')
        plt.xlabel('Time (ps)', fontsize='x-large')
        ax1.set_ylabel('Sum of chain lengths', fontsize='x-large')
        ax2.set_ylabel('Mean of chain lengths', fontsize='x-large')
        ax1.grid(alpha=0.5)
        plt.subplots_adjust(left=0.16, right=0.90)
        # plt.legend()
        plt.show()

    def get_chain_length_histo(
        self,
        angle_cutoff,
        bond_cutoff,
        bonding_rules=None,
        alist_ind_list=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        alist_ind_list (list or None)
            - List of indices of atoms object in alist to load.
            - Note) Only available for cases, len(alist_ind_list) == 1.
        """

        if alist_ind_list is None:
            from ase.io import read
            alist_ind_list = np.arange(len(read(self.alist_file, ':')))

        lengths = self.get_chain_lengths(
            angle_cutoff,
            bond_cutoff,
            bonding_rules,
            alist_ind_list,
            True,
            load_bool,
            save_bool,
            )

        # Get counts
        length_histo = []
        for i in range(len(lengths)):
            length_histo.append({})
            l_unique, l_count = np.unique(lengths[i], return_counts=True)
            for l_u, l_c in zip(l_unique, l_count):
                length_histo[-1][l_u] = l_c
        return length_histo

    def plot_chain_length_histo(
        self,
        angle_cutoff,
        bond_cutoff,
        bonding_rules=None,
        alist_ind_list=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        alist_ind (int)
            - Index of image in the alist file that you wanna see.
        alist_ind_list (list or None)
            - List of indices of atoms object in alist to load.
        """
        #
        if alist_ind_list is None:
            from ase.io import read
            alist_ind_list = np.arange(len(read(self.alist_file, ':')))

        length_histo_list = self.get_chain_length_histo(
            angle_cutoff,
            bond_cutoff,
            bonding_rules,
            alist_ind_list,
            load_bool,
            save_bool,
            )

        length_histo = {}
        for l_h in length_histo_list:
            for l, n in list(l_h.items()):
                if l in length_histo.keys():
                    length_histo[l] += n
                else:
                    length_histo[l] = n
        
        # Plot
        from matplotlib import pyplot as plt
        font = {'family':'Arial'}
        plt.rc('font', **font)
        plt.bar(length_histo.keys(), length_histo.values(), color='k', alpha=0.5)
        xticks = range(0,np.max(list(length_histo.keys())),2)
        xticklabels = list(xticks)
        xticklabels[0] = 'inf'
        plt.xticks(xticks, rotation=45, labels=xticklabels)
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.title('cut={} $\AA$ & {} deg, bond={}-{}, len_alist #{}'.format(
            bond_cutoff,
            angle_cutoff,
            self.bonding_rules_str[0],
            self.bonding_rules_str[1],
            len(alist_ind_list),
            ), fontsize='x-large')
        plt.xlabel('Length of a chain', fontsize='x-large')
        plt.ylabel('Population', fontsize='x-large')
        plt.grid(alpha=0.5)
        plt.show()

    def view_chains(
        self,
        alist_ind,
        angle_cutoff,
        bond_cutoff,
        bonding_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        alist_ind (int)
            - Index of image in the alist file that you wanna see.
        """

        chain = self.get_chain_set(
            angle_cutoff,
            bond_cutoff,
            bonding_rules,
            [alist_ind],
            load_bool,
            save_bool,
            )[0]

        from ase.io import read
        from ase.visualize import view
        alist=[]
        for c in chain:
            alist.append(read(self.alist_file, alist_ind)[np.unique(c)])
        view(alist)

    def get_avg_coord_nums(
        self,
        bond_cutoff,
        bonding_rules=None,
        alist_ind_list=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        alist_ind_list (list or None)
            - List of indices of atoms object in alist to load.
        """

        indices = self.produce_neighbor_info(
            bond_cutoff,
            bonding_rules,
            alist_ind_list,
            returns=['indices'],
            load_bool=load_bool,
            save_bool=save_bool,
            )[0]

        # Count num_bonds
        num_bonds = {}
        for ty in self.types_unique:
            num_bonds[ty] = []
        for i in range(len(indices)):
            for ty in self.types_unique:
                num_bonds[ty].append(0)
            for j in range(self.len_atoms):
                num_bonds[self.types[j]][-1] += len(indices[i][j])

        return num_bonds

    def plot_avg_coord_nums(
        self,
        bond_cutoff,
        bonding_rules=None,
        alist_ind_list=None,
        time_intvl=1.,
        load_bool=True,
        save_bool=True,
        ):
        """
        alist_ind_list (list or None)
            - List of indices of atoms object in alist to load.
        time_intvl (float)
            - This is used only in plots.
            - Time interval in picoseconds between two atoms objects included in alist_ind_list.
            - Default is 1. (x-axis will be just bin index.)
        """

        num_bonds = self.get_avg_coord_nums(
            bond_cutoff,
            bonding_rules,
            alist_ind_list,
            load_bool,
            save_bool,
            )

        # Plot
        time_arr = np.arange(len(num_bonds)) * time_intvl
        from matplotlib import pyplot as plt
        font = {'family':'Arial'}
        plt.rc('font', **font)
        colors = ['r','b','k','g','m','y']
        colors = colors[::-1]
        for ty in self.types_unique:
            plt.plot(
                time_arr,
                np.array(num_bonds[ty])/np.sum(self.types_dict[ty]),
                label=ty,
                c=colors.pop(),
                )
        plt.legend()
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.title('cut={} $\AA$, bond={}-{}, t-intvl={}'.format(bond_cutoff, self.bonding_rules_str[0], self.bonding_rules_str[1], time_intvl), fontsize='x-large')
        plt.xlabel('Time (ps)', fontsize='x-large')
        plt.ylabel('Average coordination number', fontsize='x-large')
        plt.grid(alpha=0.5)
        plt.show()

