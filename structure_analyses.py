#!/usr/bin/env python

"""
    Young-Jae Choi
    Computational Nano-Physics Laboratory
    Physics Dept. POSTECH, Pohang, South Korea.
"""

import numpy as np
import datetime
import pickle as pckl
from subprocess import call

def get_neighbors(
    alist,
    bond_cutoff,
    bond_rules=None,
    dtype='float64',
    log_fname=None,
    ):
    """
    alist (List of ase atoms objects.)
        - If memory issue occurs, use this function iteratively with smaller size of atoms list.
    bond_cutoff (Float)
        - Cutoff distance for neighbor.
    bond_rules (None or list including two groups.)
        - Only the elements in the other group will be considered as neighbor.
            e.g.) [(['Ge','Sb'],['Te'])]
        - None --> Ignore chemical kinds.
    """
    len_alist = len(alist)
    len_atoms = len(alist[0])
    types = np.array(alist[0].get_chemical_symbols())
    if bond_rules is not None:
        group1_bool = np.sum([types == t for t in bond_rules[0]], axis=0).astype(bool)
        group2_bool = np.sum([types == t for t in bond_rules[1]], axis=0).astype(bool)
    box   = []
    coord = []
    for i in range(len_alist):
        box.append(np.array(alist[i].get_cell()))
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
            if bond_rules is not None:
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

def get_3body_chain_pieces(
    indices,
    lengths,
    directions,
    angle_cutoff,
    load_path=None,
    save_path=None,
    ):
    """
    indices (list)
        - Ragged tensor in shape of (len_atoms, (number of bondings))
    lengths (list)
        - Ragged tensor in shape of (len_atoms, (number of bondings))
    directions (list)
        - Ragged tensor in shape of (len_atoms, (number of bondings), 3)
    save_path (str or None)
        - Make sure that the path already exists before this function runs.
    """

    try:
        assert load_path
        print(' *  Trying to load 3body piece pckl files.')
        piece_inds    = pckl.load(open('{}/piece_inds.pckl'   .format(load_path), 'rb'))
        piece_lengths = pckl.load(open('{}/piece_lengths.pckl'.format(load_path), 'rb'))
        piece_direcs  = pckl.load(open('{}/piece_direcs.pckl' .format(load_path), 'rb'))
    except:
        print('Failed to load pckl files from {}'.format(load_path))
    else:
        print('Loaded pckl files from {}'.format(load_path))
        return piece_inds, piece_lengths, piece_direcs

    # Main
    len_atoms = len(indices)
    #
    cos_cutoff = np.cos(angle_cutoff / 180 * np.pi)
    #
    piece_inds    = []
    piece_lengths = []
    piece_direcs  = []
    for i in range(len_atoms):
        indices_i    = np.array(indices[i])
        lengths_i    = np.array(lengths[i])
        directions_i = np.array(directions[i])
        # Gather three-body chain pieces.
        # Get upper triangle due to the duplications
        #--> shape of [number of bonds for atom j]*2
        cosines = np.triu(directions_i @ directions_i.T, 1)
        bond_inds = np.transpose(np.where(cosines < cos_cutoff))
        for inds in bond_inds:
            #                    --> shape of (3,)
            piece_inds   .append(list(np.concatenate((indices_i[inds], [i]))[[0,2,1]]))
            #                    --> shape of (2,)
            piece_lengths.append(lengths_i[inds])
            #                    --> shape of (2, 3)
            piece_direcs .append(directions_i[inds])
            piece_direcs[-1][0] *= -1.
    piece_inds    = np.reshape(piece_inds,    [-1, 3]).tolist()
    piece_lengths = np.reshape(piece_lengths, [-1, 2]).tolist()
    piece_direcs  = list(np.reshape(piece_direcs,  [-1, 2, 3]))
    if save_path is not None:
        call('mkdir -p {}'.format(save_path), shell=True)
        pckl.dump(piece_inds   , open('{}/piece_inds.pckl'      .format(save_path), 'wb'))
        pckl.dump(piece_lengths, open('{}/piece_lengths.pckl'   .format(save_path), 'wb'))
        pckl.dump(piece_direcs , open('{}/piece_direcs.pckl'.format(save_path), 'wb'))
        print('Saved 3body-piece-info pckl files at {}'.format(save_path))
    return piece_inds, piece_lengths, piece_direcs

class Structure_analyses(object):
    """
    """
    def __init__(
        self,
        alist_file,
        alist_slice=':',
        dt=0.01,
        ):
        """
        alist_file (str)
            - ASE readable alist file name including path.
            - All atoms objects must have same sequencies of chemcial symbols along the atomic indices.
        alist_slice (str)
            - Slice in python style.
        dt (float)
            - Time interval between atoms images of (alist_file).
            - Unit of picosecond.
            - Note) This parameter is independent to (alist_slice).
            - Default) 0.01 ps
        """
        #
        self.alist_file = alist_file
        #
        from ase.io import read
        self.alist = read(alist_file, alist_slice)
        if not isinstance(self.alist, list):
            self.alist = [self.alist]
        self.len_alist = len(self.alist)
        self.len_atoms = len(self.alist[0])
        # Get slice
        from ss_util import str_slice_to_list
        slice_list = str_slice_to_list(alist_slice)
        if slice_list[1] == None:
            slice_list[1] = slice_list[0] +slice_list[2] *self.len_alist
        self.alist_slice = '{}:{}:{}'.format(*slice_list)
        self.dt = dt *slice_list[2]
        # Make alist_ind_list
        self.alist_ind_list = np.arange(slice_list[1], dtype=np.int32)[slice(*slice_list)]

        # Types
        self.types = np.array(self.alist[0].get_chemical_symbols())
        self.types_unique = np.unique(self.types)
        self.types_dict = {}
        for ty in self.types_unique:
            self.types_dict[ty] = self.types == ty
        #

    def produce_neighbor_info(
        self,
        bond_cutoff,
        bond_rules=None,
        returns=['indices', 'lengths', 'directions'],
        load_bool=True,
        save_bool=True,
        ):
        """
        """
        #
        self.bond_cutoff = bond_cutoff
        if bond_rules is not None:
            self.bond_rules_str = ['', '']
            for i in range(2):
                for j in range(len(bond_rules[i])):
                    self.bond_rules_str[i] += bond_rules[i][j]
        else:
            self.bond_rules_str = ['all', 'all']
        #
        info = {}
        for ret in ['indices', 'lengths', 'directions']:
            info[ret] = []
        #
        from ase.io import read
        for i in range(self.len_alist):
            # File
            alist_ind = self.alist_ind_list[i]
            atoms = self.alist[i]

            # Calc
            #
            path = 'neigh_saved/{}.d/{}-{}-{}/{}'.format(
                self.alist_file,
                bond_cutoff,
                self.bond_rules_str[0],
                self.bond_rules_str[1],
                alist_ind,
                )
            try:
                assert load_bool == True
                print(' *  Trying to load neighbor-info pckl files.')
                indices_i    = pckl.load(open('{}/indices.pckl'   .format(path), 'rb'))
                lengths_i    = pckl.load(open('{}/lengths.pckl'   .format(path), 'rb'))
                directions_i = pckl.load(open('{}/directions.pckl'.format(path), 'rb'))
            except:
                print('Failed to load pckl files from {}'.format(path))
                indices_i, lengths_i, directions_i = get_neighbors([atoms], bond_cutoff, bond_rules, dtype='float32')
                indices_i    = indices_i[0]
                lengths_i    = lengths_i[0]
                directions_i = directions_i[0]
                #save
                if save_bool:
                    call('mkdir -p {}'.format(path), shell=True)
                    pckl.dump(indices_i,    open('{}/indices.pckl'   .format(path), 'wb'))
                    pckl.dump(lengths_i,    open('{}/lengths.pckl'   .format(path), 'wb'))
                    pckl.dump(directions_i, open('{}/directions.pckl'.format(path), 'wb'))
                    print('Saved neighbor-info pckl files at {}.'.format(path))
            else:
                print('Loaded pckl files from {}'.format(path))
            # Gather
            info['indices']   .append(indices_i)
            info['lengths']   .append(lengths_i)
            info['directions'].append(directions_i)
        return [info[ret] for ret in returns]

    def get_chain_set(
        self,
        angle_cutoff,
        bond_cutoff,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        angle_cutoff (float)
            - In degree.
        """

        #
        cos_cutoff = np.cos(angle_cutoff / 180 * np.pi)
        #
        indices, lengths, directions = self.produce_neighbor_info(
            bond_cutoff,
            bond_rules,
            returns=['indices', 'lengths', 'directions'],
            load_bool=load_bool,
            save_bool=save_bool,
            )

        ind_set      = []
        bond_direcs  = []
        bond_lengths = []
        chain_vec    = []
        for i in range(len(self.alist_ind_list)):
            alist_ind = self.alist_ind_list[i]
            path = 'neigh_saved/{}.d/{}-{}-{}/{}/{}'.format(
                self.alist_file,
                bond_cutoff,
                self.bond_rules_str[0],
                self.bond_rules_str[1],
                alist_ind,
                angle_cutoff,
                )
            try:
                assert load_bool == True
                print(' *  Trying to load chain-info pckl files.')
                ind_set_i      = pckl.load(open('{}/ind_set.pckl'     .format(path), 'rb'))
                bond_direcs_i  = pckl.load(open('{}/bond_direcs.pckl' .format(path), 'rb'))
                bond_lengths_i = pckl.load(open('{}/bond_lengths.pckl'.format(path), 'rb'))
                chain_vec_i    = pckl.load(open('{}/chain_vec.pckl'   .format(path), 'rb'))
            except:
                print('Failed to load pckl files from {}'.format(path))
            else:
                print('Loaded pckl files from {}'.format(path))
                ind_set     .append(ind_set_i)
                bond_direcs .append(bond_direcs_i)
                bond_lengths.append(bond_lengths_i)
                chain_vec   .append(chain_vec_i)
                continue

            # Gather three-body chain pieces.

            piece_inds, piece_lengths, piece_direcs = get_3body_chain_pieces(
                indices[i],
                lengths[i],
                directions[i],
                angle_cutoff,
                path,
                path,
                )
            # Classify chain pieces.
            #--> Ragged tensor in shape of ((number of chains in an image), (length of a chain))
            ind_set_i      = []
            bond_direcs_i  = []
            bond_lengths_i = []
            chain_vec_i    = []
            while piece_inds:
                ind_set_tmp      = list(piece_inds   .pop())
                bond_direcs_tmp  = list(piece_direcs .pop())
                bond_lengths_tmp = list(piece_lengths.pop())
                chain_vec_tmp    = np.sum(
                    bond_direcs_tmp * np.expand_dims(bond_lengths_tmp, axis=1),
                    axis=0,
                    )
                # Extend head or tail.
                while 1:
                    changed = False
                    chain_direc = chain_vec_tmp / np.linalg.norm(chain_vec_tmp)
                    head = ind_set_tmp[:2]
                    tail = ind_set_tmp[-2:]
                    for p in range(len(piece_inds)):
                        if piece_inds[p][:2] == tail and (piece_direcs[p][1] @ chain_direc) > -cos_cutoff:
                            #
                            ind_set_tmp.append(piece_inds[p][2])
                            del(piece_inds[p])
                            #
                            bond_direcs_tmp.append(piece_direcs[p][1])
                            del(piece_direcs[p])
                            #
                            bond_lengths_tmp.append(piece_lengths[p][1])
                            del(piece_lengths[p])
                            #
                            chain_vec_tmp += bond_direcs_tmp[-1] * bond_lengths_tmp[-1]
                            #
                            changed = True
                            break
                        elif piece_inds[p][:2] == head[::-1] and (piece_direcs[p][1] @ (-chain_direc)) > -cos_cutoff:
                            #
                            ind_set_tmp = ind_set_tmp[::-1]
                            ind_set_tmp.append(piece_inds[p][2])
                            del(piece_inds[p])
                            #
                            bond_direcs_tmp = bond_direcs_tmp[::-1]
                            bond_direcs_tmp.append(piece_direcs[p][1])
                            del(piece_direcs[p])
                            #
                            bond_lengths_tmp = bond_lengths_tmp[::-1]
                            bond_lengths_tmp.append(piece_lengths[p][1])
                            del(piece_lengths[p])
                            #
                            chain_vec_tmp *= -1.
                            chain_vec_tmp += bond_direcs_tmp[-1] * bond_lengths_tmp[-1]
                            #
                            changed = True
                            break
                        elif piece_inds[p][-2:] == head and (piece_direcs[p][0] @ chain_direc) > -cos_cutoff:
                            #
                            ind_set_tmp = ind_set_tmp[::-1]
                            ind_set_tmp.append(piece_inds[p][0])
                            del(piece_inds[p])
                            #
                            bond_direcs_tmp = bond_direcs_tmp[::-1]
                            bond_direcs_tmp.append(-piece_direcs[p][0])
                            del(piece_direcs[p])
                            #
                            bond_lengths_tmp = bond_lengths_tmp[::-1]
                            bond_lengths_tmp.append(piece_lengths[p][0])
                            del(piece_lengths[p])
                            #
                            chain_vec_tmp *= -1.
                            chain_vec_tmp += bond_direcs_tmp[-1] * bond_lengths_tmp[-1]
                            #
                            changed = True
                            break
                        elif piece_inds[p][-2:] == tail[::-1] and (-piece_direcs[p][0] @ chain_direc) > -cos_cutoff:
                            #
                            ind_set_tmp.append(piece_inds[p][0])
                            del(piece_inds[p])
                            #
                            bond_direcs_tmp.append(-piece_direcs[p][0])
                            del(piece_direcs[p])
                            #
                            bond_lengths_tmp.append(piece_lengths[p][0])
                            del(piece_lengths[p])
                            #
                            chain_vec_tmp += bond_direcs_tmp[-1] * bond_lengths_tmp[-1]
                            changed = True
                            break
                    if not changed:
                        ind_set_i     .append(ind_set_tmp)
                        bond_direcs_i .append(bond_direcs_tmp)
                        bond_lengths_i.append(bond_lengths_tmp)
                        chain_vec_i   .append(chain_vec_tmp)
                        break

            #--> Ragged tensor in shape of (len_alist, (number of chains in an image), (length of a chain))
            ind_set     .append(ind_set_i)
            #--> Ragged tensor in shape of (len_alist, (number of chains in an image), (length of a chain), 3)
            bond_direcs .append(bond_direcs_i)
            #--> Ragged tensor in shape of (len_alist, (number of chains in an image), (length of a chain))
            bond_lengths.append(bond_lengths_i)
            #--> Ragged tensor in shape of (len_alist, (number of chains in an image), 3)
            chain_vec   .append(chain_vec_i)
            # Save
            if save_bool:
                call('mkdir {}'.format(path), shell=True)
                pckl.dump(ind_set_i,      open('{}/ind_set.pckl'     .format(path), 'wb'))
                pckl.dump(bond_direcs_i,  open('{}/bond_direcs.pckl' .format(path), 'wb'))
                pckl.dump(bond_lengths_i, open('{}/bond_lengths.pckl'.format(path), 'wb'))
                pckl.dump(chain_vec_i   , open('{}/chain_vec.pckl'   .format(path), 'wb'))
                print('Saved chain-info pckl files at {}'.format(path))
        # ind_set --> Ragged tensor in shape of (len_alist, (number of chains in an image), (length of a chain))
        return ind_set, bond_direcs, bond_lengths, chain_vec

    def get_chain_lengths(
        self,
        angle_cutoff,
        bond_cutoff,
        bond_rules=None,
        inf_as_zero=False,
        load_bool=True,
        save_bool=True,
        ):
        """
        inf_as_zero (bool)
            - If true, return the length of infinite chain as zero.
        return (list)
            - List of chain lengths.
            - Note) Chain length of 0 means infinite chain (due to periodicity). 
                    They have a tail same as its head.
        """
        ind_set, bond_direcs, bond_lengths, chain_vec = self.get_chain_set(
            angle_cutoff,
            bond_cutoff,
            bond_rules,
            load_bool,
            save_bool,
            )

        # Get chain lenghts.
        lengths = []
        for i in range(len(ind_set)):
            lengths_i = []
            for j in range(len(ind_set[i])):
                if ind_set[i][j][:2] == ind_set[i][j][-2:] and inf_as_zero:
                    lengths_i.append(0)
                else:
                    lengths_i.append(len(ind_set[i][j])-1)
            lengths.append(lengths_i)
        return lengths

    def plot_chain_length_stat(
        self,
        angle_cutoff,
        bond_cutoff,
        bond_rules=None,
        inf_as_zero=False,
        therm_corr=None,
        deriv_sigma=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        inf_as_zero (bool)
            - If true, let the length of infinite chain as zero.
        therm_corr (function)
            - If a function is provided, 
        deriv_sigma (int)
            - Standard deviation (sigma) of Gaussian smearing of derivative plot.
              Unit is 'step' of final plot (dt / slice interval).
            - Set to zero for no smearing but derivative plot.
            - If 'None', don't plot the derivative plot.
        """

        lengths = self.get_chain_lengths(
            angle_cutoff,
            bond_cutoff,
            bond_rules,
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
        avg_chain = sum_chain / num_chain
        # max_chain = np.array(max_chain, dtype='int')

        # Derivatives
        if deriv_sigma is not None and not False:
            avg_chain_l = (avg_chain[1:-1] + avg_chain[:-2])/2
            avg_chain_r = (avg_chain[2:] + avg_chain[1:-1])/2
            dCdt = (avg_chain_r - avg_chain_l) / self.dt
            if deriv_sigma != 0:
                from scipy.ndimage import gaussian_filter
                dCdt = gaussian_filter(dCdt, sigma=deriv_sigma)

        # Plot
        time_arr = np.arange(len(lengths)) * self.dt
        from matplotlib import pyplot as plt
        font = {'family':'Arial'}
        plt.rc('font', **font)
        fig, ax1 = plt.subplots()
        ax1.plot(
            time_arr,
            avg_chain,
            # label='Mean of lengths',
            c='k',
            )
        # # ax1.plot(
            # # time_arr,
            # # max_chain,
            # # label='Max. of lengths',
            # # c='k',
            # # )
        ax1.tick_params(axis="y",direction="in", labelsize='x-large', labelcolor='k')
        ax1.tick_params(axis="x",direction="in", labelsize='x-large')
        ax1.set_xlabel('Time (ps)', fontsize='x-large')
        ax1.set_ylabel('Mean of chain lengths', fontsize='x-large')
        if deriv_sigma is not None and not False:
            ax2 = ax1.twinx()
            ax2.plot(
                time_arr[1:-1],
                dCdt*1000,
                label='$\Delta t$={}(ps)\nGaussian smearing $\sigma$={} (ps)'.format(self.dt, deriv_sigma *self.dt),
                c='r',
                )
            ax2.tick_params(axis="y",direction="in", labelsize='x-large', colors='r', labelcolor='r')
            ax2.set_ylabel('$\Delta C$/$\Delta t$ (ns$^{-1}$)', fontsize='x-large', c='r')
            plt.legend(loc=(0.00, 1.02), fontsize='x-large')
        # # ax2.plot(
            # # time_arr,
            # # num_chain,
            # # label='Num. of chains',
            # # c='b',
            # # )
        # ax2.plot(
            # time_arr,
            # sum_chain,
            # label='Sum of lengths',
            # c='k',
            # )
        plt.title('cut={} $\AA$ & {} deg, bond={}-{}, t-intvl={}'.format(
            bond_cutoff,
            angle_cutoff,
            self.bond_rules_str[0],
            self.bond_rules_str[1],
            self.dt,
            ), fontsize='x-large', y=1.25)
        plt.xlabel('Time (ps)', fontsize='x-large')
        ax1.grid(alpha=0.5)
        plt.subplots_adjust(left=0.10, bottom=0.12, right=0.85, top=0.77)
        plt.show()

    def get_chain_length_histo(
        self,
        angle_cutoff,
        bond_cutoff,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        """

        lengths = self.get_chain_lengths(
            angle_cutoff,
            bond_cutoff,
            bond_rules,
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
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        """
        #

        length_histo_list = self.get_chain_length_histo(
            angle_cutoff,
            bond_cutoff,
            bond_rules,
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
        l = np.array(list(length_histo.keys()))
        n = np.array(list(length_histo.values()))
        
        # Plot
        from matplotlib import pyplot as plt
        font = {'family':'Arial'}
        plt.rc('font', **font)
        fig, ax1 = plt.subplots()

        ax1.bar(l, n/np.sum(n)*100., color='k', alpha=0.5)

        xmax = np.max(l)
        plt.xlim(-1, xmax+1)
        xticks = range(0, xmax+1, 2)
        xticklabels = list(xticks)
        xticklabels[0] = 'inf'
        plt.xticks(xticks, rotation=45, labels=xticklabels)
        plt.xlabel('Length of chain', fontsize='x-large')

        title = 'cut={} $\AA$ & {} deg, bond={}-{}'.format(
            bond_cutoff,
            angle_cutoff,
            self.bond_rules_str[0],
            self.bond_rules_str[1],
            )
        if len(self.alist_ind_list) == 1:
            title += ', atoms_ind #{}'.format(self.alist_ind_list[0])
        else:
            title += ', len_alist #{}'.format(len(self.alist_ind_list))
        plt.title(title, fontsize='x-large')

        # Different scale on the right axis.
        # ax2 = ax1.twinx()
        # ax2.bar(l, n, color='k', alpha=0.5)
        ax1.tick_params(axis="both",direction="in", labelsize='x-large')
        # ax2.tick_params(axis="both",direction="in", labelsize='x-large')
        ax1.set_ylabel('Population (%)', fontsize='x-large')
        # ax2.set_ylabel('Population', fontsize='x-large')
        plt.subplots_adjust(bottom=0.14, left=0.10, right=0.88)
        ax1.grid(alpha=0.5)
        plt.show()

    def classify_3body_pieces_by_type(
        self,
        angle_cutoff,
        bond_cutoff,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        Return dictionary of classified 3-body pieces' info.
        """

        #
        cos_cutoff = np.cos(angle_cutoff / 180 * np.pi)
        #
        indices, lengths, directions = self.produce_neighbor_info(
            bond_cutoff,
            bond_rules,
            returns=['indices', 'lengths', 'directions'],
            load_bool=load_bool,
            save_bool=save_bool,
            )
        #
        piece_inds    = []
        piece_lengths = []
        # piece_direcs  = []
        for i in range(len(indices)):
            #
            path = 'neigh_saved/{}.d/{}-{}-{}/{}/{}'.format(
                self.alist_file,
                bond_cutoff,
                self.bond_rules_str[0],
                self.bond_rules_str[1],
                self.alist_ind_list[i],
                angle_cutoff,
                )
            piece_inds_i, piece_lengths_i, piece_direcs_i = get_3body_chain_pieces(
                indices[i],
                lengths[i],
                directions[i],
                angle_cutoff,
                path,
                path,
                )
            piece_inds   .append(piece_inds_i)
            piece_lengths.append(piece_lengths_i)
            # piece_direcs .append(piece_direcs_i)
        del(indices, lengths, directions)

        # Classify by species.
        length_dict = {}
        # direc_dict = {}
        for ty in self.types_unique:
            length_dict[ty] = []
            # direc_dict[ty] = []
        for i in range(len(piece_inds)):
            #
            length_i_dict = {}
            # direc_i_dict = {}
            for ty in self.types_unique:
                length_i_dict[ty] = []
                # direc_i_dict[ty] = []
            for j in range(len(piece_inds[i])):
                #
                ty = self.types[piece_inds[i][j][1]]
                length_i_dict[ty].append(piece_lengths[i][j])
                # direc_i_dict[ty].append(piece_direcs[i][j])
            #
            for ty in self.types_unique:
                #                      --> shape of (number of pieces of a type-kind in an image, 2)
                length_dict[ty].append(np.reshape(length_i_dict[ty], [-1, 2]))
                # direc_dict[ty].append(np.reshape(direc_i_dict[ty], [-1, 2, 3]))
        #      --> shape of (len_alist, (number of pieces of a type-kind in an image), 2)
        return length_dict #, direc_dict
        
    def plot_3body_pieces_stat(
        self,
        angle_cutoff,
        bond_cutoff,
        bond_rules = None,
        num_bins   = 100,
        cmap       = 'jet',
        load_bool  = True,
        save_bool  = True,
        ):
        """
        Plot bond-length correlation.
        """
        length_dict = self.classify_3body_pieces_by_type(
            angle_cutoff,
            bond_cutoff,
            bond_rules,
            load_bool,
            save_bool,
            )
        # Concate
        flat_length_dict = {}
        avgs = {}
        for ty in self.types_unique:
            concat = np.concatenate(length_dict[ty], axis=0)
            avgs[ty] = np.mean(concat)
            flat_length_dict[ty] = np.concatenate([concat, concat[:,::-1]], axis=0)

        # Plot
        from matplotlib import pyplot as plt
        font = {'family':'Arial'}
        plt.rc('font', **font)
        fig, [ax_1d_list, ax_2d_list] = plt.subplots(2,len(self.types_unique))
        ax_1d_list = list(ax_1d_list[::-1])
        ax_2d_list = list(ax_2d_list[::-1])
        for ty in self.types_unique:
            ax_1d = ax_1d_list.pop()
            ax_2d = ax_2d_list.pop()
            ax_1d.hist(
                flat_length_dict[ty][:,0],
                density=True,
                bins=num_bins,
                facecolor='k',
                alpha=0.8,
                )
            histo = ax_2d.hist2d(
                flat_length_dict[ty][:,0],
                flat_length_dict[ty][:,1],
                normed=True,
                bins=num_bins,
                cmap=cmap,
                )
            # y = x line
            bmin = np.min(histo[1])
            bmax = np.max(histo[1])
            ax_2d.plot([bmin, bmax], [bmin, bmax], c='k')
            # Style
            ax_2d.set_title('{}-center, {} pieces'.format(ty, len(flat_length_dict[ty])//2), fontsize='xx-large')
            # ax_2d.set_xlabel('mean={:.4f}'.format(avgs[ty]), fontsize='xx-large')
            ax_2d.set_aspect(1)
            ax_1d.tick_params(axis="both",direction="in", labelsize='xx-large')
            ax_2d.tick_params(axis="both",direction="in", labelsize='xx-large')
            ax_1d.set_ylabel('Normalized Population', fontsize='xx-large')
            cb = plt.colorbar(histo[3], ax=ax_2d, fraction=0.04, pad=0.03)
            if not ax_2d_list:
                cb.set_label('Normalized Density', fontsize='xx-large', labelpad=10)
            cb.ax.tick_params(axis="both",direction="in", labelsize='xx-large')
        plt.suptitle('{}, slice={}, AC={}, BC={}, BR={}'.format(
            self.alist_file,
            self.alist_slice,
            angle_cutoff,
            bond_cutoff,
            bond_rules,
            ), fontsize='xx-large')
        plt.subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, wspace=0.30)
        plt.show()

    def view_chains(
        self,
        angle_cutoff,
        bond_cutoff,
        bond_rules=None,
        chain_length=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        chain_length (int or 'inf' or None)
            - Chain length of that you wanna see.
            - 'inf' means infinitely long chains.
            - 0 (zero) also means infinitely long chains.
            - If None, show all.
        """

        if chain_length == 'inf':
            chain_length = 0

        ind_set, bond_direcs, bond_lengths, chain_vec = self.get_chain_set(
            angle_cutoff,
            bond_cutoff,
            bond_rules,
            load_bool,
            save_bool,
            )

        if chain_length is not None:
            chain_lengths = self.get_chain_lengths(
                angle_cutoff,
                bond_cutoff,
                bond_rules,
                True,
                load_bool,
                save_bool,
                )

        new_alist=[]
        for i in range(len(ind_set)):
            for j in range(len(ind_set[i])):
                if chain_length is not None:
                    if chain_lengths[i][j] == chain_length:
                        new_alist.append(self.alist[i][np.unique(ind_set[i][j])])
                    else:
                        pass
                else:
                    new_alist.append(self.alist[i][np.unique(ind_set[i][j])])
        from ase.visualize import view
        view(new_alist)

    def get_avg_coord_nums(
        self,
        bond_cutoff,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        """

        indices = self.produce_neighbor_info(
            bond_cutoff,
            bond_rules,
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
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        """

        num_bonds = self.get_avg_coord_nums(
            bond_cutoff,
            bond_rules,
            load_bool,
            save_bool,
            )

        # Plot
        time_arr = np.arange(len(list(num_bonds.values())[0])) *self.dt 
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
        plt.title('cut={} $\AA$, bond={}-{}, t-intvl={}'.format(bond_cutoff, self.bond_rules_str[0], self.bond_rules_str[1], self.dt), fontsize='x-large')
        plt.xlabel('Time (ps)', fontsize='x-large')
        plt.ylabel('Average coordination number', fontsize='x-large')
        plt.grid(alpha=0.5)
        plt.show()

    def get_positional_deviation(
        self,
        in_num_avg=10,
        out_num_avg=10,
        out_avg_dn=10,
        return_intvl=None,
        ):
        """
        in_num_avg (int)
            - For each time, the position will be defined as the average of positions of images in alist[i: i+in_num_avg].
        out_num_avg (int)
            - For each time, the positional deviation will be averaged for (out_num_avg) images.
        out_avg_dn (int)
            - Image sampling interval for averaging the positional-deviation.
        return_intvl (int or None)
            - Return the positional deviation every (return_intvl) steps.
            - Must be an integer multiple of the (out_avg_dn).
            - If None, Set to be same as (out_avg_dn).
        """
        #
        if not return_intvl:
            return_intvl = out_avg_dn
        # Gather data.
        positions = []
        cells     = []
        temps     = []
        for i in range(self.len_alist):
            positions.append(self.alist[i].get_scaled_positions())
            cells    .append(self.alist[i].get_cell())
            temps    .append(self.alist[i].get_temperature())
        positions = np.array(positions)
        cells     = np.array(cells)
        temps     = np.array(temps)

        # Inner average.
        print('Getting inner average.')
        avg_positions = []
        avg_cells     = []
        avg_temps     = []
        for i in range(0, self.len_alist -in_num_avg, out_avg_dn):
            avg_positions.append(np.mean(positions[i: i+in_num_avg], axis=0))
            avg_cells    .append(np.mean(cells    [i: i+in_num_avg], axis=0))
            avg_temps    .append(np.mean(temps    [i: i+in_num_avg], axis=0))
        avg_positions = np.array(avg_positions)
        avg_cells     = np.array(avg_cells    )
        avg_temps     = np.array(avg_temps    )
        len_bin = len(avg_positions)
        return_intvl //= out_avg_dn

        # Outer average.
        print('Getting outer average.')
        # avg_positions --> shape of (len_bin, len_atoms, 3)
        # avg_cells     --> shape of (len_bin, 3, 3)
        # avg_temps     --> shape of (len_bin)
        #
        mean_disp = []
        for i in range(0, len_bin -out_num_avg, return_intvl):
            mean_disp.append(
                np.mean(
                    np.linalg.norm(
                        np.reshape(((avg_positions[i+1: i+1 +out_num_avg] \
                            -np.expand_dims(avg_positions[i], axis=0) +np.array([0.5]*3) % 1.0) \
                            -np.array([0.5]*3)) @ cells[i], [-1, self.len_atoms, 3]),
                        axis=2,
                        ),
                    axis=0,
                    ),
                )

        # mean_disp --> shape of (len_bin, len_atoms)
        mean_disp = np.array(mean_disp)
        mean_disp_temp = mean_disp / np.sqrt(np.expand_dims(temps[:len(mean_disp)], axis=1))
        return mean_disp, mean_disp_temp

    def plot_positional_deviation(
        self,
        in_num_avg=10,
        out_num_avg=10,
        out_avg_dn=10,
        return_intvl=None,
        ):
        """
        """
        #
        if not return_intvl:
            return_intvl = out_avg_dn
        mean_disp, mean_disp_temp = self.get_positional_deviation(
            in_num_avg,
            out_num_avg,
            out_avg_dn,
            return_intvl,
            )

        dt = self.dt *return_intvl
        from matplotlib import pyplot as plt
        font = {'family':'Arial'}
        plt.rc('font', **font)
        plt.plot(np.arange(len(mean_disp)) *dt, np.mean(mean_disp, axis=1), c='r')
        # plt.plot(np.arange(len(mean_disp))*dt, np.mean(mean_disp_temp, axis=1), c='r')
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.title('slice={}, INA={}, ONA={}, OAD={}'.format(
            self.alist_slice,
            in_num_avg,
            out_num_avg,
            out_avg_dn,
            ), fontsize='x-large')
        plt.xlabel('Time (ps)', fontsize='x-large')
        plt.ylabel('Averaged positional deviation ($\AA$)', fontsize='x-large')
        plt.grid(alpha=0.5)
        plt.show()

