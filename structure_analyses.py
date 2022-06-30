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
        self.t = np.arange(self.len_alist, dtype=float) *self.dt
        # Make alist_ind_list
        self.alist_ind_list = np.arange(slice_list[1], dtype=np.int32)[slice(*slice_list)]

        # Types
        self.types = np.array(self.alist[0].get_chemical_symbols())
        self.types_unique = np.unique(self.types)
        self.types_dict = {}
        for ty in self.types_unique:
            self.types_dict[ty] = self.types == ty

        # Entries
        self.neighbor_info      = None
        self.chain_pieces       = None
        self.chain_sets         = None
        self.terminal_pieces    = None
        self.all_pieces         = None
        self.potential_energies = None

    def get_neighbor_info(
        self,
        bond_cutoff,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        if self.neighbor_info is None:
            self.produce_neighbor_info(bond_cutoff, bond_rules, load_bool, save_bool)
        indices    = self.neighbor_info['indices']
        lengths    = self.neighbor_info['lengths']
        directions = self.neighbor_info['directions']
        return indices, lengths, directions

    def produce_neighbor_info(
        self,
        bond_cutoff,
        bond_rules=None,
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
        self.neighbor_info = info

    def get_chain_pieces(
        self,
        bond_cutoff,
        angle_cutoff,
        bond_rules,
        load_bool=True,
        save_bool=True,
        ):
        if self.chain_pieces is None:
            self.produce_chain_pieces(bond_cutoff, angle_cutoff, bond_rules, load_bool, save_bool)
        piece_inds    = self.chain_pieces['piece_inds']
        piece_lengths = self.chain_pieces['piece_lengths']
        piece_direcs  = self.chain_pieces['piece_direcs']
        return piece_inds, piece_lengths, piece_direcs

    def produce_chain_pieces(
        self,
        bond_cutoff,
        angle_cutoff,
        bond_rules,
        load_bool=True,
        save_bool=True,
        ):
        """
        """
        indices, lengths, directions = self.get_neighbor_info(bond_cutoff, bond_rules, load_bool, save_bool)

        piece_inds    = []
        piece_lengths = []
        piece_direcs  = []
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
                assert load_bool
                print(' *  Trying to load chain piece pckl files.')
                piece_inds_i    = pckl.load(open('{}/piece_inds.pckl'   .format(path), 'rb'))
                piece_lengths_i = pckl.load(open('{}/piece_lengths.pckl'.format(path), 'rb'))
                piece_direcs_i  = pckl.load(open('{}/piece_direcs.pckl' .format(path), 'rb'))
            except:
                print('Failed to load pckl files from {}'.format(path))

            else:
                print('Loaded pckl files from {}'.format(path))
                piece_inds.append(piece_inds_i)
                piece_lengths.append(piece_lengths_i)
                piece_direcs.append(piece_direcs_i)
                continue

            # Main
            len_atoms = len(indices[i])
            #
            cos_cutoff = np.cos(angle_cutoff / 180 * np.pi)
            #
            piece_inds_i    = []
            piece_lengths_i = []
            piece_direcs_i  = []
            for j in range(len_atoms):
                indices_j    = np.array(indices[i][j])
                lengths_j    = np.array(lengths[i][j])
                directions_j = np.array(directions[i][j])
                # Gather three-body chain pieces.
                # Get upper triangle due to the duplications
                #--> shape of [number of bonds for atom j]*2
                cosines = np.triu(directions_j @ directions_j.T, 1)
                bond_inds = np.transpose(np.where(cosines < cos_cutoff))
                for inds in bond_inds:
                    #                    --> shape of (3,)
                    piece_inds_i   .append(list(np.concatenate((indices_j[inds], [j]))[[0,2,1]]))
                    #                    --> shape of (2,)
                    piece_lengths_i.append(lengths_j[inds])
                    #                    --> shape of (2, 3)
                    piece_direcs_i .append(directions_j[inds])
                    piece_direcs_i[-1][0] *= -1.
            piece_inds_i    = np.reshape(piece_inds_i,    [-1, 3]).tolist()
            piece_lengths_i = np.reshape(piece_lengths_i, [-1, 2]).tolist()
            piece_direcs_i  = list(np.reshape(piece_direcs_i,  [-1, 2, 3]))
            if save_bool:
                call('mkdir -p {}'.format(path), shell=True)
                pckl.dump(piece_inds_i   , open('{}/piece_inds.pckl'      .format(path), 'wb'))
                pckl.dump(piece_lengths_i, open('{}/piece_lengths.pckl'   .format(path), 'wb'))
                pckl.dump(piece_direcs_i , open('{}/piece_direcs.pckl'.format(path), 'wb'))
                print('Saved chain-piece-info pckl files at {}'.format(path))
            piece_inds.append(piece_inds_i)
            piece_lengths.append(piece_lengths_i)
            piece_direcs.append(piece_direcs_i)

        self.chain_pieces = {
            'piece_inds': piece_inds,
            'piece_lengths': piece_lengths,
            'piece_direcs': piece_direcs,
            }

    def get_chain_sets(
        self,
        bond_cutoff,
        angle_cutoff,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        if self.chain_sets is None:
            self.produce_chain_sets(bond_cutoff, angle_cutoff, bond_rules, load_bool, save_bool)
        ind_set      = self.chain_sets['ind_set']
        bond_direcs  = self.chain_sets['bond_direcs']
        bond_lengths = self.chain_sets['bond_lengths']
        chain_vec    = self.chain_sets['chain_vec']
        return ind_set, bond_direcs, bond_lengths, chain_vec

    def produce_chain_sets(
        self,
        bond_cutoff,
        angle_cutoff,
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

        # Gather three-body chain pieces.
        piece_inds, piece_lengths, piece_direcs = self.get_chain_pieces(
            bond_cutoff,
            angle_cutoff,
            bond_rules,
            load_bool,
            save_bool,
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

            #
            piece_inds_i    = piece_inds[i]
            piece_lengths_i = piece_lengths[i]
            piece_direcs_i  = piece_direcs[i]
            # Classify chain pieces.
            #--> Ragged tensor in shape of ((number of chains in an image), (length of a chain))
            ind_set_i      = []
            bond_direcs_i  = []
            bond_lengths_i = []
            chain_vec_i    = []
            while piece_inds_i:
                ind_set_tmp      = list(piece_inds_i   .pop())
                bond_lengths_tmp = list(piece_lengths_i.pop())
                bond_direcs_tmp  = list(piece_direcs_i .pop())
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
                    for p in range(len(piece_inds_i)):
                        if piece_inds_i[p][:2] == tail and (piece_direcs_i[p][1] @ chain_direc) > -cos_cutoff:
                            #
                            ind_set_tmp.append(piece_inds_i[p][2])
                            del(piece_inds_i[p])
                            #
                            bond_direcs_tmp.append(piece_direcs_i[p][1])
                            del(piece_direcs_i[p])
                            #
                            bond_lengths_tmp.append(piece_lengths_i[p][1])
                            del(piece_lengths_i[p])
                            #
                            chain_vec_tmp += bond_direcs_tmp[-1] * bond_lengths_tmp[-1]
                            #
                            changed = True
                            break
                        elif piece_inds_i[p][:2] == head[::-1] and (piece_direcs_i[p][1] @ (-chain_direc)) > -cos_cutoff:
                            #
                            ind_set_tmp = ind_set_tmp[::-1]
                            ind_set_tmp.append(piece_inds_i[p][2])
                            del(piece_inds_i[p])
                            #
                            bond_direcs_tmp = bond_direcs_tmp[::-1]
                            bond_direcs_tmp.append(piece_direcs_i[p][1])
                            del(piece_direcs_i[p])
                            #
                            bond_lengths_tmp = bond_lengths_tmp[::-1]
                            bond_lengths_tmp.append(piece_lengths_i[p][1])
                            del(piece_lengths_i[p])
                            #
                            chain_vec_tmp *= -1.
                            chain_vec_tmp += bond_direcs_tmp[-1] * bond_lengths_tmp[-1]
                            #
                            changed = True
                            break
                        elif piece_inds_i[p][-2:] == head and (piece_direcs_i[p][0] @ chain_direc) > -cos_cutoff:
                            #
                            ind_set_tmp = ind_set_tmp[::-1]
                            ind_set_tmp.append(piece_inds_i[p][0])
                            del(piece_inds_i[p])
                            #
                            bond_direcs_tmp = bond_direcs_tmp[::-1]
                            bond_direcs_tmp.append(-piece_direcs_i[p][0])
                            del(piece_direcs_i[p])
                            #
                            bond_lengths_tmp = bond_lengths_tmp[::-1]
                            bond_lengths_tmp.append(piece_lengths_i[p][0])
                            del(piece_lengths_i[p])
                            #
                            chain_vec_tmp *= -1.
                            chain_vec_tmp += bond_direcs_tmp[-1] * bond_lengths_tmp[-1]
                            #
                            changed = True
                            break
                        elif piece_inds_i[p][-2:] == tail[::-1] and (-piece_direcs_i[p][0] @ chain_direc) > -cos_cutoff:
                            #
                            ind_set_tmp.append(piece_inds_i[p][0])
                            del(piece_inds_i[p])
                            #
                            bond_direcs_tmp.append(-piece_direcs_i[p][0])
                            del(piece_direcs_i[p])
                            #
                            bond_lengths_tmp.append(piece_lengths_i[p][0])
                            del(piece_lengths_i[p])
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

            #--> Ragged tensor in shape of (len_alist, (number of chains in an image), (length of a chain)+1)
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
        self.chain_sets = {
            'ind_set': ind_set,
            'bond_direcs': bond_direcs,
            'bond_lengths': bond_lengths,
            'chain_vec': chain_vec,
            }

    def get_terminal_pieces(
        self,
        bond_cutoff,
        angle_cutoff,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        """
        if self.terminal_pieces is None:
            self.produce_terminal_pieces(bond_cutoff, angle_cutoff, bond_rules, load_bool, save_bool)
        term_piece_inds    = self.terminal_pieces['term_piece_inds']
        term_piece_lengths = self.terminal_pieces['term_piece_lengths']
        term_piece_direcs  = self.terminal_pieces['term_piece_direcs']
        return term_piece_inds, term_piece_lengths, term_piece_direcs

    def produce_terminal_pieces(
        self,
        bond_cutoff,
        angle_cutoff,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        """

        ind_set, bond_direcs, bond_lengths, chain_vec = self.get_chain_sets(
            bond_cutoff,
            angle_cutoff,
            bond_rules,
            load_bool,
            save_bool,
            )

        term_piece_inds    = []
        term_piece_direcs  = []
        term_piece_lengths = []
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
                assert load_bool
                print(' *  Trying to load terminal piece pckl files.')
                term_piece_inds_i    = pckl.load(open('{}/terminal_piece_inds.pckl'   .format(path), 'rb'))
                term_piece_direcs_i  = pckl.load(open('{}/terminal_piece_direcs.pckl' .format(path), 'rb'))
                term_piece_lengths_i = pckl.load(open('{}/terminal_piece_lengths.pckl'.format(path), 'rb'))
            except:
                print('Failed to load pckl files from {}'.format(path))

            else:
                print('Loaded pckl files from {}'.format(path))
                term_piece_inds.append(term_piece_inds_i)
                term_piece_direcs.append(term_piece_direcs_i)
                term_piece_lengths.append(term_piece_lengths_i)
                continue
            
            term_piece_inds_i    = []
            term_piece_lengths_i = []
            term_piece_direcs_i  = []
            for j in range(len(ind_set[i])):
                term_piece_inds_i.append([-1, ind_set[i][j][0], ind_set[i][j][1]])
                term_piece_direcs_i.append(np.array([np.zeros(3, dtype=float), bond_direcs[i][j][0]]))
                term_piece_lengths_i.append([0., bond_lengths[i][j][0]])

                term_piece_inds_i.append([ind_set[i][j][-2], ind_set[i][j][-1], -1])
                term_piece_direcs_i.append(np.array([bond_direcs[i][j][-1], np.zeros(3, dtype=float)]))
                term_piece_lengths_i.append([bond_lengths[i][j][-1], 0.])

            if save_bool:
                call('mkdir -p {}'.format(path), shell=True)
                pckl.dump(term_piece_inds_i   , open('{}/terminal_piece_inds.pckl'      .format(path), 'wb'))
                pckl.dump(term_piece_lengths_i, open('{}/terminal_piece_lengths.pckl'   .format(path), 'wb'))
                pckl.dump(term_piece_direcs_i , open('{}/terminal_piece_direcs.pckl'.format(path), 'wb'))
                print('Saved terminal-piece-info pckl files at {}'.format(path))
            term_piece_inds.append(term_piece_inds_i)
            term_piece_lengths.append(term_piece_lengths_i)
            term_piece_direcs.append(term_piece_direcs_i)

        self.terminal_pieces = {
            'term_piece_inds': term_piece_inds,
            'term_piece_lengths': term_piece_lengths,
            'term_piece_direcs': term_piece_direcs,
            }

    def get_all_pieces(
        self,
        bond_cutoff,
        angle_cutoff,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        """

        if self.all_pieces is None:
            piece_inds, piece_lengths, piece_direcs = self.get_chain_pieces(
                bond_cutoff,
                angle_cutoff,
                bond_rules,
                load_bool,
                save_bool,
                )

            term_piece_inds, term_piece_lengths, term_piece_direcs = self.get_terminal_pieces(
                bond_cutoff,
                angle_cutoff,
                bond_rules,
                load_bool,
                save_bool,
                )

            all_piece_inds    = []
            all_piece_lengths = []
            all_piece_direcs  = []
            for i in range(len(self.alist_ind_list)):
                all_piece_inds   .append(piece_inds   [i] + term_piece_inds   [i])
                all_piece_lengths.append(piece_lengths[i] + term_piece_lengths[i])
                all_piece_direcs .append(piece_direcs [i] + term_piece_direcs [i])

            self.all_pieces = {
                'all_piece_inds': all_piece_inds,
                'all_piece_lengths': all_piece_lengths,
                'all_piece_direcs': all_piece_direcs,
                }

        else:
            all_piece_inds    = self.all_pieces['all_piece_inds']
            all_piece_lengths = self.all_pieces['all_piece_lengths']
            all_piece_direcs  = self.all_pieces['all_piece_direcs']

        return all_piece_inds, all_piece_lengths, all_piece_direcs

    def classify_chain_pieces(
        self,
        bond_cutoff,
        angle_cutoff,
        bond_rules=None,
        unique_seq=None,
        include_terminal=True,
        merge_mirror=True,
        load_bool=True,
        save_bool=True,
        ):
        """
        unique_seq (dict)
            - You can specify it. If None, automatic.
        include_terminal (bool)
            - Include terminal pieces (including vacancy) for classification.
        merge_mirror (bool)
            - Merge mirror-image pieces in to a same set. e.g., A-B-C and C-B-A.
        """

        if include_terminal:
            piece_inds, piece_lengths, piece_direcs = self.get_all_pieces(
                bond_cutoff,
                angle_cutoff,
                bond_rules,
                load_bool,
                save_bool,
                )
            chem = np.array(self.types.tolist() + ['V'])
        else:
            piece_inds, piece_lengths, piece_direcs = self.get_chain_pieces(
                bond_cutoff,
                angle_cutoff,
                bond_rules,
                load_bool,
                save_bool,
                )
            chem = self.types

        # unique_seq.shape = (len(self.type_unique), number of unique sequences)
        # piece_addr.shape = (len(self.type_unique), number of unique sequences, number of pieces in sequence, 2)
        if unique_seq is None:
            unique_seq = dict()
        piece_addr = dict()
        for ty in self.types_unique:
            if ty in unique_seq:
                piece_addr[ty] = [[] for _ in range(len(unique_seq[ty]))]
            else:
                unique_seq[ty] = []
                piece_addr[ty] = []

        for i in range(len(self.alist_ind_list)):
            for j in range(len(piece_inds[i])):
                piece_chem = chem[piece_inds[i][j]].tolist()
                center_chem = piece_chem[1]
                if piece_chem in unique_seq[center_chem]:
                    seq_i = unique_seq[center_chem].index(piece_chem)
                    piece_addr[center_chem][seq_i].append([i,j])
                elif merge_mirror and piece_chem[::-1] in unique_seq[center_chem]:
                    seq_i = unique_seq[center_chem].index(piece_chem[::-1])
                    piece_addr[center_chem][seq_i].append([i,j])
                else:
                    unique_seq[center_chem].append(piece_chem)
                    piece_addr[center_chem].append([])
                    piece_addr[center_chem][-1].append([i,j])

        return unique_seq, piece_addr

    def get_chain_lengths(
        self,
        bond_cutoff,
        angle_cutoff,
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
        ind_set, bond_direcs, bond_lengths, chain_vec = self.get_chain_sets(
            bond_cutoff,
            angle_cutoff,
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
        bond_cutoff,
        angle_cutoff,
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
            bond_cutoff,
            angle_cutoff,
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
        plt.subplots_adjust(left=0.25, bottom=0.25, right=0.75, top=0.75, wspace=0.2, hspace=0.2)
        plt.show()

    def get_chain_length_histo(
        self,
        bond_cutoff,
        angle_cutoff,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        """

        lengths = self.get_chain_lengths(
            bond_cutoff,
            angle_cutoff,
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
        bond_cutoff,
        angle_cutoff,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        """
        #

        length_histo_list = self.get_chain_length_histo(
            bond_cutoff,
            angle_cutoff,
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
        plt.subplots_adjust(left=0.25, right=0.75, bottom=0.25, top=0.75)
        ax1.grid(alpha=0.5)
        plt.show()

    def classify_chain_length_by_type(
        self,
        bond_cutoff,
        angle_cutoff,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        Return dictionary of classified 3-body pieces' info.
        """

        #
        piece_inds, piece_lengths, piece_direcs = self.get_chain_pieces(
            bond_cutoff,
            angle_cutoff,
            bond_rules,
            load_bool,
            save_bool,
            )

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
        
    def plot_chain_pieces_stat(
        self,
        bond_cutoff,
        angle_cutoff,
        bond_rules = None,
        num_bins   = 100,
        cmap       = 'jet',
        load_bool  = True,
        save_bool  = True,
        ):
        """
        Plot bond-length correlation.
        """
        length_dict = self.classify_chain_length_by_type(
            bond_cutoff,
            angle_cutoff,
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
        bond_cutoff,
        angle_cutoff,
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

        ind_set, bond_direcs, bond_lengths, chain_vec = self.get_chain_sets(
            bond_cutoff,
            angle_cutoff,
            bond_rules,
            load_bool,
            save_bool,
            )

        if chain_length is not None:
            chain_lengths = self.get_chain_lengths(
                bond_cutoff,
                angle_cutoff,
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

        indices, lengths, directions = self.get_neighbor_info(
            bond_cutoff,
            bond_rules,
            load_bool=load_bool,
            save_bool=save_bool,
            )

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

    def get_terminal_list(
        self,
        bond_cutoff,
        angle_cutoff,
        nth_term,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        nth_term (int): Specify which n-th terminal to be counted.
            For example, nth_term=2, in a chain of "Te-(Ge)-Te-Sb-Te-Sb-Te-(Ge)-Te", two Ge atoms in parentheses are counted.
        """

        ind_set, bond_direcs, bond_lengths, chain_vec = self.get_chain_sets(
            bond_cutoff,
            angle_cutoff,
            bond_rules,
            load_bool,
            save_bool,
            )

        nth_term_py = nth_term -1

        terminals = dict()
        for ty in self.types_unique:
            terminals[ty] = [[] for i in range(len(ind_set))]
        for i in range(len(ind_set)):
            # Iter for chains
            for j in range(len(ind_set[i])):
                type_chain = self.types[ind_set[i][j]]
                terminals[type_chain[nth_term_py]][i].append(ind_set[i][j][nth_term_py])
                if nth_term_py < (len(type_chain) -nth_term_py -1):
                    terminals[type_chain[len(type_chain) -nth_term_py -1]][i].append(ind_set[i][j][len(type_chain) -nth_term_py -1])
        return terminals

    def plot_terminal_stat(
        self,
        bond_cutoff,
        angle_cutoff,
        nth_term,
        bond_rules=None,
        color_list=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        """

        terminals = self.get_terminal_list(bond_cutoff, angle_cutoff, nth_term, bond_rules, load_bool, save_bool)

        # Iter for images
        term_hist = np.zeros((len(self.types_unique), self.len_alist))
        for i in range(self.len_alist):
            for j in range(len(self.types_unique)):
                term_hist[j, i] = len(terminals[self.types_unique[j]][i])
        
        # Plot
        if color_list is None:
            color_list = [None] *len(self.types_unique)

        #
        from matplotlib import pyplot as plt
        for type_i in range(len(self.types_unique)):
            plt.plot(self.t, term_hist[type_i], label=self.types_unique[type_i], c=color_list[type_i])
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.xlim(self.t[0], self.t[-1])
        plt.ylim(0, None)
        plt.title('{}-th terminal histogram'.format(nth_term), fontsize='x-large')
        plt.xlabel('Time (ps)', fontsize='x-large')
        plt.ylabel('Population', fontsize='x-large')
        # plt.subplots_adjust(left=0.35, right=0.65)
        plt.grid(alpha=0.4)
        plt.legend(fontsize='large').set_draggable(True)
        plt.subplots_adjust(left=0.25, right=0.75, bottom=0.25, top=0.75)

        # Normalized population
        term_hist_norm = term_hist /np.sum(term_hist, axis=0)
        plt.figure()
        for type_i in range(len(self.types_unique)):
            plt.plot(self.t, term_hist_norm[type_i], label=self.types_unique[type_i], c=color_list[type_i])
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.xlim(self.t[0], self.t[-1])
        plt.ylim(0, None)
        plt.title('{}-th terminal histogram'.format(nth_term), fontsize='x-large')
        plt.xlabel('Time (ps)', fontsize='x-large')
        plt.ylabel('Population ratio', fontsize='x-large')
        # plt.subplots_adjust(left=0.35, right=0.65)
        plt.grid(alpha=0.4)
        plt.legend(fontsize='large').set_draggable(True)
        plt.subplots_adjust(left=0.25, right=0.75, bottom=0.25, top=0.75)
        plt.show()

    def view_terminal_atoms(
        self,
        bond_cutoff,
        angle_cutoff,
        nth_term,
        plot_chem=None,
        bond_rules=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        plot_chem (list of str or None)
            Specify chemical symbols to view. e.g., ['Ge', 'Sb']
        """

        #
        if plot_chem is None:
            plot_chem = self.types_unique
        else:
            plot_chem = np.array(plot_chem)

        #
        terminals = self.get_terminal_list(bond_cutoff, angle_cutoff, nth_term, bond_rules, load_bool, save_bool)

        new_alist=[]
        for i in range(self.len_alist):
            mask = []
            for ty in plot_chem:
                mask = mask + terminals[ty][i]
            new_alist.append(self.alist[i][mask])
        from ase.visualize import view
        view(new_alist)

    def get_potential_energies(
        self,
        ):
        """
        """
        if self.potential_energies is None:
            ae = []
            for i in range(len(self.alist)):
                ae.append(self.alist[i].get_potential_energies())
            self.potential_energies = np.array(ae, dtype=float)
        return self.potential_energies
            
    def plot_atomic_energy_histo(
        self,
        bond_cutoff,
        angle_cutoff,
        bond_rules=None,
        unique_seq=None,
        include_terminal=True,
        num_bins=50,
        color_list=None,
        load_bool=True,
        save_bool=True,
        ):
        """
        """

        if include_terminal:
            piece_inds, piece_lengths, piece_direcs = self.get_all_pieces(
                bond_cutoff,
                angle_cutoff,
                bond_rules,
                load_bool,
                save_bool,
                )
            chem = np.array(self.types.tolist() + ['V'])
        else:
            piece_inds, piece_lengths, piece_direcs = self.get_chain_pieces(
                bond_cutoff,
                angle_cutoff,
                bond_rules,
                load_bool,
                save_bool,
                )
            chem = self.types

        unique_seq, piece_addr = self.classify_chain_pieces(
            bond_cutoff,
            angle_cutoff,
            bond_rules,
            unique_seq,
            include_terminal,
            load_bool,
            save_bool,
            )

        pe = self.get_potential_energies()

        # pe_cls.shape = (len(self.types_unique), number of unique sequences, number of pieces in sequence)
        pe_cls = dict()
        for ty in self.types_unique:
            pe_cls[ty] = []
            for i in range(len(unique_seq[ty])):
                pe_cls[ty].append([])
                for j in range(len(piece_addr[ty][i])):
                    pa = piece_addr[ty][i][j]
                    pe_cls[ty][i].append(pe[pa[0], piece_inds[pa[0]][pa[1]][1]])

        pe_hist = dict()
        pe_bins = dict()
        for ty in self.types_unique:
            pe_hist[ty] = []
            amin = np.amin(np.concatenate(pe_cls[ty]))
            amax = np.amax(np.concatenate(pe_cls[ty]))
            for i in range(len(pe_cls[ty])):
                hist, bin_edges = np.histogram(pe_cls[ty][i], bins=num_bins, range=(amin, amax))
                pe_hist[ty].append(hist)
            pe_bins[ty] = (bin_edges[:-1] + bin_edges[1:])/2.

        # Plot energy hierarchy
        if color_list is None:
            cl = ['b', 'g', 'r', 'c', 'm', 'y'] *10
        else:
            cl = color_list
        from matplotlib import pyplot as plt
        pe_means = {}
        for ty in self.types_unique:
            pe_means[ty] = []
            print(' > Center atom: {}'.format(ty))
            plt.figure()
            for i in range(len(unique_seq[ty])):
                pe_means[ty].append(np.mean(pe_cls[ty][i]))
                print('   > Average atomic energy of {} : {:.4f} eV  {:.4f} meV'.format(
                    unique_seq[ty][i],
                    pe_means[ty][i],
                    (pe_means[ty][i] -np.mean(np.concatenate(pe_cls[ty]))) *1000,
                    ))
                plt.plot([0,1], [pe_means[ty][i]]*2, c=cl[i])
            plt.title('Atomic energy hierarchy of {} atoms'.format(ty), fontsize='x-large')
            plt.xticks([])
            plt.yticks(pe_means[ty],
                labels = ['{}-{}-{} {:.3f} eV'.format(*unique_seq[ty][i], pe_means[ty][i]) for i in range(len(unique_seq[ty]))])
            plt.tick_params(axis="both",direction="in", labelsize='large')
            plt.subplots_adjust(left=0.48, bottom=0.10, right=0.52, top=0.90, wspace=0.2, hspace=0.2)
            plt.grid(alpha=0.4)

        # Plot histogram
        for ty in self.types_unique:
            plt.figure()
            for i in range(len(pe_hist[ty])):
                plt.plot(pe_bins[ty], pe_hist[ty][i], label='{}-{}-{}'.format(*unique_seq[ty][i]), c=cl[i])
                plt.axvline(pe_means[ty][i], c=cl[i], lw=1)
            plt.xlim((pe_bins[ty][0], pe_bins[ty][-1]))
            plt.ylim((0, None))
            plt.title('Atomic energies of {} atoms'.format(ty), fontsize='x-large')
            plt.xlabel('Atomic energy (eV)', fontsize='x-large')
            plt.ylabel('Population', fontsize='x-large')
            plt.legend(fontsize='large').set_draggable(True)
            plt.tick_params(axis="both",direction="in", labelsize='x-large')
            plt.subplots_adjust(left=0.25, bottom=0.25, right=0.75, top=0.75, wspace=0.2, hspace=0.2)
            plt.grid(alpha=0.4)
        plt.show()

