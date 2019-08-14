#!/usr/bin/env python
##### CODE BY YOUNG JAE CHOI #####

from ase import Atoms, Atom
import random
from ase.build import make_supercell
from numpy import ndarray
import numpy as np

def get_chem_ind_arr(chemical_symbols):
    chem_arr = np.array(chemical_symbols)
    unique_chem = np.unique(chem_arr)
    ind = np.arange(len(chem_arr))
    chem_ind_arr = []
    for chem in unique_chem:
        chem_ind_arr.append(ind[chem_arr==chem])
    return unique_chem, np.array(chem_ind_arr)

def read_txt_matrix(file_name, start_line=0, end_line=-1):
    mat = []
    with open(file_name, 'r') as f:
        for i in range(start_line):
            f.readline()
        if end_line != -1:
            for i in range(end_line - start_line + 1):
                line = f.readline()
                seg = line.split()
                mat.append(seg)
        else:
            while True:
                line = f.readline()
                if not line: break
                seg = line.split()
                mat.append(seg)
    return mat

def nospace(string): return string.replace(' ', '_')
def E_fromAlist(alist):
    energies = []
    for atoms in alist:
        energies.append(atoms.get_potential_energy())
    return np.array(energies)

def F_fromAlist(alist):
    forces = []
    for atoms in alist:
        forces.append(atoms.get_forces())
    return np.array(forces)

def column2np(
    txtfile,
    column,
    start_line = 1,
    end_line   = -1,
    interval   = 1,
    ):
    txt = open(txtfile, "r")
    lines = txt.readlines()
    if end_line == -1:
        end_line = len(lines)
    nums = []
    for l in range(start_line, end_line, interval):
        nums.append(float(lines[l].split()[column-1]))
    array = np.array(nums)
    return array

def RanPoAtoms(cut_off_radius,
               symbols=None,
               positions=None, numbers=None,
               tags=None, momenta=None, masses=None,
               magmoms=None, charges=None,
               scaled_positions=None,
               cell=None, pbc=None, celldisp=None,
               constraint=None,
               calculator=None,
               info=None):
    if positions is not None:
        print("\npositions must not be given\n")
        exit(1)
    if scaled_positions is not None:
        print("\nscaled_positions must not be given\n")
        exit(1)
    else:
        atoms = Atoms(symbols=symbols,
                      positions=positions, numbers=numbers,
                      tags=tags, momenta=momenta, masses=masses,
                      magmoms=magmoms, charges=charges,
                      scaled_positions=None,
                      cell=cell, pbc=pbc, celldisp=celldisp,
                      constraint=constraint,
                      calculator=calculator,
                      info=info)    
        l = 0
        while True:
            l+=1
            print("trying step :: "+str(l))
            scaled_posis = []
            for i in range(len(atoms)):
                scaled_posi = []
                for j in range(3):
                    scaled_posi.append(random.random())
                scaled_posis.append(scaled_posi)
            atoms.set_scaled_positions(scaled_posis)
            supercell = make_supercell(atoms,[[2,0,0],[0,2,0],[0,0,2]])
            dist = supercell.get_all_distances()
            coll = []
            for i in range(len(supercell)):
                for j in range(len(supercell)):
                    if i is not j:
                        coll.append(dist[i][j])
            if min(coll) > cut_off_radius:
                break
        return atoms

def RanPoAtoms_2(cut_off_radius,
                random_degree,
                symbols=None,
                positions=None, numbers=None,
                tags=None, momenta=None, masses=None,
                magmoms=None, charges=None,
                scaled_positions=None,
                cell=None, pbc=None, celldisp=None,
                constraint=None,
                calculator=None,
                info=None):
    if positions is None:
        if scaled_positions is None:
            print("\nNo backbone structure is given.\n")
            exit(1)
    else:
        atoms = Atoms(symbols=symbols,
                      positions=positions, numbers=numbers,
                      tags=tags, momenta=momenta, masses=masses,
                      magmoms=magmoms, charges=charges,
                      scaled_positions=scaled_positions,
                      cell=cell, pbc=pbc, celldisp=celldisp,
                      constraint=constraint,
                      calculator=calculator,
                      info=info)    

        ############### shuffle positions ################
        array_scaled_positions = atoms.get_scaled_positions()
        shuffled_scaled_posis = array_scaled_positions.tolist()
        random.shuffle(shuffled_scaled_posis)
        atoms.set_scaled_positions(shuffled_scaled_posis)
        
        ############### get local random distribution radius ################
        supercell = make_supercell(atoms,[[2,0,0],[0,2,0],[0,0,2]])
        dist = supercell.get_all_distances()
        coll = []
        for i in range(len(supercell)):
            for j in range(len(supercell)):
                if i is not j:
                    coll.append(dist[i][j])
        ran_radi = min(coll)

        ############### shuffled position list ################
        array_shuffled_posis = atoms.get_positions()

        l = 0
        while True:
            l+=1
            print("trying step :: "+str(l))
            shuffled_posis = array_shuffled_posis.tolist()
            for i in range(len(atoms)):
                for j in range(3):
                    shuffled_posis[i][j]+=((random.random()-0.5)*2*ran_radi*random_degree)
            tmp = atoms
            tmp.set_positions(shuffled_posis)
            supercell = make_supercell(tmp,[[2,0,0],[0,2,0],[0,0,2]])
            dist = supercell.get_all_distances()
            coll = []
            for i in range(len(supercell)):
                for j in range(len(supercell)):
                    if i is not j:
                        coll.append(dist[i][j])
            if min(coll) > cut_off_radius:
                break
        atoms = tmp
        return atoms

def count2list(dict_in):
    list_out = []
    for i in range(len(dict_in.keys())):
        key = list(dict_in.keys())[i]
        num = 0
        itera = dict_in[key]
        for j in range(itera):
            list_out.append(key)
            num += 1
    return list_out

def list2count(list_inp):
    keys = list(set(list_inp))
    dict_out = dict()
    for i in keys:
        dict_out[i] = 0
    for i in list_inp:
        dict_out[i] += 1
    return dict_out

def list2numlist(list_inp):
    numlist = []
    keys=dict()
    i=0
    for ele in list_inp:
        if ele in keys:
            numlist.append(keys[ele])
        else:
            keys[ele]=i
            numlist.append(keys[ele])
            i+=1
    return numlist

def ordered_unique(in_list):
    tmp = set()
    return [x for x in in_list if not (x in tmp or tmp.add(x))]

def covalent_expect(input):
    """ Returns covalent bond expectation value of system 
    input : dict or list
    e.g. {'Si':1, 'H':2} or ['Si', 'H', 'H']

    """

    from ase.atoms import symbols2numbers as s2n
    from ase.data import covalent_radii
    
    if isinstance(input, list):
        species_dict = list2count(input)
    elif isinstance(input, dict):
        species_dict = input
    else:
        raise TypeError("input is not list nor dict")
    
    tot_num = sum(species_dict.values())
    r_sum = 0
    for key, value in species_dict.items():
        r_sum += covalent_radii[s2n([key])] * value
    expect_value = r_sum / tot_num
    return expect_value

def random_position_generator(
    backbone,
    species_kinds = None,
    species_spec  = None,
    cutoff_radi   = None,
    cutoff_frac   = None,
    random_degree = 0.9,
    strain        = None,
    strain_ratio  = [1.,1.,1.],
    vacuum        = None,
    vacuum_ratio  = None,
    ):
    """ 
    generate randomly positioned image base on backbone structure

    backbone : An ASE atoms object
        Which will be backbone position.
    species_kinds : list or None
        List of species kinds
        Identical to that of backbone object if None is provided.
    species_spec : list of int or None
        Number of each species. Will be numbered sequently.
        Identical to that of backbone object if None is provided.
    cutoff_radi : Float
        Cutoff radius for minimal distance between atoms.
        If provided with cutoff_frac simultaneously, occurs error.
    cutoff_frac : Float
        Set cutoff_radi as length scaled as 
        expectation value of covalent bonds of every atomic species.
        If provided with cutoff_radi simultaneously, occurs error.
    random_degree : Float
        Value of how much fraction of half of RDF nearest neighbor distance of backbone
        will be used as radius (from backbone) of generating new candidates.
    strain : List of three floats e.g. [0,0,5]
        Values specify how much you magnify the provided backbone cell.
        Cell gets longer along lattice vectors.
        positions will stay scaled positions of backbone atoms.
    strain_ratio : List of three floats e.g. [1,1.1,1]
        Values specify how much you magnify the provided backbone cell.
        Cell gets longer along lattice vectors.
        positions will stay scaled positions of backbone atoms.
    vacuum : List of three floats e.g. [0,0,5]
        Values specify how much you magnify the provided backbone cell with vacuum.
        Cell gets longer along lattice vectors.
        positions will stay absolute positions of backbone atoms.
        insert vacuum after strain (if provided)
    vacuum_ratio : List of three floats e.g. [1,1.1,1]
        Values specify how much you magnify the provided backbone cell with vacuum.
        Cell gets longer along lattice vectors.
        positions will stay absolute positions of backbone atoms.
        insert vacuum after strain (if provided)

    """

    backbone = backbone.copy()
    ############### collect species_spec 
    if species_kinds is None and species_spec is None:
        species = backbone.get_chemical_symbols()
    elif isinstance(species_kinds, list):
        species = []
        for i in range(len(species_kinds)):
            for j in range(species_spec[i]):
                species.append(species_kinds[i])
    else:
        raise ValueError('Somethings wrong with species_kinds and species_spec variables.'
            ' Please Check')

    ############## covalent bond length expectation value
    coval_expect = covalent_expect(species)
                
    ############# cell strain adjust
    if strain_ratio is not None and strain is not None:
        raise ValueError("strain_ratio & strain parameters provided simultaneously. \
            Just provide one.")
    if strain is not None:
        strain = np.array(strain)
        if strain.shape != (3,):
            raise ValueError("Somethings wrong with strain parameter. Please check.")
        norm = np.linalg.norm(backbone.cell, axis=1)
        strain_ratio = strain / norm + 1
    if strain_ratio is not None:
        strain_ratio = np.array(strain_ratio)
        if strain_ratio.shape != (3,):
            raise ValueError("Somethings wrong with strain_ratio parameter. Please check.")
        backbone.set_cell(
            backbone.cell * np.expand_dims(strain_ratio, axis=1),
            scale_atoms = True,
            )
    if strain_ratio is None and strain is None:
        strain_ratio = [1.,1.,1.]
        backbone.set_cell(
            backbone.cell * np.expand_dims(strain_ratio, axis=1),
            scale_atoms = True,
            )

    ############# cell vacuum adjust
    if vacuum_ratio is not None and vacuum is not None:
        raise ValueError("vacuum_ratio & vacuum parameters provided simultaneously. \
            Just provide one.")
    if vacuum is not None:
        vacuum = np.array(vacuum)
        if vacuum.shape != (3,):
            raise ValueError("Somethings wrong with vacuum parameter. Please check.")
        norm = np.linalg.norm(backbone.cell, axis=1)
        vacuum_ratio = vacuum / norm + 1
    if vacuum_ratio is not None:
        vacuum_ratio = np.array(vacuum_ratio)
        if vacuum_ratio.shape != (3,):
            raise ValueError("Somethings wrong with vacuum_ratio parameter. Please check.")
        backbone.set_cell(
            backbone.cell * np.expand_dims(vacuum_ratio, axis=1),
            scale_atoms = False,
            )
    if vacuum_ratio is None and vacuum is None:
        vacuum_ratio = [1.,1.,1.]
        backbone.set_cell(
            backbone.cell * np.expand_dims(vacuum_ratio, axis=1),
            scale_atoms = True,
            )

    ############## determine cutoff radius
    if cutoff_radi is not None and cutoff_frac is not None:
        raise ValueError("cutoff_radi & cutoff_frac parameters provided simultaneously. \
            Just provide one.")
    if cutoff_radi is not None:
        cutoff_r = cutoff_radi
    elif cutoff_frac is not None:
        cutoff_r = coval_expect * 2 * cutoff_frac
    else:
        cutoff_r = coval_expect * 2 * 0.9

    ############### get random adjust radius
    supercell = make_supercell(backbone,[[2,0,0],[0,2,0],[0,0,2]])
    from ase.optimize.precon.neighbors import estimate_nearest_neighbour_distance as rNN
    rdf_1st_peak = rNN(supercell)
    ran_radi = rdf_1st_peak / 2 * random_degree
    print("")
    print("********* Please check carefully !!!! ***********".center(80))
    print(("RDF 1st peak / 2 == %.2f" %(rdf_1st_peak/2)).center(80))
    print(("random radius degree == %.2f" %(random_degree)).center(80))
    print(("==> random adjust radius == %.2f" %(ran_radi)).center(80))
    print(("it is %.2f %% of covalent bond length expectation value." %
        (ran_radi / coval_expect * 100)).center(80))
    print("")
    print(("cf ) covalent bond length expectation value == %.2f" % coval_expect).center(80))
    print(("cf ) cutoff radius == %.2f" % cutoff_r).center(80))
    print(("cf ) cutoff radius / covalent bond expectation*2 == %.2f %%" % 
        (cutoff_r / coval_expect / 2 * 100)).center(80))
    print("")

    # ############### shuffle positions
    # backbone.set_positions(
        # np.random.permutation(backbone.get_positions()), apply_constraint = False)

    ############### Main
    new_atoms = backbone.copy()
    natoms = len(backbone)
    from time import time
    for i in range(natoms):
        new_atoms.pop()
    empty_atoms = new_atoms.copy()

    while len(new_atoms) < natoms:
        time_i = time()
        while True:
            new_atoms.append(backbone.copy()[len(new_atoms)])
            posi = new_atoms.get_positions()
            posi[-1] += (np.random.rand(3)-0.5) * 2 * ran_radi
            new_atoms.set_positions(posi, apply_constraint = False)
            dist = new_atoms.get_all_distances(mic = True)
            time_f = time()
            time_d = time_f - time_i
            if np.amin(dist + np.eye(len(dist))*100) > cutoff_r:
                print("( %d th / %d ) new atom position found" % (len(new_atoms), natoms))
                break
            elif time_d > 5:
                break
            else:
                new_atoms.pop()
        if time_d > 5:
            new_atoms = empty_atoms.copy()

    ############### shuffle positions
    new_atoms.set_positions(
        np.random.permutation(new_atoms.get_positions()), apply_constraint = False)
    ############### correct chemical symbols
    new_atoms.set_chemical_symbols(species)

    return new_atoms

class Logger(object):

    """
    **** this class is totally copied from AMP's ****
    **** ref) https://amp.readthedocs.io/en/latest/ ****
    Logger that can also deliver timing information.

    Parameters
    ----------
    file : str
        File object or path to the file to write to.  Or set to None for
        a logger that does nothing.
    """
    def __init__(self, file):
        if file is None:
            self.file = None
            return
        if isinstance(file, str):
            self.filename = file
            file = open(file, 'a')
        self.file = file
        self.tics = {}

    def tic(self, label=None):
        """Start a timer.

        Parameters
        ----------
        label : str
            Label for managing multiple timers.
        """
        import time
        if self.file is None:
            return
        if label:
            self.tics[label] = time.time()
        else:
            self._tic = time.time()

    def __call__(self, message, toc=None, tic=False, no_space=None):
        """Writes message to the log file.

        Parameters
        ---------
        message : str
            Message to be written.
        toc : bool or str
            If toc=True or toc=label, it will append timing information in
            minutes to the timer.
        tic : bool or str
            If tic=True or tic=label, will start the generic timer or a timer
            associated with label. Equivalent to self.tic(label).
        """
        import time
        if self.file is None:
            return
        dt = ''
        if toc:
            if toc is True:
                tic = self._tic
            else:
                tic = self.tics[toc]
            dt = (time.time() - tic) # ssrokyz start
            dt = ' %.1f sec.' % dt # ssrokyz end
        if self.file.closed:
            self.file = open(self.filename, 'a')
        if no_space or (no_space == None and toc):
            self.file.write(nospace(message + dt) + '\n')
        else:
            self.file.write(message + dt + '\n')
        self.file.flush()
        if tic:
            if tic is True:
                self.tic()
            else:
                self.tic(label=tic)
