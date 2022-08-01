#!/usr/bin/env python
##### CODE BY YOUNG JAE CHOI #####
import numpy as np

kB = 8.61733034e-5 # [eV/K] Boltzmann

def pick_folder_from_path(path):
    """
    path (str)
    """
    path_split = path.split('/')[:-1]
    folder = ''
    for p in path_split:
        folder += p + '/'
    return folder

def screen_species(alist, symbols):
    new_alist = []
    for i in range(len(alist)):
        chem = np.array(alist[i].get_chemical_symbols())
        mask = np.full(len(chem), False)
        for sym in symbols:
            mask += (chem == sym)
        new_alist.append(alist[i][mask])
    return(new_alist)

class bcolors:
    header = '\033[95m'
    okblue = '\033[94m'
    okcyan = '\033[96m'
    okgreen = '\033[92m'
    warning = '\033[93m'
    fail = '\033[91m'
    endc = '\033[0m'
    bold = '\033[1m'
    underline = '\033[4m'

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

def slice2str(s):
    string = ''
    if s.start is not None:
        string += str(s.start)
    string += ':'
    if s.stop is not None:
        string += str(s.stop)
    string += ':'
    if s.step is not None:
        string += str(s.step)
    return string

def parse_slice(s):
    # Slice format
    if ':' in s:
        a = [int(e) if e.strip() else None for e in s.split(":")]
    # Int format
    else:
        if int(s) == -1:
            a = [-1, None]
        else:
            a = [int(s), int(s)+1]
    return slice(*a)

def get_number_of_lines(f_obj):
    f_obj.seek(0)
    for i, l in enumerate(f_obj):
        pass
    f_obj.seek(0)
    return i + 1

def str_slice_to_list(
    str_slice,
    obj_length=None,
    ):
    """
    str_slice (str)  = String type slice. e.g. 3:-1:2
    obj_length (int) = Length of object that will use this slice object. Set to None if unknown.
    """

    #
    str_slice = str(str_slice)
    if obj_length:
        obj_length = int(obj_length)

    ## Main
    if ':' in str_slice:
        slice_list = [int(e) if e.strip() else None for e in str_slice.split(":")]
        if len(slice_list) == 2:
            slice_list.append(None)
        elif len(slice_list) != 3:
            raise ValueError('String slice option is unreadable. --> {}'.format(str_slice))
    elif int(str_slice) == -1:
        slice_list = [-1, None, None]
    else:
        slice_list = [int(str_slice), int(str_slice)+1, None]

    ## Post-process : To achieve unity.
    if obj_length:
        slice_list = slice(*slice_list).indices(obj_length)
    else:
        # 0
        if slice_list[0] == None:
            slice_list[0] = 0
        # 2
        if slice_list[2] == None:
            slice_list[2] = 1
    return slice_list

def rectify_curve(
    curve,
    rectify_cut,
    ):
    iter = True
    while True:
        test = curve[1:] - curve[:-1]
        peak_bool = np.array(list(test[:,1] > (-1 * rectify_cut)) + [True], dtype=np.bool)
        if False not in peak_bool:
            break
        curve = curve[peak_bool]
    return curve

def get_chem_ind_arr(
    chemical_symbols,
    ):
    chem_arr = np.array(chemical_symbols)
    unique_chem = np.unique(chem_arr)
    ind = np.arange(len(chem_arr))
    chem_ind_arr = []
    for chem in unique_chem:
        chem_ind_arr.append(ind[chem_arr==chem])
    return unique_chem, np.array(chem_ind_arr)

def read_txt_matrix(
    file_name,
    start_line=0,
    end_line=-1,
    ):
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
    """ Return potential energy array """
    energies = []
    for atoms in alist:
        energies.append(atoms.get_potential_energy())
    return energies
def Eperatom_fromAlist(alist):
    """ Return potential energy per atom array """
    energies = []
    for atoms in alist:
        energies.append(atoms.get_potential_energy()/len(atoms))
    return energies
def F_fromAlist(alist):
    """ Return force array """
    forces = []
    for atoms in alist:
        forces.append(atoms.get_forces())
    return forces
def v_fromAlist(alist):
    """ Return velocity array """
    velocities = []
    for atoms in alist:
        velocities.append(atoms.get_velocities())
    return velocities
def P_fromAlist(alist):
    """ Return momentum array """
    momenta = []
    for atoms in alist:
        momenta.append(atoms.get_momenta())
    return momenta
def S_fromAlist(alist):
    """ Return stress array (xx, yy, zz, yz, zx, xy) array """
    stresses = []
    for atoms in alist:
        stresses.append(atoms.get_stress())
    return stresses
def AtomicK_fromAlist(alist):
    """ Return atomic kinetic energy array """
    kinetic_e=[]
    for atoms in alist:
        kinetic_e.append(0.5 * atoms.get_masses() * np.square(np.linalg.norm(atoms.get_velocities(), axis=-1)))
    return kinetic_e
def AtomicP_fromAlist(alist):
    """ Return atomic potential energy array """
    atomic_e=[]
    for atoms in alist:
        atomic_e.append(atoms.get_potential_energies())
    return atomic_e

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
    """
    input: list.
    output: dict of counts of elements.
    """
    keys = list(set(list_inp))
    dict_out = dict()
    for i in keys:
        dict_out[i] = 0
    for i in list_inp:
        dict_out[i] += 1
    return dict_out

def list2numlist(list_inp):
    """
    input: list of any element types.
    output: list of integers. Numbers are assigned sequently to elements of first appearing.
    ex)
    inp: ['b', 'a', 'a', 'd', 'a', 'c']
    out: [ 0 ,  1 ,  1 ,  2 ,  1 ,  3 ]
    """
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

def ind_dict(list_inp):
    items = np.unique(list_inp)
    ind_dict = {}
    for item in items:
        ind_dict[item] = np.where(np.array(list_inp) == item)[0]
    return ind_dict

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
        num_spec_dict = list2count(input)
    elif isinstance(input, dict):
        num_spec_dict = input
    else:
        raise TypeError("input is not list nor dict")
    
    tot_num = np.sum(list(num_spec_dict.values()))
    r_sum = 0
    for key, value in num_spec_dict.items():
        r_sum += covalent_radii[s2n([key])[0]] * value
    expect_value = r_sum / tot_num
    return expect_value

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
            dt = ' {:.1f} sec.'.format(dt) # ssrokyz end
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
