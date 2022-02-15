import numpy as np
from ss_util import Logger

def get_new_axis(axis1, axis2):
    """
    Return three dimensional axis vectors with unit size.
    Axis1 will be new x1 coord.
    """
    axis3 = np.cross(axis1, axis2)
    axis1 = np.array(axis1) / np.linalg.norm(axis1)
    axis3 = np.array(axis3) / np.linalg.norm(axis3)
    axis2 = np.cross(axis3, axis1)
    return np.array([axis1, axis2, axis3], dtype='float')

class vector(object):
    """ 
    cutoff_num (list)       : List of maximum number of each atom-species used for describe local environment. Must have same order as type file numbering. Note) Number include itself.
    multipole_order (list)  : List that specifies which descriptor basis you will use. If 'r' is in list, (1/r) will be used as a basis. If zero or positive integers are provided, (x_i/r^n) (x_i = x,y,z) will be used as basis
        
    """
    def __init__(
        self,
        cutoff_radi,
        cutoff_num,
        multipole_order = [2,3,4],
        one_box         = True,
        logfile_name    = 'log.txt',
        ):

        ###### defining self params
        self.log             = Logger(logfile_name)
        self.cutoff_radi     = np.array(cutoff_radi, dtype = 'float64')
        self.cutoff_num      = cutoff_num
        self.multipole_order = multipole_order
        self.one_box         = one_box

    def read_npy(self, npy_path):
        self.box         = np.load(npy_path+'/box.npy')
        self.img_num     = len(self.box)
        self.coord       = np.load(npy_path+'/coord.npy')
        self.atom_num    = self.coord.shape[1]
        type_file        = open(npy_path+'/type.txt')
        type_line        = type_file.readline()
        self.types       = np.array(type_line.split(), dtype='int32')
        type_chem_line   = type_file.readline()
        self.types_chem  = type_line.split()
        
    def read_alist(self, alist):
        from traj2npy import alist2numpy as a2n
        box, coord, energy, force, stress = a2n(alist)
        
        self.box         = box
        self.img_num     = len(self.box)
        self.coord       = coord
        self.atom_num    = self.coord.shape[1]
        self.types_chem  = alist[0].get_chemical_symbols()
        if alist[-1].get_chemical_symbols() != self.types_chem:
            raise ValueError("Chemical symbols seem to be not consistent btw images. Please check")
        from ss_util import list2numlist as l2nl
        self.types       = np.array(l2nl(list(self.types_chem)), dtype='int32')

    def make_x3_supercell_coord(self):
        self.log('making x3 supercell...', tic='make_super')
        x3_supercell_set = []
        for i in range(self.img_num):
            x3_supercell = []
            for G_z in [-1,0,1]:
                for G_y in [-1,0,1]:
                    for G_x in [-1,0,1]:
                        G = np.array([[G_x], [G_y], [G_z]], dtype='int32')
                        d_r = self.box[i] * G
                        d_xyz = np.array([np.sum(d_r[:,0]),np.sum(d_r[:,1]),np.sum(d_r[:,2])], dtype='float64')
                        x3_supercell.append(self.coord[i] + d_xyz)
            x3_supercell_set.append(np.array(x3_supercell, dtype='float64'))
        self.log('x3 supercell made...', toc='make_super')
        return np.reshape(x3_supercell_set, (self.img_num,-1,3))
                
    def get_maximum_cutoff(self):
        self.log('getting maximum cutoff...', tic='max_cutoff')
        cutoff_list = []
        for i in range(self.img_num if not self.one_box else 1):
            lattice = self.box[i]
            vol = np.linalg.det(lattice)
            d1 = vol / np.linalg.norm(np.cross(lattice[1],lattice[2]))
            d2 = vol / np.linalg.norm(np.cross(lattice[2],lattice[0]))
            d3 = vol / np.linalg.norm(np.cross(lattice[0],lattice[1]))
            cutoff_list.append(np.amin([d1, d2, d3]))
        self.log('got maximum cutoff...', toc='max_cutoff')
        return np.amin(cutoff_list)

    def get_multipole_fgpt(self, r_vec, r, n_list):
        n_nlist = n_list[:]
        if 'r' in n_list:
            n_nlist.remove('r')
            return np.concatenate([1./r]+[r_vec / r**n for n in n_nlist], axis=-1)
        else:
            return np.concatenate([r_vec / r**n for n in n_nlist], axis=-1)

    # def get_multipole_fgpt_deriv(self, r_vec, r, n_list):
        # r_true = False
        # n_nlist = n_list[:]
        # if 'r' in n_list:
            # n_nlist.remove('r')
            # tmp_r = np.expand_dims(r_vec * (-1.) / r**3, axis=-1)
            # r_true = True
        # tmp = np.concatenate([np.repeat([np.eye(3, dtype='float64')], len(r), axis=0) / np.expand_dims(r**n, axis=2) - 
            # np.expand_dims(n / r**(n+2), axis=-1) * np.matmul(np.expand_dims(r_vec, axis=2), np.expand_dims(r_vec, axis=1)) for n in n_nlist], axis=-1)
        # if r_true:
            # tmp = np.concatenate([tmp_r, tmp], axis=-1)
        # return np.array([tmp[:,i,:] for i in range(3)])

    def gen_fgpts(self, inp_type, path, rotational_variation=True):
        ######## make some values
        if inp_type == 'npy':
            self.read_npy(path)
        elif inp_type == 'alist':
            self.read_alist(path)

        self.type_unique, self.type_num = np.unique(self.types, return_counts=True)
        self.x3_coord = self.make_x3_supercell_coord()
        max_cut = self.get_maximum_cutoff()
        self.log(' >>> Requirement) cutoff_radi < (maximum cutoff) <<<')
        self.log('###################################################')
        self.log('*** cutoff_radi = '+str(self.cutoff_radi))
        self.log('*** (maximum cutoff) = '+str(max_cut))
        self.log('###################################################\n')
        if self.cutoff_radi > max_cut:
            raise ValueError('cutoff radius is too big compared to cell size'
                '\n*** cutoff_radi = '+str(self.cutoff_radi)+' '
                '\n*** (maximum cutoff) = '+str(max_cut))

        self.log('getting fgpts...', tic='gen_fgpts')
        ### get distance vector
        fgpt_vectors = []
        # fgpt_deriv   = []
        Euler_angles = []
        max_radi     = []
        if rotational_variation:
            from euler_rotation import get_Euler_angles as gEa, Euler_rotation as Er
        #### fingerprint for all images
        for i in range(self.img_num):
            rel_x3_coord_i = []
            # fgpt_deriv_i   = []
            Euler_angles_i = []
            #### fingerprint for all atoms in an image
            for origin_atom in range(self.atom_num):
                #### concatenate
                rel_x3_coord_tmp = self.x3_coord[i] - self.x3_coord[i][13*self.atom_num + origin_atom]
                dist_vec = np.array(np.expand_dims(np.linalg.norm(rel_x3_coord_tmp, axis=1), axis=1), dtype='float64')
                Rx3CT_cctnted = np.concatenate((rel_x3_coord_tmp, dist_vec, np.expand_dims(
                    np.array(np.tile(self.types, 27), dtype='float64'),axis=1)), axis=1)
                # Argument sorting (ascending)
                Rx3CT_cctnted = Rx3CT_cctnted[Rx3CT_cctnted[:,3].argsort()]
                #### Rotational invariance
                if rotational_variation:
                    axis1 = np.random.rand(3)-0.5
                    axis2 = np.random.rand(3)-0.5
                    new_axis = get_new_axis(axis1, axis2)
                    Euler_angles_i.append(gEa(new_axis))
                    Rx3CT_cctnted[:,:3] = Er(Rx3CT_cctnted[:,:3],Euler_angles_i[-1])
                #### fingerprint for one (origin) atom
                origin_atom_spec = self.types[origin_atom]
                Rx3Ci_jAtom = []
                # fgpt_deriv_ij = []
                for spec in self.type_unique:
                    tmp_arr4one = Rx3CT_cctnted[Rx3CT_cctnted[:,4] == spec]
                    if origin_atom_spec == spec:
                        tmp_arr4one = tmp_arr4one[range(1,self.cutoff_num[spec]+1)]
                    else:
                        tmp_arr4one = tmp_arr4one[range(self.cutoff_num[spec])]
                    max_radi.append(tmp_arr4one[-1,3])
                    #### making multipole fgpt & it's derivative
                    Rx3Ci_jAtom.append(self.get_multipole_fgpt(tmp_arr4one[:,0:3], np.expand_dims(tmp_arr4one[:,3], axis=1), self.multipole_order))
                    # fgpt_deriv_ij.append(self.get_multipole_fgpt_deriv(tmp_arr4one[:,0:3], np.expand_dims(tmp_arr4one[:,3], axis=1), self.multipole_order))
                Rx3Ci_jAtom   = np.concatenate(tuple(Rx3Ci_jAtom  ), axis = 0)
                # fgpt_deriv_ij = np.concatenate(tuple(fgpt_deriv_ij), axis = 1)
                rel_x3_coord_i.append(Rx3Ci_jAtom)
                # fgpt_deriv_i.append(fgpt_deriv_ij)
            fgpt_vectors.append(rel_x3_coord_i)
            if rotational_variation:
                Euler_angles.append(Euler_angles_i)
            # fgpt_deriv.append(fgpt_deriv_i)
        del(self.x3_coord)
        min_max_radi = np.amin(max_radi)
        self.log('got fgpts...', toc='gen_fgpts')
        self.log(' >>> Requirement) min(maximum radius) > cutoff_radi <<<')
        self.log('####################################################')
        self.log('*** min(maximum radius) (== '+str(min_max_radi)+')')
        self.log('*** cutoff_radi == ('+str(self.cutoff_radi)+')')
        self.log('####################################################\n')
        if min_max_radi < self.cutoff_radi:
            raise ValueError('There is some local environment that cutoff_num is not sufficiently high about cutoff_radi you gave'
                '\n i.e. min(maximum radius) (== '+str(min_max_radi)+') < cutoff_radi == ('+str(self.cutoff_radi)+')')
        if rotational_variation:
            return np.array(fgpt_vectors, dtype='float64'), np.array(Euler_angles)#, np.array(fgpt_deriv, dtype='float64')
        else:
            return np.array(fgpt_vectors, dtype='float64')#, np.array(fgpt_deriv, dtype='float64')

