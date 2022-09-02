#!/usr/bin/env python
import numpy as np
from subprocess import call

class VibSolver():
    """
    ### Reference: 
    ### Code by Young Jae Choi @ POSTECH, Republic of Korea
    """

    def __init__(
        self,
        atoms_file,
        calc,
        displacement=0.03,
        plus_minus=False,
        cp_files = [],
        ):
        """
        atoms_file (ASE readable structure)
            - Only a single structure is read. If file include multiple images, it will load the final structure only.
        calc (str)
            - Choose one from (lmp,)
        displacement (float)
            -  Displacement magnitude for harmonic force constant estimation.
        plus_minus (bool)
            - True: Calculate +displacement and -displacement for a given direction, and average them.
            - False: Estimate FC by only +displacement calculation.
        cp_files (list of str)
            - Files to copy into the calculation directory.
        """
        self.atoms_file = atoms_file
        from ase.io import read
        self.atoms = read(self.atoms_file, -1)
        self.calc = calc
        self.disp = displacement
        self.pm = plus_minus
        self.cp_files = cp_files
        self.cp_files_concat = ''
        for i in range(len(self.cp_files)):
            self.cp_files_concat += self.cp_files[i] + ' '

        self.job_name = '{}-{}-d{}-pm{}'.format(calc, atoms_file, displacement, plus_minus)
        self.calc_path = './calcs/{}'.format(self.job_name)

    def _produce_structures(self):
        path='{}/poscars'.format(self.calc_path)
        call('rm -rf {}'.format(path), shell=True)
        call('mkdir -p {}'.format(path), shell=True)

        from ase.io import read, write
        posi = self.atoms.get_positions()       
        new_alist = []
        for i in range(len(self.atoms)):
            for j in range(3):
                for k in [-1, 1]:
                    if k == -1 and not self.pm:
                        pass
                    else:
                        new_posi = posi.copy()
                        new_posi[i,j] += self.disp *k
                        new_atoms = self.atoms.copy()
                        new_atoms.set_positions(new_posi)
                        new_alist.append(new_atoms)
                        # write('{}/a{}-x{}-d{}.vasp'.format(path, i, j, k), new_alist[-1])
        write('{}/poscars.traj'.format(path), new_alist)

    def _get_forces(self):
        """
        """
        try:
            forces = np.load('forces-{}.npy'.format(self.job_name))
        except:
            do_calc = True
        else:
            do_calc = False

        if do_calc:
            #
            self._produce_structures()

            #
            if self.calc == 'lmp':
                forces = self._calc_lmp()
            else:
                raise NotImplementedError('Calculator "{}" is not implemented yet.'.format(self.calc))

            #
            np.save('forces-{}.npy'.format(self.job_name), np.array(forces))

        return np.array(forces)

    def _calc_lmp(self):
        """
        """
        # Produce structure files
        lmp_pos_path = '{}/poscars/lmp'.format(self.calc_path)
        call('rm -rf {}'.format(lmp_pos_path), shell=True)
        call('mkdir {}'.format(lmp_pos_path), shell=True)
        from ase.io import read, write       
        alist = read('{}/poscars/poscars.traj'.format(self.calc_path), ':')
        for i in range(len(alist)):
            write('{}/disp-{}.lmp'.format(lmp_pos_path, i+1), alist[i], format='lammps-data')

        # Calc
        call('cp {} input-vs.in {}'.format(self.cp_files_concat, self.calc_path), shell=True)
        with open('{}/input-vs.in'.format(self.calc_path), 'a') as f:
            f.write('\nvariable        i loop {}'.format(len(alist)))
            f.write('\nlabel           loophead')
            f.write('\ndelete_atoms    group all')
            f.write('\nread_data       poscars/lmp/disp-$i.lmp add append')
            f.write('\nreset_timestep  $i')
            f.write('\nrun             0')
            f.write('\nif "$i == {}" then "jump  SELF loopout"'.format(len(alist)))
            f.write('\nnext            i')
            f.write('\njump            SELF loophead')
            f.write('\nlabel           loopout')

        call('lmp_mpi -in input-vs.in -sf intel > out', shell=True, cwd=self.calc_path)

        # Gather results
        alist = read('{}/out.dump'.format(self.calc_path), ':')
        forces = []
        for i in range(len(alist)):
            forces.append(alist[i].get_forces())
        return np.array(forces)

    def get_force_constants(self):
        """
        Returns force constants in shape (len(atoms), 3, len(atoms), 3)
        In unit of eV/Anstrom^2
        """
        forces = self._get_forces()
        if self.pm:
            forces = (forces[::2] - forces[1::2]) /2.
        forces = forces.reshape((len(self.atoms), 3, len(self.atoms), 3))
        fc = -forces /self.disp
        return fc

    def get_dynamical_matrix(self):
        fc = self.get_force_constants()
        m_vec = self.atoms.get_masses()
        m_mat = np.sqrt(np.outer(m_vec, m_vec))
        dm = fc /np.expand_dims(m_mat, (1,3))
        # ASE unit to THz
        from ase.units import fs
        dm /= fs**2
        dm = dm.reshape(3*len(self.atoms), 3*len(self.atoms))
        return dm

    def get_eigen_sets(self):
        dm = self.get_dynamical_matrix()
        # w2, eps = np.linalg.eigh((dm + dm.T)/ 2.)
        w2, eps = np.linalg.eig(dm)
        return w2, eps

    def plot_VDOS(
        self,
        nbins=200,
        ):
        w2, eps = self.get_eigen_sets()
        eigval = np.sqrt(w2)
        eigval_plot = np.real(eigval) # -np.imag(eigval)
        # eigval_plot = -np.imag(eigval)

        hist, edges = np.histogram(eigval_plot, bins=nbins)
        mids = (edges[0:-1] + edges[1:]) /2.
        from matplotlib import pyplot as plt
        plt.plot(hist, mids, c='k')
        plt.show()
