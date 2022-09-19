#!/usr/bin/env python
import numpy as np
from subprocess import call
from ase.io import read, write       

def _write_lmp(lmp_pos_path, alist, ran):
    for i in range(len(ran)):
        write('{}/disp-{}.lmp'.format(lmp_pos_path, ran[i]+1), alist[i], format='lammps-data')

class VibSolver():
    def __init__(
        self,
        atoms_file,
        calc,
        displacement=0.03,
        force_cutoff=0.,
        plus_minus=False,
        set_hermitian=True,
        cp_files = [],
        ):
        """
        atoms_file (ASE readable structure)
            - Only a single structure is read. If file include multiple images, it will load the final structure only.
        calc (str)
            - Choose one from (lmp,)
        displacement (float)
            -  Displacement magnitude for harmonic force constant estimation.
        force_cutoff (float)
            - Cutoff for force magnitude. Forces with smaller magnitude than this value will be set zero.
        plus_minus (bool)
            - True: Calculate +displacement and -displacement for a given direction, and average them.
            - False: Estimate FC by only +displacement calculation.
        set_hermitian (bool)
            - True: Make dynamical matrix Hermitian.
            - False: Not Hermitian-ized.
        cp_files (list of str)
            - Files to copy into the calculation directory.
        """
        print("""
        ### Reference: 
        ### Code by Young Jae Choi @ POSTECH, Republic of Korea
        """)
        self.atoms_file = atoms_file
        from ase.io import read
        self.atoms = read(atoms_file, -1)
        self.calc = calc
        self.disp = displacement
        self.f_cut = force_cutoff
        self.pm = plus_minus
        self.hermit = set_hermitian
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

    def _calc_lmp(self, omp=True):
        """
        """
        # Produce structure files
        lmp_pos_path = '{}/poscars/lmp'.format(self.calc_path)
        call('rm -rf {}'.format(lmp_pos_path), shell=True)
        call('mkdir {}'.format(lmp_pos_path), shell=True)
        alist = read('{}/poscars/poscars.traj'.format(self.calc_path), ':')

        if omp:
            from multiprocessing import Pool
            from os import getenv
            nproc = int(getenv('OMP_NUM_THREADS'))
            pool = Pool(nproc)
            tasks = [pool.apply_async(_write_lmp, (lmp_pos_path, alist[i:len(alist):nproc], range(i, len(alist), nproc))) for i in range(nproc)]
            for t in tasks:
                t.get()
            
        else:
            _write_lmp(lmp_pos_path, alist, range(len(alist)))

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
            forces = (forces[1::2] - forces[0::2]) /2.
        if self.f_cut not in (0., False, None):
            forces = np.where(np.expand_dims(np.linalg.norm(forces, axis=2) < self.f_cut, axis=2), 0, forces)
        forces = forces.reshape((len(self.atoms), 3, len(self.atoms), 3))
        fc = -forces /self.disp
        return fc

    def get_dynamical_matrix(self):
        fc = self.get_force_constants()
        m_vec = self.atoms.get_masses()
        m_mat = np.sqrt(np.outer(m_vec, m_vec))
        dm = fc /np.expand_dims(m_mat, (1,3))

        dm = dm.reshape(3*len(self.atoms), 3*len(self.atoms))
        return dm

    def get_eigen_sets(
        self,
        ):
        try:
            w2 = np.load('w2-{}-H{}.npy'.format(self.job_name, self.hermit))
            eps = np.load('eps-{}-H{}.npy'.format(self.job_name, self.hermit))
        except:
            print('w2 and eps files are NOT read.', flush=True)
            do_calc = True
        else:
            print('w2 and eps files are read.', flush=True)
            do_calc = False

        if do_calc:
            dm = self.get_dynamical_matrix()
            if self.hermit:
                w2, eps = np.linalg.eigh((dm + dm.T)/ 2.)
            else:
                w2, eps = np.linalg.eig(dm)
            eps = eps.T
            np.save('w2-{}-H{}.npy'.format(self.job_name, self.hermit), w2)
            np.save('eps-{}-H{}.npy'.format(self.job_name, self.hermit), eps)
        return w2, eps

    def plot_VDOS(
        self,
        nbins=200,
        nan_to_zero=True,
        plot_imaginary=True,
        gsmear_std=0.,
        ):
        """
        gsmear_std (float)
            - Sigma (std) for Gaussian smearing.
        """
        
        w2, eps = self.get_eigen_sets()
        w = np.sqrt(w2)

        # ASE unit to THz (below three scalings are equivalent).
        # from phonopy.units import VaspToTHz
        # w *= (VaspToTHz)
        from ase.units import fs
        w *= (fs*1e+3/2/np.pi)
        # from ase import units
        # w *= np.sqrt(units._e *1e20 /units._amu /1e24 /(2*np.pi)**2)

        if nan_to_zero:
            w = np.where(np.isnan(w), 0, w)
        eigval_plot = np.real(w)
        if plot_imaginary:
            eigval_plot -= np.where(np.imag(w) < 0, np.imag(w), 0)

        hist, edges = np.histogram(eigval_plot, bins=nbins, density=True)
        mids = (edges[0:-1] + edges[1:]) /2.
        dw = mids[1] - mids[0]
        if gsmear_std not in (0., False, None):
            print(' Gaussian smearing...', flush=True)
            from scipy.ndimage.filters import gaussian_filter1d
            hist = gaussian_filter1d(hist, gsmear_std /dw)

        from matplotlib import pyplot as plt
        for i in range(2):
            plt.figure()
            if i==0:
                plt.plot(hist, mids, c='k')
                plt.xlabel('$g(\omega)$', fontsize='x-large')
            else:
                plt.plot(hist/mids**2, mids, c='k')
                plt.xlabel('$g(\omega)/\omega^2$', fontsize='x-large')
            plt.ylabel('Frequency (THz)', fontsize='x-large')
            plt.subplots_adjust(left=0.40, bottom=0.25, right=0.60, top=0.752, wspace=0.2, hspace=0.2)
            plt.tick_params(axis="both",direction="in", labelsize='x-large')
            plt.grid(alpha=0.5)
        plt.show()
