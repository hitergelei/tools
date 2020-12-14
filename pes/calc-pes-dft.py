#!/usr/bin/env python
import numpy as np

# Params
disp_min = -0.5
disp_max = 0.5
disp_nbins = 8
direc = [1, 1, 1]

# Main
dx = (disp_max-disp_min) /disp_nbins
disp = np.array(range(disp_nbins+1)) /(disp_nbins) *(disp_max-disp_min) +disp_min
direc /= np.linalg.norm(direc)
from ase.io import read, write
atoms = read('supercell_2x2x2_si-conv.traj')

from subprocess import call
call('rm -rf calc-dft', shell=True)
call('mkdir calc-dft', shell=True)
for i in range(len(disp)):
    folder = 'calc-dft/{}_{}'.format(i, disp[i])
    call('mkdir {}'.format(folder), shell=True)
    new_atoms = atoms.copy()
    posi = new_atoms.get_positions()
    posi[0] += disp[i] *direc
    new_atoms.set_positions(posi)
    write('{}/POSCAR'.format(folder), new_atoms)
    call('cp INCAR POTCAR {}'.format(folder), shell=True)
    call('mpiexec.hydra -np $OMP_NUM_THREADS vasp_std > out', shell=True, cwd=folder)

    



