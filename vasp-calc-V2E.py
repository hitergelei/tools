#!/usr/bin/env python 
import numpy as np

# params
var = 0.001
poscar = 'POSCAR'
# calc = True
calc = False

#
from ase.io import read, write
atoms = read(poscar)
cell = atoms.get_cell().copy()

from subprocess import call
if calc:
    call('rm -rf calc && mkdir calc', shell=True)
E = []
Sxx = []
for i in range(11):
    if calc:
        new_atoms = atoms.copy()
        new_cell = cell * (1+(i-5)*var)
        new_cell[2] = [0., 0., 20.]
        new_atoms.set_cell(new_cell, scale_atoms=True)
        call('mkdir calc/{}'.format(i), shell=True)
        write('calc/{}/POSCAR'.format(i), new_atoms)
        call('cp INCAR POTCAR KPOINTS calc/{}'.format(i), shell=True)
        call('mpirun -np 32 vasp_std > out'.format(i), shell=True, cwd='calc/{}'.format(i))
        print('Calculation done: calc-{}'.format(i))
    result = read('calc/{}/vasprun.xml'.format(i), -1)
    Sxx.append(result.get_stress()[0])
    E.append(result.get_potential_energy())

from matplotlib import pyplot as plt
plt.plot(E)
plt.title('variation={}'.format(var), fontsize='x-large')
plt.tick_params(axis="both",direction="in", labelsize='x-large')
plt.xlabel('Step', fontsize='x-large')
plt.ylabel('E (eV)', fontsize='x-large')
plt.grid(alpha=0.5)
plt.subplots_adjust(left=0.20, bottom=0.15, right=0.95, top=0.90)

plt.figure()
from ase import units
plt.plot(np.array(Sxx) / 1e5 / units.Pascal /1e3)
plt.title('variation={}'.format(var), fontsize='x-large')
plt.tick_params(axis="both",direction="in", labelsize='x-large')
plt.xlabel('Step', fontsize='x-large')
plt.ylabel(r'$S_{{xx}}$ (kBar)', fontsize='x-large')
plt.grid(alpha=0.5)
plt.subplots_adjust(left=0.20, bottom=0.15, right=0.95, top=0.90)
plt.show()
