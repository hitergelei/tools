#!/usr/bin/env python

# > PARAMS
prim_cell    = 'gete-alpha-prim-dpmd.vasp'
N            = 9
inp_file     = 'input-md.in'
iter_range   = range(0,30,1)
temp         = 300
dt           = 0.010
sample_intvl = 1
corr_len     = 20000
simul_len    = 100000

# > MAIN
from subprocess import call
from os import environ
environ['CUDA_VISIBLE_DEVICES'] = ''
NNN = ((N,0,0),(0,N,0),(0,0,N))

# Make cells
from phonon_distrib import get_phonon_distrib_alist
alist = get_phonon_distrib_alist(
    prim_cell,
    NNN,
    temp,
    seed_range = iter_range,
    calc       = 'lmp',
    cp_files   = ('frozen_model.pb',),
    plus_minus = True,
    )

# Do calc
from ase.io import write
for i in iter_range:
    # Input
    lines = [
        '# Sample LAMMPS input script for thermal conductivity of solid Ar\n',
        '\n',
        'units       metal\n',
        'variable    T equal {}\n'.format(temp),
        'variable    T2 equal $T*1.\n',
        'variable    dt equal {}\n'.format(dt),
        'variable    s equal {}    # sample interval\n'.format(sample_intvl),
        'variable    p equal {}    # correlation length\n'.format(corr_len),
        'variable    d equal {}    # dump interval\n'.format(simul_len),
        '\n',
        '# @ convert from LAMMPS metal units to SI\n',
        'variable    kB equal 8.61733034e-5    # [eV/K] Boltzmann\n',
        'variable    e equal 1.60218e-19\n',
        'variable    eV2J equal $e\n',
        'variable    ps2s equal 1.0e-12\n',
        'variable    A2m equal 1.0e-10\n',
        'variable    convert equal ${eV2J}/${ps2s}/${A2m} # [W/mK] = [J/smK]\n',
        '\n',
        '# setup problem\n',
        '\n',
        # 'package      omp {}\n'.format(environ['OMP_NUM_THREADS']),
        'package      omp 4\n',
        'dimension    3\n',
        'boundary     p p p\n',
        'box          tilt large\n',
        'read_data    structure.in\n',
        'mass         1 72.64\n',
        'mass         2 127.60\n',
        'pair_style   deepmd frozen_model.pb\n',
        'pair_coeff   \n',
        'timestep     ${dt}\n',
        'thermo_style custom step etotal ke temp pe press vol\n',
        'thermo       $s\n',
        '\n',
        '# equilibration and thermalization\n',
        '\n',
        'fix          NVT all nvt/omp temp $T $T 1.\n',
        'run          5000\n',
        # '\n',
        # 'unfix        NVT\n',
        # 'fix          NVT all nvt/omp temp $T $T 1.\n',
        # 'run          {}\n'.format(int(50/dt)),
        '\n',
        '# thermal conductivity calculation, switch to NVE if desired\n',
        '\n',
        # 'unfix       NVT\n',
        # 'fix         NVE all nve/omp\n',
        # '\n',
        'reset_timestep 0\n',
        'compute      myKE all ke/atom\n',
        'compute      myPE all pe/atom\n',
        'compute      1 all pe/atom\n',
        'compute      myStress all stress/atom NULL virial\n',
        'compute      flux all heat/flux myKE myPE myStress\n',
        'fix          JJ all ave/correlate $s $p $d &\n',
        '             c_flux[1] c_flux[2] c_flux[3] type auto file J0Jt.dat ave running\n',
        'thermo_style custom step etotal ke temp pe press vol c_flux[1] c_flux[2] c_flux[3]\n',
        'thermo       $s\n',
        'dump         1 all custom 100 out.dump id element mass type x y z fx fy fz vx vy vz c_1\n',
        'dump_modify  1 element Ge Te\n',
        'run          $d\n',
        ]

    # 
    call('rm -rf job-{} && mkdir job-{}'.format(i, i), shell=True)
    call('cp run.sh frozen_model.pb job-{}/'.format(i), shell=True)
    # Write input file
    write('job-{}/init.traj'.format(i), alist[i])
    write(
        'job-{}/structure.in'.format(i),
        alist[i],
        format='lammps-data',
        velocities=True,
        )
    with open('job-{}/{}'.format(i, inp_file), 'w') as f:
        for j in range(len(lines)):
            f.write(lines[j])

    call('sh run.sh', cwd='./job-{}'.format(i), shell=True)
