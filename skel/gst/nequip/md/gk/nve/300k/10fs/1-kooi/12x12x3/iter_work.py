#!/usr/bin/env python

# > PARAMS
inp_file     = 'input-md.in'
iter_range   = range( 0, 1,1)
temp         = 300
dt           = 0.01
sample_intvl = 1
corr_len     = 1000000
simul_len    = 1000000

# > MAIN
from subprocess import call
from os import environ

# Do calc
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
        'dimension    3\n',
        'newton       off\n',
        'boundary     p p p\n',
        # 'box          tilt large\n',
        'read_data    structure.in\n',
        'mass         1 72.64\n',
        'mass         2 121.76\n',
        'mass         3 127.60\n',
        'pair_style	  nequip\n',
        'pair_coeff   * * ../../../../../../../../frozen_model.pb Ge Sb Te\n',
        'timestep     0.01\n',
        'thermo_style custom step temp etotal ke pe press pxx pyy pzz pyz pxz pxy vol\n',
        'thermo_modify   format float %.15g\n',
        'thermo       1000\n',
        'log          pre-log.lammps\n',
        '\n',
        '# equilibration and thermalization\n',
        '\n',
        'velocity     all create ${T2}'+' {} mom yes dist gaussian\n'.format(i+1),
        'fix          NPT all npt temp ${T} ${T} 1. aniso 0 0 100.0\n',
        'compute      pe all pe/atom\n',
        'dump         0 all custom 1000 pre-out.dump id element mass type x y z fx fy fz vx vy vz c_pe\n',
        'dump_modify  0 element Ge Sb Te\n',
        'run          5000\n',
        'undump       0\n',
        # '\n',
        # 'unfix        NVT\n',
        # 'fix          NVT all nvt temp $T $T 1.\n',
        # 'run          {}\n'.format(int(50/dt)),
        '\n',
        '# thermal conductivity calculation, switch to NVE if desired\n',
        '\n',
        'unfix       NPT\n',
        # 'fix         NPT all npt temp ${T} ${T} 100. aniso 0 0 100.0\n',
        'fix         NVE all nve\n',
        # '\n',
        'reset_timestep 0\n',
        'timestep     ${dt}\n',
        'compute      myKE all ke/atom\n',
        'compute      myPE all pe/atom\n',
        'compute      myStress all stress/atom NULL virial\n',
        'compute      flux all heat/flux myKE myPE myStress\n',
        'fix          JJ all ave/correlate $s $p $d &\n',
        '             c_flux[1] c_flux[2] c_flux[3] type auto file J0Jt.dat ave running\n',
        'thermo_style custom step temp etotal ke pe press pxx pyy pzz pyz pxz pxy vol c_flux[1] c_flux[2] c_flux[3]\n',
        'thermo_modify   format float %.15g\n',
        'thermo       1000\n',
        'log          log.lammps\n',
        'dump         1 all custom 1000 out.dump id element mass type x y z fx fy fz vx vy vz c_pe\n',
        'dump_modify  1 element Ge Sb Te\n',
        'run          $d\n',
        ]

    # 
    call('rm -rf job-{} && mkdir job-{}'.format(i, i), shell=True)
    call('cp run.sh structure.in job-{}/'.format(i), shell=True)
    # Write input file
    with open('job-{}/{}'.format(i, inp_file), 'w') as f:
        for j in range(len(lines)):
            f.write(lines[j])

    call('sh run.sh', cwd='./job-{}'.format(i), shell=True)
