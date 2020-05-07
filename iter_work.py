#!/usr/bin/env python

# > PARAMS
inp_file  = 'input.in'
step_init  = 40000 # steps
step_multi = 8
num_iter = 4 # step_final will be step_init * step_multi**(num_iter-1)

# > MAIN
from subprocess import call

# Input
lines = [
    'units           metal\n',
    'atom_style      atomic\n',
    'boundary        p p p\n',
    '\n',
    'box             tilt large\n',
    'read_data       structure.in\n',
    '\n',
    'mass            1 72.64\n',
    'mass            2 121.76\n',
    'mass            3 127.60\n',
    '\n',
    'pair_style      deepmd frozen_model.pb\n',
    'pair_coeff      \n',
    '\n',
    'compute         1 all pressure thermo_temp\n',
    '\n',
    'read_dump       init.dump 20000 x y z vx vy vz\n',
    'fix             1 all npt temp 300 900 1.0 iso 0 0 1.0\n',
    'thermo_style    multi\n',
    'dump            1 all custom 100 out.dump id element mass type x y z fx fy fz vx vy vz\n',
    'dump_modify     1 element Ge Sb Te\n',
    '\n',
    'timestep        0.01\n',
    'thermo          100\n',
    # 'run             ',
    ]

# Do calc
for i in range(num_iter):
    step = step_init * step_multi**i
    call('rm -rf {} && mkdir {}'.format(step, step), shell=True)
    call('cp frozen_model.pb init.dump structure.in run.sh {}/'.format(step), shell=True)
    # Write input file
    new_lines = lines[:17]
    # new_lines += ['fix             1 all npt temp 300 900 {} iso 0 0 1.0\n'.format(step/5000*1.0)]
    new_lines.extend(lines[17:])
    new_lines += ['run             {}'.format(step)]
    with open('{}/{}'.format(step, inp_file), 'w') as f:
        for j in range(len(new_lines)):
            f.write(new_lines[j])

    call('sh run.sh', cwd='./{}'.format(step), shell=True)
