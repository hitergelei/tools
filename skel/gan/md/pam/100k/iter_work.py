#!/usr/bin/env python

# > PARAMS
inp_file   = 'input-md.in'
iter_range = range( 0, 1,1)
T          = 100
dt         = 0.001
steps      = 100000

# > MAIN
from subprocess import call
from os import environ

# Do calc
for i in iter_range:
    # Input
    lines = """
        # GaN NEMD
        
        variable        T equal {}
        variable        dt equal {}
        variable        seed equal {}+1
        variable        steps equal {}
        """.format(T, dt, i, steps) + """
        units           metal
        atom_style      atomic
        boundary        f f f
        
        newton          off
        
        # box             tilt large
        read_data       structure.in
        
        mass            1 69.723
        mass            2 14.0067
        
        pair_style      nequip
        pair_coeff      * * ../../../../frozen_model.pb Ga N
        
        variable        T2 equal ${T}*2
        # velocity        all create ${T2} ${seed} mom yes dist gaussian
        compute         pe all pe/atom
        # compute         2 all pressure thermo_temp
        
        # read_dump       init.dump 5000 x y z vx vy vz box yes replace yes
        region          source block   49.50  50.50  INF INF  INF INF
        region          sink   block  125.44 126.44  INF INF  INF INF
        group           source region source
        group           sink   region sink
        
        # Equilibrium

        # fix             1 all nve
        fix             1 all nvt temp ${T} ${T} 0.01
        # fix             1 all npt temp 200 200 10. aniso 0 0 1.
        # fix             1 all box/relax aniso 0.0 vmax 1e-4

        thermo_style    custom step temp etotal pe ke
        thermo_modify   format float %.15g
        thermo          100
        log             pre-log.lammps
        
        dump            1 all custom 100 pre-out.dump id element mass type x y z fx fy fz vx vy vz c_pe
        dump_modify     1 element Ga N
        
        timestep        ${dt}
        run             5000

        # Non-equilibrium

        unfix           1
        fix             1 all nve
        fix             2 sink   heat 1 -0.5
        fix             3 source heat 1  0.5
        
        thermo_style    custom step temp etotal pe ke
        thermo_modify   format float %.15g
        thermo          100
        log             log.lammps

        undump          1
        dump            1 all custom 100 out.dump id element mass type x y z fx fy fz vx vy vz c_pe
        dump_modify     1 element Ga N
        run             ${steps} 
        """

    # 
    call('rm -rf job-{} && mkdir job-{}'.format(i, i), shell=True)
    call('cp run.sh structure.in job-{}/'.format(i), shell=True)
    # Write input file
    with open('job-{}/{}'.format(i, inp_file), 'w') as f:
        f.write(lines)

    call('sh run.sh', cwd='./job-{}'.format(i), shell=True)
