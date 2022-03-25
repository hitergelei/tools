#!/usr/bin/env python

from subprocess import call

# liquid part
job_list = ['1-melt', '2-quench']

from ase.io import read, write
len_atoms = len(read('init-liquid.traj'))

for job in job_list:
    call('rm -rf {}/'.format(job), shell=True)
    call('mkdir {}/'.format(job), shell=True)
    call('cp {}.in run.sh frozen_model.pb {}/'.format(job, job), shell=True)

call('cp init-liquid.traj {}/init.traj'.format(job_list[0]), shell=True)

for i in range(len(job_list)):
    call('mv {}.in input.in'.format(job_list[i]), cwd=job_list[i], shell=True)
    call('ase2poscar init.traj', cwd=job_list[i], shell=True)
    call('lmp-pos2lmp.awk vasp_images_init.traj.d/POSCAR_000000 > structure.in', cwd=job_list[i], shell=True)
    call('sh run.sh', cwd=job_list[i], shell=True)
    if i != len(job_list) -1:
        # call('ase convert -n -1 lmp-result.traj ../{}/init.traj'.format(job_list[i+1]), cwd=job_list[i], shell=True)
        call('tail -n {} out.dump > ../{}/init.dump'.format(len_atoms+9, job_list[i+1]), cwd=job_list[i], shell=True)
        call('ase convert init.dump init.traj', cwd=job_list[i+1], shell=True)

# Concatenate
backbone = read('backbone-whole.vasp')
mask = backbone.get_scaled_positions()[:,1] >= 0.5
liquid_posi = read('2-quench/out.dump', -1).get_scaled_positions()
liquid_posi[:,1] = liquid_posi[:,1] /2 + 0.5
backbone_posi = backbone.get_scaled_positions()
backbone_posi[mask] = liquid_posi[:]
backbone.set_scaled_positions(backbone_posi)
write('2-quench/init-whole.traj', backbone)

# whole part
job_list = ['3-heat', '4-cryst', '6-cool']

len_atoms = len(read('2-quench/init-whole.traj'))

for job in job_list:
    call('rm -rf {}/'.format(job), shell=True)
    call('mkdir {}/'.format(job), shell=True)
    call('cp {}.in run.sh frozen_model.pb {}/'.format(job, job), shell=True)

call('cp 2-quench/init-whole.traj {}/init.traj'.format(job_list[0]), shell=True)

for i in range(len(job_list)):
    call('mv {}.in input.in'.format(job_list[i]), cwd=job_list[i], shell=True)
    call('ase2poscar init.traj', cwd=job_list[i], shell=True)
    call('lmp-pos2lmp.awk vasp_images_init.traj.d/POSCAR_000000 > structure.in', cwd=job_list[i], shell=True)
    call('sh run.sh', cwd=job_list[i], shell=True)
    if i != len(job_list) -1:
        # call('ase convert -n -1 lmp-result.traj ../{}/init.traj'.format(job_list[i+1]), cwd=job_list[i], shell=True)
        call('tail -n {} out.dump > ../{}/init.dump'.format(len_atoms+9, job_list[i+1]), cwd=job_list[i], shell=True)
        call('ase convert init.dump init.traj', cwd=job_list[i+1], shell=True)
