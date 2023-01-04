#!/usr/bin/env python

from subprocess import call

job_list = [
    '1-melt',
    '2-quench',
    '3-heat',
    '4-cryst',
    '6-cool',
    '7-vacf',
    ]

from ase.io import read
len_atoms = len(read('init.traj'))

for job in job_list:
    call('rm -rf {}/'.format(job), shell=True)
    call('mkdir {}/'.format(job), shell=True)
    call('cp {}.in run.sh {}/'.format(job, job), shell=True)

call('cp init.traj {}/'.format(job_list[0]), shell=True)

for i in range(len(job_list)):
    call('mv {}.in input.in'.format(job_list[i]), cwd=job_list[i], shell=True)
    call('ase2poscar init.traj', cwd=job_list[i], shell=True)
    call('lmp-pos2lmp.awk vasp_images_init.traj.d/POSCAR_000000 > structure.in', cwd=job_list[i], shell=True)
    call('sh run.sh', cwd=job_list[i], shell=True)
    if i != len(job_list) -1:
        # call('ase convert -n -1 lmp-result.traj ../{}/init.traj'.format(job_list[i+1]), cwd=job_list[i], shell=True)
        call('tail -n {} out.dump > ../{}/init.dump'.format(len_atoms+9, job_list[i+1]), cwd=job_list[i], shell=True)
        call('ase convert init.dump init.traj', cwd=job_list[i+1], shell=True)

