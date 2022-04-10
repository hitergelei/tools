#!/usr/bin/env python
from subprocess import call
for i in range( 0, 1,1):
    # call('lmp2traj.py -l &', shell=True, cwd='job-{}'.format(i))
    # call('unfold_positions.py lmp-results.traj &', shell=True, cwd='job-{}'.format(i))
    # call('therm_conduct.py lmp-results.traj 0.01 100 -a 10 -p -l -t &', shell=True, cwd='job-{}'.format(i))
    # call('therm_conduct_new.py lmp-results.traj 0.01 100 -a 10 -p -t &', shell=True, cwd='job-{}'.format(i))
    # call('therm_conduct_new.py unfold_lmp-results.traj 0.01 20 -a 10 -p -t &', shell=True, cwd='job-{}'.format(i))
    # call('therm_conduct_new.py unfold_lmp-results.traj 0.01 20 -a 10 -p -t -e &', shell=True, cwd='job-{}'.format(i))
    # call('einstein_relation.py unfold_lmp-results.traj 0.01 20 -a 10 &', shell=True, cwd='job-{}'.format(i))
    call('plot-j0jt.py 300 0.01 -p -l &', shell=True, cwd='job-{}'.format(i))
    # call('therm_conduct_mask.py lmp-results.traj 0.01 20 -a 10 -p &', shell=True, cwd='job-{}'.format(i))
