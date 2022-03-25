
from subprocess import call
for i in range(30):
    call('vel2dos.py lmp-results.traj 0.01 -p &', shell=True, cwd='job-{}'.format(i))
