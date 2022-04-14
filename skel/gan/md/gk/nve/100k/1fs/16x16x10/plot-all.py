#!/usr/bin/env python
from subprocess import call

for i in range(1):
    call('plot-j0jt.py 100 0.01 &', shell=True, cwd='job-{}'.format(i))
