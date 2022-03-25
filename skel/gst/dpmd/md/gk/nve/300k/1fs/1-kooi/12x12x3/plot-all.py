#!/usr/bin/env python
from subprocess import call

for i in range(1):
    call('plot-j0jt.py 300 0.001 &', shell=True, cwd='job-{}'.format(i))
