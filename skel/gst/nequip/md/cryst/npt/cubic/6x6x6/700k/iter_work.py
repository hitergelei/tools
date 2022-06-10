#!/usr/bin/env python

# Multi-threading
start_ind = 0
end_ind = '' 
num_threads = 20

# > MAIN
from subprocess import call
for i in range(num_threads):
    alist_slice = '{}:{}:{}'.format(start_ind+i, end_ind, num_threads)
    call('python do-structure-anal.py {} &'.format(alist_slice), shell=True)
