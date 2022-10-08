#!/usr/bin/env python

#
from sys import argv
if len(argv) != 2:
    raise ValueError('Wrong number of variables')

#
from os import environ
from subprocess import call
for i in range(2,13):
    call('ssh pcs_gpu{} "{}"'.format(i, argv[1]), shell=True)
for i in range(1,11):
    call('ssh pcs_cpu{} "{}"'.format(str(i).zfill(2), argv[1]), shell=True)
