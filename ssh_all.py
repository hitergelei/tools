#!/usr/bin/env python

#
from sys import argv
if len(argv) != 2:
    raise ValueError('Wrong number of variables')

#
from os import environ
from subprocess import call
for i in range(2,12):
    call('ssh pcs_gpu{} "{}"'.format(i, argv[1]), shell=True)
for i in range(1,7):
    call('ssh pcs_cpu0{} "{}"'.format(i, argv[1]), shell=True)
