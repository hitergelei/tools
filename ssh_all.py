#!/usr/bin/env python

#
from sys import argv
if len(argv) != 2:
    raise ValueError('Wrong number of variables')

#
from os import environ
from subprocess import call
if environ['HOSTNAME'] == 'pcs_gpu1':
    for i in range(2,12):
        call('ssh pcs_gpu{} "{}"'.format(i, argv[1]), shell=True)
