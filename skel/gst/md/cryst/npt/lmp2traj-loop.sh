#!/usr/bin/env sh

lmp2traj.py -n=:20000
lmp2traj.py -n=20000:40000
lmp2traj.py -n=40000:60000
lmp2traj.py -n=60000:80000
lmp2traj.py -n=80000:
