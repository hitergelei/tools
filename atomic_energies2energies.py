#!/usr/bin/env python

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    -
    """)
    # Positional arguments
    parser.add_argument('traj_file', type=str, help='ASE readable atoms list file name.')
    return parser.parse_args()

if __name__ == '__main__':
    ## Intro
    import datetime
    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ Phys. Dep. of POSTECH in Korea <<<<<'.center(120))
    print(('Code runtime : '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('-'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args = argparse()

from ase.io import read, write, Trajectory
alist = read(args.traj_file, ':')
traj = Trajectory('new-'+args.traj_file, 'w')

for atoms in alist:
    atoms._calc.results['energies'] = atoms._calc.results['atomic_energies'].copy()
    del(atoms._calc.results['atomic_energies'])
    traj.write(atoms)

from subprocess import call
call('mv {} old-{}'.format(args.traj_file, args.traj_file), shell=True)
