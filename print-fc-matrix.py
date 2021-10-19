#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Print force constant matrix in float form.
    """)
    # Positional arguments
    parser.add_argument('pho_pckl', type=str, help='Provide the name of a pickle file of phonopy.')
    # # Optional arguments
    parser.add_argument('-r', '--rot_angle', nargs=3, type=float, help='3 angles for Euler rotation in degree unit.')

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
    print('Print force constant matrix in float form.'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()

    # Load FCs
    import pickle as pckl
    with open(args.pho_pckl, 'rb') as f:
        pho = pckl.load(f)
    fc = pho.get_force_constants()

    if args.rot_angle is not None:
        # Rotation
        rot_ang = np.array(args.rot_angle) /180 *np.pi
        # Rot 1
        R = np.array([[ np.cos(rot_ang[0]), np.sin(rot_ang[0]),                  0],
                      [-np.sin(rot_ang[0]), np.cos(rot_ang[0]),                  0],
                      [                  0,                  0,                  1]])
        R = np.tile(R, (fc.shape[0], fc.shape[1], 1, 1))
        fc = np.matmul(
            np.matmul(R, fc),
            np.transpose(R, (0,1,3,2)),
            )
        # Rot 2
        R = np.array([[                  1,                  0,                  0],
                      [                  0, np.cos(rot_ang[1]), np.sin(rot_ang[1])],
                      [                  0,-np.sin(rot_ang[1]), np.cos(rot_ang[1])]])
        R = np.tile(R, (fc.shape[0], fc.shape[1], 1, 1))
        fc = np.matmul(
            np.matmul(R, fc),
            np.transpose(R, (0,1,3,2)),
            )
        # Rot 3
        R = np.array([[ np.cos(rot_ang[2]), np.sin(rot_ang[2]),                  0],
                      [-np.sin(rot_ang[2]), np.cos(rot_ang[2]),                  0],
                      [                  0,                  0,                  1]])
        R = np.tile(R, (fc.shape[0], fc.shape[1], 1, 1))
        fc = np.matmul(
            np.matmul(R, fc),
            np.transpose(R, (0,1,3,2)),
            )

    from ss_util import bcolors
    fc = np.transpose(fc, [0,2,1,3]).reshape(fc.shape[0]*fc.shape[-1], -1)
    line = '        '
    for i in range(len(fc)):
        if i%3 == 0:
            alpha = 'x'
        elif i%3 == 1:
            alpha = 'y'
        else:
            alpha = 'z'
        line += ' ({},{})'.format(i//3, alpha)
        if i % 3 == 2:
            line +='  '
    print(bcolors.okblue + line + bcolors.endc)
    for i in range(len(fc)):
        if i%3 == 0:
            alpha = 'x'
        elif i%3 == 1:
            alpha = 'y'
        else:
            alpha = 'z'
        line = bcolors.okblue + '  ({},{})'.format(i//3, alpha) + bcolors.endc
        for j in range(len(fc)):
            num = '{:9.4f}'.format(fc[i,j])
            if num == '  -0.0000':
                num = '   0.0000'
            if num != '   0.0000':
                num = bcolors.okgreen + num + bcolors.endc
            line += num
            if j % 3 == 2:
                line +='  '
        print(line)
        if i % 3 == 2:
            print()
    print()
