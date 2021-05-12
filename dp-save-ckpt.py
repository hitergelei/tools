#!/usr/bin/env python

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    This code will save checkpoint files at desired training steps.
    Only works for DeePMD.
    """)
    # Positional arguments
    parser.add_argument('save_step_intvl', type=int, help='Interval of steps needed to be saved (int type).')
    # Optional arguments
    parser.add_argument('-f', '--lcurve_file', type=str, default='lcurve.out', help='lcurve.out file name.')
    parser.add_argument('-t', '--time_intvl', type=float, default='5.', help='Time interval between step checks.')
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
    print('This code will save checkpoint files at desired training steps.'.center(120))
    print('Only works for DeePMD.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args = argparse()
    s_intvl = args.save_step_intvl
    l_file = args.lcurve_file
    t_intvl = args.time_intvl

    # Wait for calculation to start
    from time import sleep
    from os.path import isfile
    from subprocess import check_output, call
    while True:
        if isfile(l_file):
            try:
                tmp = int(str(check_output('tail -n 1 {}'.format(l_file), shell=True))[2:].split()[0])
            except:
                sleep(t_intvl)
            else:
                call("echo ' >>_dp-save-ckpt.py:_Start_saving_checkpoint_files.' >> {}".format(l_file), shell=True)
                break
        else:
            sleep(t_intvl)

    # Main
    call('mkdir bu', shell=True)
    last_save = None
    while True:
        words = str(check_output('tail -n 1 {}'.format(l_file), shell=True))[2:].split()
        if len(words) > 2 and words[0].isnumeric():
            step = int(words[0])
            if last_save is None:
                last_save = (step // s_intvl) * s_intvl
            elif last_save + s_intvl <= step:
                call('rm -rf bu/{}; mkdir bu/{}'.format(last_save + s_intvl, last_save + s_intvl), shell=True)
                sleep(1)
                call('cp * bu/{}/'.format(last_save + s_intvl), shell=True)
                call("echo ' >>_dp-save-ckpt.py:_Step_{}_saved!' >> {}".format(last_save + s_intvl, l_file), shell=True)
                last_save += s_intvl
        sleep(t_intvl)
    call("echo ' >>_dp-save-ckpt.py:_terminated.' >> {}".format(l_file), shell=True)
