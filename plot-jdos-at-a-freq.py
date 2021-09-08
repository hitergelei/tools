#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Plot JDOS at a frequency point.
    """)
    # Positional arguments
    parser.add_argument('mesh', type=int, nargs=3, help='Provide three q-mesh numbers along x, y, and z.')
    parser.add_argument('freq_point', type=int, help='Provide the index of the freqeuncy point.')
    # # Optional arguments

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
    print('Plot JDOS at a frequency point.'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()
    mesh = args.mesh
    w_p = args.freq_point

    ir_points = np.load('ir_grid_points-{}{}{}.npy'.format(*mesh))
    ir_weights = np.load('ir_grid_weights-{}{}{}.npy'.format(*mesh))

    jdos1 = []
    jdos2 = []
    for g in ir_points:
        dat = np.loadtxt('jdos-m{}{}{}-g{}.dat'.format(*mesh, g))
        if g == 0:
            dat = np.zeros([201,3])
        w = dat[:, 0]
        jdos1.append(dat[:, 1])
        jdos2.append(dat[:, 2])
    # jdos.shape = (2, len(ir_points), len(w))
    jdos = np.array([jdos1, jdos2])

    # Get one Freq. point
    # jdos.shape = (2, len(ir_points)-1)
    jdos = jdos[:, 1:, w_p]
    
    from matplotlib import pyplot as plt
    plt.title(r'q=({}x{}x{}), $\omega$#{}'.format(*mesh, w_p), fontsize='x-large')
    plt.plot(range(len(ir_points)-1), jdos[0], label='Diff')
    plt.plot(range(len(ir_points)-1), jdos[1], label='Summ')
    plt.xlabel('Grid point index', fontsize='x-large')
    plt.ylabel(r'JDOS (THz$^{-1}$)', fontsize='x-large')
    plt.legend(fontsize='large').set_draggable(True)
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.grid(alpha=0.5)
    plt.show()
