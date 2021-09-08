#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Plot weighted phonon joint DOS.
    """)
    # Positional arguments
    parser.add_argument('mesh', type=int, nargs=3, help='Provide three q-mesh numbers along x, y, and z.')
    # # Optional arguments
    parser.add_argument('-t', '--temperature', type=int, default=None, nargs='+',
        help='Set temperature manually. Multiple temperature can be set. Only int type. **IF NOT SET, calculate JDOS**')
    parser.add_argument('-g', '--grid_point', type=int, help='Plot at single grid point. (Not sum of all)')
    parser.add_argument('-l', '--plot_half', action='store_true', help='Plot only half range of frequency.')

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
    print('Plot weighted phonon joint DOS.'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()
    mesh = args.mesh
    temp = args.temperature
    if temp is None:
        temp = [None]

    ir_points = np.load('ir_grid_points-{}{}{}.npy'.format(*mesh))
    ir_weights = np.load('ir_grid_weights-{}{}{}.npy'.format(*mesh))
    if args.grid_point is not None:
        ir_points = [args.grid_point]
        ir_weights = [1.]

    jdos1 = []
    jdos2 = []
    for t in temp:
        jdos1.append([])
        jdos2.append([])
        for g in ir_points:
            if temp[0] is None:
                dat = np.loadtxt('jdos-m{}{}{}-g{}.dat'.format(*mesh, g))
                if g == 0:
                    dat[0,1]=0.
            else:
                dat = np.loadtxt('jdos-m{}{}{}-g{}-t{}.dat'.format(*mesh, g, t))
            w = dat[:, 0]
            jdos1[-1].append(dat[:, 1])
            jdos2[-1].append(dat[:, 2])
    # jdos.shape = (2, len(temp), len(ir_points), len(w))
    jdos = np.array([jdos1, jdos2])
    # jdos.shape = (len(temp), 2, len(ir_points), len(w))
    jdos = np.transpose(jdos, [1, 0, 2, 3])

    #
    for i in range(len(ir_weights)):
        jdos[:, :, i] *= ir_weights[i]
    # jdos.shape = (len(temp), 2, len(w))
    jdos = np.sum(jdos, axis=2)
    # jdos = np.amin(jdos[:,:,1:], axis=2)
    if not args.grid_point is not None:
        jdos /= mesh[0] * mesh[1] * mesh[2]

    from matplotlib import pyplot as plt
    for i in range(len(temp)):
        plt.figure()
        title = 'q=({}x{}x{})'.format(*mesh)
        if temp[0] is not None:
            title += ', {}K'.format(temp[i])
        if args.grid_point is not None:
            title += ', grid point: {}'.format(ir_points[0])
        plt.title(title, fontsize='x-large')
        plt.plot(w, jdos[i, 0], label='Diff', c='b')
        plt.plot(w, jdos[i, 1], label='Summ', c='r')
        plt.plot(w, np.sum(jdos[i], axis=0), label='Total', c='k')
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.legend(fontsize='large').set_draggable(True)
        plt.xlabel('Frequency (THz)', fontsize='x-large')
        if temp[0] is None:
            plt.ylabel(r'JDOS (THz$^{-1}$)', fontsize='x-large')
        else:
            plt.ylabel(r'Weighted JDOS (THz$^{-1}$)', fontsize='x-large')
        if args.plot_half:
            plt.xlim([0, w[-1]/2.])
        else:
            plt.xlim([0, w[-1]])
        plt.grid(alpha=0.5)
    plt.show()
