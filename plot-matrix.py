#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    This code plots matrix with colors.
    """)
    # Positional arguments
    parser.add_argument('mat_npy', type=str, help='Matrix saved as npy format.')
    # Optional arguments
    parser.add_argument('-c', '--cmin', type=float, help='Minimum value for color plot.')
    parser.add_argument('-d', '--cmax', type=float, help='Maximum value for color plot.')
    parser.add_argument('-e', '--hermit', action='store_true', help='Hemitianize the matrix')
    return parser.parse_args()

def plot_matrix(
    mat,
    cmin=None,
    cmax=None,
    ):
    from matplotlib import pyplot as plt
    for i in range(2):
        plt.figure()
        if i == 0:
            img = plt.imshow(np.real(mat), cmap='hot', vmin=cmin, vmax=cmax)
            plt.title('Real part', fontsize='x-large')
        else:
            img = plt.imshow(np.imag(mat), cmap='hot', vmin=cmin, vmax=cmax)
            plt.title('Imaginary part', fontsize='x-large')
        plt.colorbar(img, cmap='hot')
        plt.xlim(0, None)
        plt.ylim(None, 0)
    plt.show()

if __name__ == '__main__':
    ## Intro
    import datetime
    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ POSTECH <<<<<'.center(120))
    print(('Code runtime : '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('This code plots matrix with colors.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    mat = np.load(args.mat_npy)
    if args.hermit:
        mat = (mat.T + np.conj(mat)) /2.
    plot_matrix(mat, args.cmin, args.cmax)


