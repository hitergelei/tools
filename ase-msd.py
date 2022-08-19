#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    Plot partial mean square displacement.
    NOTE!!!! YOU MUST PROVIDE UNWRAPPED IMAGES !!!!!
    """)
    # Positional arguments
    parser.add_argument('alist_file', type=str, help='ASE readable atoms list file name.')
    parser.add_argument('n_sample', type=int, help='The number of samples to get reference position. Set the resference structure as the structure at t=0, if n_sample=0')
    parser.add_argument('dt', type=float, help='dt of the provided images in ps unit. The value must change by the third part of the img_slice parameter.')
    # Optional arguments
    parser.add_argument('-n', '--img_slice', type=str, default=':', help='Image range following python convention. default=":" (e.g.) -n :1000:10')
    parser.add_argument('-g', '--gsmear', type=int, default=0, help='Positions will be Gaussian smeared with this sigma before the MSD calculation [default: No smear, 0]')
    parser.add_argument('-c', '--color_list', type=str, nargs='+', help='Color list for the plot')
    return parser.parse_args()

def get_MSD(
    alist,
    n_sample = 10,
    gsmear   = None,
    ):
    """
    MSD means mean square displacement

    n_sample (int)
        - The number of samples to get reference position
        - Must be larger than 2
    gsmear (int)
        - Positions will be Gaussian smeared before the MSD calculation.
    """
    #
    if n_sample < 3 and n_sample != 0:
        raise ValueError('n_sample must be zero or an integer larger than 2.')

    # Gather data.
    posi = []
    for i in range(len(alist)):
        posi.append(alist[i].get_positions())
    posi = np.array(posi)

    # Inner average.
    if gsmear not in [None, False, 0.]:
        from scipy.ndimage import gaussian_filter
        posi = gaussian_filter(posi, sigma=gsmear, axis=0, mode='nearest')
    
    # Main
    # MSD --> shape of (len(alist), len_atoms)
    MSD = np.zeros(posi.shape[:-1], dtype=float)
    for i in range(len(alist)):
        if n_sample == 0:
            MSD[i] = np.square(np.linalg.norm(posi[i] - posi[0], axis=-1))
        else:
            start = i - n_sample //2
            end = i + n_sample //2 + n_sample %2
            if start < 0:
                start = 0
            if end > len(alist):
                end = len(alist)
            MSD[i] = np.square(np.linalg.norm(posi[i] - np.mean(posi[start:end], axis=0), axis=-1))

    return MSD

def plot_MSD(
    alist,
    n_sample,
    dt,
    gsmear     = None,
    color_list = None,
    ):
    """
    """
    # dt to ns unit
    # dt /= 1000.

    #
    MSD = get_MSD(
        alist,
        n_sample,
        gsmear,
        )

    chem = np.array(alist[0].get_chemical_symbols())
    chem_unique = np.unique(chem)
    mask = []
    avg_MSD = []
    for i in range(len(chem_unique)):
        mask.append(chem == chem_unique[i])
        avg_MSD.append(np.mean(MSD[:, mask[i]], axis=1))
    avg_MSD = np.array(avg_MSD)
    DW = avg_MSD *8 *np.pi**2 /3

    tit_msd = ''
    tit_dw = ''
    print('    Species       Mean-square displacement (Angst.^2)        Debye-Waller factor (Angst.^2)')
    for i in range(len(chem_unique)):
        print('      {}  {:30.8f}     {:35.8f}'.format(chem_unique[i], np.mean(avg_MSD[i]), np.mean(DW[i])))
        tit_msd += '{}:{:.4f} '.format(chem_unique[i], np.mean(avg_MSD[i]))
        tit_dw += '{}:{:.2f} '.format(chem_unique[i], np.mean(DW[i]))

    # Gaussian smearing for plot
    if color_list in [None, False, []]:
        color_list = [None] *len(chem_unique)
    from matplotlib import pyplot as plt
    for j in range(2):
        plt.figure()
        for i in range(len(chem_unique)):
            if j==0:
                plt.plot(np.arange(len(avg_MSD[i]))[1:] *dt, avg_MSD[i][1:], c=color_list[i], label=chem_unique[i])
                plt.ylabel('MSD ($\AA ^2$)', fontsize='x-large')
                plt.title('NS={}, dt={}, GS={}\n{}'.format(n_sample, dt, gsmear, tit_msd), fontsize='small', pad=5)
            if j==1:
                plt.plot(np.arange(len(DW[i]))[1:] *dt, DW[i][1:], c=color_list[i], label=chem_unique[i])
                plt.ylabel('Debye-Waller factor ($\AA ^2$)', fontsize='x-large')
                plt.title('NS={}, dt={}, GS={}\n{}'.format(n_sample, dt, gsmear, tit_dw), fontsize='small', pad=5)
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.xlabel('Time (ps)', fontsize='x-large')
        plt.subplots_adjust(left=0.30, bottom=0.25, right=0.70, top=0.75, wspace=0.2, hspace=0.2)
        plt.legend(fontsize='large').set_draggable(True)
        plt.grid(alpha=0.5)
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
    print('Plot partial mean square displacement'.center(120))
    print('NOTE!!!! YOU MUST PROVIDE UNWRAPPED IMAGES !!!!!'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()
    if args.alist_file[:7] != 'unwrap_':
        raise RuntimeError('You must provide unwrapped images. File name must starts with "unwrap_".')

    from ase.io import read
    alist = read(args.alist_file, args.img_slice)
    if not isinstance(alist, list):
        alist = [alist]

    plot_MSD(
        alist,
        args.n_sample,
        args.dt,
        args.gsmear,
        args.color_list,
        )
