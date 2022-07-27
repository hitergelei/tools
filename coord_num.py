#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    DESCRIPTION
    """)
    # Positional arguments
    parser.add_argument('chem', type=str, nargs=2, help='Two chemical symbols consisting bonds. e.g., Ge Te | "a": any symbols')
    parser.add_argument('alist_file', type=str, help='ASE readable atoms list file name.')
    parser.add_argument('bond_range', type=float, nargs=2, help='Range of the bonding to count. e.g., 3.0 5.2')
    # Optional arguments
    parser.add_argument('-n', '--img_slice', type=str, default=':', help='Image range following python convention. default=":" (e.g.) -n :1000:10')
    parser.add_argument('-t', '--dt', type=float, default=None, help='Time interval for plot in ps unit. [default: None]')
    return parser.parse_args()

def get_coord_num(rdf, bond_range, num_den):
    """
    rdf (array of shape (number of distance points, 2))
        - rdf[:, 0] are distance and rdf[:, 1] are RDF value.
    bond_range (tuple of shape (2,))
        - bond_range = (minimum distance, maximum distance)
    num_den (float)
        - number density of the species of interest.
    """
    mask = bond_range[0] <= rdf[:, 0]
    mask = mask * (rdf[:, 0] < bond_range[1])
    rdf_range = rdf[mask]
    delta_d = rdf[1,0] -rdf[0,0]

    coord_num = 0.
    for i in range(len(rdf_range)):
        coord_num += 4.*np.pi*rdf_range[i,0]**2 *rdf_range[i,1] *num_den *delta_d
    return coord_num

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
    print('DESCRIPTION'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()
    chem = args.chem
    bond_range = args.bond_range
    if args.dt is not None:
        dt = args.dt /1000.
    else:
        dt = None

    from ase.io import read
    alist = read(args.alist_file, args.img_slice)
    if not isinstance(alist, list):
        alist = [alist]

    t = np.arange(len(alist)) +1
    if dt is not None:
        t = t *dt

    coord_num = []
    from ase_rdf import get_RDF
    for i in range(len(alist)):
        rdf = get_RDF([alist[i]], bond_range[1], 300, (chem[0], chem[1]), log=True)
        num_den1 = np.sum(np.array(alist[i].get_chemical_symbols()) == chem[0]) /alist[i].get_volume()
        num_den2 = np.sum(np.array(alist[i].get_chemical_symbols()) == chem[1]) /alist[i].get_volume()
        coord_num.append([
            get_coord_num(rdf, bond_range, num_den1),
            get_coord_num(rdf, bond_range, num_den2),
            ])
    print(' Coordination number of bonding range ({}, {}):'.format(*bond_range))
    print('  {}    {}'.format(chem[0], chem[1]))
    for i in range(len(coord_num)):
        print('{:.4f}, {:.4f}'.format(*coord_num[i]))
    coord_num = np.array(coord_num)
    
    #
    from scipy.stats import linregress as lr
    linfit = []
    for i in range(len(chem)):
        linfit.append(lr(np.log10(t), coord_num[:,0]))

    from matplotlib import pyplot as plt
    for i in range(len(chem)):
        plt.figure()
        plt.plot(t, coord_num[:, i], c='k')
        plt.plot(
            t,
            linfit[i].slope *np.log10(t) +linfit[i].intercept,
            c='r',
            label=r'$y=a\rm{log}(x)+b$'+'\n'+r'a={:.2f}, b={:.2f}'.format(linfit[i].slope, linfit[i].intercept),
            )
        plt.xscale('log')
        # plt.yscale('log')
        plt.title(r'chem:{}, bond:{}-{} ({}~{})$\rm \AA$'.format(chem[i], chem[0], chem[1], *bond_range), fontsize='large', pad=10.)
        plt.ylabel('Coordination number ({})'.format(chem[i]), fontsize='x-large')
        if dt is not None:
            plt.xlabel('Time (ns)', fontsize='x-large')
        else:
            plt.xlabel('Frame index', fontsize='x-large')
        plt.tick_params(axis="both",direction="in", labelsize='x-large')
        plt.subplots_adjust(left=0.30, bottom=0.25, right=0.70, top=0.75, wspace=0.2, hspace=0.2)
        plt.grid(alpha=0.5)
        plt.legend(fontsize='large').set_draggable(True)
    plt.show()
