#!/usr/bin/env python
import numpy as np

## Functions
def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    Will run nappy ADF plot.
    """)
    # Positional arguments
    parser.add_argument('symbol1', type=str, help='symbol1,2,3, are chemical symbols consisting bonds')
    parser.add_argument('symbol2', type=str, help='around the angle like symbol1-2-3. "a": any symbols.')
    parser.add_argument('symbol3', type=str, help='e.g.1. Ge Te Ge | e.g.2. a a a | e.g.3 a Te a')
    parser.add_argument('file_list', type=str, nargs='+', help='ASE readable atoms list file name. Multiple input files can be read.')
    # Optional arguments
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='Image slice following python convention. default=":" (e.g.) -n :1000:10')
    parser.add_argument('-w', '--dDeg', type=float, default=1., help='Width of the angular degree. [default: 1.0]')
    parser.add_argument('-r', '--rcut', type=float, default=3., help='Cutoff radius of the bonding pair. [default: 3.0]')
    parser.add_argument('-f', '--f_format', type=str, default=None, help='Input file format. [default: Automatic if ASE readable.]')
    parser.add_argument('-g', '--gsmear', type=float, default=0., help='Width(simga, STD) of Gaussian smearing in degree unit. Zero means no smearing. [default: 0]')
    parser.add_argument('-a', '--no_average', dest='avg_bool', action='store_false', help='Not to take average over files. [default: take average]')
    parser.add_argument('-s', '--dont_save', dest='save_bool', action='store_false', help='If provided, ADOS arrays will not be saved. Default: Save array')
    parser.add_argument('-o', '--dont_load', dest='load_bool', action='store_false', help='If provided, ADOS arrays will not be loaded. Default: Load if possible')
    parser.add_argument('-m', '--Nprocs', type=int, default=1, help='Number of process for multiprocessing. [Default: serial compute]')
    parser.add_argument('-u', '--adf_upper', type=float, default=None, help='Upper bound for ADF plot [Default: automatic]')
    parser.add_argument('-l', '--adf_lower', type=float, default=0, help='Lower bound for ADF plot [Default: 0]')
    return parser.parse_args()

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
    print('Will run nappy ADF plot.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()
    from nappy.adf import adf_atom, adf, adf_gather
    from chemical_symbol_number_inverter import invert_chem_sym_num as ics

    ## def
    infiles  = args.file_list
    symbol1  = args.symbol1
    symbol2  = args.symbol2
    symbol3  = args.symbol3
    dang     = float(args.dDeg)
    drad     = np.pi *dang/180.0
    rcut     = float(args.rcut)
    sigma    = int(args.gsmear)
    avg_bool = args.avg_bool
    ffmt     = args.f_format
    na       = int(180.0/dang) +1
    ## Slice process
    if ':' in args.image_slice:
        slice_list = [int(e) if e.strip() else None for e in args.image_slice.split(":")]
        if len(slice_list) == 2:
            slice_list.append(None)
    else:
        slice_list = [int(args.image_slice), int(args.image_slice)+1, None]
    if slice_list[0] == None:
        slice_list[0] = 0
    if slice_list[2] == None:
        slice_list[2] = 1

    fname = 'adf_saved/{}-slice_{}_{}_{}-_{}-{}-{}_dDeg-{}_rcut-{}_avg-{}.npz'.format(
        infiles[0], slice_list[0], slice_list[1], slice_list[2], symbol1, symbol2, symbol3, dang, rcut, avg_bool)
    try:
        assert args.load_bool == True
        assert len(infiles) == 1
        npz = np.load('{}'.format(fname))
    except:
        if args.load_bool == True:
            print('Failed to load saved npz file. Calculation will be carried out')
            print('Case 1) Number of input file must be 1 to load npz. len(infiles)=={}'.format(len(infiles)))
            print('Case 2) Failed to load npz file "{}"'.format(fname))

        alist = []
        from ase.io import read
        for infname in infiles:
            read_obj = read(infname, args.image_slice)
            if isinstance(read_obj, list):
                alist.extend(read_obj)
            else:
                alist.append(read_obj)
        if len(alist) == 0: raise ValueError('No an image provided.')

        ## Multiprocessing
        from time import time
        time_init = time()
        print('Caculation started.')
        from multiprocessing import Pool
        pool = Pool(args.Nprocs)
        tasks = [pool.apply_async(adf_gather, (alist[_::args.Nprocs], dang, rcut, symbol1, symbol2, symbol3)) for _ in range(args.Nprocs)]

        gather = []
        for task in tasks:
            task.wait(10)
            gather.append(task.get())
        gather = np.array(gather)

        for items in gather:
            angd = gather[0,0]
            agr  = np.sum(gather[:,1], axis=0)
            nsum = np.sum(gather[:,2], axis=0)
        agr /= dang ## Normalized to continueous version. Integral = average angle-pair numbers
        if args.avg_bool:
            agr /= nsum
        print('Caculation ended. Elapse time = {} (s)'.format(time()-time_init))
        print('Totally {} atoms in {} images have been calculated'.format(nsum, len(alist)))

        if args.save_bool and len(infiles) == 1:
            from subprocess import call
            call('mkdir adf_saved', shell=True)
            np.savez(fname, angd=angd, agr=agr)
            print('{} file has saved'.format(fname))
    else:
        print('File "{}" has been loaded successfully.'.format(fname))
        angd, agr = npz['angd'], npz['agr']

    if not sigma == 0:
        print(' Gaussian smearing...')
        # from gaussian_smear import gsmear
        # agr= gsmear(angd,agr,sigma)
        from scipy.ndimage.filters import gaussian_filter1d
        agr = gaussian_filter1d(agr, sigma /dang)
    ## Debug option
    print('Average number of normalized angle-pairs.={}'.format(np.trapz(agr, angd)))

    ## Plot
    import matplotlib.pyplot as plt
    plt.plot(angd,agr,'-')
    if (symbol1, symbol2, symbol3) == ('a','a','a'):
        plt.ylabel('Total ADF (deg$^{-1}$)', fontsize='x-large')
    else:
        print('Partial ADF (deg){}'.format(chem[0]+chem[1]+chem[2]))
        plt.ylabel('Partial ADF$_{{{}}}$'.format(chem[0]+chem[1]+chem[2])+' (deg$^{-1}$)', fontsize='x-large')
    plt.xlabel('Bond angle (deg)', fontsize='x-large')
    plt.subplots_adjust(left=0.16, bottom=0.28, right=0.94, top=0.75, wspace=0.20, hspace=0.20)
    plt.xticks(range(0,181,20),fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xlim(0., 180.)
    plt.ylim(args.adf_lower, args.adf_upper)
    
    plt.show()
