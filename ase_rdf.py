#!/usr/bin/env python
import numpy as np

def get_RDF(
    alist,
    rcut,
    nBin=500,
    symbol_tuple=None,
    log=False,
    ):
    from asap3.analysis.rdf import RadialDistributionFunction as RDF
    RDFobj = RDF(
        atoms=alist[0],
        rMax=rcut,
        nBins=nBin,
        )
    for i in range(1,len(alist)):
        RDFobj.atoms = alist[i]
        RDFobj.update()
        if log and i % 1000 == 999:
            print('\t Updating '+str(i+1)+" th image's RDF")

    ## Total RDF
    if symbol_tuple == ('a', 'a'):
        rdf = RDFobj.get_rdf()
    ## Partial RDF
    else:
        # Get normalize constant
        (unique, counts) = np.unique(alist[0].get_chemical_symbols(), return_counts=True)
        norm_const = counts[list(unique).index(symbol_tuple[1])] / np.sum(counts, dtype=np.float)
        #
        from chemical_symbol_number_inverter import invert_chem_sym_num
        spec_inds = invert_chem_sym_num(symbol_tuple)
        #
        rdf = RDFobj.get_rdf(elements=tuple(spec_inds)) / norm_const
    x = np.arange(nBin) / float(nBin) * rcut
    ## Return curve
    return np.transpose(np.concatenate(([x], [rdf])))

def get_s_factor(
    r,
    RDF,
    ):
    """
                              inf                sin(kr) 
    S(k) = 1 + 4 \pi \rho dr (sum) r^2 {g(r)-1} ---------
                              r=0                  kr    
    where \rho: density
          g(r): RDF
    """
    dr = r[1] - r[0]
    k = np.fft.fftfreq(len(r)) / dr
    kr_matrix = k.reshape(-1,1) *r.reshape(-1,1).T
    S = 1. +4*np.pi *0.030 *dr *np.sum(
        np.reshape(r**2 *(RDF-1), (1,-1)) *np.sinc(kr_matrix/np.pi),
        axis=1,
        )
    # print(np.reshape(r**2 *(RDF-1), (1,-1)) *np.sinc(kr_matrix/np.pi))
    # print(np.sum(np.reshape(r**2 *(RDF-1), (1,-1)) *np.sinc(kr_matrix/np.pi), axis=1))
    realpart = k >= 0.
    return k[realpart], S[realpart]

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    This code will give you the (Total/Partial) Raidial Distribution Function.
    Return npy file.
    """)
    # Positional arguments
    parser.add_argument('symbol1', type=str, help='symbol1,2, are chemical symbols consisting bonds.')
    parser.add_argument('symbol2', type=str, help='e.g. Ge-Te == Te-Ge | "a": any symbols.')
    parser.add_argument('file_list', type=str, nargs='+', help='ASE readable atoms list file name. Multiple files can be read.')
    # Optional arguments
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='Image slice following python convention. default=":" (e.g.) -n :1000:10')
    parser.add_argument('-r', '--rcut', type = float, default=8.5, help='Maximum radius for RDF. Default: 8.5')
    parser.add_argument('-b', '--nBin', type=int, default=500, help='Number of bins. Default: 500')
    parser.add_argument('-g', '--gsmear', type=float, default=0., help='Width(simga, STD) of Gaussian smearing in Angstrom unit. Zero means no smearing. [default: 0]')
    parser.add_argument('-e', '--rectify-cut', type=float, default=None, help='All of drastic kink higher than this will be omitted. [Default: no rectify]')
    parser.add_argument('-s', '--dont_save', dest='save_bool', action='store_false', help='If provided, npy will not be saved. Default: Save array')
    parser.add_argument('-o', '--dont_load', dest='load_bool', action='store_false', help='If provided, npy will not be loaded. Default: Load if possible')
    parser.add_argument('-u', '--rdf_upper', type=float, default=None, help='Upper bound for RDF plot [Default: automatic]')
    parser.add_argument('-l', '--rdf_lower', type=float, default=0, help='Lower bound for RDF plot [Default: 0]')
    parser.add_argument('-p', '--s_upper', type=float, default=None, help='Upper bound for S(Q) plot [Default: automatic]')
    parser.add_argument('-q', '--s_lower', type=float, default=0, help='Lower bound for S(Q) plot [Default: 0]')
    parser.add_argument('-x', '--xtick_list', type=float, nargs='+', default=None, help='Specify x ticks of RDF. [Default: automatic]')
    parser.add_argument('-y', '--ytick_list', type=float, nargs='+', default=None, help='Specify y ticks of RDF. [Default: automatic]')
    parser.add_argument('-v', '--s_xtick_list', type=float, nargs='+', default=None, help='Specify x ticks of S(Q). [Default: automatic]')
    parser.add_argument('-w', '--s_ytick_list', type=float, nargs='+', default=None, help='Specify y ticks of S(Q). [Default: automatic]')
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
    print('This code will give you the (Total/Partial) Raidial Distribution Function'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    ## Read input params
    # params
    symbol1     = args.symbol1
    symbol2     = args.symbol2
    rcut        = args.rcut
    nBin        = args.nBin
    dr          = rcut/ nBin
    gsmear_std  = args.gsmear
    rectify_cut = args.rectify_cut
    # Slice process
    from ss_util import str_slice_to_list
    slice_list = str_slice_to_list(args.image_slice)
    # out file
    file_list = args.file_list
    out_fname  = 'rdf-saved/{}_slice-{}-{}-{}_sym-{}-{}_nBin-{}_rcut-{}_.npy'.format(
        file_list[0], slice_list[0], slice_list[1], slice_list[2], symbol1, symbol2, nBin, rcut)
    out_fname2 = 'rdf-saved/{}_slice-{}-{}-{}_sym-{}-{}_nBin-{}_rcut-{}_.npy'.format(
        file_list[0], slice_list[0], slice_list[1], slice_list[2], symbol2, symbol1, nBin, rcut)

    ## Main
    try:
        assert args.load_bool == True
        assert len(file_list) == 1
        curve = np.load(out_fname)
    except:
        try:
            assert args.load_bool == True
            assert len(file_list) == 1
            curve = np.load(out_fname2)
        except:
            do_calc = True
            if args.load_bool:
                print('Failed to load saved npy file. Calculation will be carried out')
                print('Case 1) Number of input file must be 1 to load npy. len(file_list)=={}'.format(len(file_list)))
                print('Case 2) Failed to load npy file "{}"'.format(out_fname))
                print('             or equivalent data "{}"'.format(out_fname2))
        else:
            print('File "{}" has been loaded.'.format(out_fname2))
            do_calc = False
        if do_calc:
            ## Read inputs
            alist = []
            from ase.io import read
            for infname in file_list:
                read_obj = read(infname, args.image_slice)
                if isinstance(read_obj, list):
                    alist.extend(read_obj)
                else:
                    alist.append(read_obj)
            if len(alist) == 0: raise ValueError('No an image provided.')
            curve = get_RDF(alist, rcut, nBin, (symbol1, symbol2), log=True)
            if args.save_bool and len(file_list) == 1:
                from subprocess import call
                call('mkdir rdf-saved', shell=True)
                np.save(out_fname, curve)
                print('=================================================================================================='.center(120))
                print('RDF saved! ----------> {}'.format(out_fname).center(120))
                print('=================================================================================================='.center(120))
    else:
        print('File "{}" has been loaded.'.format(out_fname))

    # @ Rectify curve
    if rectify_cut:
        from ss_util import rectify_curve
        curve = rectify_curve(curve, rectify_cut)

    if not gsmear_std == 0:
        print(' Gaussian smearing...')
        # from gaussian_smear import gsmear
        # agr= gsmear(angd,agr,gsmear_std)
        from scipy.ndimage.filters import gaussian_filter1d
        curve[:,1] = gaussian_filter1d(curve[:,1], gsmear_std /dr)
    # Debug option
    print('Integration of RDF.={}'.format(np.trapz(curve[:,1], curve[:,0])))

    # @ Get structure factor
    k, S = get_s_factor(curve[:,0], curve[:,1])

    # @ Plot
    import matplotlib.pyplot as plt
    font = {'family':'Arial'}
    plt.rc('font', **font)
    # Plot RDF
    plt.plot(curve[:,0], curve[:,1], c='r')
    if (symbol1, symbol2) == ('a', 'a'):
        plt.ylabel('Total RDF $(\AA^{-1})$', fontsize='x-large')
    else:
        plt.ylabel('Partial RDF$_{{{}}}$'.format(symbol1+symbol2)+' $(\AA^{-1})$', fontsize='x-large')
    plt.xlabel('Distance $(\AA)$', fontsize='x-large')
    plt.subplots_adjust(left=0.09, bottom=0.30, right=0.97, top=0.60, wspace=0.2, hspace=0.2)
    if args.xtick_list is not None:
        plt.xticks(args.xtick_list ,fontsize='x-large')
    else:
        from math import ceil
        plt.xticks(range(1,int(ceil(rcut))), fontsize='x-large')
    if args.ytick_list is not None:
        plt.yticks(args.ytick_list ,fontsize='x-large')
    else:
        plt.yticks(fontsize='x-large')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xlim(0., rcut)
    plt.ylim(args.rdf_lower, args.rdf_upper)
    plt.axhline(1., linestyle='dashed', linewidth=1, c='k')
    plt.title(out_fname[11:-4], pad=10)
    plt.grid(alpha=0.2)

    # Plot S(Q)
    plt.figure()
    plt.plot(k, S, c='r')
    if (symbol1, symbol2) == ('a', 'a'):
        plt.ylabel('Structure factor, $S(Q)$', fontsize='x-large')
    else:
        plt.ylabel('Partial structure factor, $S_{{{}}}(Q)$'.format(symbol1+symbol2), fontsize='x-large')
    plt.xlabel('Q $(\AA^{-1})$', fontsize='x-large')
    # plt.subplots_adjust(left=0.09, bottom=0.30, right=0.97, top=0.60, wspace=0.2, hspace=0.2)
    if args.s_xtick_list is not None:
        plt.xticks(args.s_xtick_list ,fontsize='x-large')
    else:
        from math import ceil
        plt.xticks(range(1,int(ceil(np.max(k)))), fontsize='x-large')
    if args.s_ytick_list is not None:
        plt.yticks(args.s_ytick_list ,fontsize='x-large')
    else:
        plt.yticks(fontsize='x-large')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xlim(0., np.max(k))
    plt.ylim(args.s_lower, args.s_upper)
    plt.axhline(1., linestyle='dashed', linewidth=1, c='k')
    plt.title(out_fname[11:-4], pad=10)
    plt.grid(alpha=0.2)
    plt.show()
