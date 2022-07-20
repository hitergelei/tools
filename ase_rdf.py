#!/usr/bin/env python
import numpy as np
from subprocess import call
from ase.io import read

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
    rho,
    ):
    """
                              inf                sin(kr) 
    S(k) = 1 + 4 \pi \rho dr (sum) r^2 {g(r)-1} ---------
                              r=0                  kr    
    where \rho: Number density
          g(r): RDF
    """
    dr = r[1] - r[0]
    k = np.fft.fftfreq(len(r)) / dr
    kr_matrix = k.reshape(-1,1) *r.reshape(-1,1).T
    S = 1. +4*np.pi *rho *dr *np.sum(
        np.reshape(r**2 *(RDF-1), (1,-1)) *np.sinc(kr_matrix/np.pi),
        axis=1,
        )
    # print(np.reshape(r**2 *(RDF-1), (1,-1)) *np.sinc(kr_matrix/np.pi))
    # print(np.sum(np.reshape(r**2 *(RDF-1), (1,-1)) *np.sinc(kr_matrix/np.pi), axis=1))
    realpart = k >= 0.
    return k[realpart], S[realpart]

def get_curve(
    alist,
    image_slice,
    alist_file,
    symb1,
    symb2,
    nBin,
    rcut,
    load_bool,
    save_bool,
    rectify_cut,
    gsmear_std,
    dr,
    ):
    # Slice process
    from ss_util import str_slice_to_list
    slice_list = str_slice_to_list(image_slice)
    # out file
    out_fname  = 'rdf-saved/{}_slice-{}-{}-{}_sym-{}-{}_nBin-{}_rcut-{}_.npy'.format(
        alist_file, *slice_list, symb1, symb2, nBin, rcut)
    out_fname2 = 'rdf-saved/{}_slice-{}-{}-{}_sym-{}-{}_nBin-{}_rcut-{}_.npy'.format(
        alist_file, *slice_list, symb2, symb1, nBin, rcut)

    ## Main
    try:
        assert load_bool == True
        curve = np.load(out_fname)
    except:
        try:
            assert load_bool == True
            curve = np.load(out_fname2)
        except:
            do_calc = True
            if load_bool:
                print('Failed to load saved npy file. Calculation will be carried out')
                print(' Failed to load npy file "{}"'.format(out_fname))
                print('      or equivalent data "{}"'.format(out_fname2))
        else:
            print('File "{}" has been loaded.'.format(out_fname2))
            do_calc = False
        if do_calc:
            curve = get_RDF(alist, rcut, nBin, (symb1, symb2), log=True)
            if save_bool:
                from ss_util import pick_folder_from_path as pffp
                folder = pffp(out_fname)
                call('mkdir -p {}'.format(folder), shell=True)
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

    return curve

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    This code will give you the (Total/Partial) Raidial Distribution Function.
    Return npy file.
    """)
    # Positional arguments
    parser.add_argument('symbol1', type=str, help='symbol1,2, are chemical symbols consisting bonds.')
    parser.add_argument('symbol2', type=str, help='e.g. Ge Te | "a": any symbols, "x": do all partial.')
    parser.add_argument('alist_file', type=str, help='ASE readable atoms list file name.')
    # Optional arguments
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='Image slice following python convention. default=":" (e.g.) -n :1000:10')
    parser.add_argument('-r', '--rcut', type = float, default=8.5, help='Maximum radius for RDF. Default: 8.5')
    parser.add_argument('-b', '--nBin', type=int, default=500, help='Number of bins. Default: 500')
    parser.add_argument('-g', '--gsmear', type=float, default=0., help='Width(simga, STD) of Gaussian smearing in Angstrom unit. Zero means no smearing. [default: 0]')
    parser.add_argument('-e', '--rectify-cut', type=float, default=None, help='All of drastic kink higher than this will be omitted. [Default: no rectify]')
    parser.add_argument('-m', '--multiply', type=float, default=1., help='Multiply this value to RDF (re-scale). [default: 1.]')
    parser.add_argument('-s', '--dont_save', dest='save_bool', action='store_false', help='If provided, npy will not be saved. Default: Save array')
    parser.add_argument('-o', '--dont_load', dest='load_bool', action='store_false', help='If provided, npy will not be loaded. Default: Load if possible')
    parser.add_argument('-t', '--dont_share_y', action='store_true', help='Subplots will not share y-axes if provided.')
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

    #
    den_list = []
    ## Read inputs
    alist = read(args.alist_file, args.image_slice)
    if not isinstance(alist, list):
        alist = [alist]

    den_list = []
    for atoms in alist:
        den_list.append(len(atoms) / atoms.get_volume())
    num_den = np.mean(den_list)

    # In case symbol is 'x'
    chem_list = np.unique(alist[0].get_chemical_symbols()).tolist()
    if symbol1 == 'x':
        symb1_list = chem_list[:]
    else:
        symb1_list = [symbol1]
    if symbol2 == 'x':
        symb2_list = chem_list[:]
    else:
        symb2_list = [symbol2]
    
    # Make symbol_sets
    symbol_sets = []
    if len(symb1_list) == 1 or len(symb2_list) == 1:
        for s1 in symb1_list:
            for s2 in symb2_list:
                symbol_sets.append([s1, s2])
    else:
        for i in range(len(chem_list)):
            for j in range(i,len(chem_list)):
                symbol_sets.append([chem_list[i], chem_list[j]])

    # Main
    curve_list = []
    for symb_set in symbol_sets:
        cv = get_curve(
            alist,
            args.image_slice,
            args.alist_file,
            symb_set[0],
            symb_set[1],
            nBin,
            rcut,
            args.load_bool,
            args.save_bool,
            rectify_cut,
            gsmear_std,
            dr,
            )
        cv[:,1] *= args.multiply,
        curve_list.append(cv)

    # @ Get structure factor
    k_list = []
    S_list = []
    for curve in curve_list:
        k, S = get_s_factor(curve[:,0], curve[:,1], num_den)
        k_list.append(k)
        S_list.append(S)

    # @ Plot
    title  = '{} slice-{} symb-{},{} nBin-{} rcut-{}'.format(
        args.alist_file, args.image_slice, symbol1, symbol2, nBin, rcut)
    import matplotlib.pyplot as plt
    if args.dont_share_y:
        fig, axs = plt.subplots(len(curve_list), sharex=True)
    else:
        fig, axs = plt.subplots(len(curve_list), sharex=True, sharey=True)
    if not isinstance(axs, np.ndarray):
        axs = [axs]

    # Plot RDF
    if args.rdf_upper is not None:
        rdf_upper = args.rdf_upper
    else:
        rdf_upper = np.max(np.array(curve_list)[:,:,1]) * 1.10
    for i in range(len(curve_list)):
        #
        axs[i].plot(curve_list[i][:,0], curve_list[i][:,1], c='k', lw=2)
        #
        if (symbol_sets[i][0], symbol_sets[i][1]) == ('a', 'a'):
            axs[i].set_ylabel(r'$g_{tot}(r)$', fontsize='x-large')
        else:
            axs[i].set_ylabel(r'$g_{{{}}}$(r)'.format(symbol_sets[i][0]+symbol_sets[i][1]), fontsize='x-large')
        #
        if args.xtick_list is not None:
            axs[i].set_xticks(args.xtick_list)
        else:
            intvl = int(rcut // 10 + 1)
            axs[i].set_xticks(range(0, int(rcut)+1, intvl))
        #
        if args.ytick_list is not None:
            axs[i].set_yticks(args.ytick_list)
        else:
            intvl = int(rdf_upper // 4 + 1)
            axs[i].set_yticks(range(0, int(rdf_upper)+1, intvl))
        #
        axs[i].tick_params(axis="both",direction="in", labelsize='x-large', labelbottom=False)
        axs[i].set_xlim(0., rcut)
        axs[i].set_ylim(args.rdf_lower, rdf_upper)
        axs[i].axhline(1., linestyle='dashed', linewidth=1, c='k')
        axs[i].grid(alpha=0.4)
    axs[-1].tick_params(axis="both",direction="in", labelsize='x-large', labelbottom=True)
    axs[-1].set_xlabel('Distance $(\AA)$', fontsize='x-large')
    axs[0].set_title(title, pad=10)
    bottom = (1.-len(axs)*0.1) /2.
    plt.subplots_adjust(left=0.20, bottom=bottom, right=0.80, top=1-bottom, wspace=0.20, hspace=0.20)

    # Plot S(Q)
    if args.s_upper is not None:
        s_upper = args.s_upper
    else:
        s_upper = np.max(np.array(S_list)) * 1.10
    if args.dont_share_y:
        fig, axs = plt.subplots(len(curve_list), sharex=True)
    else:
        fig, axs = plt.subplots(len(curve_list), sharex=True, sharey=True)
    if not isinstance(axs, np.ndarray):
        axs = [axs]
    for i in range(len(curve_list)):
        #
        axs[i].plot(k_list[i], S_list[i], c='k', lw=2)
        #
        if (symbol_sets[i][0], symbol_sets[i][0]) == ('a', 'a'):
            axs[i].set_ylabel('$S_{tot}(Q)$', fontsize='x-large')
        else:
            axs[i].set_ylabel('$S_{{{}}}(Q)$'.format(symbol_sets[i][0]+symbol_sets[i][1]), fontsize='x-large')
        #
        if args.s_xtick_list is not None:
            axs[i].set_xticks(args.s_xtick_list)
        else:
            intvl = int(np.max(k_list[i]) // 10 + 1)
            axs[i].set_xticks(range(0, int(np.max(k_list[i]))+1, intvl))
        #
        if args.s_ytick_list is not None:
            axs[i].set_yticks(args.s_ytick_list)
        else:
            intvl = int(s_upper // 4 + 1)
            axs[i].set_yticks(range(0, int(s_upper)+1, intvl))
        #
        axs[i].tick_params(axis="both",direction="in", labelsize='x-large', labelbottom=False)
        axs[i].set_xlim(0., np.max(k_list[i]))
        axs[i].set_ylim(args.s_lower, s_upper)
        axs[i].axhline(1., linestyle='dashed', linewidth=1, c='k')
        axs[i].grid(alpha=0.4)
    axs[-1].tick_params(axis="both",direction="in", labelsize='x-large', labelbottom=True)
    axs[-1].set_xlabel('Q $(\AA^{-1})$', fontsize='x-large')
    axs[0].set_title(title, pad=10)
    bottom = (1.-len(axs)*0.1) /2.
    plt.subplots_adjust(left=0.20, bottom=bottom, right=0.80, top=1-bottom, wspace=0.20, hspace=0.20)
    plt.show()
