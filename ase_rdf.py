#!/usr/bin/env python

import numpy as np

def get_RDF(alist, rMax, nBin=500, chem=None, log=False):
    from asap3.analysis.rdf import RadialDistributionFunction as RDF
    RDFobj = RDF(atoms=alist[0],
                 rMax=rMax,
                 nBins=nBin)
    for i in range(1,len(alist)):
        RDFobj.atoms = alist[i]
        RDFobj.update()
        if log and i % 1000 == 999:
            print('\t Updating '+str(i+1)+" th image's RDF")

    ## Total RDF
    if chem == None:
        rdf = RDFobj.get_rdf()
    ## Partial RDF
    else:
        # Get normalize constant
        (unique, counts) = np.unique(alist[0].get_chemical_symbols(), return_counts=True)
        norm_const = counts[list(unique).index(chem[1])] / np.sum(counts, dtype=np.float)
        #
        from chem_num_inverter import invert_chem_num
        spec_inds = invert_chem_num(chem)
        #
        rdf = RDFobj.get_rdf(elements=tuple(spec_inds)) / norm_const
    x = np.arange(nBin) * rMax / nBin
    ## Return curve
    return np.transpose(np.concatenate(([x], [rdf])))

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    This code will give you the (Total/Partial) Raidial Distribution Function.
    Return npy file.
    """)
    # Positional arguments
    parser.add_argument('inp_file', type=str, help='ASE readable atoms list file name')
    # Optional arguments
    parser.add_argument('-p', '--partial', type=str, default='total-total', help='If you need partial RDF. Default: total-total RDF. E.g.: -p Ge-Te')
    parser.add_argument('-r', '--rmax', type = float, default=12., help='Maximum radius for RDF. Default: 12.')
    parser.add_argument('-n', '--nbin', type=int, default=500, help='Number of bins. Default: 500')
    parser.add_argument('-c', '--rectify-cut', type=float, default=False, help='All of drastic variation higher than this value will be omitted in output. Default: no rectify')
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
    # --partial
    if args.partial == 'total-total':
        chem = None
    elif len(args.partial.split('-')) == 2:
        chem = args.partial.split('-')
    else:
        raise ValueError('--partial (or -p) argument (={}) is wrong.'.format(args.partial))
    # --rmax, --nbin
    rMax = args.rmax
    nBin = args.nbin
    # --rectify_cut
    rectify_cut = args.rectify_cut
    if rectify_cut == False:
        rectify_cut = False
    # --inp_file
    inp_fname = args.inp_file
    out_fname = 'rdf-{}-{}-{:.2f}-{}.npy'.format(inp_fname, args.partial, rMax, nBin)
    if chem is not None:
        out_fname2 = 'rdf-{}-{}-{:.2f}-{}.npy'.format(inp_fname, chem[1]+'-'+chem[0], rMax, nBin)

    ## Load saved file
    try:
        curve = np.load(out_fname)
    except:
        print('File "{}" was not found.'.format(out_fname).center(120))
        do_calc = True
        if chem is not None:
            try:
                curve = np.load(out_fname2)
            except:
                print('Also, file "{}" was not found. Calculation will be carried out.'.format(out_fname2).center(120))
                do_calc = True
            else:
                out_fname = out_fname2
                print('File "{}" is loaded.'.format(out_fname).center(120))
                do_calc = False
        if do_calc:
            ## Read inputs
            from ase.io import read
            alist = read(inp_fname, ':')
            curve = get_RDF(alist, rMax, nBin, chem, log=True)
            np.save(out_fname, curve)
            print('=================================================================================================='.center(120))
            print('Return ----------> {}'.format(out_fname).center(120))
            print('=================================================================================================='.center(120))
    else:
        print('File "{}" is loaded.'.format(out_fname).center(120))

    ## Rectify curve
    if rectify_cut:
        from ss_util import rectify_curve
        curve = rectify_curve(curve, rectify_cut)
    ## Plot
    import matplotlib.pyplot as plt
    # spline (optional)
    # from scipy.interpolate import make_interp_spline, BSpline
    # spl = make_interp_spline(curve[:-1, 0][::20], curve[:-1, 1][::20])
    # curve[:-1, 1] = spl(curve[:-1, 0])
    plt.plot(curve[:-1, 0], curve[:-1, 1], linewidth=3)
    if chem == None:
        plt.ylabel(r'$g_{tot}(r)$', fontsize='xx-large')
    else:
        plt.ylabel(r'$g_{{{}{}}}(r)$'.format(chem[0], chem[1]), fontsize='xx-large')
    plt.xlabel(r'$r\ /\ \AA$', fontsize='xx-large')
    plt.subplots_adjust(left=0.15, bottom=0.37, right=0.95, top=0.63, wspace=0.2, hspace=0.2)
    plt.xticks(fontsize='xx-large')
    plt.yticks(fontsize='xx-large')
    plt.xlim(0., rMax)
    plt.hlines(1., 0., rMax, linestyle='dashed', linewidth=2)
    plt.show()
