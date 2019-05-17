#!/usr/bin/env python

import numpy as np

def get_RDF(alist, rMax, nBins=500, chem=None, log=False):
    from asap3.analysis.rdf import RadialDistributionFunction as RDF
    RDFobj = RDF(atoms=alist[0],
                 rMax=rMax,
                 nBins=nBins)
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
        norm_const = counts[list(unique).index(chem2)] / np.sum(counts)
        #
        from chem_num_inverter import invert_chem_num
        spec_inds = invert_chem_num(chem)
        #
        rdf = RDFobj.get_rdf(elements=tuple(spec_inds)) / norm_const
    x = np.arange(nBins) * rMax / nBins
    ## Return curve
    return np.transpose(np.concatenate(([x], [rdf])))

def rectify_RDF(curve, rectify_cut):
    iter = True
    while True:
        test = curve[1:] - curve[:-1]
        peak_bool = np.array(list(test[:,1] > (-1 * rectify_cut)) + [True], dtype=np.bool)
        if False not in peak_bool:
            break
        curve = curve[peak_bool]
    return curve

if __name__ == '__main__':
    import sys
    import datetime

    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ POSTECH <<<<<'.center(120))
    print(('code started time: '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('This code will give you (Total/Partial) Raidial Distribution Function'.center(120))
    print('')
    print('Useage    ==> ./ase_rdf.py >1st_spec< >2nd_spec< >atoms list file< >rMax< >nBins< (>rectify_cut<)'.center(120))
    print('Example 1 ==> ./ase_rdf.py     all        all         si.traj        12     500         0.2      '.center(120))
    print('Example 2 ==> ./ase_rdf.py      Ge         Te         si.traj        12     500         0.2      '.center(120))
    print('')
    print('Return ----------> rdf-(original name)-(rMax)-(nBins).npy'.center(120))
    print('')
    print('=================================================================================================='.center(120))

    ## Read input params
    chem1 = sys.argv[1]
    chem2 = sys.argv[2]
    inp_fname = sys.argv[3]
    rMax = float(sys.argv[4])
    nBins = int(sys.argv[5])
    if len(sys.argv) == 7:
        rectify_cut = float(sys.argv[6])
    else:
        rectify_cut = None
    fname = 'rdf-{}-{}{}-{:.2f}-{}.npy'.format(inp_fname, chem1, chem2, rMax, nBins)
    if chem1 != 'all' and chem2 != 'all':
        chem = [chem1, chem2]
    elif chem1 == 'all' and chem2 == 'all':
        chem = None
    else:
        raise ValueError('chem1 ({}) or chem2 ({}) argument is wrong. (Or both)'.format(chem1, chem2))

    ## Load saved file
    try:
        curve = np.load(fname)
    except:
        print('File "{}" not found. Calculation will be carried out.'.format(fname).center(120))
        ## Read inputs
        from ase.io import read
        alist = read(inp_fname, ':')
        curve = get_RDF(alist, rMax, nBins, chem, log=True)
        np.save(fname, curve)
        print('=================================================================================================='.center(120))
        print('Return ----------> {}'.format(fname).center(120))
        print('=================================================================================================='.center(120))
    else:
        print('File "{}" is loaded.'.format(fname).center(120))


    ## Rectify curve
    if rectify_cut:
        curve = rectify_RDF(curve, rectify_cut)
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
