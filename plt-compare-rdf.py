#!/usr/bin/env python

import numpy as np

## Hyper params
legend_bool = False
ylim_low    = None
ylim_up     = 16.5

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
    print('This code will give you comparison graph of two RDFs'.center(120))
    print('')
    print('Useage  ==> ./plt-compare-rdf.py >DFT RDF file< >ML RDF file< (>rectify_cut<)'.center(120))
    print('Example ==> ./plt-compare-rdf.py   ....DFT.npy   ....DPMD.npy       0.2      '.center(120))
    print('')
    print('=================================================================================================='.center(120))

    ## Read Inputs
    rdf1_f = sys.argv[1]
    rdf2_f = sys.argv[2]
    rdf1 = np.load(rdf1_f)
    rdf2 = np.load(rdf2_f)
    if len(sys.argv) == 4:
        rectify_cut = np.float(sys.argv[3])
    else:
        rectify_cut = None
    ## Split
    [chem11, chem12, rMax1] = rdf1_f.split('-')[-4:-1]
    [chem21, chem22, rMax2] = rdf2_f.split('-')[-4:-1]
    chem1 = chem11+chem12
    chem2 = chem21+chem22
    if chem1 != chem2 or rMax1 != rMax2:
        raise ValueError('Two files seem not proper to compare. ({} != {} or {} != {})'.format(chem1, chem2, rMax1, rMax2))
    elif chem1 != 'totaltotal':
        chem = chem1
    else:
        chem = None
    rMax = np.float(rMax1)

    ## Rectify curve
    if rectify_cut:
        from ase_rdf import rectify_RDF
        rdf1 = rectify_RDF(rdf1, rectify_cut)
        rdf2 = rectify_RDF(rdf2, rectify_cut)
    ## Plot
    import matplotlib.pyplot as plt
    # spline (optional)
    # from scipy.interpolate import make_interp_spline, BSpline
    # spl = make_interp_spline(rdf1[:-1, 0][::20], rdf1[:-1, 1][::20])
    # rdf1[:-1, 1] = spl(rdf1[:-1, 0])
    font = {'family':'Arial'}
    plt.rc('font', **font)
    plt.plot(rdf1[:-1, 0], rdf1[:-1, 1], linewidth=3, label='DFT')
    plt.plot(rdf2[:-1, 0], rdf2[:-1, 1], linewidth=3, label='ML-Pot', linestyle='dashed', color='r')
    if chem == None:
        plt.ylabel(r'$g_{tot}(r)$', fontsize='xx-large')
    else:
        plt.ylabel(r'$g_{{{}}}(r)$'.format(chem), fontsize='xx-large')
    plt.xlabel(r'$r\ (\AA)$', fontsize='xx-large')
    plt.subplots_adjust(left=0.15, bottom=0.37, right=0.95, top=0.63, wspace=0.2, hspace=0.2)
    plt.xticks(fontsize='xx-large')
    plt.yticks(fontsize='xx-large')
    plt.tick_params(axis="y",direction="in")
    plt.tick_params(axis="x",direction="in")
    plt.xlim(0., rMax)
    plt.ylim(ylim_low, ylim_up)
    plt.hlines(1., 0., rMax, linestyle='dashed', linewidth=2)
    ## Legend options
    if legend_bool:
        plt.legend(loc='upper right', fontsize='large')
    
    plt.show()

