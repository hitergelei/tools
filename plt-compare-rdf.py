#!/usr/bin/env python
import numpy as np

## Functions
def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    This code will give you comparison graph of two RDFs
    """)
    # Positional arguments
    parser.add_argument('rdf_npy_list', type=str, nargs='+', help='RDF npy files.')
    # Optional arguments
    parser.add_argument('-e', '--rectify_cut', type=float, default=0, help='Throw away data of kink larger than cutoff. [Default: no rectify]')
    parser.add_argument('-g', '--gsmear', type=float, default=0., help='Width(simga, STD) of Gaussian smearing in Angstrom unit. [default: 0]')
    parser.add_argument('-u', '--rdf_upper', type=float, default=None, help='Upper bound for RDF plot [Default: automatic]')
    parser.add_argument('-l', '--rdf_lower', type=float, default=0, help='Lower bound for RDF plot [Default: 0]')
    parser.add_argument('-b', '--label_list', type=str, nargs='+', default=None, help='Legend label list. len(label_list) == len(rdf_npy_list). [default: no legend].')
    parser.add_argument('-s', '--lstyle_list', type=str, nargs='+', default=None, help='Line style list. len(lstyle_list) == len(rdf_npy_list). [default: solid line].')
    parser.add_argument('-c', '--lcolor_list', type=str, nargs='+', default=None, help='Line color list. len(lcolor_list) == len(rdf_npy_list). [default: automatic].')
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
    print('This code will give you comparison graph of two RDFs'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    ## Read Inputs
    fname_list   = args.rdf_npy_list
    rdf_npy_list = [np.load(fname_list[i]) for i in range(len(fname_list))]
    radii_arr    = rdf_npy_list[0][:,0]
    rdf_arr_list = [rdf_npy_list[i][:,1] for i in range(len(fname_list))]
    rectify_cut  = args.rectify_cut
    label_list   = args.label_list
    lstyle_list  = args.lstyle_list
    lcolor_list   = args.lcolor_list
    if label_list and len(label_list) != len(rdf_npy_list):
        raise ValueError('len(label_list)=={}  != len(rdf_npy_list)=={}'.format(len(label_list), len(rdf_npy_list)))
    if not label_list:
        label_list = [None]*len(rdf_npy_list)
    if lstyle_list and len(lstyle_list) != len(rdf_npy_list):
        raise ValueError('len(lstyle_list)=={}  != len(rdf_npy_list)=={}'.format(len(label_list), len(rdf_npy_list)))
    if not lstyle_list:
        lstyle_list = [None]*len(rdf_npy_list)
    if lcolor_list and len(lcolor_list) != len(rdf_npy_list):
        raise ValueError('len(lcolor_list)=={}  != len(rdf_npy_list)=={}'.format(len(label_list), len(rdf_npy_list)))
    if not lcolor_list:
        lcolor_list = [None]*len(rdf_npy_list)

    ## Prepare data
    symbol_list, nBin, rcut = [],[],[]
    for fname in fname_list:
        symbol_tmp, nBin_tmp, rcut_tmp = fname.split('_')[-4:-1]
        symbol_tmp = symbol_tmp.split('-')[1:]
        symbol_list.append(symbol_tmp); nBin.append(nBin_tmp); rcut.append(rcut_tmp)
        if nBin_tmp != nBin[0] or rcut_tmp != rcut[0]:
            raise ValueError('Files are not proper to compare. \nMust have same "nBin, cutoff radii". {} != {} or {} != {}'.format(nBin_tmp, nBin[0], rcut_tmp, rcut[0]))
        if symbol_tmp != symbol_list[0] and symbol_tmp[::-1] != symbol_list[0]:
            raise ValueError('Files are not proper to compare. \nMust have same "Chemical symbols". {} != {}'.format(symbol_tmp, symbol_list[0]))
    symbol_list = symbol_list[0]
    nBin        = float(nBin[0].split('-')[1])
    rcut        = float(rcut[0].split('-')[1])
    dr   = rcut/ nBin

    ## Rectify curve
    if rectify_cut:
        from ss_util import rectify_curve
        rdf_arr_list = [rectify_curve(rdf_arr_list[i], rectify_cut) for i in range(len(rdf_arr_list))]
    ## Gaussian smearing
    if args.gsmear:
        from scipy.ndimage.filters import gaussian_filter1d
        rdf_arr_list = [gaussian_filter1d(rdf_arr_list[i], args.gsmear /dr) for i in range(len(rdf_arr_list))]

    ## Plot
    import matplotlib.pyplot as plt
    font = {'family':'Arial'}
    plt.rc('font', **font)
    [plt.plot(radii_arr, rdf_arr_list[i], label=label_list[i], linestyle=lstyle_list[i], c=lcolor_list[i]) for i in range(len(rdf_arr_list))]
    if symbol_list == ['a','a']:
        plt.ylabel('Total RDF $(\AA^{-1})$', fontsize='x-large')
    else:
        plt.ylabel('Partial RDF$_{{{}}}$'.format(symbol_list[0]+symbol_list[1])+' $(\AA^{-1})$', fontsize='x-large')
    plt.xlabel('Distance $(\AA)$', fontsize='x-large')
    plt.subplots_adjust(left=0.09, bottom=0.40, right=0.97, top=0.60, wspace=0.2, hspace=0.2)
    from math import ceil, floor
    plt.xticks(range(1,int(ceil(rcut))), fontsize='x-large')
    plt.yticks(range(0,int(ceil(np.amax(rdf_arr_list))),2), fontsize='x-large')
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xlim(0., rcut)
    plt.ylim(args.rdf_lower, args.rdf_upper)
    plt.axhline(1., linestyle='dashed', linewidth=1, c='k')
    plt.title('tmp', pad=10)
    plt.grid(alpha=0.2)
    ## Legend options
    if label_list != [None]*len(label_list):
        plt.legend(loc='upper left', fontsize='large')
    else:
        plt.legend().set_visible(False)
    plt.show()

