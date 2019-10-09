#!/usr/bin/env python

import numpy as np

## Functions
def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    This code will give you comparison graph of two ADFs
    """)
    # Positional arguments
    parser.add_argument('adf_npz_list', type=str, nargs='+', help='ADF npz files. First file will be drawn by black solid line.')
    # Optional arguments
    parser.add_argument('-e', '--rectify_cut', type=float, default=0, help='Throw away data of kink larger than cutoff. [Default: 0]')
    parser.add_argument('-g', '--gsmear', type=float, default=0., help='Width(simga, STD) of Gaussian smearing in degree unit. Zero means no smearing. [default: 0]')
    parser.add_argument('-u', '--adf_upper', type=float, default=None, help='Upper bound for ADF plot [Default: automatic]')
    parser.add_argument('-l', '--adf_lower', type=float, default=0, help='Lower bound for ADF plot [Default: 0]')
    parser.add_argument('-b', '--label_list', type=str, nargs='+', default=None, help='Legend label list. len(label_list) == len(adf_npz_list). [default: no legend].')
    parser.add_argument('-s', '--lstyle_list', type=str, nargs='+', default=None, help='Line style list. len(lstyle_list) == len(adf_npz_list). [default: solid line].')
    parser.add_argument('-c', '--lcolor_list', type=str, nargs='+', default=None, help='Line color list. len(lcolor_list) == len(adf_npz_list). [default: automatic].')
    # parser.add_argument('-n', '--image_slice', type=str, default=':', help='Image slice following python convention. default=":" (e.g.) -n :1000:10')
    # parser.add_argument('-r', '--rcut', type=float, default=3., help='Cutoff radius of the bonding pair. [default: 3.0]')
    # parser.add_argument('-f', '--f_format', type=str, default=None, help='Input file format. [default: Automatic if ASE readable.]')
    # parser.add_argument('-g', '--gsmear', type=float, default=0., help='Width(simga, STD) of Gaussian smearing in degree unit. Zero means no smearing. [default: 0]')
    # parser.add_argument('-a', '--no_average', dest='avg_bool', action='store_false', help='Not to take average over files. [default: take average]')
    # parser.add_argument('-s', '--dont_save', dest='save_bool', action='store_false', help='If provided, ADOS arrays will not be saved. Default: Save array')
    # parser.add_argument('-l', '--dont_load', dest='load_bool', action='store_false', help='If provided, ADOS arrays will not be loaded. Default: Load if possible')
    # parser.add_argument('-m', '--Nprocs', type=int, default=1, help='Number of process for multiprocessing. [Default: serial compute]')
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
    print('This code will give you comparison graph of two ADFs'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()

    ## Read Inputs
    fname_list   = args.adf_npz_list
    adf_npz_list = [np.load(fname_list[i]) for i in range(len(fname_list))]
    angle_arr    = adf_npz_list[0]['angd']
    adf_arr_list = [adf_npz_list[i]['agr'] for i in range(len(fname_list))]
    rectify_cut  = args.rectify_cut
    label_list   = args.label_list
    lstyle_list  = args.lstyle_list
    lcolor_list   = args.lcolor_list
    if label_list and len(label_list) != len(adf_npz_list):
        raise ValueError('len(label_list)=={}  != len(adf_npz_list)=={}'.format(len(label_list), len(adf_npz_list)))
    if not label_list:
        label_list = [None]*len(adf_npz_list)
    if lstyle_list and len(lstyle_list) != len(adf_npz_list):
        raise ValueError('len(lstyle_list)=={}  != len(adf_npz_list)=={}'.format(len(label_list), len(adf_npz_list)))
    if not lstyle_list:
        lstyle_list = [None]*len(adf_npz_list)
    if lcolor_list and len(lcolor_list) != len(adf_npz_list):
        raise ValueError('len(lcolor_list)=={}  != len(adf_npz_list)=={}'.format(len(label_list), len(adf_npz_list)))
    if not lcolor_list:
        lcolor_list = [None]*len(adf_npz_list)
    ## Prepare data
    chem, dDeg, rcut = [],[],[]
    for fname in fname_list:
        chem_tmp, dDeg_tmp, rcut_tmp = fname.split('_')[-4:-1]
        chem.append(chem_tmp); dDeg.append(dDeg_tmp); rcut.append(rcut_tmp)
        if chem_tmp != chem[0] or dDeg_tmp != dDeg[0] or rcut_tmp != rcut[0]:
            raise ValueError('Files are not proper to compare. Must have same "Chemical symbols, dDeg, cutoff radii".')
    chem = chem[0].split('-')
    dDeg = dDeg[0].split('-')[1]
    rcut = rcut[0].split('-')[1]

    ## Rectify curve
    if rectify_cut:
        from ss_util import rectify_curve
        adf_arr_list = [rectify_curve(adf_arr_list[i], rectify_cut) for i in range(len(adf_arr_list))]
    ## Gaussian smearing
    if args.gsmear:
        from scipy.ndimage.filters import gaussian_filter1d
        adf_arr_list = [gaussian_filter1d(adf_arr_list[i], args.gsmear /dDeg) for i in range(len(adf_arr_list))]

    ## Plot
    import matplotlib.pyplot as plt
    font = {'family':'Arial'}
    plt.rc('font', **font)
    if label_list:
        [plt.plot(angle_arr, adf_arr_list[i], label=label_list[i], linestyle=lstyle_list[i], c=lcolor_list[i]) for i in range(len(adf_arr_list))]
    else:
        plt.plot(angle_arr, adf_arr_list[0], c='k')
        [plt.plot(angle_arr, adf_arr_list[i], linestyle='dashed') for i in range(1,len(adf_arr_list))]
    if chem == ['a','a','a']:
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
    ## Legend options
    if label_list:
        plt.legend(loc='best', fontsize='large')
    
    plt.show()

