#!/usr/bin/env python
import numpy as np

def get_atom_inten_hist(
    alist,
    direc,
    nbins=200,
    ):
    """
    direc (array of shape (3,)): Directional vector along the axis of plot.
    """
    # Normalize
    direc = np.array(direc) /np.linalg.norm(direc)

    # pos.shape = (len(alist), len(atoms), 3)
    pos = []
    for i in range(len(alist)):
        pos.append(alist[i].get_positions())
    
    # dist.shape = (len(alist), len(atoms))
    dist = np.dot(pos, direc)
    dist_min = np.amin(dist)
    dist_max = np.amax(dist)
    
    # inten_hist.shape = (len(atoms), nbins)
    inten_hist = []
    for i in range(len(alist[0])):
        hist, bin_edges = np.histogram(
            dist[:,i],
            nbins,
            (dist_min, dist_max),
            )
        inten_hist.append(hist)

    return np.array(inten_hist), bin_edges

def get_peaks(
    hist,
    bin_centers,
    threshold,
    ):
    max_hei = np.amax(hist)
    #
    peak_pos = []
    peak_hei = []
    local_pos = bin_centers[0]
    local_hei = 0
    on_peak = False
    for i in range(len(hist)):
        if hist[i] /max_hei >= threshold:
            on_peak = True
            if hist[i] > local_hei:
                local_pos = bin_centers[i]
                local_hei = hist[i]
        else:
            if on_peak:
                peak_pos.append(local_pos)
                peak_hei.append(local_hei)
            on_peak = False
            local_hei = 0
    if on_peak:
        peak_pos.append(local_pos)
        peak_hei.append(local_hei)

    return np.array(peak_pos), np.array(peak_hei)

def argparse():
    import argparse
    parser = argparse.ArgumentParser(description = """
    Plot atomic intensity along a direction.
    """)
    # Positional arguments
    parser.add_argument('inp_file', type=str, help='ASE readable atoms list file name.')
    parser.add_argument('direc', type=float, nargs=3, help='Directional vector along the axis of plot.')
    # Optional arguments
    parser.add_argument('-n', '--image_slice', type=str, default=':', help='Image slice following python convention. [default=":"] e.g. -n :1000:10')
    parser.add_argument('-b', '--num_bins', type=int, default=200, help='Number of bins for histogram. [default=200]')
    parser.add_argument('-s', '--symbols', type=str, nargs='+', default=None, help='Specify chemical symbols for counting. [default: all species]')
    # parser.add_argument('-g', '--gsmear', type=float, default=0., help='Width(simga, STD) of Gaussian smearing in Angstrom unit. Zero means no smearing. [default: 0]')
    parser.add_argument('-t', '--peak_thres', type=float, default=0.1, help='Threshold for peak print/save. [default: 0.1 times maximum of partial histogram]')
    parser.add_argument('-c', '--color_list', type=str, nargs='+', default=None, help='Specify colors for partial histogram plots [default: auto].')
    parser.add_argument('-l', '--no_legend', action='store_false', dest='plot_legend', help='No legend plot [default: plot]')
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
    print('Plot atomic intensity along a direction.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    args = argparse()
    
    from ase.io import read
    alist = read(args.inp_file, args.image_slice)
    if not isinstance(alist, list):
        alist = [alist]

    inten_hist, bin_edges = get_atom_inten_hist(
        alist,
        args.direc,
        args.num_bins,
        )
    bin_centers = (bin_edges[0:-1] + bin_edges[1:]) /2.

    # Chemical sorting
    chem = np.array(alist[0].get_chemical_symbols())
    unique_chem = np.unique(chem)
    chem_mask = dict()
    for c in unique_chem:
        chem_mask[c] = chem == c

    #
    inten_hist_chem = dict()
    for c in unique_chem:
        inten_hist_chem[c] = np.sum(inten_hist[chem_mask[c]], axis=0)


    #
    if args.symbols:
        plot_chem = np.array(args.symbols)
    else:
        plot_chem = unique_chem

    #
    inten_hist_tot = np.sum([inten_hist_chem[c] for c in unique_chem], axis=0)

    # Peaks
    peak_pos_dict = dict()
    peak_hei_dict = dict()
    # Total
    peak_pos_dict['tot'], peak_hei_dict['tot'] = get_peaks(
        inten_hist_tot,
        bin_centers,
        args.peak_thres,
        )
    print("\n")
    print(" > Peaks for total histogram:") 
    print("  Ind  |  Peak position |  Interval |     Peak height   |  Height ratio to max")
    peak_pos_prev = 0.
    for i in range(len(peak_pos_dict['tot'])):
        print(' #{:3d}{:18.2f}{:12.2f}{:18d}{:18.2f} %'.format(
            i+1,
            peak_pos_dict['tot'][i],
            peak_pos_dict['tot'][i] -peak_pos_prev,
            peak_hei_dict['tot'][i],
            peak_hei_dict['tot'][i] /np.amax(peak_hei_dict['tot']) *100.,
            ))
        peak_pos_prev = peak_pos_dict['tot'][i]
    np.savez('peak_total.npz', pos=peak_pos_dict['tot'], hei=peak_hei_dict['tot'])
    # Partial
    for c in unique_chem:
        peak_pos_dict[c], peak_hei_dict[c] = get_peaks(
            inten_hist_chem[c],
            bin_centers,
            args.peak_thres,
            )
        print("\n")
        print(" > Peaks for {} histogram:".format(c)) 
        print("  Ind  |  Peak position |  Interval |     Peak height   |  Height ratio to max")
        peak_pos_prev = 0.
        for i in range(len(peak_pos_dict[c])):
            print(' #{:3d}{:18.2f}{:12.2f}{:18d}{:18.2f} %'.format(
                i+1,
                peak_pos_dict[c][i],
                peak_pos_dict[c][i] -peak_pos_prev,
                peak_hei_dict[c][i],
                peak_hei_dict[c][i] /np.amax(peak_hei_dict[c]) *100.,
                ))
            peak_pos_prev = peak_pos_dict[c][i]
        np.savez('peak_{}.npz'.format(c), pos=peak_pos_dict[c], hei=peak_hei_dict[c])

    #
    from matplotlib import pyplot as plt
    if not args.color_list:
        color_list = [None]*len(plot_chem)
    else:
        color_list = args.color_list
    for i in range(len(plot_chem)):
        plt.plot(inten_hist_chem[plot_chem[i]], bin_centers, c=color_list[i], label=plot_chem[i])
    if not args.symbols:
        plt.plot(inten_hist_tot, bin_centers, c='k')
    if args.plot_legend:
        plt.legend(fontsize='large').set_draggable(True)
    plt.tick_params(axis="both",direction="in", labelsize='x-large')
    plt.xticks([])
    plt.ylim(bin_edges[0], bin_edges[-1])
    plt.xlim(1e-4, None)
    plt.title('({}, {}, {}) direction'.format(*args.direc))
    plt.xlabel('Hist.', fontsize='x-large')
    plt.ylabel(r'Position ($\AA$)', fontsize='x-large')
    plt.subplots_adjust(left=0.35, right=0.65)
    plt.grid(alpha=0.4)
    plt.show()

