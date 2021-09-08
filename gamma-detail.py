#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Check phonopy gamma_detail file.
    """)
    # Positional arguments
    parser.add_argument('pho3_pckl', type=str, help='Phono3py pickle file returned with gamma_detail hdf5 file.')
    parser.add_argument('grid_address', type=int, nargs=3, help='Three integers of grid address.')
    # # Optional arguments

    return parser.parse_args()

if __name__ == '__main__':
    ## Intro
    import datetime
    now = datetime.datetime.now()
    time = now.strftime('%Y-%m-%d %H:%M:%S')
    print('')
    print('>>>>> Code by Young Jae Choi @ Phys. Dep. of POSTECH in Korea <<<<<'.center(120))
    print(('Code runtime : '+time).center(120))
    print('')
    print('=================================================================================================='.center(120))
    print('Check phonopy gamma_detail file.'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()
    ga = args.grid_address
    pho3_pckl = args.pho3_pckl

    # Main

    import pickle as pckl
    with open(pho3_pckl, 'rb') as f:
        tc = pckl.load(f).thermal_conductivity

    grid_address = tc._grid_address
    ir_grid_points = tc.get_grid_points()
    grid_map = tc._grid_mapping_table
    # gp is an index for all grid_points and gi is for ir_grid_points.
    gp = grid_address.tolist().index(ga)
    gi = ir_grid_points.tolist().index(gp)
    temp = tc.get_temperatures()
    mesh = tc.get_mesh_numbers()
    # freq.shape = (len(q_points), len(band))
    freq = tc.get_frequencies()
    # q_points.shape = (len(grid_address))
    q_points = grid_address / mesh

    # gam_i = tc.gamma[0,0,gp]


    import h5py
    gam_hdf5 = 'gamma_detail-m{}-g{}.hdf5'.format(pho3_pckl.split('-')[-3][2:], gp)
    gam_db = h5py.File(gam_hdf5, 'r')
    # gam.shape = (len(T), len(triplet), len(band), len(band), len(band))
    gam_ir = gam_db['gamma_detail'][()]
    # gam_wei.shape = (len(triplet))
    gam_wei = gam_db['weight'][()]
    # trip.shape = (len(triplet), 3)
    trip = gam_db['triplet'][()]

    # gam_i.shape = (len(T), len(triplet), len(band), len(band), len(band))
    gam_i = gam_ir * np.expand_dims(gam_wei, [0, 2, 3, 4])

    # Print
    from ss_util import bcolors
    print("""\n=======================================================================\n
        > Mesh: {}x{}x{}
        > Grid address: ({}, {}, {}) -> q = ({:.5f}, {:.5f}, {:.5f})
        > Frequencies: {}
        """.format(
        *mesh,
        *ga, *q_points[gp],
        freq[gi],
        ))
    for i in range(len(temp)):
        print(
            """
        > T = {}K
        > Linewidths: {}
            """.format(
            temp[i],
            np.sum(gam_i, (1, 3, 4))[i],
            ))

    band_ind = int(input('\n***Enter band index to search: \n(NOTE! python numbering! i.e. starts from 0, 1, 2,...)\n> '))
    crit_width = float(input('\n***Enter the critical linewidth to search: \n(Info shows if linewidth is larger than this.\n> '))

    print("""
        > Selected Band: #{}   ({}-th band)
        > Critical linewidth: {:e}
        """.format(
        band_ind, band_ind+1,
        crit_width,
        ))

    # inds.shape = (len(search results), 4)
    # inds[i] components = (T_ind, triplet_ind, band2_ind, band3_ind)
    inds = np.transpose(np.where(gam_i[:, :, band_ind] > crit_width))
    gam_sel = []
    for i in range(len(inds)):
        gam_sel.append(gam_i[inds[i,0], inds[i,1], band_ind, inds[i,2], inds[i,3]])
    gam_sel = np.array(gam_sel)
    decending_order = np.argsort(gam_sel)[::-1]

    for i in decending_order:
        trip_i = trip[inds[i,1]]
        trip_ir = [0, 0, 0]
        for j in range(3):
            trip_ir[j] = ir_grid_points.tolist().index(grid_map[trip_i[j]])
        print("""      #{}: {}gamma_lambda_1({}K) = {:.5e} THz{}
            {}Phonon 1{}: band #{}, q_1=({:.5f}, {:.5f}, {:.5f}), w={:.5f}THz
            {}Phonon 2{}: band #{}, q_2=({:.5f}, {:.5f}, {:.5f}), w={:.5f}THz
            {}Phonon 3{}: band #{}, q_3=({:.5f}, {:.5f}, {:.5f}), w={:.5f}THz
            """.format(
            i, bcolors.okgreen, temp[inds[i][0]], gam_sel[i], bcolors.endc,
            bcolors.okblue, bcolors.endc, band_ind,  *q_points[trip_i[0]], freq[trip_ir[0], band_ind ],
            bcolors.okblue, bcolors.endc, inds[i,2], *q_points[trip_i[1]], freq[trip_ir[1], inds[i,2]],
            bcolors.okblue, bcolors.endc, inds[i,3], *q_points[trip_i[2]], freq[trip_ir[2], inds[i,3]],
            ))
