#!/usr/bin/env python

from asap3.analysis.rdf import RadialDistributionFunction as RDF
import sys
import matplotlib.pyplot as plt
import numpy as np

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
    print('This code will give you Raidial Distribution Function'.center(120))
    print('')
    print('Useage  ==> ./ase-rdf.py >atoms list file< >rMax< >nBins< (>rectify_cut<)'.center(120))
    print('Example ==> ./ase-rdf.py si.traj 12 500 0.2                              '.center(120))
    print('')
    print('Return ----------> rdf-(original name)-(rMax)-(nBins).npy'.center(120))
    print('')
    print('=================================================================================================='.center(120))

    inp_fname = sys.argv[1]
    rMax = float(sys.argv[2])
    nBins = int(sys.argv[3])
    if len(sys.argv) == 5:
        rectify_cut = float(sys.argv[4])
    else:
        rectify_cut = None
    fname = 'rdf-{}-{:.2f}-{}.npy'.format(inp_fname, rMax, nBins)
    try:
        curve = np.load(fname)
        print('File {} is loaded.'.format(fname).center(120))
    except:
        print('File {} not found. Calculation will be carried out.'.format(fname).center(120))
        ## Read inputs
        from ase.io import read
        alist = read(inp_fname, ':')
        RDFobj = RDF(atoms=alist[0],
                     rMax=rMax,
                     nBins=nBins)

        for i in range(1,len(alist)):
            RDFobj.atoms = alist[i]  # Fool RDFobj to use the new atoms
            RDFobj.update()           # Collect data
            if i % 1000 == 999:
                print('\t Updating '+str(i+1)+" th image's RDF")

        rdf = RDFobj.get_rdf()
        x = np.arange(nBins) * rMax / nBins

        ## Writing output
        curve = np.transpose(np.concatenate(([x], [rdf])))
        np.save(fname, curve)
        print('=================================================================================================='.center(120))
        print('Return ----------> {}'.format(fname).center(120))
        print('=================================================================================================='.center(120))

    ## Rectify curve
    if rectify_cut:
        iter = True
        while True:
            test = curve[1:] - curve[:-1]
            peak_bool = np.array(list(test[:,1] > (-1 * rectify_cut)) + [True], dtype=np.bool)
            if False not in peak_bool:
                break
            curve = curve[peak_bool]

    ## Plot
    plt.plot(curve[:-1, 0], curve[:-1, 1])
    plt.show()
