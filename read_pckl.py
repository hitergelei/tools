#!/usr/bin/env python

import cPickle as pckl
import sys
import os
import numpy as np

def read_cew(fname):
    """read cew file and return as list"""
    load = pckl.load(open(fname, "rb"))
    return load

if __name__ == '__main__':

    print("\n\n#######################################################################################\n")
    print("      %%%%%%%%%%% This code will give you Raidial Distribution Function %%%%%%%%%")
    print("useage ==> ./ase-rdf.py 'trajectory file' 'rMax' 'nBins'")
    print("           EXAMPLE) ./ase-rdf.py si.traj 12 500")
    print("                    OUTPUT file is 'rdf-si.traj.log'")
    print("#######################################################################################")
    print("The Number of cew-pckl files is "+str(len(sys.argv)-1)+".\n\n")
    if len(sys.argv) <= 1:
        print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
        sys.exit(1)
    np.set_printoptions(threshold=np.nan)
    for fname in sys.argv[1:]:
        print("\n\nfor file "+fname)
        output = read_cew(fname)
        print(output)



