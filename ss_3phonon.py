#!/usr/bin/env python

        #####    Code by YJ Choi of CNPL, Dep. of Phys. POSTECH, Korea      #####
        #####                    ssrokyz@postech.ac.kr                      #####
        ##### Some codes are copied from phonopy. Code for Phono3py v1.12.0 #####

import numpy as np
from subprocess import call

def bu_and_mkdir(calc_dir, ndir):
    call(['rm -rf '+calc_dir+'/bu-'+ndir], shell=True)
    call(['mv '+calc_dir+'/'+ndir+' '+calc_dir+'/bu-'+ndir], shell=True)
    call(['mkdir -p '+calc_dir+'/'+ndir], shell=True)


