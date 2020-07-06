import pyprocar

pyprocar.bandsplot(
    'PROCAR',
    outcar   = 'OUTCAR',
    elimit   = [-14,8],
    kticks   = list(range(0,400,100))+[399],
    # knames   = [u'$L$',u'$\\Gamma$',u'$X$',u'$U|K$',u'$\\Gamma$'],
    cmap     = 'Reds',
    # mode     = 'scatter',
    # mode     = 'atomic',
    mode     = 'parametric',
    # markersize = 200,
    atoms    = [1],
    # orbitals = [0],
    # orbitals = [1],
    # orbitals = [2],
    # orbitals = [3],
    # orbitals = [1,2,3],
    # orbitals = [4,5,6,7,8],
    orbitals = [0,1,2,3,4,5,6,7,8],
    kpointsfile  = 'KPOINTS',
    )
