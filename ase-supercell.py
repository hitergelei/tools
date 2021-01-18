#!/usr/bin/env python
##### CODE BY YOUNG JAE CHOI #####
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Make supercell of ASE readable.
    useage ==> ./ase-supercell.py 'trajactory file' 'supercell integers'
               EXAMPLE) ./supercell.py structure.traj 2 2 1
    """)
    # Positional arguments
    parser.add_argument('traj_file', type=str, help='ASE readable traj file.')
    parser.add_argument('s_a1', type=int, help='Number of duplication for a1 direction.')
    parser.add_argument('s_a2', type=int, help='Number of duplication for a2 direction.')
    parser.add_argument('s_a3', type=int, help='Number of duplication for a3 direction.')
    #
    parser.add_argument('-n', '--traj_slice', type=str, default=':', help='ASE understandable slice format. default=":" (e.g.) -n :1000:10')
    parser.add_argument('-s', '--sequencial', dest='order_chem', action='store_false', help='If called, reordering of chemical symbols is skipped.')
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
    print('Make supercell of ASE readable'.center(120))
    print("useage ==> ./ase-supercell.py 'trajactory file' 'supercell integers' ".center(120))
    print("   e.g.) ./ase-supercell.py structure.traj 2 2 1 ".center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args       = argparse()
    traj_file  = args.traj_file
    s_a1       = args.s_a1
    s_a2       = args.s_a2
    s_a3       = args.s_a3
    traj_slice = args.traj_slice
    order_chem = args.order_chem

    from ase.io import read, write
    alist = read(traj_file, traj_slice)
    frames = len(alist)

    ncell = s_a1*s_a2*s_a3
    matrix = [[s_a1,0,0],
              [0,s_a2,0],
              [0,0,s_a3]]

    from ase.io.trajectory import Trajectory
    supertraj = Trajectory("supercell_"+str(s_a1)+"x"+str(s_a2)+"x"+str(s_a3)+"_"+traj_file, "w")

    print("I'll make supercell of '"+traj_file+"' file with "+str(frames)+" frames.")

    #
    from copy import deepcopy
    from ase.build import make_supercell
    from permute_sequence import permute_sequence
    supercells = []
    for i in range(frames):
        if i % 1000 == 999:
            print("Making supercell of "+str(i+1)+"th frame")
        atoms = alist[i]
        super = make_supercell(atoms, matrix)
        calc = atoms._calc
        if calc != None:
            scalc = deepcopy(calc)
            super.set_calculator(scalc)
            scalc.atoms = super
        
            try:
                energy=calc.results['energy']
            except KeyError:
                print("No energy info")
                energy=None
            except:
                print("Something goes wrong with energy info")
                raise ValueError
            else:
                scalc.results['energy']=energy*ncell
                
            try:
                forces=calc.results['forces']
            except KeyError:
                print("No forces info")
                forces=None
            except:
                print("Something goes wrong with forces info")
                raise ValueError
            else:
                newf = forces
                for j in range(ncell-1):
                    newf=np.concatenate((newf, forces))
                scalc.results['forces']=newf
            
            # # Not implemented yet.
            #try:
            #    stress=calc.results['stress']
            #except KeyError:
            #    print("No stress info")
            #    stress=None
            #except:
            #    print("Something goes wrong with stress info")
            #    raise ValueError
            #else:
            #    pass
        #print(calc.__dict__)
        #print(scalc.__dict__)
        supercells.append(super)
            
    #
    if order_chem:
        #
        from ss_util import ordered_unique
        chem = alist[0].get_chemical_symbols()
        unique_chem = ordered_unique(chem)
        #
        from ss_util import ind_dict
        chem_inds = ind_dict(supercells[0].get_chemical_symbols())
        chem_order = []
        for c in unique_chem:
            chem_order.append(chem_inds[c])
        chem_order = np.concatenate(chem_order)
        print(chem_order)
        #
        for i in range(len(supercells)):
            supercells[i] = permute_sequence(supercells[i], chem_order)
            supertraj.write(supercells[i])

    print("\n#############  Primitive cell  #############")
    print("Total number of atom = "+str(len(atoms)))
    print("Cell =")
    print(atoms.get_cell())

    print("\n#############  Supercell cell  #############")
    print("Transformation matrix, P ="),;print(matrix)
    print("")
    print("Total number of atom = "+str(len(super)))
    print("Cell =")
    print(super.get_cell())
    print("")
    print("pbc ="),;print(super.get_pbc())
    print("")

    supertraj.close()
    print("Supercell trajectory file :: 'supercell_"+str(s_a1)+"x"+str(s_a2)+"x"+str(s_a3)+"_"+traj_file+"'\n")

