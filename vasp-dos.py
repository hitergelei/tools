#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    This code will plot (total/atomic/chemical) DOS figure.
    """)
    # Positional arguments
        # -Blank-
    # Optional arguments
    parser.add_argument('-d', '--doscar_path', type=str, default='DOSCAR', help='DOSCAR file path. Default: ./DOSCAR')
    parser.add_argument('-p', '--poscar_path', type=str, default='POSCAR', help='POSCAR file path. Default: ./POSCAR')
    parser.add_argument('-c', '--cdos_bool', action='store_true', help='If provided, Plot chemical-DOS. Even when not provided, will be "True" if cdos_sym is provided.')
    parser.add_argument('-m', '--orbit_mom', action='store_true', help='If provided, DOS of different momentum are split.')
    parser.add_argument('-s', '--cdos_sym', type=str, nargs='+', default=None, help='Chemical symbols for chemical DOS plot. Multiple symbols are readable. Default: Plot all')
    parser.add_argument('-a', '--ados_bool', action='store_true', help='If provided, Plot atomic-DOS. Even when not provided, will be "True" if ados_ind is provided.')
    parser.add_argument('-n', '--ados_ind', type=int, nargs='+', default=None, help='Atomic indices for atomic DOS plot. Multiple integers are readable. Default: Plot all')
    parser.add_argument('-o', '--pdos_orbit', type=str, nargs='+', default=('s','p','d'), help='Orbitals for chemical DOS plot. Multiple orbitals are readable. Available choices: s, p, d or sum. Default: s,p,d')
    parser.add_argument('-r', '--no_pdos_residue', dest='pdos_residue', action='store_false', help='If provided, PDOS residue will not be plotted')
    parser.add_argument('-f', '--dont_flip_xy', dest='flip_xy', action='store_false', help='If provided, DOS will not be flipped.')
    parser.add_argument('-e', '--no_legend', dest='legend_bool', action='store_false', help='If provided, PDOS legend will not be plotted.')
    parser.add_argument('-u', '--Elim_up', type=float, default=None, help='Upper limit of energy range. Default: automatic')
    parser.add_argument('-l', '--Elim_low', type=float, default=None, help='Lower limit of energy range. Default: automatic')
    parser.add_argument('-v', '--DOSlim_up', type=float, default=None, help='Upper limit of DOS range. Default: automatic')
    parser.add_argument('-i', '--DOSlim_low', type=float, default=None, help='Lower limit of DOS range. Default: automatic')
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
    print('This code will plot (total/atomic/chemical) DOS figure.'.center(120))
    print('=================================================================================================='.center(120))
    print('')
    ## Argparse
    args         = argparse()
    doscar_path  = args.doscar_path
    poscar_path  = args.poscar_path
    cdos_sym     = args.cdos_sym
    cdos_bool    = True if cdos_sym else args.cdos_bool
    ados_ind     = args.ados_ind
    ados_bool    = True if ados_ind else args.ados_bool
    pdos_orbit   = args.pdos_orbit
    pdos_residue = args.pdos_residue
    flip_xy      = args.flip_xy
    legend_bool  = args.legend_bool
    Elim_up      = args.Elim_up
    Elim_low     = args.Elim_low
    DOSlim_up    = args.DOSlim_up
    DOSlim_low   = args.DOSlim_low
    orbit_mom    = args.orbit_mom

# Main
if __name__ == '__main__':
    with open(doscar_path, 'r') as f:
        # First line
        words = f.readline().split()
        natoms = int(words[0])
        if cdos_bool:
            if int(words[2]) is not 1:
                raise ValueError('DOSCAR file is not including Partial-DOS information. Please check.')
        for i in range(4):
            f.readline()
        # line 6
        words = f.readline().split()
        Ef   = float(words[2])
        Emax = float(words[0]) - Ef
        Emin = float(words[1]) - Ef
        if not Elim_up:
            Elim_up = Emax
        if not Elim_low:
            Elim_low = Emin

        # Get nE
        nE=0
        while 1:
            words = f.readline().split()
            if len(words) == 3:
                nE+=1
            else:
                break
        dE = (Emax-Emin) / (nE-1)
        E_arr = np.arange(Emin, Emax+1e-12, dE)

        # Check if it is non-collinear calc.
        words = f.readline().split()
        if len(words) == 37:
            collinear_bool = False
        elif len(words) == 10:
            collinear_bool = True
        else:
            raise NotImplementedError

        # Return to line 6
        f.seek(0)
        for i in range(6):
            f.readline()

        ## Total-DOS
        dos_arr = []
        for i in range(nE):
            dos_arr.append(float(f.readline().split()[1]))
        dos_arr = np.array(dos_arr)
        # DOS range
        if not DOSlim_up:
            DOSlim_up = np.amax(dos_arr[(E_arr <= Elim_up) * (E_arr >= Elim_low)])*1.1
        if not DOSlim_low:
            DOSlim_low = 0.

        # Plot
        from matplotlib import pyplot as plt
        font = {'family':'sans-serif', 'sans-serif':'Arial'}
        plt.rc('font', **font)
        if not cdos_bool and not ados_bool:
            # Plot
            if flip_xy:
                plt.plot(dos_arr, E_arr)
            else:
                plt.plot(E_arr, dos_arr)
        elif cdos_bool and ados_bool:
            raise ValueError('Atomic-PDOS and chemical-PDOS both enabled. Please switch one at a time.')
        else:
            # Plot total-DOS
            if flip_xy:
                # plt.fill_between(dos_arr, E_arr, color='k', alpha=0.3)
                plt.plot(dos_arr, E_arr, color='k')
            else:
                # plt.fill_betweenx(dos_arr, E_arr, color='k', alpha=0.3)
                plt.plot(E_arr, dos_arr, color='k')

            # atomic-PDOS
            atomic_pdos_arr = []
            for i in range(natoms):
                f.readline() # Skip first line
                atomic_pdos_arr_i = []
                for j in range(nE):
                    atomic_pdos_arr_i.append(np.array(f.readline().split()[1:], dtype=float))
                atomic_pdos_arr.append(atomic_pdos_arr_i)
            atomic_pdos_arr = np.array(atomic_pdos_arr, dtype=float)
            if not collinear_bool:
                atomic_pdos_arr = atomic_pdos_arr[:,:,range(0,36,4)]
            # Transpose (atom, E, orbit) to be (atom, orbit, E)
            atomic_pdos_arr = np.transpose(atomic_pdos_arr, (0,2,1))
            ## Plot atomic-PDOS
            if ados_bool:
                if not ados_ind:
                    ados_ind = np.arange(natoms)
                for atom_ind in ados_ind:
                    if 'sum' in pdos_orbit:
                        if flip_xy:
                            plt.plot(np.sum(atomic_pdos_arr[atom_ind], axis=0), E_arr, label='atom{}-(s+p+d)'.format(atom_ind))
                        else:
                            plt.plot(E_arr, np.sum(atomic_pdos_arr[atom_ind], axis=0), label='atom{}-(s+p+d)'.format(atom_ind))
                    else:
                        # s-orbital
                        if 's' in pdos_orbit:
                            if flip_xy:
                                plt.plot(atomic_pdos_arr[atom_ind, 0], E_arr, label='atom{}-s'.format(atom_ind))
                            else:
                                plt.plot(E_arr, atomic_pdos_arr[atom_ind, 0], label='atom{}-s'.format(atom_ind))
                        # p-orbital
                        if 'p' in pdos_orbit:
                            if flip_xy:
                                if orbit_mom:
                                    plt.plot(atomic_pdos_arr[atom_ind, 1], E_arr, label='atom{}-py'.format(atom_ind))
                                    plt.plot(atomic_pdos_arr[atom_ind, 2], E_arr, label='atom{}-pz'.format(atom_ind))
                                    plt.plot(atomic_pdos_arr[atom_ind, 3], E_arr, label='atom{}-px'.format(atom_ind))
                                else:
                                    plt.plot(np.sum(atomic_pdos_arr[atom_ind, 1:4], axis=0), E_arr, label='atom{}-p'.format(atom_ind))
                            else:
                                if orbit_mom:
                                    plt.plot(E_arr, atomic_pdos_arr[atom_ind, 1], label='atom{}-py'.format(atom_ind))
                                    plt.plot(E_arr, atomic_pdos_arr[atom_ind, 2], label='atom{}-pz'.format(atom_ind))
                                    plt.plot(E_arr, atomic_pdos_arr[atom_ind, 3], label='atom{}-px'.format(atom_ind))
                                else:
                                    plt.plot(E_arr, np.sum(atomic_pdos_arr[atom_ind, 1:4], axis=0), label='atom{}-p'.format(atom_ind))
                        # d-orbital
                        if 'd' in pdos_orbit:
                            if flip_xy:
                                if orbit_mom:
                                    plt.plot(atomic_pdos_arr[atom_ind, 4], E_arr, label='atom{}-dxy'.format(atom_ind))
                                    plt.plot(atomic_pdos_arr[atom_ind, 5], E_arr, label='atom{}-dyz'.format(atom_ind))
                                    plt.plot(atomic_pdos_arr[atom_ind, 6], E_arr, label='atom{}-d$z^2$'.format(atom_ind))
                                    plt.plot(atomic_pdos_arr[atom_ind, 7], E_arr, label='atom{}-dxz'.format(atom_ind))
                                    plt.plot(atomic_pdos_arr[atom_ind, 8], E_arr, label='atom{}-d$x^2$-$y^2$'.format(atom_ind))
                                else:
                                    plt.plot(np.sum(atomic_pdos_arr[atom_ind, 4:9], axis=0), E_arr, label='atom{}-d'.format(atom_ind))
                            else:
                                if orbit_mom:
                                    plt.plot(E_arr, atomic_pdos_arr[atom_ind, 4], label='atom{}-dxy'.format(atom_ind))
                                    plt.plot(E_arr, atomic_pdos_arr[atom_ind, 5], label='atom{}-dyz'.format(atom_ind))
                                    plt.plot(E_arr, atomic_pdos_arr[atom_ind, 6], label='atom{}-d$z^2$'.format(atom_ind))
                                    plt.plot(E_arr, atomic_pdos_arr[atom_ind, 7], label='atom{}-dxz'.format(atom_ind))
                                    plt.plot(E_arr, atomic_pdos_arr[atom_ind, 8], label='atom{}-d$x^2$-$y^2$'.format(atom_ind))
                                else:
                                    plt.plot(E_arr, np.sum(atomic_pdos_arr[atom_ind, 4:9], axis=0), label='atom{}-d'.format(atom_ind))
            ## Plot chemical-PDOS
            elif cdos_bool:
                # Get species sequence
                from ase.io import read
                atoms = read(poscar_path)
                chemical_symbols = atoms.get_chemical_symbols()
                if len(chemical_symbols) != natoms:
                    raise ValueError('Numbers of atoms for POSCAR(->{:d}) and DOSCAR(->{:d}) are not consistent.'.format(len(chemical_symbols), natoms))
                from ss_util import get_chem_ind_arr
                unique_chem, pdos_indices = get_chem_ind_arr(chemical_symbols)
                # Gather along chem
                chem_pdos_arr = []
                for spec_i in range(len(unique_chem)):
                    chem_pdos_arr.append(np.sum(atomic_pdos_arr[pdos_indices[spec_i]], axis=0))
                chem_pdos_arr = np.array(chem_pdos_arr, dtype=float)
                # Plot
                if not cdos_sym:
                    cdos_sym = unique_chem
                for spec_i in range(len(cdos_sym)):
                    spec_i = list(unique_chem).index(cdos_sym[spec_i])
                    if 'sum' in pdos_orbit:
                        if flip_xy:
                            plt.plot(np.sum(chem_pdos_arr[spec_i], axis=0), E_arr, label=unique_chem[spec_i]+'-(s+p+d)')
                        else:
                            plt.plot(E_arr, np.sum(chem_pdos_arr[spec_i], axis=0), label=unique_chem[spec_i]+'-(s+p+d)')
                    else:
                        # s-orbital
                        if 's' in pdos_orbit:
                            if flip_xy:
                                plt.plot(chem_pdos_arr[spec_i, 0], E_arr, label=unique_chem[spec_i]+'-s')
                            else:
                                plt.plot(E_arr, chem_pdos_arr[spec_i, 0], label=unique_chem[spec_i]+'-s')
                        # p-orbital
                        if 'p' in pdos_orbit:
                            if flip_xy:
                                if orbit_mom:
                                    plt.plot(chem_pdos_arr[spec_i, 1], E_arr, label=unique_chem[spec_i]+'-py')
                                    plt.plot(chem_pdos_arr[spec_i, 2], E_arr, label=unique_chem[spec_i]+'-pz')
                                    plt.plot(chem_pdos_arr[spec_i, 3], E_arr, label=unique_chem[spec_i]+'-px')
                                else:
                                    plt.plot(np.sum(chem_pdos_arr[spec_i, 1:4], axis=0), E_arr, label=unique_chem[spec_i]+'-p')
                            else:
                                if orbit_mom:
                                    plt.plot(E_arr, chem_pdos_arr[spec_i, 1], label=unique_chem[spec_i]+'-py')
                                    plt.plot(E_arr, chem_pdos_arr[spec_i, 2], label=unique_chem[spec_i]+'-pz')
                                    plt.plot(E_arr, chem_pdos_arr[spec_i, 3], label=unique_chem[spec_i]+'-px')
                                else:
                                    plt.plot(E_arr, np.sum(chem_pdos_arr[spec_i, 1:4], axis=0), label=unique_chem[spec_i]+'-p')
                        # d-orbital
                        if 'd' in pdos_orbit:
                            if flip_xy:
                                if orbit_mom:
                                    plt.plot(chem_pdos_arr[spec_i, 4], E_arr, label=unique_chem[spec_i]+'-dxy')
                                    plt.plot(chem_pdos_arr[spec_i, 5], E_arr, label=unique_chem[spec_i]+'-dyz')
                                    plt.plot(chem_pdos_arr[spec_i, 6], E_arr, label=unique_chem[spec_i]+'-d$z^2$')
                                    plt.plot(chem_pdos_arr[spec_i, 7], E_arr, label=unique_chem[spec_i]+'-dxz')
                                    plt.plot(chem_pdos_arr[spec_i, 8], E_arr, label=unique_chem[spec_i]+'-d$x^2$-$y^2$')
                                else:
                                    plt.plot(np.sum(chem_pdos_arr[spec_i, 4:9], axis=0), E_arr, label=unique_chem[spec_i]+'-d')
                            else:
                                if orbit_mom:
                                    plt.plot(E_arr, chem_pdos_arr[spec_i, 4], label=unique_chem[spec_i]+'-dxy')
                                    plt.plot(E_arr, chem_pdos_arr[spec_i, 5], label=unique_chem[spec_i]+'-dyz')
                                    plt.plot(E_arr, chem_pdos_arr[spec_i, 6], label=unique_chem[spec_i]+'-d$z^2$')
                                    plt.plot(E_arr, chem_pdos_arr[spec_i, 7], label=unique_chem[spec_i]+'-dxz')
                                    plt.plot(E_arr, chem_pdos_arr[spec_i, 8], label=unique_chem[spec_i]+'-d$x^2$-$y^2$')
                                else:
                                    plt.plot(E_arr, np.sum(chem_pdos_arr[spec_i, 4:9], axis=0), label=unique_chem[spec_i]+'-d')
            if pdos_residue:
                if flip_xy:
                    plt.plot(dos_arr - np.sum(np.sum(atomic_pdos_arr, axis=1), axis=0), E_arr, label='residue')
                else:
                    plt.plot(E_arr, dos_arr - np.sum(np.sum(atomic_pdos_arr, axis=1), axis=0), label='residue')
            if legend_bool:
                plt.legend(loc=(1.05, 0.), fontsize='x-large')
    #
    if flip_xy:
        plt.xlim((DOSlim_low, DOSlim_up))
        plt.ylim((Elim_low, Elim_up))
        plt.xlabel('DOS(arb. units)', fontsize='xx-large')
        plt.ylabel('Energy(eV)', fontsize='xx-large')
    else:
        plt.xlim((Elim_low, Elim_up))
        plt.ylim((DOSlim_low, DOSlim_up))
        plt.xlabel('Energy(eV)', fontsize='xx-large')
        plt.ylabel('DOS(arb. units)', fontsize='xx-large')
    plt.tick_params(axis="both",direction="in", labelsize='xx-large')
    plt.subplots_adjust(left=0.20, bottom=0.15, right=0.50, top=0.95)
    plt.grid(alpha=0.2)
    plt.show()

