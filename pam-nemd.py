#!/usr/bin/env python
import numpy as np

def argparse():
    import argparse
    from ss_util import parse_slice
    parser = argparse.ArgumentParser(description = """
    Young-Jae Choi, POSTECH, Korea, Rep. of.
    Calculate PAM of selected region.
    """)
    # Positional arguments
    parser.add_argument('alist_ase', type=str,
        help='ASE readable atoms-object list. NOTE) Positions should not be folded by boundaries.')
    parser.add_argument('r0_ase', type=str,
        help='Groud-state positions of the alist_ase file. ASE readable atoms-object. Should be equally unfolded with alist_ase file.')
    # # Optional arguments
    parser.add_argument('-b', '--bc', default=['n']*6, nargs=6,
        help='Boundary condition to calculate PAM. form: (x low  x up  y low  y up  z low  z up). "n" means no limit. [Default: Include all range.]')
    parser.add_argument('-n', '--img_slice', type=str, default=':', help='ASE understandable slice [Default: ":" ]')
    parser.add_argument('-t', '--T_grad', action='store_true', help='Print gradient T if provided.')

    return parser.parse_args()

def calc_pam(disp, mom, vol):
    """
         1      N        ___             
    J = --- * (sum) (r_i-r_i) X p_i 
         V      i                        

    N = len(atoms)

    INPUT
    disp [array, shape=(len(alist), len(atoms), 3)]: Displacements.
    mom [array, shape=(len(alist), len(atoms), 3)]: Momenta.
    vol [float]: volume of selected region.

    RETURN
    PAM per volume. Unit=( eV fs / A^3 )
    """

    # disp = posi - np.mean(posi, axis=0)
    return np.mean(np.sum(np.cross(disp, mom), axis=1), axis=0) /vol

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
    print('Calculate PAM of selected region.'.center(120))
    print('=================================================================================================='.center(120))
    print('')

    ## Argparse
    args = argparse()

    from ase.io import read
    alist = read(args.alist_ase, args.img_slice)
    if not isinstance(alist, list):
        alist = [alist]
    #
    posi = []
    mom = []
    # vel = []
    for atoms in alist:
        posi.append(atoms.get_positions())
        mom.append(atoms.get_momenta())
        # vel.append(atoms.get_velocities())
    posi = np.array(posi)
    mom = np.array(mom)
    # vel = np.array(vel)

    r0 = read(args.r0_ase).get_positions()
    disp = posi - r0

    print(
        # np.std(vel),
        # np.mean(vel, axis=0),
        "np.std(mom)=",
        np.std(mom),
        "\nnp.mean(np.mean(mom, axis=0), axis=0)=",
        np.mean(np.mean(mom, axis=0), axis=0),
        "\nnp.mean(mom, axis=0)=\n",
        np.mean(mom, axis=0),
        "\n\nnp.std(disp)=",
        np.std(disp),
        "\nnp.mean(np.mean(disp, axis=0), axis=0)=",
        np.mean(np.mean(disp, axis=0), axis=0),
        "\nnp.mean(disp, axis=0)=\n",
        np.mean(disp, axis=0),
        "\n",
        )

    # Masking
    mask = np.ones(r0.shape, dtype=int)
    if args.bc[0] != 'n':
        mask[:, 0] *= r0[:, 0] >= float(args.bc[0])
    if args.bc[1] != 'n':
        mask[:, 0] *= r0[:, 0] <= float(args.bc[1])
    if args.bc[2] != 'n':
        mask[:, 1] *= r0[:, 1] >= float(args.bc[2])
    if args.bc[3] != 'n':
        mask[:, 1] *= r0[:, 1] <= float(args.bc[3])
    if args.bc[4] != 'n':
        mask[:, 2] *= r0[:, 2] >= float(args.bc[4])
    if args.bc[5] != 'n':
        mask[:, 2] *= r0[:, 2] <= float(args.bc[5])
    mask = np.array(np.prod(mask, axis=1), dtype=bool)
    natoms = int(np.sum(mask))
    vol = alist[0].get_volume() /len(alist[0]) *natoms

    disp = disp[:, mask]
    mom = mom[:, mask]
    # from ase.visualize  import view
    # view(atoms[mask])

    J = calc_pam(disp, mom, vol)

    # Unit conversion ( eV fs / A^3 ) to ( J s / m^3 )
    from ase import units
    J *= units._e * 1e-15 / 1e-30
    
    print('PAM per volume = {} ( J s / m^3 )'.format(J))
