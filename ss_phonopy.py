def calc_vasp(phonon, verbose = False):
    """ Calculate Force Constant with Vasp """
    ################### all same from now 
    import subprocess as sp
    import numpy as np
    import pickle
    np.set_printoptions(threshold=np.nan)
        
    delta = np.linalg.norm(phonon.get_displacements()[0][1:4])
    directions = phonon.get_displacement_directions()
    N = phonon.get_supercell_matrix()[0][0]
    image_num = len(directions)
    tf = {0:'nonsym', 1:'sym'}
    sym = tf[phonon._is_symmetry]
    job_name = "x" + str(N) + "_d" + str(delta) + "_" + sym
    import sys
    if (sys.version_info > (3,0)):
        pwd = str((sp.check_output("pwd"))[:-1])[2:-1]
    else:
        pwd = str((sp.check_output("pwd"))[:-1])
    calc_dir = pwd + "/calcs/" + job_name
    sp.call(["rm -rf " + calc_dir + "/poscars"], shell = True)
    sp.call(["mkdir -p " + calc_dir + "/poscars"], shell = True)
    sp.call(["cp POSCAR-* SPOSCAR " + calc_dir + "/poscars"], shell = True)
    sp.call(["rm -rf POSCAR-* SPOSCAR"], shell = True)

    try:
        phonon = pickle.load(open("pickle-"+job_name+".p", "rb"))
        if phonon.get_force_constants() is None:
            raise ValueError
        ################### all same until now 
    except:
        print("******* There is no saved pickle file. ********".center(80))
        print("******* Vasp calc will be carried. ********".center(80))

        forces = []
        for i in range(image_num):
            ndir = "pos" + str(i+1).zfill(4) + "_atom" + str(directions[i][0]).zfill(4) + "_direc" + \
                   str(directions[i][1]) + str(directions[i][2]) + str(directions[i][3])
            sp.call(["rm", "-rf", calc_dir+"/BU-"+ndir])
            sp.call(["mv", calc_dir+"/"+ndir, calc_dir+"/BU-"+ndir])
            sp.call(["mkdir", "-p", calc_dir+"/"+ndir])
            sp.call(["cp", "INCAR", "POTCAR", "KPOINTS", calc_dir+"/poscars/POSCAR-"+str(i+1).zfill(3), calc_dir+"/"+ndir])
            sp.call(["cp POSCAR-"+str(i+1).zfill(3)+" POSCAR"], cwd = calc_dir+"/"+ndir, shell = True)
            if (sys.version_info > (3,0)):
                sp.call(["mpiexec.hydra -np $NSLOTS vasp_gpu > "+ pwd+"/calcs/"+ndir+"/out"], \
                         cwd = pwd+"/calcs/"+ndir, shell = True)
            else:
                sp.call(["mpiexec.hydra -np $NSLOTS vasp_std > "+ pwd+"/calcs/"+ndir+"/out"], \
                         cwd = pwd+"/calcs/"+ndir, shell = True)
            import io
            with io.open(calc_dir+"/"+ndir+"/vasprun.xml", "rb") as xmlfile:
                from phonopy.interface.vasp import VasprunxmlExpat
                vasprun = VasprunxmlExpat(xmlfile)
                if vasprun.parse():
                    force_now = (vasprun.get_forces()).tolist()
#                   epsilon = vasprun.get_epsilon()
#                   borns = vasprun.get_born()
#                   lattice = vasprun.get_lattice()[-1]
#                   points = vasprun.get_points()[-1]
#                   symbols = vasprun.get_symbols()
#                   unit_cell = PhonopyAtoms(symbols=symbols,
#                                            scaled_positions=points,
#                                            cell=lattice)
                    forces.extend(force_now)
                    if verbose:
                        print("force_now"+str(i+1)+"\n"+str(np.asarray(force_now)))
        forces_arr = np.asarray(forces)
        if verbose:
            print("\n\n"+"forces"+"\n"+str(forces_arr))
        phonon.set_forces(forces_arr)
        phonon.produce_force_constants()
    
        if verbose:
            print("\n\ndisplacement_dataset =>\n\n")
            print(phonon.get_displacement_dataset())
            print("\n\nforce_constants =>\n\n")
            print(str(2*N**3)+" x "+str(2*N**3)+" x "+str(3)+" x "+str(3)+" matrix\n\n")
            print(phonon.get_force_constants())
    
        pickle.dump(phonon, open("pickle-"+job_name+".p", "wb"))
    else:
        print("******* Pickle file has been loaded. ********".center(80))

    return phonon

#def get_fc(directions, delta):
#    """ Get forces, FCs """
#    import subprocess as sp
#
#    image_num = len(directions)
#    pwd = (sp.check_output("pwd"))[:-1]
#    for i in range(image_num):
#        ndir = str(delta) + "-" + str(directions[i][0]) + "-" + \
#               str(directions[i][1]) + str(directions[i][2]) + str(directions[i][3])

def calc_dpmd(phonon, verbose = False):
    """ Calculate Force Constant with DPMD """
    ################### all same from now 
    import subprocess as sp
    import numpy as np
    import pickle
    np.set_printoptions(threshold=np.nan)
        
    delta = np.linalg.norm(phonon.get_displacements()[0][1:4])
    directions = phonon.get_displacement_directions()
    N = phonon.get_supercell_matrix()[0][0]
    image_num = len(directions)
    tf = {0:'nonsym', 1:'sym'} 
    sym = tf[phonon._is_symmetry]
    job_name = "x" + str(N) + "_d" + str(delta) + "_" + sym
    import sys
    if (sys.version_info > (3,0)):
        pwd = str((sp.check_output("pwd"))[:-1])[2:-1]
    else:
        pwd = str((sp.check_output("pwd"))[:-1])
    calc_dir = pwd + "/calcs/" + job_name
    sp.call(["rm -rf " + calc_dir + "/poscars"], shell = True)
    sp.call(["mkdir -p " + calc_dir + "/poscars"], shell = True)
    sp.call(["cp POSCAR-* SPOSCAR " + calc_dir + "/poscars"], shell = True)
    sp.call(["rm -rf POSCAR-* SPOSCAR"], shell = True)

    try:
        phonon = pickle.load(open("pickle-"+job_name+".p", "rb"))
        if phonon.get_force_constants() is None:
            raise ValueError
        ################### all same until now 
    except:
        print("******* There is no saved pickle file. ********".center(80))
        print("******* DPMD calc will be carried. ********".center(80))

        forces = []
        for i in range(image_num):
            ndir = "pos" + str(i+1).zfill(4) + "_atom" + str(directions[i][0]).zfill(4) + "_direc" + \
                   str(directions[i][1]) + str(directions[i][2]) + str(directions[i][3])
            sp.call(["rm", "-rf", calc_dir+"/BU-"+ndir])
            sp.call(["mv", calc_dir+"/"+ndir, calc_dir+"/BU-"+ndir])
            sp.call(["mkdir", "-p", calc_dir+"/"+ndir])
            sp.call(["cp","frozen_model.pb", "input.in", calc_dir+"/poscars/POSCAR-"+str(i+1).zfill(3), calc_dir+"/"+ndir])
            sp.call(["lmp-pos2lmp.awk POSCAR-"+str(i+1).zfill(3)+" > structure.in"], cwd = calc_dir+"/"+ndir, shell = True)
            sp.call(["mpiexec.hydra -np $NSLOTS lmp_mpi -in input.in > out"], cwd = calc_dir+"/"+ndir, shell = True)
            from ase.io.lammpsrun import read_lammps_dump as read_dump
            atoms = read_dump(calc_dir+"/"+ndir+"/out.dump", index=0, order=True)
            #print(atoms.__dict__)
            force_now = [atoms.get_forces(apply_constraint=False).tolist()]
            forces.extend(force_now)
            if verbose:
                print("force_now"+str(i+1)+"\n"+str(np.asarray(force_now)))
        forces_arr = np.asarray(forces)
        if verbose:
            print("\n\n"+"forces"+"\n"+str(forces_arr))
        phonon.set_forces(forces_arr)
        phonon.produce_force_constants()
 
        if verbose:
            print("\n\ndisplacement_dataset =>\n\n")
            print(phonon.get_displacement_dataset())
            print("\n\nforce_constants =>\n\n")
            print(str(2*N**3)+" x "+str(2*N**3)+" x "+str(3)+" x "+str(3)+" matrix\n\n")
            print(phonon.get_force_constants())
 
        pickle.dump(phonon, open("pickle-"+job_name+".p", "wb"))
    else:
        print("******* Pickle file has been loaded. ********".center(80))

    return phonon

def calc_amp(phonon, nn, verbose = False):
    """ Calculate Force Constant with AMP """
    ################### all same from now 
    import subprocess as sp
    import numpy as np
    import pickle
    np.set_printoptions(threshold=np.nan)
        
    delta = np.linalg.norm(phonon.get_displacements()[0][1:4])
    directions = phonon.get_displacement_directions()
    N = phonon.get_supercell_matrix()[0][0]
    image_num = len(directions)
    tf = {0:'nonsym', 1:'sym'}
    sym = tf[phonon._is_symmetry]
    job_name = "x" + str(N) + "_d" + str(delta) + "_" + sym
    import sys
    if (sys.version_info > (3,0)):
        pwd = str((sp.check_output("pwd"))[:-1])[2:-1]
    else:
        pwd = str((sp.check_output("pwd"))[:-1])
    calc_dir = pwd + "/calcs/" + job_name
    sp.call(["rm -rf " + calc_dir + "/poscars"], shell = True)
    sp.call(["mkdir -p " + calc_dir + "/poscars"], shell = True)
    sp.call(["cp POSCAR-* SPOSCAR " + calc_dir + "/poscars"], shell = True)
    sp.call(["rm -rf POSCAR-* SPOSCAR"], shell = True)

    try:
        phonon = pickle.load(open("pickle-"+job_name+".p", "rb"))
        if phonon.get_force_constants() is None:
            raise ValueError
        ################### all same until now 
    except:
        print("******* There is no saved pickle file. ********".center(80))
        print("******* AMP calc will be carried. ********".center(80))

        forces = []
        for i in range(image_num):
            ndir = "pos" + str(i+1).zfill(4) + "_atom" + str(directions[i][0]).zfill(4) + "_direc" + \
                   str(directions[i][1]) + str(directions[i][2]) + str(directions[i][3])
            sp.call(["rm", "-rf", calc_dir+"/BU-"+ndir])
            sp.call(["mv", calc_dir+"/"+ndir, calc_dir+"/BU-"+ndir])
            sp.call(["mkdir", "-p", calc_dir+"/"+ndir])
            sp.call(["cp", calc_dir+"/poscars/POSCAR-"+str(i+1).zfill(3), calc_dir+"/"+ndir])

            ########### calculate forces & atomic energies with amp ############
            from ase.io import read
            atoms = read("calcs/"+ndir+"/POSCAR-"+str(i+1).zfill(3),
                         format = "vasp")
            atoms.set_pbc(True)
            from amp import Amp
            calc = Amp.load(nn)
            atoms.set_calculator(calc)
            #print(atoms.__dict__)

            #********** numerical force must precede ***********
            force_now = [(calc.calculate_numerical_forces(atoms)).tolist()] # alternative
            #force_now = [atoms.get_forces(apply_constraint=False).tolist()] # alternative
            #stress_now = calc.calculate_numerical_stress(atoms) # just in case
            #********** energy calculation ***********
            energy_now = atoms.get_potential_energy()
            energies_now = calc.get_atomic_potentials(atoms)
            #********** force information restore ***********
            atoms._calc.results['forces'] = np.asarray(force_now[0])
            if verbose:
                print(energies_now)
                print(atoms._calc.results['forces'])
            forces.extend(force_now)
            #********** traj file gen ***********
            from ase.io.trajectory import Trajectory
            traj = Trajectory(calc_dir+"/"+ndir+"/output.traj", "w")
            traj.write(atoms)
            traj.close()
            
            if verbose:
                print("force_now"+str(i+1)+"\n"+str(np.asarray(force_now)))
        forces_arr = np.asarray(forces)
        if verbose:
            print("\n\n"+"forces"+"\n"+str(forces_arr))
        phonon.set_forces(forces_arr)
        phonon.produce_force_constants()
 
        if verbose:
            print("\n\ndisplacement_dataset =>\n\n")
            print(phonon.get_displacement_dataset())
            print("\n\nforce_constants =>\n\n")
            print(str(2*N**3)+" x "+str(2*N**3)+" x "+str(3)+" x "+str(3)+" matrix\n\n")
            print(phonon.get_force_constants())
 
        pickle.dump(phonon, open("pickle-"+job_name+".p", "wb"))
    else:
        print("******* Pickle file has been loaded. ********".center(80))

    return phonon

def calc_amp_tf(phonon, nn, verbose = False):
    """ Calculate Force Constant with AMP with tensorflow """
    ################### all same from now 
    import subprocess as sp
    import numpy as np
    import pickle
    np.set_printoptions(threshold=np.nan)
        
    delta = np.linalg.norm(phonon.get_displacements()[0][1:4])
    directions = phonon.get_displacement_directions()
    N = phonon.get_supercell_matrix()[0][0]
    image_num = len(directions)
    tf = {0:'nonsym', 1:'sym'}
    sym = tf[phonon._is_symmetry]
    job_name = "x" + str(N) + "_d" + str(delta) + "_" + sym
    import sys
    if (sys.version_info > (3,0)):
        pwd = str((sp.check_output("pwd"))[:-1])[2:-1]
    else:
        pwd = str((sp.check_output("pwd"))[:-1])
    calc_dir = pwd + "/calcs/" + job_name
    sp.call(["rm -rf " + calc_dir + "/poscars"], shell = True)
    sp.call(["mkdir -p " + calc_dir + "/poscars"], shell = True)
    sp.call(["cp POSCAR-* SPOSCAR " + calc_dir + "/poscars"], shell = True)
    sp.call(["rm -rf POSCAR-* SPOSCAR"], shell = True)

    try:
        phonon = pickle.load(open("pickle-"+job_name+".p", "rb"))
        if phonon.get_force_constants() is None:
            raise ValueError
        ################### all same until now 
    except:
        print("******* There is no saved pickle file. ********".center(80))
        print("******* AMP calc will be carried. ********".center(80))

        forces = []
        for i in range(image_num):
            ndir = "pos" + str(i+1).zfill(4) + "_atom" + str(directions[i][0]).zfill(4) + "_direc" + \
                   str(directions[i][1]) + str(directions[i][2]) + str(directions[i][3])
            sp.call(["rm", "-rf", calc_dir+"/BU-"+ndir])
            sp.call(["mv", calc_dir+"/"+ndir, calc_dir+"/BU-"+ndir])
            sp.call(["mkdir", "-p", calc_dir+"/"+ndir])
            sp.call(["cp", calc_dir+"/poscars/POSCAR-"+str(i+1).zfill(3), calc_dir+"/"+ndir])

            ########### calculate forces & atomic energies with amp ############
            from ase.io import read
            atoms = read(calc_dir+"/"+ndir+"/POSCAR-"+str(i+1).zfill(3),
                         format = "vasp")
            atoms.set_pbc(True)
            calc = nn
            atoms.set_calculator(calc)
            #print(atoms.__dict__)

            #********** numerical force must precede ***********
            force_now = [(np.squeeze(calc.calculate_numerical_forces(atoms), axis=2)).tolist()] # alternative
            #force_now = [atoms.get_forces(apply_constraint=False).tolist()] # alternative
            #stress_now = calc.calculate_numerical_stress(atoms) # just in case
            #********** energy calculation ***********
            energy_now = atoms.get_potential_energy()
            energies_now = calc.get_atomic_potentials(atoms)
            #********** force information restore ***********
            atoms._calc.results['forces'] = np.asarray(force_now[0])
            if verbose:
                print(energies_now)
                print(atoms._calc.results['forces'])
            forces.extend(force_now)
            #********** traj file gen ***********
            from ase.io.trajectory import Trajectory
            traj = Trajectory(calc_dir+"/"+ndir+"/output.traj", "w")
            traj.write(atoms)
            traj.close()
            
            if verbose:
                print("force_now"+str(i+1)+"\n"+str(np.asarray(force_now)))
        forces_arr = np.asarray(forces)
        if verbose:
            print("\n\n"+"forces"+"\n"+str(forces_arr))
        phonon.set_forces(forces_arr)
        phonon.produce_force_constants()
 
        if verbose:
            print("\n\ndisplacement_dataset =>\n\n")
            print(phonon.get_displacement_dataset())
            print("\n\nforce_constants =>\n\n")
            print(str(2*N**3)+" x "+str(2*N**3)+" x "+str(3)+" x "+str(3)+" matrix\n\n")
            print(phonon.get_force_constants())
 
        pickle.dump(phonon, open("pickle-"+job_name+".p", "wb"))
    else:
        print("******* Pickle file has been loaded. ********".center(80))

    return phonon

def make_band(path, N_q):
    import numpy as np
    bands = []
    for i in range(len(path)):
        q_start  = np.array(path[i][0])
        q_end    = np.array(path[i][1])
        band = []
        for j in range(N_q+1):
            band.append(q_start + (q_end - q_start) / N_q * j)
        bands.append(band)
    return bands

def plot_band(phonon, labels=None):
    import numpy as np
    import matplotlib.pyplot as plt
    if labels:
        from matplotlib import rc
        rc('font',**{'family':'serif','sans-serif':['Times']})
        rc('text', usetex=False)

    fig, ax = plt.subplots()
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_tick_params(which='both', direction='in')
    ax.yaxis.set_tick_params(which='both', direction='in')

    phonon._band_structure.plot(plt, labels=labels)
    return plt

def plot_band_and_dos(phonon, pdos_indices=None, labels=None):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    if labels:
        from matplotlib import rc
        rc('font',**{'family':'serif','sans-serif':['Times']})
        rc('text', usetex=False)


    plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1])
    ax1 = plt.subplot(gs[0, 0])
    ax1.xaxis.set_ticks_position('both')
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_tick_params(which='both', direction='in')
    ax1.yaxis.set_tick_params(which='both', direction='in')
    phonon._band_structure.plot(plt, labels=labels)
    plt.ylabel('Frequency(meV)', fontsize=18)
    plt.xlabel('')
    plt.grid(True)
    plt.title('Phonon dispersion', fontsize=20)
    plt.yticks(np.arange(0, 100, 10), fontsize=16)

    ax2 = plt.subplot(gs[0, 1], sharey=ax1)
    ax2.xaxis.set_ticks_position('both')
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_tick_params(which='both', direction='in')
    ax2.yaxis.set_tick_params(which='both', direction='in')
    plt.subplots_adjust(wspace=0.08)
    plt.setp(ax2.get_yticklabels(), visible=False)

    if pdos_indices is None:
        phonon._total_dos.plot(plt,
                               ylabel="",
                               draw_grid=True,
                               flip_xy=True)
    else:
        phonon._pdos.plot(plt,
                          indices=pdos_indices,
                          ylabel="",
                          draw_grid=True,
                          flip_xy=True)

    ax2.set_xlim((0, None))
    plt.title('DOS', fontsize=20)
    plt.xlabel('')
    plt.xticks([])

    return plt
