        ##### Code by YJ Choi of CNPL, Dep. of Phys. POSTECH, Korea #####
        #####                 ssrokyz@postech.ac.kr                 #####
        ##### Some codes are copied from phonopy

import numpy as np

def calc_vasp(phonon, verbose = False, acoustic_sum_rule = False):
    """ Calculate Force Constant with Vasp """
    #>>>>>>>>>>>>>>>>>>> all same from now <<<<<<<<<<<<<<<<<<
    import sys
    import subprocess as sp
    import pickle
    np.set_printoptions(threshold=np.nan)
        
    delta      = np.linalg.norm(phonon.get_displacements()[0][1:4])
    directions = phonon.get_displacement_directions()
    N1         = phonon.get_supercell_matrix()[0][0]
    N2         = phonon.get_supercell_matrix()[1][1]
    N3         = phonon.get_supercell_matrix()[2][2]
    image_num  = len(directions)
    sym_dict   = {0:'nonsym', 1:'sym'}
    sym        = sym_dict[phonon._is_symmetry]
    job_name   = "x" + str(N1) + str(N2) + str(N3) + "_d" + str("%.3f" % delta) + "_" + sym
    pckl_name  = 'pickle-'+job_name+'.p'
    import os
    pwd = os.getcwd()
    calc_dir = pwd + "/calcs/" + job_name
    sp.call(["rm -rf " + calc_dir + "/POSCARS"], shell = True)
    sp.call(["mkdir -p " + calc_dir + "/POSCARS"], shell = True)
    sp.call(['cp SPOSCAR POSCAR-000'], shell=True)
    sp.call(["mv POSCAR-* SPOSCAR " + calc_dir + "/POSCARS"], shell = True)

    try:
        phonon = pickle.load(open(pckl_name, 'rb'))
        if phonon.get_force_constants() is None:
            raise ValueError
    except:
        print("******* There is no saved pickle file. ********".center(120))
        do_calc = True
    else:
        print("******* Pickle file has been loaded. ********".center(120))
        do_calc = False
    #>>>>>>>>>>>>>>>>>>> all same until now <<<<<<<<<<<<<<<<<<
    if do_calc:
        print("******* VASP calc will be carried out. ********".center(120))
        forces = []
        for i in range(image_num):
            print(' >>> Starting {:d}-th image calculation <<< '.format(i+1).center(120))
            ndir = "pos" + str(i+1).zfill(4) + "_atom" + str(directions[i][0]).zfill(4) + \
                "_direc" + str(directions[i][1]) + str(directions[i][2]) + str(directions[i][3])
            sp.call(["rm", "-rf", calc_dir+"/BU-"+ndir])
            sp.call(["mv", calc_dir+"/"+ndir, calc_dir+"/BU-"+ndir])
            sp.call(["mkdir", "-p", calc_dir+"/"+ndir])
            sp.call(["cp", "INCAR", "POTCAR", "KPOINTS", 
                calc_dir+"/POSCARS/POSCAR-"+str(i+1).zfill(3), calc_dir+"/"+ndir])
            sp.call(["cp POSCAR-"+str(i+1).zfill(3)+" POSCAR"], 
                cwd = calc_dir+"/"+ndir, shell = True)
            if (sys.version_info > (3,0)):
                sp.call(["mpiexec.hydra -np $NSLOTS vasp_gpu > "+ calc_dir+"/"+ndir+"/out"], \
                         cwd = calc_dir+"/"+ndir, shell = True)
            else:
                sp.call(["mpiexec.hydra -np $NSLOTS vasp_std > "+ calc_dir+"/"+ndir+"/out"], \
                         cwd = calc_dir+"/"+ndir, shell = True)
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
        if acoustic_sum_rule:
            phonon.symmetrize_force_constants()
    
        if verbose:
            print("\n\ndisplacement_dataset =>\n\n")
            print(phonon.get_displacement_dataset())
            print("\n\nforce_constants =>\n\n")
            print(phonon.get_force_constants())
    
        pickle.dump(phonon, open(pckl_name, "wb"))
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

def calc_dpmd(phonon, acoustic_sum_rule=False, F_0_correction=True, verbose=False):
    """ Calculate Force Constant with DPMD """
    #>>>>>>>>>>>>>>>>>>> all same from now <<<<<<<<<<<<<<<<<<
    import sys
    import subprocess as sp
    import pickle
    np.set_printoptions(threshold=np.nan)
        
    delta      = np.linalg.norm(phonon.get_displacements()[0][1:4])
    directions = phonon.get_displacement_directions()
    N1         = phonon.get_supercell_matrix()[0][0]
    N2         = phonon.get_supercell_matrix()[1][1]
    N3         = phonon.get_supercell_matrix()[2][2]
    image_num  = len(directions)
    sym_dict   = {0:'nonsym', 1:'sym'}
    sym        = sym_dict[phonon._is_symmetry]
    job_name   = "x" + str(N1) + str(N2) + str(N3) + "_d" + str("%.3f" % delta) + "_" + sym
    pckl_name  = 'pickle-'+job_name+'.p'
    import os
    pwd = os.getcwd()
    calc_dir = pwd + "/calcs/" + job_name
    sp.call(["rm -rf " + calc_dir + "/POSCARS"], shell = True)
    sp.call(["mkdir -p " + calc_dir + "/POSCARS"], shell = True)
    sp.call(['cp SPOSCAR POSCAR-000'], shell=True)
    sp.call(["mv POSCAR-* SPOSCAR " + calc_dir + "/POSCARS"], shell = True)

    try:
        phonon = pickle.load(open(pckl_name, 'rb'))
        if phonon.get_force_constants() is None:
            raise ValueError
    except:
        print("******* There is no saved pickle file. ********".center(120))
        do_calc = True
    else:
        print("******* Pickle file has been loaded. ********".center(120))
        do_calc = False
    #>>>>>>>>>>>>>>>>>>> all same until now <<<<<<<<<<<<<<<<<<
    if do_calc:
        print("******* DPMD calc will be carried out. ********".center(120))
        forces = []

        from ase.io.lammpsrun import read_lammps_dump as read_dump
        from ase.io import write
        #### To get F_0 for F' ( F' = F - F_0 )
        directions = [[0,0,0,0]] + directions
        
        for i in range(image_num+1):
            print(' >>> Starting {:d}-th image calculation <<< '.format(i).center(120))
            ndir = "pos" + str(i).zfill(4) + "_atom" + str(directions[i][0]).zfill(4) + \
                "_direc" + str(directions[i][1]) + str(directions[i][2]) + str(directions[i][3])
            sp.call(["rm", "-rf", calc_dir+"/BU-"+ndir])
            sp.call(["mv", calc_dir+"/"+ndir, calc_dir+"/BU-"+ndir])
            sp.call(["mkdir", "-p", calc_dir+"/"+ndir])
            sp.call(["cp","frozen_model.pb", "input.in", 
                calc_dir+"/POSCARS/POSCAR-"+str(i).zfill(3), calc_dir+"/"+ndir])
            sp.call(["lmp-pos2lmp.awk POSCAR-"+str(i).zfill(3)+" > structure.in"],
                cwd = calc_dir+"/"+ndir, shell = True)
            sp.call(["mpiexec.hydra -np $NSLOTS lmp_mpi -in input.in > out"],
                cwd = calc_dir+"/"+ndir, shell = True)
            atoms = read_dump(calc_dir+"/"+ndir+"/out.dump", index=0, order=True)
            if i == 0:
                F_0 = atoms.get_forces(apply_constraint=False)
            else:
                if F_0_correction:
                    atoms._calc.results['forces'] -= F_0
                force_now = [atoms.get_forces(apply_constraint=False).tolist()]
                forces.extend(force_now)
                if verbose:
                    print("force_now"+str(i)+"\n"+str(np.asarray(force_now)))
            write(calc_dir+'/'+ndir+'/'+'atoms.traj', atoms)

        forces_arr = np.asarray(forces)
        if verbose:
            print("\n\n"+"forces"+"\n"+str(forces_arr))
        phonon.set_forces(forces_arr)
        phonon.produce_force_constants(forces_arr)
        if acoustic_sum_rule:
            phonon.symmetrize_force_constants()
 
        if verbose:
            print("\n\ndisplacement_dataset =>\n\n")
            print(phonon.get_displacement_dataset())
            print("\n\nforce_constants =>\n\n")
            print(phonon.get_force_constants())
 
        pickle.dump(phonon, open(pckl_name, "wb"))
    return phonon


def calc_amp(phonon, calc, verbose = False, numeric_F_dx=0.001, parallel = True, acoustic_sum_rule = False):
    """ Calculate Force Constant with AMP """
    #>>>>>>>>>>>>>>>>>>> all same from now <<<<<<<<<<<<<<<<<<
    import sys
    import subprocess as sp
    import pickle
    np.set_printoptions(threshold=np.nan)
        
    delta      = np.linalg.norm(phonon.get_displacements()[0][1:4])
    directions = phonon.get_displacement_directions()
    N1         = phonon.get_supercell_matrix()[0][0]
    N2         = phonon.get_supercell_matrix()[1][1]
    N3         = phonon.get_supercell_matrix()[2][2]
    image_num  = len(directions)
    sym_dict   = {0:'nonsym', 1:'sym'}
    sym        = sym_dict[phonon._is_symmetry]
    job_name   = "x" + str(N1) + str(N2) + str(N3) + "_d" + str("%.3f" % delta) + "_" + sym
    pckl_name  = 'pickle-'+job_name+'.p'
    import os
    pwd = os.getcwd()
    calc_dir = pwd + "/calcs/" + job_name
    sp.call(["rm -rf " + calc_dir + "/POSCARS"], shell = True)
    sp.call(["mkdir -p " + calc_dir + "/POSCARS"], shell = True)
    sp.call(['cp SPOSCAR POSCAR-000'], shell=True)
    sp.call(["mv POSCAR-* SPOSCAR " + calc_dir + "/POSCARS"], shell = True)

    try:
        phonon = pickle.load(open(pckl_name, 'rb'))
        if phonon.get_force_constants() is None:
            raise ValueError
    except:
        print("******* There is no saved pickle file. ********".center(120))
        do_calc = True
    else:
        print("******* Pickle file has been loaded. ********".center(120))
        do_calc = False
    #>>>>>>>>>>>>>>>>>>> all same until now <<<<<<<<<<<<<<<<<<
    if do_calc:
        print("******* AMP calc will be carried out. ********".center(120))
        forces = []
        for i in range(image_num):
            print(' >>> Starting {:d}-th image calculation <<< '.format(i+1).center(120))
            ndir = "pos" + str(i+1).zfill(4) + "_atom" + str(directions[i][0]).zfill(4) + \
                "_direc" + str(directions[i][1]) + str(directions[i][2]) + str(directions[i][3])
            sp.call(["rm", "-rf", calc_dir+"/BU-"+ndir])
            sp.call(["mv", calc_dir+"/"+ndir, calc_dir+"/BU-"+ndir])
            sp.call(["mkdir", "-p", calc_dir+"/"+ndir])
            sp.call(["cp", calc_dir+"/POSCARS/POSCAR-"+str(i+1).zfill(3), calc_dir+"/"+ndir])

            ########### calculate forces & atomic energies with amp ############
            from ase.io import read
            atoms = read(calc_dir+"/"+ndir+"/POSCAR-"+str(i+1).zfill(3),
                         format = "vasp")
            atoms.set_pbc(True)
            from amp import Amp
            atoms.set_calculator(calc)
            from amp.utilities import Logger
            log = atoms._calc._log
            log("**********************************************")
            log(str(i+1)+" th image calculation start", tic = 'image')
            log("**********************************************")

            #********** numerical force must precede ***********
            # force_now = [calc.calculate_numerical_forces(
                # atoms,
                # d = numeric_F_dx,
                # parallel = parallel,
                # )] # alternative
            force_now = [atoms.get_forces(apply_constraint=False).tolist()] # alternative
            #stress_now = calc.calculate_numerical_stress(atoms) # just in case
            #********** energy calculation ***********
            energy_now = atoms.get_potential_energy()
            energies_now = atoms.get_potential_energies()
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
            log("**********************************************")
            log(str(i+1)+" th image calculated", toc = 'image')
            log("**********************************************")
            
        forces_arr = np.asarray(forces)
        if verbose:
            print("\n\n"+"forces"+"\n"+str(forces_arr))
        phonon.set_forces(forces_arr)
        phonon.produce_force_constants()
        if acoustic_sum_rule:
            phonon.symmetrize_force_constants()
 
        if verbose:
            print("\n\ndisplacement_dataset =>\n\n")
            print(phonon.get_displacement_dataset())
            print("\n\nforce_constants =>\n\n")
            print(phonon.get_force_constants())
 
        pickle.dump(phonon, open(pckl_name, "wb"))
    return phonon


def calc_amp_tf(phonon, calc, verbose = False, numeric_F_dx=0.001, parallel = True, acoustic_sum_rule = False):
    """ Calculate Force Constant with AMP with tensorflow """
    #>>>>>>>>>>>>>>>>>>> all same from now <<<<<<<<<<<<<<<<<<
    import sys
    import subprocess as sp
    import pickle
    np.set_printoptions(threshold=np.nan)
        
    delta      = np.linalg.norm(phonon.get_displacements()[0][1:4])
    directions = phonon.get_displacement_directions()
    N1         = phonon.get_supercell_matrix()[0][0]
    N2         = phonon.get_supercell_matrix()[1][1]
    N3         = phonon.get_supercell_matrix()[2][2]
    image_num  = len(directions)
    sym_dict   = {0:'nonsym', 1:'sym'}
    sym        = sym_dict[phonon._is_symmetry]
    job_name   = "x" + str(N1) + str(N2) + str(N3) + "_d" + str("%.3f" % delta) + "_" + sym
    pckl_name  = 'pickle-'+job_name+'.p'
    import os
    pwd = os.getcwd()
    calc_dir = pwd + "/calcs/" + job_name
    sp.call(["rm -rf " + calc_dir + "/POSCARS"], shell = True)
    sp.call(["mkdir -p " + calc_dir + "/POSCARS"], shell = True)
    sp.call(['cp SPOSCAR POSCAR-000'], shell=True)
    sp.call(["mv POSCAR-* SPOSCAR " + calc_dir + "/POSCARS"], shell = True)

    try:
        phonon = pickle.load(open(pckl_name, 'rb'))
        if phonon.get_force_constants() is None:
            raise ValueError
    except:
        print("******* There is no saved pickle file. ********".center(120))
        do_calc = True
    else:
        print("******* Pickle file has been loaded. ********".center(120))
        do_calc = False
    #>>>>>>>>>>>>>>>>>>> all same until now <<<<<<<<<<<<<<<<<<
    if do_calc:
        print("******* AMP-TF calc will be carried out. ********".center(120))
        forces = []
        for i in range(image_num):
            print(' >>> Starting {:d}-th image calculation <<< '.format(i+1).center(120))
            ndir = "pos" + str(i+1).zfill(4) + "_atom" + str(directions[i][0]).zfill(4) + \
                "_direc" + str(directions[i][1]) + str(directions[i][2]) + str(directions[i][3])
            sp.call(["rm", "-rf", calc_dir+"/BU-"+ndir])
            sp.call(["mv", calc_dir+"/"+ndir, calc_dir+"/BU-"+ndir])
            sp.call(["mkdir", "-p", calc_dir+"/"+ndir])
            sp.call(["cp", calc_dir+"/POSCARS/POSCAR-"+str(i+1).zfill(3), calc_dir+"/"+ndir])

            ########### calculate forces & atomic energies with amp ############
            from ase.io import read
            atoms = read(calc_dir+"/"+ndir+"/POSCAR-"+str(i+1).zfill(3),
                         format = "vasp")
            atoms.set_pbc(True)
            atoms.set_calculator(calc)
            from amp.utilities import Logger
            log = atoms._calc._log
            log("**********************************************")
            log(str(i+1)+" th image calculation start", tic = 'image')
            log("**********************************************")

            #********** numerical force must precede ***********
            force_now = [calc.calculate_numerical_forces(
                atoms,
                d = numeric_F_dx,
                parallel = parallel,
                )] # alternative
            #force_now = [atoms.get_forces(apply_constraint=False).tolist()] # alternative
            #stress_now = calc.calculate_numerical_stress(atoms) # just in case
            #********** energy calculation ***********
            energy_now = atoms.get_potential_energy()
            energies_now = atoms.get_potential_energies()
            #********** force information restore ***********
            atoms._calc.results['forces'] = force_now[0]
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
            log("**********************************************")
            log(str(i+1)+" th image calculated", toc = 'image')
            log("**********************************************")

        forces_arr = np.asarray(forces)
        if verbose:
            print("\n\n"+"forces"+"\n"+str(forces_arr))
        phonon.set_forces(forces_arr)
        phonon.produce_force_constants()
        if acoustic_sum_rule:
            phonon.symmetrize_force_constants()
 
        if verbose:
            print("\n\ndisplacement_dataset =>\n\n")
            print(phonon.get_displacement_dataset())
            print("\n\nforce_constants =>\n\n")
            print(phonon.get_force_constants())
 
        pickle.dump(phonon, open(pckl_name, "wb"))
    return phonon

def calc_amp_tf_bunch(phonon, calc, verbose = False, numeric_F_dx=0.001, parallel = True, acoustic_sum_rule = False):
    """ Calculate Force Constant with AMP with tensorflow (fast version) """
    #>>>>>>>>>>>>>>>>>>> all same from now <<<<<<<<<<<<<<<<<<
    import sys
    import subprocess as sp
    import pickle
    np.set_printoptions(threshold=np.nan)
        
    delta      = np.linalg.norm(phonon.get_displacements()[0][1:4])
    directions = phonon.get_displacement_directions()
    N1         = phonon.get_supercell_matrix()[0][0]
    N2         = phonon.get_supercell_matrix()[1][1]
    N3         = phonon.get_supercell_matrix()[2][2]
    image_num  = len(directions)
    sym_dict   = {0:'nonsym', 1:'sym'}
    sym        = sym_dict[phonon._is_symmetry]
    job_name   = "x" + str(N1) + str(N2) + str(N3) + "_d" + str("%.3f" % delta) + "_" + sym
    pckl_name  = 'pickle-'+job_name+'.p'
    import os
    pwd = os.getcwd()
    calc_dir = pwd + "/calcs/" + job_name
    sp.call(["rm -rf " + calc_dir + "/POSCARS"], shell = True)
    sp.call(["mkdir -p " + calc_dir + "/POSCARS"], shell = True)
    sp.call(['cp SPOSCAR POSCAR-000'], shell=True)
    sp.call(["mv POSCAR-* SPOSCAR " + calc_dir + "/POSCARS"], shell = True)

    try:
        phonon = pickle.load(open(pckl_name, 'rb'))
        if phonon.get_force_constants() is None:
            raise ValueError
    except:
        print("******* There is no saved pickle file. ********".center(120))
        do_calc = True
    else:
        print("******* Pickle file has been loaded. ********".center(120))
        do_calc = False
    #>>>>>>>>>>>>>>>>>>> all same until now <<<<<<<<<<<<<<<<<<
    if do_calc:
        print("******* AMP-TF-bunch calc will be carried out. ********".center(120))
        forces = []
        for i in range(image_num):
            print(' >>> Starting {:d}-th image calculation <<< '.format(i+1).center(120))
            ndir = "pos" + str(i+1).zfill(4) + "_atom" + str(directions[i][0]).zfill(4) + \
                "_direc" + str(directions[i][1]) + str(directions[i][2]) + str(directions[i][3])
            sp.call(["rm", "-rf", calc_dir+"/BU-"+ndir])
            sp.call(["mv", calc_dir+"/"+ndir, calc_dir+"/BU-"+ndir])
            sp.call(["mkdir", "-p", calc_dir+"/"+ndir])
            sp.call(["cp", calc_dir+"/POSCARS/POSCAR-"+str(i+1).zfill(3), calc_dir+"/"+ndir])

            ########### calculate forces & atomic energies with amp ############
            from ase.io import read
            atoms = read(calc_dir+"/"+ndir+"/POSCAR-"+str(i+1).zfill(3),
                         format = "vasp")
            atoms.set_pbc(True)
            atoms.set_calculator(calc)
            from amp.utilities import Logger
            log = atoms._calc._log
            log("**********************************************")
            log(str(i+1)+" th image calculation start", tic = 'image')
            log("**********************************************")

            #********** numerical force must precede ***********
            force_now = [calc.calculate_numerical_forces(
                atoms,
                d = numeric_F_dx,
                parallel = parallel,
                )] # alternative
            #force_now = [atoms.get_forces(apply_constraint=False).tolist()] # alternative
            #stress_now = calc.calculate_numerical_stress(atoms) # just in case
            #********** energy calculation ***********
            energy_now = atoms.get_potential_energy()
            energies_now = atoms.get_potential_energies()
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
            log("**********************************************")
            log(str(i+1)+" th image calculated", toc = 'image')
            log("**********************************************")
            
        forces_arr = np.asarray(forces)
        if verbose:
            print("\n\n"+"forces"+"\n"+str(forces_arr))
        phonon.set_forces(forces_arr)
        phonon.produce_force_constants()
        if acoustic_sum_rule:
            phonon.symmetrize_force_constants()
 
        if verbose:
            print("\n\ndisplacement_dataset =>\n\n")
            print(phonon.get_displacement_dataset())
            print("\n\nforce_constants =>\n\n")
            print(phonon.get_force_constants())
 
        pickle.dump(phonon, open(pckl_name, "wb"))
    return phonon

def make_band(path, N_q):
    bands = []
    for i in range(len(path)):
        q_start  = np.array(path[i][0])
        q_end    = np.array(path[i][1])
        band = []
        for j in range(N_q+1):
            band.append(q_start + (q_end - q_start) / N_q * j)
        bands.append(band)
    return bands

# def plot_band(phonon, labels=None):
    # import matplotlib.pyplot as plt
    # if labels:
        # from matplotlib import rc
        # rc('font',**{'family':'serif','sans-serif':['Times']})
        # rc('text', usetex=False)

    # fig, ax = plt.subplots()
    # ax.xaxis.set_ticks_position('both')
    # ax.yaxis.set_ticks_position('both')
    # ax.xaxis.set_tick_params(which='both', direction='in')
    # ax.yaxis.set_tick_params(which='both', direction='in')

    # phonon._band_structure.plot(plt, labels=labels)
    # return plt

def set_projection(phonon, proj_eigvec):
    self = phonon._band_structure
    self._projections = {}
    self._proj_freq = {}

    for key in proj_eigvec.keys():
        proj_tmp = []
        freq_tmp = []
        for _path in self._paths:
            proj_tmp2 = []
            freq_tmp2 = []
            for _q in _path:
                freq, eigvec = phonon.get_frequencies_with_eigenvectors(_q)
                eigvec = eigvec.T   # numpy.linalg.eigh returns transposed eigen vector
                #### Gather
                proj_tmp2.append(eigvec)
                freq_tmp2.append(freq)
            #### Gather
            proj_tmp.append(proj_tmp2)
            freq_tmp.append(freq_tmp2)
        #### Transform to array
        proj_tmp = np.array(proj_tmp)
        freq_tmp = np.array(freq_tmp)
        #### calculate similarity
        proj_tmp = np.square(np.absolute(np.squeeze(np.matmul(np.conj(proj_tmp), np.expand_dims(proj_eigvec[key], axis=-1)))))
        self._projections[key] = proj_tmp
        self._proj_freq[key]   = freq_tmp

def bs_plot(self, plt, ax, proj_size_factor, proj_colors, proj_alpha, reverse_seq, labels=None):
    if self._projections is not None:
        #### Define key list
        key_list = list(self._projections.keys())
        if reverse_seq:
            key_list.reverse()
        #### Pick colors
        proj_colors = proj_colors[len(proj_colors)-len(key_list):]
        if reverse_seq:
            proj_colors.reverse()
        #### 
        legend = []
        #### Iter for projector eigenvectors
        for key in key_list:
            #### Iter for q_path fragments
            # for distances, frequencies, projections, proj_freq in zip(self._distances,
            for distances, projections, freq in zip(self._distances,
                                                         self._projections[key],
                                                         self._proj_freq[key]):
                #### Iter for band lines
                for i in range(len(freq.T)):
                    plt.plot(distances, freq.T[i], 'k-')
                    legend_tmp = plt.scatter(
                        distances,
                        freq.T[i],
                        proj_size_factor * projections.T[i],
                        proj_colors[-1],
                        alpha=proj_alpha,
                        edgecolors='none',
                        label=key,
                        )
            #### Gather just the one sample legend
            legend.append(legend_tmp)
            #### Throw away used color
            proj_colors.pop()
        #### Legend plot
        if reverse_seq:
            legend.reverse(); key_list.reverse()
        plt.legend(legend, key_list, scatterpoints = 1, fontsize='x-large')
    else:
        for distances, frequencies in zip(self._distances,
                                          self._frequencies):
            for freqs in frequencies.T:
                plt.plot(distances, freqs, 'k-')

    plt.ylabel('Frequency')
    plt.xlabel('Wave vector')

    if labels and len(labels) == len(self._special_points):
        plt.xticks(self._special_points, labels, fontsize = 20) # ssrokyz
    else:
        plt.xticks(self._special_points, [''] * len(self._special_points))
    plt.xlim(0, self._distance)
    plt.axhline(y=0, linestyle=':', linewidth=0.5, color='k')

def plot_band_and_dos(
    phonon,
    pdos_indices     = None,
    labels           = None,
    unit             = 'THz',
    proj_eigvec      = None,
    proj_size_factor = 400.,
    proj_colors      = ['r', 'g', 'b', 'c', 'm', 'y'],
    proj_alpha       = 0.5,
    ylim_lower       = None,
    ylim_upper       = None,
    reverse_seq      = False,
    ):
    """
    proj_eigvec = (dict) = Eigenvectors that will be used for projection. Keys of dict will be used as label for pyplot.
    proj_colors = (list) = Use colors sequently. Duplication is possible.
    reverse_seq = (bool) = Change sequence of plotting scatters. (Decides which one goes up. Only esthetic purpose.)
    """
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    if labels:
        from matplotlib import rc
        rc('font',**{'family':'serif','sans-serif':['Times']})
        rc('text', usetex=False)

    #### Variable setting
    proj_colors.reverse()

    plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1])
    ax1 = plt.subplot(gs[0, 0])
    ax1.xaxis.set_ticks_position('both')
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_tick_params(which='both', direction='in')
    ax1.yaxis.set_tick_params(which='both', direction='in')

    #### projection
    phonon._band_structure._projections = None
    if proj_eigvec is not None:
        # dtype managing
        for key in proj_eigvec.keys():
            proj_eigvec[key] = np.array(proj_eigvec[key], dtype=np.complex128)
        set_projection(phonon, proj_eigvec)

    bs_plot(phonon._band_structure, plt, ax1, proj_size_factor, proj_colors, proj_alpha, reverse_seq, labels=labels)
    if unit == 'meV':
        plt.ylabel('Frequency(meV)', fontsize=22)
    elif unit == 'THz':
        plt.ylabel('Frequency(THz)', fontsize=22)
    plt.xlabel('')
    plt.grid(True)
    plt.title('Phonon dispersion', fontsize=24)
    plt.yticks(fontsize=20)

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
    plt.title('DOS', fontsize=24)
    plt.xlabel('')
    plt.xticks([])
    plt.ylim(ylim_lower, ylim_upper)

    return plt

