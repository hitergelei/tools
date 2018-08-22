#!/usr/bin/env python2
##### CODE BY YOUNGJAE CHOI #####

from ase import Atoms, Atom
import random
from ase.build import make_supercell
from numpy import ndarray

def RanPoAtoms(cut_off_radius,
			   symbols=None,
			   positions=None, numbers=None,
               tags=None, momenta=None, masses=None,
               magmoms=None, charges=None,
               scaled_positions=None,
               cell=None, pbc=None, celldisp=None,
               constraint=None,
               calculator=None,
               info=None):
	if positions is not None:
		print("\npositions must not be given\n")
		exit(1)
	if scaled_positions is not None:
		print("\nscaled_positions must not be given\n")
		exit(1)
	else:
		atoms = Atoms(symbols=symbols,
                      positions=positions, numbers=numbers,
                      tags=tags, momenta=momenta, masses=masses,
                      magmoms=magmoms, charges=charges,
                      scaled_positions=None,
                      cell=cell, pbc=pbc, celldisp=celldisp,
                      constraint=constraint,
                      calculator=calculator,
                      info=info)	
		l = 0
		while True:
			l+=1
			print("trying step :: "+str(l))
			scaled_posis = []
			for i in range(len(atoms)):
				scaled_posi = []
				for j in range(3):
					scaled_posi.append(random.random())
				scaled_posis.append(scaled_posi)
			atoms.set_scaled_positions(scaled_posis)
			supercell = make_supercell(atoms,[[2,0,0],[0,2,0],[0,0,2]])
			dist = supercell.get_all_distances()
			coll = []
			for i in range(len(supercell)):
				for j in range(len(supercell)):
					if i is not j:
						coll.append(dist[i][j])
			if min(coll) > cut_off_radius:
				break
		return atoms

def RanPoAtoms_2(cut_off_radius,
				random_degree,
			    symbols=None,
			    positions=None, numbers=None,
                tags=None, momenta=None, masses=None,
                magmoms=None, charges=None,
                scaled_positions=None,
                cell=None, pbc=None, celldisp=None,
                constraint=None,
                calculator=None,
                info=None):
	if positions is None:
		if scaled_positions is None:
			print("\nNo backbone structure is given.\n")
			exit(1)
	else:
		atoms = Atoms(symbols=symbols,
                      positions=positions, numbers=numbers,
                      tags=tags, momenta=momenta, masses=masses,
                      magmoms=magmoms, charges=charges,
                      scaled_positions=scaled_positions,
                      cell=cell, pbc=pbc, celldisp=celldisp,
                      constraint=constraint,
                      calculator=calculator,
                      info=info)	
		
############### shuffle positions ################
		array_scaled_positions = atoms.get_scaled_positions()
		shuffled_scaled_posis = array_scaled_positions.tolist()
		random.shuffle(shuffled_scaled_posis)
		atoms.set_scaled_positions(shuffled_scaled_posis)
		
############### get local random distribution diameter ################
		supercell = make_supercell(atoms,[[2,0,0],[0,2,0],[0,0,2]])
		dist = supercell.get_all_distances()
		coll = []
		for i in range(len(supercell)):
			for j in range(len(supercell)):
				if i is not j:
					coll.append(dist[i][j])
		ran_diameter = min(coll)

############### shuffled position list ################
		array_shuffled_posis = atoms.get_positions()

		l = 0
		while True:
			l+=1
			print("trying step :: "+str(l))
			shuffled_posis = array_shuffled_posis.tolist()
			for i in range(len(atoms)):
				for j in range(3):
					shuffled_posis[i][j]+=((random.random()-0.5)*ran_diameter*random_degree)
			tmp = atoms
			tmp.set_positions(shuffled_posis)
			supercell = make_supercell(tmp,[[2,0,0],[0,2,0],[0,0,2]])
			dist = supercell.get_all_distances()
			coll = []
			for i in range(len(supercell)):
				for j in range(len(supercell)):
					if i is not j:
						coll.append(dist[i][j])
			if min(coll) > cut_off_radius:
				break
		atoms = tmp
		return atoms
