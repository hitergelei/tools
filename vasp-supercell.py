#!/usr/bin/env python
import sys
import ase.io.vasp

if len(sys.argv) < 3:
   sys.exit("Usage: %s filename N1xN2xN3" % sys.argv[0]) 

cell = ase.io.vasp.read_vasp(sys.argv[1])

repeats= [int(x) for x in sys.argv[2].split("x")]
if len(repeats) != 3:
    sys.exit("Erroring parsing supercell repeat. It should be three positive numbers separated by \"x\"")

ase.io.vasp.write_vasp(sys.argv[1]+"-"+sys.argv[2]+".vasp",cell*(repeats[0],repeats[1],repeats[2]),label="supercell-"+sys.argv[2],direct=True,sort=True)
