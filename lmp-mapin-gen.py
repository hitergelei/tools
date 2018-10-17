#!/usr/bin/env python

import sys

print("\n\n#######################################################################################")
print("              %%%%%%%%%%% This code will make you a map.in file %%%%%%%%%")
print("                 useage ==> ./mapin-gen.py 'lx' 'ly' 'lz' 'basis#'")
print("                        e.g.) ./mapin-gen.py 8 8 8 2")
print("                       The result file name will be map.in")
print("#######################################################################################")
if len(sys.argv) is 5:
    print("                The Number of arguments is correct.\n\n")
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
    sys.exit(1)

lx = int(sys.argv[1])
ly = int(sys.argv[2])
lz = int(sys.argv[3])
nbasis = int(sys.argv[4])

mapin = open("map.in", "w")

mapin.write(str(lx)+" "+str(ly)+" "+str(lz)+" "+str(nbasis)+"\n#l1 l2 l3 k atom_id\n")

i=0
for x in range(lx):
    for y in range(ly):
        for z in range(lz):
            for b in range(nbasis):
                i += 1
                mapin.write(str(x)+" "+str(y)+" "+str(z)+" "+str(b)+" "+str(i)+"\n")

