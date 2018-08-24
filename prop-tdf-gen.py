#!/usr/bin/env python

import sys

print("\n\n#######################################################################################\n")
print("      %%%%%%%%%%% This code will make PROPhet 'training_data' file  %%%%%%%%%")
print("   useage ==> ./prop-tdf-gen.py 'location of data files directory' 'first' 'last'")
print("           EXAMPLE) ./prop-tdf-gen.py ../data/ 0 10")
print("                    OUTPUT file is 'training_data'")
print("#######################################################################################")
if len(sys.argv) is 4:
	print("          The Number of arguments is correct.\n\n")
else:
	print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
	sys.exit(1)

dd = sys.argv[1]
ni = int(sys.argv[2])
nf = int(sys.argv[3])

tdf = open("training_data", "w")
tdf.write('# start with a "data" statement. Data for this example were generated\n# with VASP. The "data" statement affects all data examples that come\n# after it, until the next "data" statement.\n\n# Then provide the base directory for each training example, each \n# on its own line.  Often it is best to specify absolute paths,\n# so that different training runs can be conducted on the same data. \n# Here, relative paths are used.')
tdf.write("\n\ndata code=vasp density=Chg_density")

for i in range(ni, nf + 1):
	tdf.write("\n"+dd+"/"+str(i))

tdf.close()
