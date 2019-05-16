#!/usr/bin/env python

from asap3.analysis.rdf import RadialDistributionFunction as RDF
from ase.io.trajectory import Trajectory
import sys
import matplotlib.pyplot as plt
import numpy as np

print("\n\n#######################################################################################\n")
print("      %%%%%%%%%%% This code will give you Raidial Distribution Function %%%%%%%%%")
print("useage ==> ./ase-rdf.py 'trajectory file' 'rMax' 'nBins'")
print("           EXAMPLE) ./ase-rdf.py si.traj 12 500")
print("                    OUTPUT file is 'rdf-si.traj.log'")
print("#######################################################################################")
if len(sys.argv) is 4:
    print("          The Number of arguments is correct.\n\n")
else:
    print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
    sys.exit(1)

inp_fname = sys.argv[1]
rMax = sys.argv[2]
rMax = float(rMax)
nBins = sys.argv[3]
nBins = int(nBins)

traj = Trajectory(inp_fname, "r")
RDFobj = RDF(atoms=traj[0],
             rMax=rMax,
             nBins=nBins)

for i in range(1,len(traj)):
    RDFobj.atoms = traj[i]  # Fool RDFobj to use the new atoms
    RDFobj.update()           # Collect data
    if i % 1000 == 999:
        print("\t Updating "+str(i+1)+" th image's RDF")

rdf = RDFobj.get_rdf()
x = np.arange(nBins) * rMax / nBins

out = open("rdf-"+inp_fname+".log", "w")
for i in range(len(x)):
    out.write(str("%.7f" %x[i])+"\t"+str(rdf[i])+"\n")

print("\n\n#######################################################################################\n")
print("      %%%%%%%%%%% This code will give you Raidial Distribution Function %%%%%%%%%")
print("useage ==> ./ase-rdf.py 'trajectory file' 'rMax' 'nBins'")
print("           EXAMPLE) ./ase-rdf.py si.traj 12 500")
print("                    OUTPUT file is 'rdf-si.traj.log'")
print("#######################################################################################")

plt.plot(x, rdf)
plt.show()

