#!/usr/bin/env python
import sys


print("\nuseage ==> ./vasp_gap.py 'fermi energy'\n")
if len(sys.argv) is 2:
	print(" ")
else:
	print("*****ERROR***** The number of arguments is not correct *****ERROR*****\n\n")
	sys.exit(1)

fe = float(sys.argv[1])

input = open("EIGENVAL", "r")
output = open("gap_EIGEN.txt", "w")

def isN(i):
  try:
    float(i) and i is int(i)
    return True
  except:
    return False

lines = input.readlines()
for line in lines:
	output.write(line[:-1]+"	")
	column = None
	if len(line.split()) > 1:
		column = line.split()[1]
		if isN(line.split()[0]):
			output.write(str(float(column)-fe)+"\n")
		else:
			###output.write("\n")
			output.write("\n")
	else:
		output.write("\n")
			 
