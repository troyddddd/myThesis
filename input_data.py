
#file includes: input data functin

#transfer data from csv to list

import csv

#input data for matrix s,l,r
def data_to_srl(filename):
	import numpy as np
	result = np.zeros([r+1,n])
	with open(filename,'r') as file:
		for line in file:
			for cha in line.split(" "):
				result[i,j] = cha
				i+=1
			j+=1

	return result

#input data for BPF