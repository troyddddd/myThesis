
#file includes: input data functin

#transfer data from csv to list

import csv
import numpy as np
#input data for matrix s,l,r
def data_to_srl(filename):
	import numpy as np
	result = []
	with open(filename,'r') as f:
		result = [line.strip().split(" ") for line in f.readlines()]
	f.close()
	result = np.asarray(result)
	result = result.astype(np.float64)
	return result

#input data for BPF
def data_to_BPF(filename):
	import numpy as np
	result = []
	with open (filename,'r') as f:
		for line in f:
			line = line.strip()
			for ele in line.split(" "):
				result.append(int(ele))
	f.close()
	result = np.asarray(result)
	return result

def data_to_ptup(filename):
	result = []
	with open(filename,'r') as f:
		for line in f:
			line = line.strip()
			result.append((line.split(" ")[0],line.split(" ")[1]))
	f.close()
	return result
