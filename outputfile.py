
def writeOriginalmatrix(S,edit,filename):
	import numpy as np
	with open (filename,edit) as openfile:
		for i in range(S.shape[0]):
			for j in range(S.shape[1]):
			    openfile.write(str(S[i][j])+' ')
			openfile.write('\n')

def writeBPF(S,end,filename):
	with open(filename,'a') as openfile:
		if end == False:
			openfile.write(str(S)+' ')
		else:
			openfile.write(str(S)+'\n')