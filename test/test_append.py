
import os
import numpy as np

temp = np.random.rand(3,2)
temp1 = np.random.rand(2,3)
print temp
print temp[1,1]
print temp1

if os.path.isfile('./123.txt'):
	with open ('123.txt','a') as f:
		f.write('\n')
		for i in range(temp.shape[0]):
			for j in range(temp.shape[1]):
				f.write(str(temp[i,j])+' ')
			f.write('\n')
	f.close()
else:
	with open ('123.txt','w') as f:
		for i in range(temp1.shape[0]):
			for j in range(temp1.shape[1]):
				f.write(str(temp1[i,j])+' ')
			f.write('\n')
	f.close()

