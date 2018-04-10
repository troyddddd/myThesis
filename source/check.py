# import sys
# import os

# dir_name = ""
# test = os.listdir(dir_name)

# for item in test:
#     if item.endswith(".zip"):
#         os.remove(os.path.join(dir_name, item))


import numpy as np

class Solution:
	def __init__(self,x,y):
		self.first = x
		self.second = y
	def test(self):
		temp = np.ones((self.first,self.second))
		for i in range(self.first):
			for j in range(self.second):
				temp[i,j] = i+j
		return temp
solution = Solution(2,3)
print(solution.test())

