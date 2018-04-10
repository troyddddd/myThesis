

import numpy as np

class Solution:
	def __init__(self,x,y):
		self.first = x
		self.second = y
	def compose(self):
		result = np.zeros((self.first,self.second))
		print result
		for i in range(self.first):
			for j in range(self.second):
				result[i,j] = i+j
		return result

solution = Solution(2,2)
result = solution.compose()
print result