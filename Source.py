

import csv
import os
import numpy as np
from config import Config

class DataFile():
	def __init__(self,q,n,r,p,Ru,Ro,choice):
		self.q_val = float(q)
		self.col_num = n
		self.search_radius = r
		self.p_val = float(p)
		self.ru_val = float(Ru)
		self.ro_val = float(Ro)
		self.choice = choice

		# self.S = np.zeros([self.search_radius+1, self.col_num])
		# self.R = np.ones([self.search_radius+1, self.col_num])
		# self.L = np.ones([self.search_radius+1, self.col_num])*-1
		# self.BPF = np.ones(self.col_num)*-1
		# self.BPF = self.BPF.astype(int)
		# self.p_tup = []
	def Initialize(self):
		count = 0

		# Set up initial matrix for data initialization
		# S: Type of different Simulation Site
		# R: Resistance matrix
		# L: the site state perceived by company

		S = np.zeros((self.search_radius+1, self.col_num))
		R = np.ones((self.search_radius+1, self.col_num))
		print R
		L = np.ones((self.search_radius+1, self.col_num))*-1
		BPF = np.ones(self.col_num)*-1
		BPF = BPF.astype(int)
		p_tup = []

		if self.choice == 0:
			for i in range(self.search_radius+1):
				for j in range(self.col_num):
					rand = np.random.uniform(0,1)
					Rrand = np.random.lognormal(self.ru_val,self.ro_val)
					R[i,j] = Rrand
					if rand <= self.q_val:
						S[i,j] = 1
						if np.random.uniform(0, 1) <= self.p_val and i != 0:
							p_tup.append((i,j))
			for x in range(0,self.col_num):		#Changes half of the State 1 squares to State 2 in the bottom row
				if S[0,x]==1:
					seed=np.random.uniform(0, 1)
					if seed <= 0.5:
						S[0,x]=2
						L[0,x]=S[0,x]
						count+=1
						BPF[x]=0
			if count == 0:
				random = np.random.randint(0,self.col_num)
				S[0,random] = 2
				L[0,random] = S[0,random]
				BPF[random] = 0
		if self.choice == 1:
			for i in range(self.search_radius+1):
				for j in range(self.col_num):
					rand = np.random.uniform(0, 1)
					Rrand = np.random.lognormal(self.ru_val, self.ro_val)
					R[i,j] = Rrand
					if rand <= self.q_val:
						S[i,j] = 1
						if i != 0:
							if np.random.uniform(0,1) <= float(self.p_val)/np.log(9+i): #This p/i is what makes the probabiility of a revenue generator decrease as the height increases
								p_tup.append((i,j))

			for x in range(self.col_num):
				if S[0,x] == 1:
					seed = np.random.uniform(0,1)
					if seed <= 0.5:
						S[0,x] = 2
						L[0,x] = S[0,x]
						count += 1
						BPF[x] = 0
			if count == 0:
				random = np.random.randint(0,self.col_num)
				S[0,random] = 2
				L[0,random] = S[0,random]
				BPF[random] = 0
		print R
		return S,BPF,p_tup,R,L

	def Generate_Source(self,S,R,L,p_tup,BPF):

		with open('source/source_S.csv','w') as fs:
			for i in range(len(S)):
				for j in range(len(S[i])):
					fs.write(str(S[i][j]) + ' ')
				fs.write('\n')

		with open('source/source_BPF.csv','w') as fBPF:
			for i in range(len(BPF)):
				fBPF.write(str(BPF[i]) + ' ')

		with open('source/source_R.csv','w') as fr:
			for i in range(len(R)):
				for j in range(len(R[i])):
					fr.write('%.2f'%round(R[i][j]) + ' ')
				fr.write('\n')

		with open('source/source_L.csv','w') as fl:
			for i in range(len(L)):
				for j in range(len(L[i])):
					fl.write(str(L[i][j]) + ' ')
				fl.write('\n')

		with open('source/source_ptup.csv','w') as fp:
			fp.write('\n'.join('%s %s' % x for x in p_tup))

if __name__ == "__main__":
	if not os.path.isdir('/source'):
		os.system("sudo mkdir /source")
	dir_name = "source/"
	file_lst = os.listdir(dir_name)
	for item in file_lst:
		if item.endswith(".csv"):
			os.remove(os.path.join(dir_name,item))
	config = Config()
	source = DataFile(config.q_val,config.num_col,config.search_radius,config.p_val,config.ru_val,config.ro_val,config.choice)
	S,BPF,p_tup,R,L=source.Initialize()
	source.Generate_Source(S,R,L,p_tup,BPF)






























