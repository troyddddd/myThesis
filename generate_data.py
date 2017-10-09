
import csv

def Initialize(q,n,r,p,Ru,Ro,choice):
    
    """
This function creates 3 matrices of size [r+1,n].

S - State Matrix
    Either State 0 (w.p. 1-q) or 1 (w.p. q)  ... Bottom Row Seeded with State 2 
    State 1 Sites above the bottom row are also designated as "revenue generators" with probability p
R - Resistance Matrix
    Values Based on Distribution with mean Ru and standard dev Ro
L - Lattice (As Seen by the Company)
    Sites are state -1 except for the bottom row seeds from matrix S
       
    """
    count=0
    S=np.zeros([r+1,n])
    R=np.ones([r+1,n])
    L=np.ones([r+1,n])*-1
    BPF=np.ones(n)*-1
    BPF=BPF.astype(int)
    p_tup=[]
    
    if choice==0:
        for i in range(r+1):
            for j in range(n):
                rand=np.random.uniform(0,1)
                Rrand=np.random.lognormal(Ru,Ro)
                R[i,j]=Rrand
                if rand<=q:
                    S[i,j]=1
                    if np.random.uniform(0,1)<=p and i!=0:
                        p_tup.append((i,j))
        for x in range(n):       #Changes half of the State 1 squares to State 2 in the bottom row
            if S[0,x]==1:
                seed=np.random.uniform(0,1)
                if seed<=0.5:
                    S[0,x]=2
                    L[0,x]=S[0,x]
                    count+=1
                    BPF[x]=0
        if count==0: #ensures that at least one state 2 will be seeded
            random=np.random.randint(0,n)
            S[0,random]=2
            L[0,random]=S[0,random]
            BPF[random]=0
    if choice==1: #New to Thesis 5.0 this allows the matrices to be dependent upon the height.
        for i in range(r+1):
            for j in range(n):
                rand=np.random.uniform(0,1)
                Rrand=np.random.lognormal(Ru,Ro)
                R[i,j]=Rrand
                if rand<=q:
                    S[i,j]=1
                    if i!=0:
                        if np.random.uniform(0,1)<=float(p)/np.log(9+i): #This p/i is what makes the probabiility of a revenue generator decrease as the height increases
                            p_tup.append((i,j))
        for x in range(n):       #Changes half of the State 1 squares to State 2 in the bottom row
            if S[0,x]==1:
                seed=np.random.uniform(0,1)
                if seed<=0.5:
                    S[0,x]=2
                    L[0,x]=S[0,x]
                    count+=1
                    BPF[x]=0
        if count==0: #ensures that at least one state 2 will be seeded
            random=np.random.randint(0,n)
            S[0,random]=2
            L[0,random]=S[0,random]
            BPF[random]=0
    
    """
    S - State Matrix
    BPF - Initial Best Practice Frontier
    p_tup - Array of tuples designating the sites designated as "revenue generators."
    R - Resistance Matrix
    L- Lattic
    """
    return  S,BPF,p_tup,R,L

nums_runs = 0
q = 0.0
n = 0
r = 0
p = 0.0
Ru = 0.0
Ro = 0.0
choice=str(raw_input("Would you like to perform multiple runs across multiple variables? (Y or N) "))
        prob_choice=int(raw_input("How do you want the revenue generators within the lattice to be distributed?\n\n0.) Evenly Distributed\n1.) More Dense at the Bottom and More Valuable at the Top\n"))
        if choice.lower()=='n' or choice.lower=='no':
                names=['num_runs','q','n','r','p','Ru','Ro','pu','po','initial_b','b_percent','E_min', 't_max']
                while True:
                    try:
                        num_runs=int(raw_input("How many runs would you like to conduct?"))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter an integer greater than zero!")
                        continue
                    if num_runs<=0:
                        print ("Please input an integer greater than zero")
                        continue
                    else:
                        break
                while True:
                    try:
                        q=float(raw_input("What value for the percolation probability (0<q<1)?"))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter a real number between zero and one!")
                        continue
                    if q<0 or q>1:
                        print ("Please input a probability between zero and one!")
                        continue
                    else:
                        break
                while True:
                    try:
                        n=int(raw_input("How many columns do you want to initialize the matrix to (n)?"))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter an integer greater than zero!")
                        continue
                    if n<=0:
                        print ("Please input an integer greater than zero")
                        continue
                    else:
                        break
                while True:
                    try:
                        r=int(raw_input("What search radius (r) would you like to use? "))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter an integer greater than zero!")
                        continue
                    if r<=0:
                        print ("Please input an integer greater than zero")
                        continue
                    else:
                        break
                while True:
                    try:
                        p=float(raw_input("What value for the probability that a state 1 site is a prize (p)?"))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter a real number between zero and one!")
                        continue
                    if p<0 or p>1:
                        print ("Please input a probability between zero and one!")
                        continue
                    else:
                        break
                while True:
                    try:
                        Ru=float(raw_input("What is the mean of the values in the resistance matrix (Ru)?"))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter a real number greater than zero!")
                        continue
                    if Ru<0:
                        print ("Please input a positive number!")
                        continue
                    else:
                        break
                while True:
                    try:
                        Ro=float(raw_input("What is the standard deviation of the values in the resistance matrix (Ro)?"))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter a real number greater than zero!")
                        continue
                    if Ro<0:
                        print ("Please input a positive number!")
                        continue
                    else:
                        break
                
S,BPF,p_tup,R,L = Initialize(q,n,r,p,Ru,Ro,choice)

with open('source_s.csv','w') as fs:
	for i in range(len(S)):
		for j in range(len(S[i])):
			fs.write("%.2f"%round(S[i][j])," ")
		fs.write('\n')
with open('source_BPF.csv','w') as fBPF:
	for i in range(len(BPF)):
		fBPF.write("%.2f"%round(BPF[i]))
with open('source_ptup.csv','w') as ftup:
