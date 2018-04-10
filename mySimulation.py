# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 21:22:28 2016

@author: Hal
"""
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 20:25:01 2016

@author: Hal
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import copy
import math as math
from mpl_toolkits.mplot3d import Axes3D
import random as rand   #First used in Thesis 2.0
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from datetime import datetime
import outputfile as op #self built function module
import os
import random_generator as rg
import input_data as inp #self built function
import makegraph
import sys
from config import Config
import csv
"""Changes from Thesis 2.0:

1.) Addition of a line_plot function which allows the user to track how a metric changes over time during multiple simulation runs.


"""

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

    S = inp.data_to_srl('./source/source_S.csv')
    R = inp.data_to_srl('./source/source_R.csv')
    L = inp.data_to_srl('./source/source_L.csv')
    BPF =inp.data_to_BPF('./source/source_BPF.csv')
    p_tup = inp.data_to_ptup('./source/source_ptup.csv')

    """
    S - State Matrix
    BPF - Initial Best Practice Frontier
    p_tup - Array of tuples designating the sites designated as "revenue generators."
    R - Resistance Matrix
    L- Lattice
    """
    return  S,BPF,p_tup,R,L
    

def New(q,n,r,p,p_tup,h,Ru,Ro,choice): # h gives previous height

    """
    This function gives new rows to add to each of the 3 matrices.  It is called by the function 'Update'
    
    """
    S_add=np.zeros([r,n])
    R_add=np.zeros([r,n])
    L_add=np.ones([r,n])*-1
    for i in range(r):
        for j in range(n):
            rand=np.random.uniform(0,1)
            Rrand=np.random.lognormal(Ru,Ro)
            R_add[i,j]=Rrand
            if rand<=q:
                S_add[i,j]=1
                rand2=np.random.uniform(0,1)
                if choice==0:
                    if rand2<=p:
                       p_tup.append((h+i,j)) 
                if choice==1:
                    if rand2<=float(p)/np.log10(9+h+i): #see above comment in Initialize about p/i
                        p_tup.append((h+i,j))
    return  S_add,p_tup,R_add,L_add
    
    
def Update(S,R,L,BPF,r,q,p,p_tup,Ru,Ro,choice):
    
    """
    This function adds rows to the 3 matrices to ensure that R&D search has a lattice site that is big enough
    for the given r value.  
    
    """
    m=S.shape[0]
    n=S.shape[1]
    check=max(BPF)+r+1-m  # this is the key value where if greater than zero rows need to be added
    if check>0:   #add rows
        newrows=New(q,n,int(check),p,p_tup,m,Ru,Ro,choice)
        S_update=np.vstack([S,newrows[0]])
        p_tup_update=newrows[1]
        R_update=np.vstack([R,newrows[2]])
        L_update=np.vstack([L,newrows[3]])
        
    else:  #no need to add rows
        S_update=S
        p_tup_update=p_tup
        R_update=R
        L_update=L
    return S_update,R_update,L_update,p_tup_update

        


def BPF(S,BPF_old):
    
    """This function finds the height (row value) for each column in the matrix"""
    m=S.shape[0]
    n=S.shape[1]
    BPF_new=[]
    for j in range(n):
        BPF_new.append(BPF_old[j])
        for i in range(int(BPF_old[j]+1),m): #only need to check from one space above the current BPF to the end of the matrix
                if S[i,j]==2:
                    BPF_new[j]=i
    
                
    return BPF_new #an array of the BPF values


def Search_Index(m,n,BPF_x,BPF_y,r,L,strategy):
    
    """
    This function gives the coordinate values of the diamond search.  It is called by the 'Search' function
    
    """
    col=BPF_y
    index=BPF_x
    num_sites=0
    y_val=[]
    x_val=[]
    if strategy == 'lr':
        for i in range(max(0,col-r),min(col+r+1,n)):
            if i != BPF_y:
                if L[index,i]==-1:
                    y_val.append(i)
                    x_val.append(index)
                    num_sites+=1
    elif strategy == 'ud':
        for a in range(-r,r+1):
            if index+a >=0 and index+a<=m-1:
                if index+a!=BPF_x:
                    if L[(index+a),col] == -1:
                        y_val.append(col)
                        x_val.append(index+a)
                        num_sites+=1
    else:
        for i in range(max(0,col-r),min(col+r+1,n)):
            var=abs(abs(col-i)-r)
            for a in range(-var,var+1):
                if index+a>=0 and index+a<=m-1:
                    if index+a!=BPF_x or i!=BPF_y:
                        if L[(index+a),i]==-1: #Key to the R&D search not searching sites that have already been discovered (i.e. not state -1)
                            y_val.append(i)
                            x_val.append(index+a)
                            num_sites+=1
                            
    return x_val,y_val,num_sites
    
    
def RD(x_val,y_val,num_sites,E,S,R,L):
    
    """
    This function creates the budget per lattice site and then updates the corresponding site of the resistance matrix (R)
    by subtracting the budget amount. If the resistance site becomes negative then the L site is changed to the corresponding S site.
    If that S site is 1 then the coordinate is returned for further testing to see if it should become state 2.
    
    """
    budget=float(E)/num_sites  #ASSUMPTION that the budget is equally distributed amongst all hidden sites in the diamond search

    one_index_x=[]
    one_index_y=[]   
    for a in range(num_sites):
        R[x_val[a],y_val[a]]=R[x_val[a],y_val[a]]-budget
        if R[x_val[a],y_val[a]]<0:   #negativity check
            L[x_val[a],y_val[a]]=S[x_val[a],y_val[a]]
            if L[x_val[a],y_val[a]]==1:

                one_index_x.append(x_val[a])
                one_index_y.append(y_val[a])
                 
    return one_index_x,one_index_y  
    
def twocheck(S,L,twos_i,twos_j,p_tup,p_win,win_col,j,strategy): #this function checks to see if states should be changed from 2-3
    m=S.shape[0]
    n=S.shape[1]
    """ This function takes an i,j input of a lattice site that is known to be
    state 2.  It then looks in a cross pattern (up,down,left,right).  If these sites
    are in state 1 they become state 2.  These values are recorded so that the search can then be
    performed using them.  The sites are also checked to see if they were prizes."""
    count=len(twos_i)
    while count>0:
        i=twos_i[0] #value for the row
        j=twos_j[0] #value for the column (firm)
        for x in range(max(0,i-1),min(m,i+2)): #searches up and down
            if L[x,j]==1:
                S[x,j]=2
                L[x,j]=S[x,j]
                twos_i.append(x)
                twos_j.append(j)
                if (x,j) in p_tup and ((x,j) not in p_win):
                    p_win.append((x,j))
                    win_col.append(j)
                                             
        for y in range(max(0,j-1),min(j+2,n)): #searches left and right
            if L[i,y]==1:
                S[i,y]=2
                L[i,y]=S[i,y]
                twos_i.append(i)
                twos_j.append(y)  
                if (i,y) in p_tup and ((i,y) not in p_win):
                    p_win.append((i,y))
                    win_col.append(j)
                                    
        del twos_i[0]
        del twos_j[0]
        count=len(twos_i) #keeps track of how many sites that are now 2 still need checking
        
    return S,p_tup,p_win,win_col
    
def Search(S,R,L,r,E,fold,p_tup,individual_bankrupt,strategy,concur,check_height): #searches squares within a given radius r to do R&D on.  This R&D effort is given by E and if successful states are changed from 1 to 2
    """
    This is the key function where R&D search is performed.  It takes the coordinates from the 'Search_Index' function and if num_sites
    is greater than zero it continues on to perform the RD function.  After the RD function any states changed from -1 to 1 go onto further
    testing to see if they should become state 2.  If this state 2 happens the twocheck function is then performed. The Search function
    is performed centered around the BPF of each column in the lattice.  This means this function is run n times each period (t).
    """
    m=int(S.shape[0])
    n=int(S.shape[1])
    p_win=[]
    win_col=[]
    order=range(0,n) #gives an array of n column values
    rand.shuffle(order) #randomizes the column order for the R&D search
    
    for j in order:
        BPF_y=j
        BPF_x=fold[j]
        if BPF_x>=0 and individual_bankrupt[j]==0: #ensures that there is a BPF point around with R&D can be conducted and the column has a budget
            # op.writeBPF(BPF_y,False,filename)
            # op.writeBPF(BPF_x,True,filename)
            x_val,y_val,num_sites=Search_Index(m,n,BPF_x,BPF_y,r,L,strategy)
            if concur == True:
                temp_sequence = sorted(rg.generate_sequence(j,n,5))
                for ele in temp_sequence:
                    if fold[ele] >= 0:
                        # op.writeBPF(ele,False,filename)
                        # op.writeBPF(fold[ele],True,filename)
                        x_temp_val,y_temp_val,temp_num_sites = Search_Index(m,n,fold[ele],ele,r,L,strategy)
                        for elex in x_temp_val:
                            x_val.append(elex)
                        for eley in y_temp_val:
                            y_val.append(eley)
                        num_sites = num_sites + temp_num_sites
            '''
            only search for the max bpf
            '''
            if check_height == True:
                temp_sqeuence = sorted(rg.generate_sequence(j,n,5))
                candidate = []
                for ele in temp_sequence:
                    candidate.append((ele,fold[ele]))
                max_candidate = max(candidate, key=lambda item: item[0])
                x_val,y_val,num_sites = Search_Index(m,n,max_candidate[1],max_candidate[0],r,L)
            if num_sites>0:
                one_index=RD(x_val,y_val,num_sites,E[j],S,R,L)
                x_val=one_index[0]
                y_val=one_index[1]
                for v in range(len(one_index[0])):
                    c=x_val[v]
                    d=y_val[v]
                    count=0
                    for x in range(max(0,c-1),min(m,c+2)): #searches up and down
                        if S[x,d]==2:
                            S[c,d]=2
                            L[c,d]=S[c,d]
                            if (c,d) in p_tup and ((c,d) not in p_win):
                                p_win.append((c,d))
                                win_col.append(j)
                            count=1
                    for y in range(max(0,d-1),min(d+2,n)): #searches left and right
                        if S[c,y]==2:
                            S[c,d]=2
                            L[c,d]=S[c,d]
                            if (c,d) in p_tup and ((c,d) not in p_win):
                                p_win.append((c,d))
                                win_col.append(j)
                            count=1
                    if count==1:
                        twos_i=[c]
                        twos_j=[d]
                        Y=twocheck(S,L,twos_i,twos_j,p_tup,p_win,win_col,j,strategy)  #searches like a chain for further changes to state 2
                        S=Y[0]
                        p_tup=Y[1]
                        p_win=Y[2]
                        win_col=Y[3]

    return p_win,win_col #array of tuples corresponding to the 'prizes' discovered during this round of Search and win_col is an array of the columns that were doing search when the prizes were found
    
def win(individual_budget,p_win_old,p_win,pu,po,p_tup,win_col,choice):
    """
    This function updates the budget by  adding the amount of revenue
    gained from the prizes found."""
    
    
    for a in win_col: #a is the column coordinate of columns that found a prize this search period
        dollar=np.random.lognormal(pu,po)
        num=0
        if choice==0:
            dollar=dollar
        elif choice==1: #New to Thesis 5.0 this is what makes the revenue generators worth more the higher up in the lattice they appear.
            dollar=dollar*np.log10(9+p_win[num][0])
            num+=1
        if dollar<0: #ensures a non-negative revenue
            dollar=0
        
        individual_budget[a]+=dollar # only the column that found the prize has its budget increased
        
    for y in p_win:
        p_tup.remove(y)
    for i in p_win:
        p_win_old.append(i)
        
    return individual_budget,p_win_old, p_tup  
    
def lattice(L,names,val,p_win_old,key):
    """
    This plots the Lattice grid.
    
    """
    fig=plt.figure("Lattice")
    
    prize_array=p_win_old
    ax=plt.subplot(111)
    S_new=copy.deepcopy(L)
    for i in prize_array: # Allows for the prize lattices to appear yellow once found
        S_new[i[0],i[1]]=3

    X=np.matrix(S_new)
    m=X.shape[0]
    n=X.shape[1]
    x=[]
    y=[]
        
    for i in range (m+1):
        for j in range (n+1):
               x.append(i)
               y.append(j)
    x=np.reshape(x,(m+1,n+1))
    y=np.reshape(y,(m+1,n+1))
    
    cmap=mpl.colors.ListedColormap(['white','gray','green','blue','yellow'])
    bounds=[-1.5,-.5,.5,1.5,2.5,3.5]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    plt.ylim(0,m)
    plt.xlim(0,n)
    
    #This gives the color blocks seen in the legend below the bottom of the lattice
    
    yellow_patch = mpatches.Patch(facecolor='yellow',label='Revenue Producer',ec='black')
    blue_patch = mpatches.Patch(facecolor='blue', label='Viable',ec='black')
    green_patch = mpatches.Patch(facecolor='green', label='Feasible',ec='black')
    gray_patch=mpatches.Patch(facecolor='gray', label='Not Feasible',ec='black')
    white_patch=mpatches.Patch(facecolor='white', label='Hidden',ec='black')
    
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

    # Put a legend below current axis
    plt.legend(handles=[yellow_patch,blue_patch,green_patch,gray_patch,white_patch],loc='upper center', bbox_to_anchor=(0.5,-.01))
    # txt=''
    # for i in range(1,len(names)):
    #                 txt+=names[i]+' = '+str(val[i])
    #                 if i!=len(names)-1:
    #                     txt+='\n'
                     
    # fig.text(.01,.5,txt,bbox=dict(facecolor='white', ec='black', alpha=1.0),fontsize=15) #This gives the parameters of the run conducted
    plt.pcolormesh(y,x,np.array(X),cmap=cmap,norm=norm,edgecolor='k')
    plt.title(key)
    plt.draw()
    plt.show()
    fig.savefig(key+'.png')  
    
def plot(X,Y,Z,names,val,one,two,zlabel):
    """
    This plots the averaged metrics over the given number of runs.
    """
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X=np.array(X)
    Y=np.array(Y)
    Z=np.array(Z)
    l=np.linspace(min(Z),max(Z),6)
    bounds=np.floor(l)
    cmap=mpl.colors.ListedColormap(['b','g','y','m','c'])    
    '''
    b: blue
    g: green
    r: red
    c: cyan
    m: magenta
    y: yellow
    k: black
    w: white
    '''
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
    txt=''
    for i in range(len(names)):
        if i!=one and i!=two:
            txt+=names[i]+' = '+str(val[i])
            if i != len(names)-1:
                txt+='\n'
    
    
    ax.plot_trisurf(X,Y,Z,cmap=cmap,norm=norm)
    fig.text(0.01,.5,txt,bbox=dict(facecolor='white', ec='black', alpha=0.75))
    ax.set_xlabel(names[one])
    ax.set_ylabel(names[two])
    ax.set_zlabel(zlabel)
    
    plt.show()
    
def line_plot(X,Y,names,val,metric_num,val_metric,y_label):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    plot_num=len(Y)
    
    line_color=['b','r','g','k','y']+['m','c']
    line_marker=['o']*5+['o']*3
    line_array=[]
    x_max=max(X)
    y_max=0
    
    for j in range(plot_num):
        a=max(Y[j])
        if a>y_max:
            y_max=a
    
    txt=''
    for i in range(len(names)):
        if i!=metric_num:
            txt+=names[i]+' = '+str(val[i])
            if i != len(names)-1:
                txt+='\n'
    
    fig.text(0.01,.5,txt,bbox=dict(facecolor='white', ec='black', alpha=0.75))
    
    ax.set_xlabel('Time (t)')
    ax.set_ylabel(y_label)


    line_array=[]
    # red dashes, blue squares and green triangles
    for i in range(plot_num):
        ax.plot(X,Y[i],line_color[i]+line_marker[i],label='r',markersize=5)
        line_array.append(str(val_metric[i]))
    
    plt.legend(line_array,bbox_to_anchor=(1, .6), loc=2, borderaxespad=0.)
    print (str(datetime.now()))
    plt.show()

def Iterate(q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max,choice,strategy,concur,shared,check_height):
    """
    This is the key function that allows the simulation to run. It runs at most t_max iterations of the Search function.  It updates all
    the matrices, BPF, budgets, etc. If ever the budget goes below zero the iteration stops and a 'Bankruptcy' is called.
    """
    """
    These dictionaries are for store the initial lattice for different strategy
    Use dictionaries to store the corresponding matrices, such as S,L,R,...
    As well as necessary checking parameters, such as bankruptcy, t_bankruptcy, ...
    """
    dict_s = {}
    dict_fold = {}
    dict_p_tup = {}
    dict_r = {}
    dict_l = {}
    dict_p_win_old = {}
    dict_individual_b = {}
    dict_total_b = {}
    dict_E = {}
    dict_individual_bankrupt = {}
    dict_Bankrupt = {}
    dict_B = {}
    dict_B_percent = {}
    dict_T = {}
    dict_max_BPF = {}
    dict_t = {}
    dict_percentage = {}
    dict_t_Bankrupt = {}
    dict_B_percent = {}
    for i in range(len(strategy)):
        key_name = 'strategy: '+ strategy[i] + ',concurrence: '+str(concur[i])
        S,Fold,p_tup,R,L=Initialize(q,n,r,p,Ru,Ro,choice)
        dict_s[key_name] = S
        dict_fold[key_name] = Fold
        dict_p_tup[key_name] = p_tup
        dict_r[key_name] = R
        dict_l[key_name] = L
        dict_p_win_old[key_name] = []
        dict_individual_b[key_name] = [initial_b]*n
        total_b = float(initial_b)*n
        dict_total_b[key_name] = [total_b]
        dict_E[key_name] = [initial_b*b_percent]*n
        dict_percentage[key_name] = []
        for j in range(n):
            if dict_E[key_name][j] < E_min and dict_individual_b[key_name][j] >= E_min:
                dict_E[key_name][j] = E_min
            if dict_E[key_name][j] < E_min and dict_individual_b[key_name][j] < E_min:
                dict_E[key_name][j] = dict_individual_b[key_name][j]
        dict_individual_bankrupt[key_name] = [0]*n   # structure: strategy[i]: individual bankrupt array initialized at zero
        dict_Bankrupt[key_name] = 0
        dict_t_Bankrupt[key_name] = -1
        dict_B[key_name] = []
        dict_B[key_name].append(total_b)
        dict_B_percent[key_name] = []
        dict_max_BPF[key_name] = [0]
        dict_B_percent[key_name].append(1)
        dict_t[key_name] = 0 # time step count for different strategies
        dict_T[key_name] = 0
    count=0
    tmpResult = []
    t = 0
    T = []
    T.append(t)
    negone_count = 0
    for key,val in dict_t_Bankrupt.items():
    	if val == -1:
    		negone_count+=1
    while t<t_max and negone_count>1:
        print 'time%s:'%str(t)
        # element in strategy array is the name of the certain strategy
        for i in range(len(strategy)):
            key_name = 'strategy: '+ strategy[i] + ',concurrence: '+str(concur[i])
            if dict_t_Bankrupt[key_name] == -1:
                p_win,win_col=Search(dict_s[key_name],dict_r[key_name],
                                    dict_l[key_name],r,dict_E[key_name],
                                    dict_fold[key_name],dict_p_tup[key_name],
                                    dict_individual_bankrupt[key_name],
                                    strategy[i],concur[i],check_height)
                Frontier=BPF(dict_s[key_name],dict_fold[key_name])
                print key_name,p_win
                for j in range(n):
                    dict_individual_b[key_name][j]=dict_individual_b[key_name][j]-dict_E[key_name][j]
                dict_individual_b[key_name],dict_p_win_old[key_name],dict_p_tup[key_name] = win(dict_individual_b[key_name],
                                                                                                dict_p_win_old[key_name],
                                                                                                p_win,pu,po,dict_p_tup[key_name],win_col,choice)
                
                for c in range(n):
                    if dict_individual_b[key_name][c]<0.0001:
                        dict_individual_bankrupt[key_name][c] = 1
                        dict_individual_b[key_name][c] = 0
                total_b=sum(dict_individual_b[key_name])

                if sum(dict_individual_bankrupt[key_name])==n and count!=1:
                    Bankrupt=1
                    dict_t_Bankrupt[key_name]=t
                
                dict_s[key_name],dict_r[key_name],dict_l[key_name],dict_p_tup[key_name]=Update(dict_s[key_name],dict_r[key_name],dict_l[key_name],Frontier,r,q,p,dict_p_tup[key_name],Ru,Ro,choice)
                dict_fold[key_name]=Frontier
                if shared == True:
                    '''
                    sum the total available budget from the previous iteration
                    and put it into the differernt evenly.
                    '''
                    individual_amount = total_b/len(individial_b)
                    for z in range(len(individial_b)):
                        dict_individual_b[key_name][z] = individual_amount
                for a in range(n):
                    dict_E[key_name][a]=dict_individual_b[key_name][a]*b_percent
                    if dict_E[key_name][a]<E_min and dict_individual_b[key_name][a]>=E_min:
                        dict_E[key_name][a]=E_min
                    if dict_E[key_name][a]<E_min and dict_individual_b[key_name][a]<E_min:
                        dict_E[key_name][a]=dict_individual_b[key_name][a]

                dict_B[key_name].append(total_b)
                dict_total_b[key_name].append(total_b)
       	if os.path.isfile('budget.txt'):
       		with open ('budget.txt','a') as f:
       			f.write('\n')
       			f.write('time%s'%str(t))
       			for key,val in dict_total_b.items():
       				f.write('\n')
       				f.write(key+':'+str(val[-1]))
       	else:
       		with open ('budget.txt','w') as f:
       			f.write('time%s'%str(t))
       			for key,val in dict_total_b.items():
       				f.write('\n')
       				f.write(key+':'+str(val[-1]))
	    # calcualte the denominator of the ratio
        for k,v in dict_total_b.items():
            print type(dict_total_b[k])
        currentSum=0
        print 'cur: ',type(currentSum)
        for k,v in dict_total_b.items():
			currentSum += dict_total_b[k][-1]
        factor = 1.0
        temp_increment = {}
        for k,v in dict_total_b.items():
        	temp_increment[k] = v[-1]-v[-2]
        max_val = max(temp_increment.values())
        denominator = 0
        min_allocation_ratio = 0.1
        for k,v in temp_increment.items():
        	denominator+=math.exp(factor*(temp_increment[k]-max_val))
        for k,v in dict_percentage.items():
        	dict_percentage[k] = math.exp(factor*(temp_increment[k]-max_val))/denominator
        if t < 30:
            min_count = 0
            for k,v, in dict_percentage.items():
                if dict_percentage[k] < min_allocation_ratio:
                    dict_percentage[k] = min_allocation_ratio
                    min_count+=1
            for k,v in dict_percentage.items():
                if dict_percentage[k] > min_allocation_ratio:
                    dict_percentage[k] = (1-min_allocation_ratio*min_count)/(len(strategy)-min_count)
        print dict_percentage
        for k,v in dict_individual_b.items():
            for i in range(len(dict_individual_b[k])):
        	   dict_individual_b[k][i] = currentSum * dict_percentage[k]/n
        negone_count = 0
        for k,v in dict_t_Bankrupt.items():
        	if v == -1:
        		negone_count+=1
        t+=1
    return dict_l,dict_t_Bankrupt,dict_p_win_old
    
    
def runs(runs,q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max,choice,strategy,concur,shared,check_height):
    """This function simply performs multiple runs of the Iteration function.  It also performs some statistics, basically
    averaging the key metrics over the number of runs"""
    max_BPF=[]
    max_BPF_array_array=[]
    average_BPF=[]
    last_t=[]
    
    B_percent_final=[]
    
    Bankrupt_num=0
    t_Bankrupt_array=[]
    
    p_win_old_array=[]
    
    B_percent_array=[]
    '''
    holding the total budget left in the different time period: Yijun
    '''
    budget_x = []
    budget_y = []
    x=0
    while x<runs:
        L,T,B,B_percent,Bankrupt,t_Bankrupt,Frontier,p_win_old,max_BPF_array,budget_x,budget_y,filename=Iterate(q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max,choice,strategy,concur,shared,check_height,budget_x,budget_y)
        '''
        added by Yijun Dai. Making graph
        '''
        left_budget = budget_y[-1]
        if runs == 1:
            makegraph.makegraph(budget_x,budget_y,filename+"_money.png")

        max_BPF_val=max(Frontier)
        average_BPF_val=np.average(Frontier)
        
        B_percent_final_val=B_percent[-1]
        
        last_t_val=T[-1]
        
        Bankrupt_num+=Bankrupt
        
        
        
        #Appending to Array
        max_BPF.append(max_BPF_val)
        average_BPF.append(average_BPF_val)
        last_t.append(last_t_val)
        B_percent_final.append(B_percent_final_val)
        t_Bankrupt_array.append(t_Bankrupt)
        p_win_old_array.append(p_win_old)
        B_percent_array.append(B_percent)
        max_BPF_array_array.append(max_BPF_array)
        
        x+=1

        print "Run %s complete!" % str(x)
        
    #Average Values of the Arrays  
    max_BPF_average=np.average(max_BPF)
    average_BPF_average=np.average(average_BPF)
  
    B_percent_final_average=np.average(B_percent_final)
    last_t_average=np.average(last_t)
    total=[]
    for i in range(len(t_Bankrupt_array)):
        if t_Bankrupt_array[i]!=-1:
            total.append(t_Bankrupt_array[i])
    if len(total)>0:
        t_Bankrupt_array_average=np.average(total)
    else:
        t_Bankrupt_array_average=[]
    
    Bankrupt_percent=float(Bankrupt_num)/runs*100
    if runs == 1:
        return L,max_BPF_average, average_BPF_average, last_t_average, B_percent_final_average, Bankrupt_num, Bankrupt_percent, t_Bankrupt_array_average,p_win_old_array,B_percent_array,max_BPF_array_array,left_budget
    else:
        return L,max_BPF_average, average_BPF_average, last_t_average, B_percent_final_average, Bankrupt_num, Bankrupt_percent, t_Bankrupt_array_average,p_win_old_array,B_percent_array,max_BPF_array_array
            
""" 
    
q=.65
n=60
r=5
p=.5
Ru=250
Ro=5
pu=1000
po=5
initial_b=1500
b_percent=.2
E_min=25
t_max=20


num=1

Y=runs(num,q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max)
    
for i in range(len(Y)):
    print Y[i]
"""

def Value_String(choices,names,graph_choice):
    """
    This function is key to allowing user entry of numbers to be interpreted as corresponding to a parameter.  It takes in the 
    users 'choices' array and creates arrays of multiple values for the runs function to perform on.  The user enters a beginning and ending
    value and the number of values to create.
    
    """
    X=[0]
    for x in range(1,len(names)):
        if x in choices:
            x_start=float(raw_input("What is the starting "+ names[x]+ " value? "))
            x_end=float(raw_input("What is the ending "+ names[x]+ " value? "))
            a=False
            while a==False:
                try:
                    x_num=float(raw_input("How many values of " + names[x]+" will you be testing? "))
                except ValueError:
                    print('That is not a valid input!')
                if graph_choice==2:
                    a=True
                elif graph_choice==1 and x_num<=10 and x_num>0:
                    a=True
                else:
                    print('Please enter a number between 1-10 for the number of values to test!')
                    
                        
            if x==2 or x==3 or x==12: # These refer to the index of names that need to be integer values
                Y=np.linspace(x_start,x_end,x_num) #line that creates the array
                Y=np.floor(Y) #only because an integer value is necessary
                Y=Y.astype(int)
                X.append(Y)
            else:
                Y=np.linspace(x_start,x_end,x_num)
                X.append(Y)
        elif x==2 or x==3 or x==12:
            Y=int(raw_input("What is the constant value for "+ names[x] +" ? "))
            X.append(Y)
        else:
            Y=float(raw_input("What is the constant value for "+ names[x] +" ? "))
            X.append(Y)
    return X

def Percolation():
    """Function combining all of the previous into a nicer user interface.  It asks for user inputs for the parameters."""
    
    print ("\nWelcome to the Percolation Simulator!\n")
    # These variables are for getting values from the Config.py
    Conf = Config()
    pre_run = Conf.num_run
    pre_q = Conf.q_val
    pre_n = Conf.num_col
    pre_r = Conf.search_radius
    pre_p = Conf.p_val
    pre_Ru = Conf.ru_val
    pre_Ro = Conf.ro_val
    pre_pu = Conf.pu_val
    pre_po = Conf.po_val
    pre_initial_b = Conf.initial_b
    pre_b_percent = Conf.b_percent
    pre_E_min = Conf.e_min
    pre_t_max = Conf.t_max
    epoch = 10
    d = False

    # msg string used to warn users during program processing
    int_msg = "Sorry that input from Configuration (Config.py) was not valid!" 
    "Please enter an integer greater than zero! "
    zero_msg = "Please input an integer greater than zero"

    while d == False:
        prob_choice = int(raw_input("How do you want the revenue generators within the lattice to be distributed?\n\n0.) Evenly Distributed\n1.) More Dense at the Bottom and More Valuable at the Top\n"))
        names = ['num_runs','q','n','r','p','Ru','Ro','pu','po','initial_b','b_percent','E_min', 't_max']
        while True:
            try: 
                num_runs = int(pre_run)
            except ValueError:
                print int_msg
                continue
            if num_runs <= 0:
                print zero_msg
                continue
            else:
                break
        while True:
            try:
                q = float(pre_q)
            except ValueError:
                print "Sorry that input was not valid! Please enter a real number between "
                "zero and one. "
                continue
            if q < 0 or q > 1:
                print "Please input a probability between zero and one! "
                continue
            else:
                break
        while True:
            try:
                n = int(pre_n)
            except ValueError:
                print int_msg
                continue
            if n <= 0:
                print zero_msg
                continue
            else:
                break
        while True:
            try:
                r = int(pre_r)
            except ValueError:
                print ing_msg
                continue
            if r <= 0:
                print zero_msg
                continue
            else:
                break
        while True:
            try:
                p = float(pre_p)
            except ValueError:
                print "Sorry that input from Configuration (Config.py) wasn't valid! Please enter a real number between"
                "zero and one!"
                continue
            if p < 0 or p > 1:
                print "Please input a probability between zero and one! "
                continue
            else:
                break
        while True:
            try:
                Ru = float(pre_Ru)
            except ValueError:
                print "Sorry that input from Configuration (Config.py) was not valid! "
                "Please enter a real number greater than zero!"
                continue
            if Ru < 0:
                print "Please input a positive number! "
                continue
            else:
                break
        while True:
            try:
                Ro = float(pre_Ro)
            except ValueError:
                print "Sorry taht input from Configuration (Config.py) was not valid! "
                "Please enter a real number greater than zero! "
                continue
            if Ro < 0:
                print "Please input a positive number! "
                continue
            else:
                break
        while True:
            try:
                pu = float(pre_pu)
            except ValueError:
                print "Sorry that input from Configuration (Config.py) was not valid! "
                "Please enter a real number greater than zero! "
                continue
            if pu < 0:
                print "Please input a positive number! "
                continue
            else:
                break
        while True:
            try:
                po = float(pre_po)
            except ValueError:
                print "Sorry that input from Configuration (Config.py) was not valid! "
                "Please enter a real number greater than zero! "
                continue
            if po < 0:
                print "Please input a positive number! "
                continue
            else:
                break
        while True:
            try:
                initial_b = float(pre_initial_b)
            except ValueError:
                print "Sorry that input from Configuration (Config.py) was not valid! "
                "Please enter a real number greater than zero! "
                continue
            if initial_b < 0:
                print "Please input a positive number! "
                continue
            else:
                break
        while True:
            try:
                b_percent = float(pre_b_percent)
            except ValueError:
                print "Sorry that input from Configuration (Config.py) was not valid! "
                "Please enter a real number between zero and one! "
                continue
            if b_percent < 0 or b_percent > 1:
                print "Please input a valid percentage! "
                continue
            else:
                break
        while True:
            try:
                E_min = float(pre_E_min)
            except ValueError:
                print "Sorry that input from Configuration (Config.py) was not valid! "
                "Please enter a real number between zero and initial budget! ", 
                str(initial_b), "!"
                continue
            if E_min < 0 or E_min > initial_b:
                print "Please input a number between zero and the initial budget ",str(initial_b),"!"
                continue
            else:
                break
        while True:
            try:
                t_max = int(pre_t_max)
            except ValueError:
                print "Sorry that input from Configuration (Config.py) was not valid! "
                "Please enter a real number greater than zero! "
                continue
            if t_max <= 0:
                print "Please input an integer gerater than zero! "
            else:
                break
        strategy = ['all','ud','lr','all']
        concur = [False,False,False,True]
        # temporary use: shared = False, check_height = False
        shared = False
        check_height = False
        result_dict_l,result_bankrupt,result_p_win_old = Iterate(q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max,prob_choice,strategy,concur,shared,check_height)
        print result_p_win_old
        val = [num_runs,q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max]
        for key,val in result_dict_l.items():
            lattice(result_dict_l[key],names,val,result_p_win_old[key],key)
        d = True

Percolation()
    
