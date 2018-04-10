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

