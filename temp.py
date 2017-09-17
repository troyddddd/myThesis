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


def Search_Index(m,n,BPF_x,BPF_y,r,L):
    
    """
    This function gives the coordinate values of the diamond search.  It is called by the 'Search' function
    
    """
    col=BPF_y
    index=BPF_x
    num_sites=0
    y_val=[]
    x_val=[]
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
    
def twocheck(S,L,twos_i,twos_j,p_tup,p_win,win_col,j): #this function checks to see if states should be changed from 2-3
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
    
def Search(S,R,L,r,E,fold,p_tup,individual_bankrupt,strategy,con_distance): #searches squares within a given radius r to do R&D on.  This R&D effort is given by E and if successful states are changed from 1 to 2
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
            x_val,y_val,num_sites=Search_Index(m,n,BPF_x,BPF_y,r,L)
            if con_distance > 0 and j+con_distance < n:
                x_temp_val,y_temp_val,temp_num_sites = Search_Index(m,n,fold[j+con_distance],j,r,L)
                for elex in x_temp_val:
                    x_val.append(elex)
                for eley in y_temp_val:
                    y_val.append(eley)
                num_sites = num_sites + temp_num_sites
            if num_sites>0:
                one_index=RD(x_val,y_val,num_sites,E[j],S,R,L)
                x_val=one_index[0]
                y_val=one_index[1]
                for v in range(len(one_index[0])):
                        c=x_val[v]
                        d=y_val[v]
                        count=0
                        if strategy == 'ud' or strategy == 'all':
                            for x in range(max(0,c-1),min(m,c+2)): #searches up and down
                                if S[x,d]==2:
                                    S[c,d]=2
                                    L[c,d]=S[c,d]
                                    if (c,d) in p_tup and ((c,d) not in p_win):
                                        p_win.append((c,d))
                                        win_col.append(j)
                                    count=1
                        if strategy == 'lr' or strategy == 'all':
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
                            Y=twocheck(S,L,twos_i,twos_j,p_tup,p_win,win_col,j)  #searches like a chain for further changes to state 2
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
    
def lattice(L,names,val,p_win_old,con_distance,strategy):
    """
    This plots the Lattice grid.
    
    """
    fig=plt.figure("Lattice")
    
    prize_array=p_win_old[0]
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
    txt=''
    for i in range(1,len(names)):
                    txt+=names[i]+' = '+str(val[i])
                    if i!=len(names)-1:
                        txt+='\n'
                     
    fig.text(.01,.5,txt,bbox=dict(facecolor='white', ec='black', alpha=1.0),fontsize=15) #This gives the parameters of the run conducted
    plt.pcolormesh(y,x,np.array(X),cmap=cmap,norm=norm,edgecolor='k')
    temps= ""
    temps = temps + "Distance: " + str(con_distance) + "Direction: " + str(strategy)
    plt.title(temps)
    #plt.draw()
    plt.show()  
    
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
    
def Iterate(q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max,choice,strategy,con_distance):
    """
    This is the key function that allows the simulation to run. It runs at most t_max iterations of the Search function.  It updates all
    the matrices, BPF, budgets, etc. If ever the budget goes below zero the iteration stops and a 'Bankruptcy' is called.
    """
    
    S,Fold,p_tup,R,L=Initialize(q,n,r,p,Ru,Ro,choice)
    p_win_old=[]
  
    individual_b=[initial_b]*n #difference from thesis 3.0 ...  each column now has a budget separate from the other columns
    total_b=float(initial_b)*n  #key change from thesis 3.0
    
    total_initial_b=total_b
   
    
    E=[initial_b*b_percent]*n
    for i in range(n):
        if E[i]<E_min and individual_b[i]>=E_min:
                E[i]=E_min
        if E[i]<E_min and individual_b[i]<E_min:
                E[i]=individual_b[i]
    individual_bankrupt=[0]*n
    Bankrupt=0
    t_Bankrupt=-1
    t=0
    
    B=[]
    B.append(total_b)
    T=[]
    T.append(t)
    B_percent=[]
    max_BPF=[0]
    B_percent.append(1)
    count=0
    
    while t<t_max and t!=-1:
        
        p_win,win_col=Search(S,R,L,r,E,Fold,p_tup,individual_bankrupt,strategy,con_distance)
       
        Frontier=BPF(S,Fold)
        
        for j in range(n):
            individual_b[j]=individual_b[j]-E[j]

        individual_b,p_win_old,p_tup=win(individual_b,p_win_old,p_win,pu,po,p_tup,win_col,choice)
        
        total_b=sum(individual_b)
        
        
        for c in range(n):
            if individual_b[c]<0.0001:
                individual_bankrupt[c]=1
        
        if sum(individual_bankrupt)==n and count!=1:
            Bankrupt=1
            t_Bankrupt=t
            t=-2
            
        S,R,L,p_tup=Update(S,R,L,Frontier,r,q,p,p_tup,Ru,Ro,choice)
           
            
        Fold=Frontier
        
        for a in range(n):
            E[a]=individual_b[a]*b_percent
            if E[a]<E_min and individual_b[a]>=E_min:
                E[a]=E_min
            if E[a]<E_min and individual_b[a]<E_min:
                E[a]=individual_b[a]
            
        t+=1
        B.append(total_b)
        T.append(t)
        B_percent.append(total_b/total_initial_b)
        max_BPF.append(max(Fold))
        
    if t==-1:
            T[-1]=t_Bankrupt
            
    """
    L - Final Lattice Structure
    T - Array of t values
    B - Array of budget values
    B_percent - Array of budget proprotion values compared to the initial budget
    Bankrupt - Binary variable for showing 'Bankruptcy' (0 - Not Bankrupt)  (1 - Bankrupt)
    t_Bankrupt - gives the time when Bankruptcy occurred (if no Bankruptcy variable defaults to -1)
    Frontier - Final BPF array
    p_win_old - Array of tuples showing all of the coordinates where prizes were found
    
    """
    return L,T,B,B_percent,Bankrupt,t_Bankrupt,Frontier,p_win_old,max_BPF
    
    
def runs(runs,q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max,choice,strategy,con_distance):
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
   
    x=0
    while x<runs:
        L,T,B,B_percent,Bankrupt,t_Bankrupt,Frontier,p_win_old,max_BPF_array=Iterate(q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max,choice,strategy,con_distance)
        
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
    d=False
    while d==False:
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
                while True:
                    try:
                        pu=float(raw_input("What is the mean of the prize value (pu)?"))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter a real number greater than zero!")
                        continue
                    if pu<0:
                        print ("Please input a positive number!")
                        continue
                    else:
                        break
                while True:
                    try:
                        po=float(raw_input("What is the standard deviation of the prize value (po)?"))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter a real number greater than zero!")
                        continue
                    if po<0:
                        print ("Please input a positive number!")
                        continue
                    else:
                        break
                while True:
                    try:
                        initial_b=float(raw_input("What is the starting budget (per column) for the company (initial_b)?"))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter a real number greater than zero!")
                        continue
                    if initial_b<0:
                        print ("Please input a positive number!")
                        continue
                    else:
                        break
                while True:
                    try:
                        b_percent=float(raw_input("What percentage of the remaining budget do you want given to R&D search each period (b_percent)?"))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter a real number between zero and one!")
                        continue
                    if b_percent<0 or b_percent>1:
                        print ("Please input a percentage between 0 and 1!")
                        continue
                    else:
                        break
                while True:
                    try:
                        E_min=float(raw_input("What is the minimum amount of budget that will be given to each lattice site every R&D period (E_min)?"))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter a real number between zero and the inital budget ",initial_b,"!")
                        continue
                    if E_min<0 or E_min>initial_b:
                        print ("Please input a number between zero and the inital budget ",initial_b,"!")
                        continue
                    else:
                        break
                while True:
                    try:
                        t_max=int(raw_input("What is the maximum number of iterations (t_max) you want each run to conudct?"))
                    except ValueError:
                        print ("Sorry that input wasn't valid! Please enter an integer greater than zero!")
                        continue
                    if t_max<=0:
                        print ("Please input an integer greater than zero")
                        continue
                    else:
                        break
                val = [num_runs,q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max]
                val1 = [num_runs,q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max]
                val2 = [num_runs,q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max]
                Z = runs(num_runs,q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max,prob_choice,'all',0)
                Z1 = runs(num_runs,q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max,prob_choice,'ud',0)
                Z2 = runs(num_runs,q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max,prob_choice,'lr',0)
                Z3 = runs(num_runs,q,n,r,p,Ru,Ro,pu,po,initial_b,b_percent,E_min,t_max,prob_choice,'all',20)
                if num_runs==1:
                    o=False
                    while o==False:
                        lat_choice=raw_input("Would you like to display a graph of the final lattice structure? (Y or N)")
                        if lat_choice.lower()=='y' or lat_choice.lower()=='yes':
                            lattice(Z[0],names,val,Z[8],0,'all')
                            lattice(Z1[0],names,val,Z1[8],0,'up-down')
                            lattice(Z2[0],names,val,Z2[8],0,'left-right')
                            lattice(Z3[0],names,val,Z3[8],20,'all')
                            o=True
                        elif lat_choice.lower()=='no' or lat_choice.lower()=='n':
                            o=True
                        else:
                            print("Not a valid entry! Please enter yes (Y) or no (N). ")
                print "Statistics"
                print "Max_BPF_Average all direction = ",Z[1]
                print "Max_BPF_Average up and down = ", Z1[1]
                print "Max_BPF_Average left and right = ", Z2[1]
                print "Average Percent of Initial Budget Remaining all direction = ",Z[4]*100
                print "Average Percent of Initial Budget Remaining up and down = ",Z1[4]*100
                print "Average Percent of Initial Budget Remaining left and right = ",Z2[4]*100
                print "All Direction: Number of Times Gone Bankrupt = ",Z[5], " This means there was a ",Z[6],"% Bankruptcy rate!"
                print "Up and down: Number of Times Gone Bankrupt = ",Z1[5], " This means there was a ",Z1[6],"% Bankruptcy rate!"
                print "Left and right: Number of Times Gone Bankrupt = ",Z2[5], " This means there was a ",Z2[6],"% Bankruptcy rate!"
                d=True
        elif choice.lower()=='y' or choice.lower=='yes':
            s=False
            graph_choice=0
            while s==False:
                try:
                    graph_choice=int(raw_input('Would you like to run tests changing one variable (enter 1) or two variables (enter 2)?'))
                except ValueError:
                    print ("That's not a valid entry! Please enter 1 or 2!")
                if graph_choice==1 or graph_choice==2:
                    s=True
                else:
                    print('Please enter 1 for changing a single variable or 2 for changing two variables across multiple runs!')
            if graph_choice==2:       
                print("""
    Please Select 2 of the following metrics to change values across (Select the Number to the Right of the Arrow)!
    
    
    Percolation Probability (q) -> 1
    Number of firms (n) -> 2
    Search radius (r) -> 3
    Prize probability (p) -> 4
    Mean of the resistance matrix (Ru) -> 5
    Standard deviation of the resistance matrix (Ro) -> 6
    Mean of the prizes (pu) -> 7
    Standard deviation of the prizes (po) -> 8
    Initial budget (initial_b) -> 9
    Percent of budget for R&D each period (b_percent) -> 10
    R&D budget minimum (E_min) -> 11
    Maximum number of iterations 't' per run (t_max) -> 12
    """)
                e=False
                f=False
                while e==False:
                    try:
                        one=int(raw_input("1st Choice? (x-coordinate) "))
                    except ValueError:
                        print("That's not an integer!")
                        one=0
                    if one>=1 and one<=12:
                        e=True
                    else:
                        print("Please enter a valid number from 1-12")
                while f==False:
                    try:
                        two=int(raw_input("2nd Choice? (y-coordinate) "))
                    except ValueError:
                        print("That's not an int!")
                        two=0
                    if two>=1 and two<=12 and two!=one:
                        f=True
                    else:
                        print("Please enter a valid number from 1-12 and not %i which you already picked! "%(one))
                names=['num_runs','q','n','r','p','Ru','Ro','pu','po','initial_b','b_percent','E_min', 't_max']
                X=Value_String([one,two],names,graph_choice)
                x_num=[]
                y_num=[]
                z_num=[]
                z_num2=[]
                z_num3=[]
                
                num=False
                while num==False:
                    try:
                        num_runs=int(raw_input("How many runs would you like done with these variables? "))
                        num=True
                    except ValueError:
                        print("That is not an integer! Please enter an integer number of runs to perform ")
                    X[0]=num_runs
                    A=copy.deepcopy(X)
                    """
                    print("What would you like the variable for the plot height to be?\n")
                    
                    w=False
                    while w==False:
                            try:
                                z_choice=int(raw_input("1. Maximum Height Achieved\n2. Percent of Initial Budget Remaining\n3. Percent of Runs that were Bankruptcy"))
                            except ValueError:
                                print("That was not an integer value!")
                            if z_choice==1:
                                zlabel='Maximum Height Achieved'
                                mark=1
                                w=True
                            elif z_choice==2:
                                zlabel='% of Initial Budget Remaining'
                                mark=4
                                w=True
                            elif z_choice==3:
                                zlabel='% Bankruptcy'
                                mark=6
                                w=True
                            else:
                                print('Please enter a number between 1 and 3')  
                    """
                    for j in X[one]:
                        A[one]=j
                        for z in X[two]:
                            A[two]=z
                            print (names[one]+' = '+str(j)+'  '+names[two]+' = '+str(z))
                            x_num.append(j)
                            y_num.append(z)
                            Z=runs(A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],A[9],A[10],A[11],A[12],prob_choice)
                            z_num.append(Z[1])
                            z_num2.append(Z[4])
                            z_num3.append(Z[6])
                            
                    plot(x_num,y_num,z_num,names,X,one,two,'Maximum Height (BPF) Achieved')
                    plot(x_num,y_num,z_num2,names,X,one,two,'Proportion of Initial Budget Remaining')
                    plot(x_num,y_num,z_num3,names,X,one,two,'Bankruptcy Percent')
                    print (str(datetime.now()))
            if graph_choice==1:
                print("""
    Please Select 1 of the following metrics to change values across (Select the Number to the Right of the Arrow)!
    
    
    Percolation Probability (q) -> 1
    Number of firms (n) -> 2
    Search radius (r) -> 3
    Prize probability (p) -> 4
    Mean of the resistance matrix (Ru) -> 5
    Standard deviation of the resistance matrix (Ro) -> 6
    Mean of the prizes (pu) -> 7
    Standard deviation of the prizes (po) -> 8
    Initial budget (initial_b) -> 9
    Percent of budget for R&D each period (b_percent) -> 10
    R&D budget minimum (E_min) -> 11
    Maximum number of iterations 't' per run (t_max) -> 12
    """)
                e=False
                while e==False:
                    try:
                        one=int(raw_input("Choice? The variable over which the values will change. "))
                    except ValueError:
                        print("That's not an integer!")
                        one=0
                    if one>=1 and one<=12:
                        e=True
                    else:
                        print("Please enter a valid number from 1-12")   
                names=['num_runs','q','n','r','p','Ru','Ro','pu','po','initial_b','b_percent','E_min', 't_max']       
                X=Value_String([one],names,graph_choice)
                t_max=X[-1]
                num=False
                while num==False:
                    try:
                        num_runs=int(raw_input("How many runs would you like done with these variables? "))
                        num=True
                    except ValueError:
                        print("That is not an integer! Please enter an integer number of runs to perform ")
                    X[0]=num_runs
                    A=copy.deepcopy(X)
                    B_percent=[]
                    Max_Height=[]
                    for j in X[one]:
                        A[one]=j
                        print (names[one]+' = '+str(j))
                        
                        Z=runs(A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],A[9],A[10],A[11],A[12],prob_choice)
                        Max_Height.append(Z[-1])
                        B_percent.append(Z[-2])
                
                
                T=range(t_max+1)
                x_num=X[one]
             
                for a in range(len(x_num)):
                    for b in range(num_runs):
                        diff=t_max-len(B_percent[a][b])
                        if diff>0:
                            B_percent[a][b]=B_percent[a][b]+[0]*(diff+1)
                            Max_Height[a][b]=Max_Height[a][b]+[Max_Height[a][b][-1]]*(diff+1)
                            
                            
                B_percent_Average_array=[]
                Max_Height_Average_array=[]
                for d in range(len(x_num)):
                    B_percent_Average=[]
                    Max_Height_Average=[]
                    for c in range(t_max+1):
                        total=0
                        total_2=0
                        for f in range(num_runs):
                            total=total+B_percent[d][f][c]
                            total_2=total_2+Max_Height[d][f][c]
                        total_average=float(total)/num_runs
                        total_2_average=float(total_2)/num_runs
                        B_percent_Average.append(total_average)
                        
                        Max_Height_Average.append(total_2_average)
                        
                    B_percent_Average_array.append(B_percent_Average)
                    
                    Max_Height_Average_array.append(Max_Height_Average)
                
                
                line_plot(T,B_percent_Average_array,names,X,one,x_num,'Proportion of Initial Budget Remaining')
                line_plot(T,Max_Height_Average_array,names,X,one,x_num,'Max Height')
                
            d=True
        else:
            print("Sorry that input wasn't valid! Please enter Y or N! ")

Percolation()
    
