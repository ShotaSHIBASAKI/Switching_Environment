#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 15:52:34 2021

@author: shibasakishota
"""
import numpy as np
import scipy as sp
from scipy.integrate import solve_ivp #ode integration
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def FigA5():
    """
    Two species abundances under EFs alone
    """
    R=np.array([200, 50])  # some of resource and toxin inflows
    alpha =0.1
    Tf=400
    delta_array=np.array([0.1, 0.2,0.4,0.6, 0.8])
    nu_array=np.array([-3,-1, 1])
    iteration=8
    for d in range(np.size(delta_array)):
        delta=delta_array[d]
        fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(12,4))
        for v in range(np.size(nu_array)):
            
            if nu_array[v]<0:
                nu=np.power(0.1, -nu_array[v])
            else:
                nu=np.power(10, nu_array[v])
            S1=[]
            S2=[]
            for i in range(iteration):
                s1=[]
                s2=[]
                Ti=0
                T_switch=np.random.exponential(1/nu)
                #print(T_switch)
                init=np.array([150, 100, 10, 10])
                if i==0:
                    state=0
                elif i==1:
                    state=1
                else:
                    if np.random.rand()<=0.5:
                        state=0
                    else:
                        state=1
                    
                while T_switch<Tf:
                    # run ode
                    #print(init)
                    sol = solve_ivp(ODE, [Ti, T_switch], y0=init, method='LSODA',
                                    args=[R[state], alpha, delta]) # solving ode
                    # update T_switch
                    Ti=T_switch
                    T_switch += np.random.exponential(1/nu)
                    state=1-state
                    s1.extend(sol.y[2, :])
                    s2.extend(sol.y[3, :])
                    init=sol.y[:,-1].reshape(4)

                # run final simulation
                sol = solve_ivp(ODE, [Ti, Tf], y0=init, method='LSODA',
                                    args=[R[state], alpha, delta]) # solving ode
                s1.extend(sol.y[2, :])
                s2.extend(sol.y[3, :])
                S1.append(s1)
                S2.append(s2)
            #print((delta, nu))
            if v==0:
                for j in range(len(S1)):
                    ax1.plot(S1[j], S2[j], linewidth=2)
                ax1.set_xlim(0, 125)
                #ax1.set_xlabel('abundance of species 1', fontsize=20)
               
                ax1.set_ylim(0, 60)
                ax1.set_ylabel('abundance of species 2', fontsize=20)
                tname=str('delta=%.1f, nu= 10 ^ %d' %(delta, nu_array[v]))
                ax1.set_title(tname, fontsize=16)
            elif v==1:
                for j in range(len(S1)):
                    ax2.plot(S1[j], S2[j], linewidth=2)
                ax2.set_xlim(0, 125)
                ax2.set_xlabel('abundance of species 1', fontsize=20)
                ax2.set_ylim(0, 60)
                #ax2.set_ylabel('abundance of species 2', fontsize=20)
                tname=str('delta=%.1f, nu= 10 ^ %d' %(delta, nu_array[v]))
                ax2.set_title(tname, fontsize=16)
            else:
                for j in range(len(S1)):
                    ax3.plot(S1[j], S2[j], linewidth=2)
                ax3.set_xlim(0, 125)
                #ax3.set_xlabel('abundance of species 1', fontsize=20)
                
                ax3.set_ylim(0, 60)
                #ax3.set_ylabel('abundance of species 2', fontsize=20)
                tname=str('delta=%.1f, nu= 10 ^ %d' %(delta, nu_array[v]))
                ax3.set_title(tname, fontsize=16)
        fname=str('OnlyEF_abundances_delta%.1f' %(delta))
        plt.savefig(fname+'_nolines.pdf', )
        ax1.vlines(x=5, ymin=0, ymax=60, color='k', linestyle='--')
        ax1.vlines(x=10, ymin=0, ymax=60, color='k', linestyle='--')
        ax2.vlines(x=5, ymin=0, ymax=60,color='k', linestyle='--')
        ax2.vlines(x=10, ymin=0, ymax=60,color='k', linestyle='--')
        ax3.vlines(x=5, ymin=0, ymax=60, color='k', linestyle='--')
        ax3.vlines(x=10, ymin=0, ymax=60, color='k', linestyle='--')
        ax1.hlines(y=5, xmin=0, xmax=125, color='k', linestyle='--')
        ax1.hlines(y=10, xmin=0, xmax=125, color='k', linestyle='--')
        ax2.hlines(y=5, xmin=0, xmax=125,color='k', linestyle='--')
        ax2.hlines(y=10, xmin=0, xmax=125,color='k', linestyle='--')
        ax3.hlines(y=5, xmin=0, xmax=125, color='k', linestyle='--')
        ax3.hlines(y=10, xmin=0, xmax=125, color='k', linestyle='--')
        plt.savefig(fname+'_withlines.pdf', bbox_inches='tight',pad_inches=0.05)
        plt.show()
            


def FigA6(iteration=1000, delta=0.2):
    """
    PDMP and showing distribution of total population size n and \hat{n}
    """
    R=np.array([200, 50])  # some of resource and toxin inflows
    alpha =0.1
    Tf=200
    nu_array=np.array([-3,-1, 1])
    for v in range(np.size(nu_array)):
        if nu_array[v]<0:
            nu=np.power(0.1, -nu_array[v])
        else:
            nu=np.power(10, nu_array[v])
        result=[]
        for i in range(iteration):
            Ti=0
            T_switch=np.random.exponential(1/nu)
            #print(T_switch)
            init=np.array([150, 100, 10, 10])
            if np.random.rand()<=0.5:
                state=0
            else:
                    state=1
                
            while T_switch<Tf:
                # run ode
                #print(init)
                sol = solve_ivp(ODE, [Ti, T_switch], y0=init, method='LSODA',
                            args=[R[state], alpha, delta]) # solving ode
                # update T_switch
                Ti=T_switch
                T_switch += np.random.exponential(1/nu)
                state=1-state
                init=sol.y[:,-1].reshape(4)
                
            # run final simulation
            sol = solve_ivp(ODE, [Ti, Tf], y0=init, method='LSODA',
                            args=[R[state], alpha, delta]) # solving ode
            result.append(sol.y[:,-1])
        np.savetxt(str('PDMP_ver2_nu%d.csv' %(nu_array[v])), result,
                       delimiter=',',header='r1,t1,s1,s2', fmt='%.3f')
  
def Plot(nu_array=np.array([-3, -1, 1])):
    #ymax=[500, 40, 300]
    y_top=[0.005, 0.00045, 0.0035]
    for i in range (np.size(nu_array)):
        data=np.loadtxt(str('PDMP_ver2_nu%d.csv' %(nu_array[i])), delimiter=',',
                        skiprows=1)
        n=np.sum(data, 1)
        hat_n=data[:,0] - data[:, 1] + data[:,2] + data[:, 3]
        df=pd.DataFrame({"n": n, "hat_n":hat_n})
        sns.set(rc={'axes.labelsize':20,
            'xtick.labelsize':16,
            'ytick.labelsize':16
                })
        ax=sns.kdeplot(x=df.n, y=df.hat_n, cmap="Blues", shade=True, 
                       thresh=0, cbar=True)
        ax.set_xlabel('total population size $n$')
        ax.set_ylabel(r'$\hat{n}$')
        ax.set_xticks([125, 150, 175, 200, 225])
        ax.set_yticks([-75, -25,  25, 75])
        ax.grid(False)
        plt.savefig(str('PDMP_nu%d_contour.pdf' %(nu_array[i])), 
              bbox_inches='tight',pad_inches=0.05)
        plt.show()
        
        # histgram of n
        plt.hist(df.n,  range=(110, 220),bins=20, color='b', alpha=0.5)
        plt.xlabel('total population size $n$', fontsize=20)
        plt.ylabel('count', fontsize=20)
        plt.xticks([125, 150, 175, 200, 225])
        plt.grid(False)
        plt.savefig(str('PDMP_nu%d_n.pdf' %(nu_array[i])), 
              bbox_inches='tight',pad_inches=0.05)
       
        plt.show()
        
        # histgram of hat_n
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.hist(df.hat_n,  range=(-76, 76),bins=20, color='b', alpha=0.5)
        ax1.set_xlabel(r'$\hat{n}$', fontsize=20)
        ax1.set_ylabel(r'count', fontsize=20)
        ax1.set_xticks([-75, -25, 25, 75])
        ax1.grid(False)
        #ax.set(yscale='log')
        
        ax2 = ax1.twinx()
        hat=np.linspace(-74, 74, 3000)
        if nu_array[i]<0:
            nu=np.power(0.1, -nu_array[i])
            theo_hat_n=((200-125-hat)*(hat-50+125))**(nu/0.1-1)
        else:
            nu=np.power(10, nu_array[i])
            theo_hat_n=((200-125-hat)*(hat-50+125))**(80-1)
            # nu =10^1 is too large to couase overflow
        
        #log10_hat_n= (nu/0.1-1)* (np.log10(200-125-hat)+np.log10(hat-50+125))
        ax2.plot(hat, theo_hat_n/(sum(theo_hat_n)), color='r', 
                 linestyle='--', linewidth=2)
        #ax2.plot(hat, log10_hat_n-np.log10(((200-125-hat)*(hat-50+125))**(nu/0.1-1)), color='r')
        ax2.set_ylabel(r'probability '+r'$q(\hat{n})$', fontsize=20)
        ax2.set_ylim(bottom=0, top=y_top[i])
        ax2.grid(False)
        plt.savefig(str('PDMP_nu%d_hatn.pdf' %(nu_array[i])), 
              bbox_inches='tight',pad_inches=0.05)
        plt.show()
        
        """
        sns.set(rc={'axes.labelsize':20,
            'xtick.labelsize':16,
            'ytick.labelsize':16})
        ax=sns.jointplot(data=df,x='n', y='hat_n', kind='kde',
                         xlim=(100, 230),
                         ylim=(-150, 150))
        ax.set_axis_labels('total population size $n$', 
               r'$\hat{n}$')
        """
        """
        plt.hist(n, bins=40)
        plt.xlabel('total population', fontsize=20)
        plt.ylabel('Count', fontsize=20)
        plt.title(r'$ \nu = 10$'+ str('^ %d' %(nu_array[i])), fontsize=20)
        #plt.vlines(U[1], 0, ymax[i], 'k', linewidth=2, linestyle='--')
        #plt.vlines(0, 0, ymax[i], 'k', linewidth=2, linestyle='-')
        #plt.vlines(U[0], 0, ymax[i], 'k', linewidth=2, linestyle='--')
        #plt.savefig(str('PDMP_nu%d.pdf' %(nu_array[i])), 
        #            bbox_inches='tight',pad_inches=0.0)
        plt.show()
        plt.hist(hat_n, bins=40)
        plt.xlabel(r'$\hat{n}$', fontsize=20)
        plt.ylabel('Count', fontsize=20)
        plt.title(r'$ \nu = 10$'+ str('^ %d' %(nu_array[i])), fontsize=20)
        
        ax.ax_joint.text(110, 125, r"$\nu=10$"+str("^%d" %(nu_array[i])), fontsize=20)
        plt.savefig(str('PDMP_nu%d.pdf' %(nu_array[i])), 
              bbox_inches='tight',pad_inches=0.05)
        plt.show()
        """
        

def ODE(t, x, R, alpha, delta):
    dxdt=np.zeros([4])
    dxdt[0]= alpha*(R-x[0]) - (1.0*x[2]+0.91*x[3])*x[0]/(100+x[0])
    dxdt[1] = alpha*(125-x[1]) -delta*x[1]*(x[2]+x[3])/(100+x[1])
    dxdt[2]= (1.0*x[0]/(100+x[0]) - delta*x[1]/(100+x[1])-alpha)*x[2]
    dxdt[3]= (0.91*x[0]/(100+x[0]) - delta*x[1]/(100+x[1])-alpha)*x[3]
    return dxdt
