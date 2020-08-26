#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 12:16:41 2020
Camcluating competitive exclusion
for two species scenario. See Fig. 4A
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp #ode integration

def CompProb(s1, s2, nu):
    #s1 is expected to persist in absence of noise
    #calculate prob s2 excludes s1
    prob=np.zeros(np.size(nu))
    for i in range(np.size(nu)):
        switch=nu[i]
        data=np.loadtxt(str("Dynamics_frequency_nu%d_last.csv" %(switch)), 
                        delimiter=',')
        counter=0
        for j in range (np.size(data, 0)):
            if data[j,s1] <pow(10, -6) and data[j,s2] >=pow(10, -6):
                counter+=1
        prob[i]=counter/(np.size(data, 0))
    return prob
                    
def MeanFieldODE(k,l, beta_matrix, K_matrix):
    """
    solving ode 
    in two species 
    in absenc of noise
    """
    Tf=500
    beta_k=beta_matrix[k, :]
    K_k=K_matrix[k, :]
    beta_l=beta_matrix[l, :]
    K_l=K_matrix[l, :]
     
     
    init_c=np.ones(np.size(beta_k))*150
    init_s=np.ones(2)*10
    init=np.concatenate((init_c, init_s))
    sol = solve_ivp(ODE, [0, Tf], y0=init, method='LSODA',
                    args=[beta_k, K_k, beta_l, K_l]) # solving ode
    # determine which species should be species1 (i.e., superior in mean fieald)?
    if sol.y[-2, -1]>pow(10, -3) and sol.y[-1, 0]<pow(10,-3):
        # species k outcompetes l
        return [k,l]
    elif sol.y[-2, -1]<pow(10, -3) and sol.y[-1, 0]>pow(10, -3):
        # species l outcompetes k
        return [l, k]
    elif sol.y[-2, -1]>pow(10, -3) and sol.y[-1, 0]>pow(10, -3):
        #coexistence  depending on population size
        if sol.y[-2, -1]>sol.y[-1,-1]:
            return [k,l]
        else:
            return [l,k]
    
    else:
        # both extinction:    randomly chosen
        if np.random.rand()<0.5:
            return [k,l]
        else:
             return [l, k]
            
def ODE (t, x, beta_k, K_k, beta_l, K_l):
    alpha=0.1  # dilution rate
    Cin=125   # inflow
    num_compound=np.size(beta_k)  
    dxdt=np.zeros(num_compound+2)
    for i in range(num_compound):
        dxdt[i]=alpha*(Cin-x[i]) -(abs(beta_k[i])*x[i]/(x[i]+K_k[i])*x[-2]) -\
        (abs(beta_l[i])*x[i]/(x[i]+K_l[i])*x[-1])
        dxdt[-2]+=beta_k[i]*x[i]/(x[i]+K_k[i])
        dxdt[-1]+=beta_l[i]*x[i]/(x[i]+K_l[i])
    dxdt[-2]=(dxdt[-2]-alpha)*x[-2]
    dxdt[-1]=(dxdt[-1]-alpha)*x[-1]
    return dxdt
    
 
def main():
    """
    Calculating probability that the inferior species (species 2) in the mean feald 
    excludes the superior species in the stochastic simulations
    """
    Species=2
    death_array=np.array([0.1, 0.2, 0.4, 0.6, 1.0])#mean toxin sensitivity
    iteration=100
    nu=np.linspace(-5, 3, 9)
    for i in range(np.size(death_array)):
        comp_excl=np.zeros([0,9])
        death=death_array[i]
        path=str("./death%.1f/Species_number%d" %(death, Species))
        os.chdir(path)
        for j in range(1, iteration+1):
            print(str("iteration %d" %(j)))
            os.chdir(str("./parameter%d"%(j)))
            K_matrix=np.loadtxt(str("Prameter_K_set%d.csv" %(j)),
                                delimiter=',', skiprows=0)
            beta_matrix=np.loadtxt(str("Prameter_beta_set%d.csv" %(j)),
                                delimiter=',', skiprows=0)
            
            #combination of two species
            k=0 
            l=1
            # Check which species outcompetes the other in mean-filed approx.
            s1, s2 = MeanFieldODE(k, l, beta_matrix, K_matrix)
            #s1 excludes s2 in the absence of noise
            
            #calculating competitive exclusion
            prob=CompProb(s1, s2, nu)
            prob=prob.reshape([1, np.size(nu)])
            comp_excl=np.concatenate((comp_excl, prob), axis=0) 
            os.chdir("../")
        
        #Plot violin plot
        print(comp_excl)
        plt.violinplot(comp_excl)
        comp_mean=np.mean(comp_excl, 0)
        plt.plot(np.linspace(1,9,9), comp_mean, marker='D', color='k')
        plt.xticks([2,4,6,8],['$-4$', '$-2$', '$0$', '$2$'],
                       fontsize=24)
        plt.xlabel('$\log_{10}$'+' switching rate', fontsize=28)  
        plt.yticks(fontsize=24)
        fname=str("CompetitiveExclusion_violin_Species%d_death%.0f.pdf" %(Species, death*10))
        plt.savefig(fname,bbox_inches="tight", pad_inches=0.05)
        plt.show()
        
        os.chdir("../../")



if __name__ == "__main__":
    main()
