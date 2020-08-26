#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 10:13:22 2020
Plot Richness distribution over the switching rate (Figs.4C, E etc)
"""


import numpy as np
import matplotlib.pyplot as plt
import os

def Richness_count(data_sp, S):
    #counting a number of species survived
    counter=np.zeros(S+1) # richness counter
    for i in range(np.size(data_sp,0)):
        count=0
        for j in range(S):
            if data_sp[i,j]>0:
                count+=1 # richness
        counter[count]+=1
    return counter
        


def main(S):
    # S number of species at the beginnng of each simulatuon
    parm_set=100 # number of meta-community given mean toxin sensitivity and switching rate
    delta_array=np.array([0.1,0.2, 0.4, 0.6, 1.0]) # mean toxin sensitivities
    
    nu=np.linspace(-5, 3, 9)
    for i in range(np.size(delta_array)):
        delta=delta_array[i]
        Richness=np.zeros((S+1, np.size(nu))) # richness can be 0,1, ...S
        """
        #use the following codes to calculate species richness from the results of main_diversity.c
        path=str("./death%.1f/Species_number%d" %(delta, S))
        os.chdir(path)
        for j in range(parm_set):
            path=str("./parameter%d" %(j+1))
            os.chdir(path)
            for k in range(np.size(nu)):
                fname=str("Dynamics_nu%d_last.csv" %(nu[k]))
                data=np.loadtxt(fname, delimiter=',', skiprows=1)
                data_sp=data[:,  S+1:2*S+1]  # species data in each community given parameter set
                Richness[:, k] += Richness_count(data_sp, S)
            os.chdir("../")
        """    
        #otherwise use the following
        fname=str("Prop_richness_species%d_death%d.csv" %(S, 10*delta))
        Richness.T= np.loadtxt(fname, delimiter=',',skiprows=1)
        
        
        
        #plot data at mean toxin seinsitivity = delta
        base=np.zeros(np.size(nu))
        colo_list=['#8dd3c7','#ffffb3', '#bebada', '#fb8072', '#80b1d3', 
                   '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd','#ccebc5' ]
        #print(Richness)
        for i in range(S+1):
            plt.bar(nu, Richness[i, :], bottom=base, color=colo_list[i],
                        label=str('richness %d' %(i)))
            base+=Richness[i, :]
        plt.xticks([-4,-2,0,2],['$-4$', '$-2$', '$0$', '$2$'],
                       fontsize=24)
        plt.xlabel('$\log_{10}$'+' switching rate', fontsize=28)  
        plt.yticks(fontsize=24)
        #plt.ylabel("number of case", fontsize=18)
        #plt.legend(loc='lower center', bbox_to_anchor=(.5, 1.1), ncol=3, fontsize=18)
        tname=str("Prop_richness_species%d_death%d.pdf" %(S, 10*delta))
        plt.savefig(tname,bbox_inches="tight", pad_inches=0.05)
        plt.show()
        
        #save csv file
        #proportion of species richness x switching rate nu
        h=np.arange(S+1)
        header=','.join(str(x) for x in h) #headerfor richness
        cname=str("Prop_richness_species%d_death%d.csv" %(S, 10*delta))
        np.savetxt(cname, Richness.T, delimiter=',', header='richness'+header, fmt='%d')
        #os.chdir("../../")
        