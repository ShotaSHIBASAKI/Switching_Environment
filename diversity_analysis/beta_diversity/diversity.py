#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 14:25:09 2020
Convert abundance data into frequency
Then calculate alpha, beta, gamma diversity
follwoing Jost 2007 and Chao 2012
Corresponding figures are Figs.4B, D, etc.
"""

import numpy as np
from plotnine import *
import matplotlib.pyplot as plt
import os




def ConvertFreq(data):
    """
        data should be abundance data
        each row: each simulation
        each column: each species
    """
    data_conv=data
    
    for i in range(np.size(data,0)):
        if sum(data[i,:])==0:
            data_conv[i,:]=np.zeros(np.size(data,1))
        else:
            data_conv[i, :] =data[i,:]/sum(data[i,:])
    return data_conv

def Alpha(freq_matrix, q, weight):
    """
    calculate alpha diversity using frequency data

    Parameters
    ----------
    freq_matrix : 2dim array for species frequency
    row: simulation
    col: species
    q : effectvive number of species
    wieght: weight vector for each simulation

    Returns: alpha diversity
    ------

    """
    if q==1:
        # use approx to q ->1
        alpha=np.exp(np.dot(weight,np.sum(InfoGain(freq_matrix),1)))
    else:
        if np.sum(freq_matrix)==0:
            alpha=0
        else:
            alpha=(np.dot(weight**q, np.sum((freq_matrix**q),1))/np.sum(weight**q))**(1/(1-q))
        
    return alpha

def Gamma(freq_matrix, q, weight):
    """
    calculate gamma diversity using frequency data

    Parameters
    ----------
    freq_matrix : 2dim array for species frequency
    row: simulation
    col: species
    q : effectvive number of species
    wieght: weight vector for each simulation
    Returns: gamma diversity
    ------
    """
    freq_mean=np.dot(weight, freq_matrix)
    if q==1:
        # use approximation of q -> 1
        gamma=np.exp(np.sum(InfoGain(freq_mean)))
    else:
        gamma=(np.sum(freq_mean**q))**(1/(1-q))
    return gamma
        





def InfoGain(p):
    #if p <=0, the infromation gain should be 0
    return np.where(p<=0, 0, -p*np.log(p))
        
def main():
    q=1 #  parameter for calculating diversity
    parm_set=100  # number of parameter sets
    iteration=100 # number of iteration using a parameter set
    delta_array=np.array([0.1,0.2, 0.4, 0.6, 1.0])
    species_array=np.array([2,4,6,8,10])
    nu_exp=np.linspace(-5,3,9)
    
    for i in range(np.size(delta_array)):
        delta=delta_array[i]
        #path=str('./death%.1f' %(delta))
        #os.chdir(path)
        for j in range(np.size(species_array)):
            n=species_array[j]
            #path=str('species_number%d' %(n))
            
            #initialize diversity matrix
            alpha_div=np.zeros([iteration, np.size(nu_exp)]) #alpha diversity
            beta_div=np.zeros([iteration, np.size(nu_exp)]) #alpha diversity
            gamma_div=np.zeros([iteration, np.size(nu_exp)]) #alpha diversity
            #os.chdir(path)
            
            #use this codes when you plot the results of main_diversity.c
            """
            for k in range(parm_set):
                #os. chdir(str('./parameter%d' %(k+1)))
                for l in range(np.size(nu_exp)):
                    
                    #convert aboundance data into csv data
                    fname=str("Dynamics_nu%d_last.csv" %(nu_exp[l]))
                    data=np.loadtxt(fname, delimiter=',', skiprows=1)
                    weight=np.sum(data,1)
                    if np.sum(weight)>0:
                         weight=weight/np.sum(weight) # normalize
                    data_freq=ConvertFreq(data[:,  n+1:2*n+1])
                    sname=str("Dynamics_frequency_nu%d_last.csv" %(nu_exp[l]))
                    np.savetxt(sname, data_freq, delimiter=',')
                    
                    #calculate alpha, beta, gamma diversity
                    alpha_div[k,l]=Alpha(data_freq,q,weight)
                    gamma_div[k,l]=Gamma(data_freq,q, weight)
                    beta_div[k,l]=gamma_div[k,l]/alpha_div[k,l]
                
                #os.chdir('../')  # end analyzing parameter k data
            """
            #otherwise plot the beta diversity using beta_div_speciesN_death_delta.csv files
            #Note that you cannot plot either alpha or gamma diversities in this case
            fname=str("beta_div_species%d_death%d.csv" %(n, 10*delta))
            beta_div=np.loadtxt(fname, delimiter=',', skiprows=1)
            
                    
               
            """
            #plot distribution of diversities
            for l in range(np.size(nu_exp)):
                # alpha diversity
                plt.hist(alpha_div[:,l])
                plt.xlabel(str("alpha diversituy with q=%d" %(q)), fontsize=18)
                plt.ylabel('frequency', fontsize=18)
                gname_alpha=str("AlphaDiv_weighted_q%d_species%d_death%d_nu%d.pdf" 
                                %(q,n, 10*delta, nu_exp[l]))
                plt.savefig(gname_alpha, bbox_inches="tight", pad_inches=0.05)
                plt.show()
                 # beta diversity
                plt.hist(beta_div[:,l])
                plt.xlabel(str("beta diversituy with q=%d" %(q)), fontsize=18)
                plt.ylabel('frequency', fontsize=18)
                gname_beta=str("BetaDiv_weighted_q%d_species%d_death%d_nu%d.pdf" 
                               %(q, n, 10*delta, nu_exp[l]))
                plt.savefig(gname_beta, bbox_inches="tight", pad_inches=0.05)
                plt.show()
                 # gama diversity
                plt.hist(gamma_div[:,l])
                plt.xlabel(str("gamma diversituy with q=%d" %(q)), fontsize=18)
                plt.ylabel('frequency', fontsize=18)
                gname_gamma=str("GammaDiv_weighted_q%d_species%d_death%d_nu%d.pdf" 
                                %(q,n, 10*delta, nu_exp[l]))
                plt.savefig(gname_gamma, bbox_inches="tight", pad_inches=0.05)
                plt.show()    
            """
            """
            # violin plot over swiching rate
            fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, 
                                                figsize=(9, 4), sharey=False)
            #alpha diversity
            ax1.set_title('alpha diversity')
            ax1.violinplot(alpha_div)
            alpha_mean=np.mean(alpha_div, 0)
            ax1.plot(np.linspace(1,9,9), alpha_mean, marker='D', color='k')
            ax1.set_xticks([2,4,6,8])
            ax1.set_xticklabels(['$10^{-4}$', '$10^{-2}$', '$10^0$', '$10^2$'])
            
            #beta diversity
            ax2.set_title('beta diversity')
            ax2.set_xlabel('switching rate', fontsize=18)
            ax2.violinplot(beta_div)
            beta_mean=np.mean(beta_div, 0)
            ax2.plot(np.linspace(1,9,9), beta_mean, marker='D', color='k')
            ax2.set_xticks([2,4,6,8])
            ax2.set_xticklabels(['$10^{-4}$', '$10^{-2}$', '$10^0$', '$10^2$'])
            
            #gamma diversity
            ax3.set_title('gamma diversity')
            ax3.violinplot(gamma_div)
            gamma_mean=np.mean(gamma_div, 0)
            ax3.plot(np.linspace(1,9,9), gamma_mean, marker='D', color='k')
            ax3.set_xticks([2,4,6,8])
            ax3.set_xticklabels(['$10^{-4}$', '$10^{-2}$', '$10^0$', '$10^2$'])
            
            plt.savefig(str("ABGdiversities_q%d_species%d_death%d.pdf" %(q, n, 10*delta)),
                        bbox_inches="tight", pad_inches=0.05)
            plt.show()
            """
            # showing only beta div
            plt.violinplot(beta_div)
            beta_mean=np.mean(beta_div, 0)
            plt.plot(np.linspace(1,9,9), beta_mean, marker='D', color='k')
            plt.xticks([2,4,6,8],['$-4$', '$-2$', '$0$', '$2$'],
                       fontsize=24)
            plt.xlabel('$\log_{10}$'+' switching rate', fontsize=28)  
            plt.yticks(fontsize=24)
            #plt.ylabel('beta diversity', fontsize=18)
            plt.savefig(str("Betadiversities_q%d_species%d_death%d.pdf" %(q, n, 10*delta)),
                        bbox_inches="tight", pad_inches=0.05)
            plt.show()
            
            """
            #save csv for beta diversity
            header='switching rate: 10^-5, -4, -3, -2, -1, 0, 1, 2, 3'
            cname=str("beta_div_species%d_death%d.csv" %(n, 10*delta))
            np.savetxt(cname, beta_div, delimiter=',', header=header, fmt='%.6f')
            """
            #os.chdir("../")
        #os.chdir("../")


if __name__ == "__main__":
    main()
