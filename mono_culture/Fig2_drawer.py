import os
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import stats
import pymc3 as pm
from scipy.stats import dirichlet
import pandas as pd
import seaborn as sns

import pymc3 as pm
from scipy.stats import dirichlet

def SmaplingDIr(obs, index):
    """
    Estimating 95% Heighest posterior deinsty (HDP) in Dirichlet distribution.
    Although observed data follows a multinomial-distribution,
    we can also use the cases when the data are binary (in such cases,DIrichlet is Beta dist)
    Parameters
    ----------
    obs : 1D array
        number of obseraving each event.
    index : int
        event to analyze its occuring probability.
    Returns
    -------
    95% HPD of occuring event index
    """
    a=np.ones(np.size(obs)) # prior is uniform
    a+=obs  # posterior after observing data
    #posterior = dirichlet(a)
    sample=dirichlet.rvs(a, size=10000)#sampleing parameter values from posterior distirbution
    sample_index=sample[:, index] # focus on parameter[index] alone
    """
    plt.hist(sample_index)
    plt.xlim(0,1)
    plt.show()
    """
    # calculate 95% HDI
    return pm.stats.hpd(sample_index)

def main(model, val=0):
    # model: scenario of environmental switching:  1,2,or 3
    # val; which species 2 to comp@ete with: in the main text, choose 0. 
    # You can also chose val =1,..., 4 for more sloer grower
    """
    plotting the probabilities of
    diff in extinction
    """
    d_array=np.linspace(0.1, 1.0, 10) # toxin sensitivity
    nu_array=np.linspace(-5, 3, 9)  # log scale in 3d is problematic in matplot
    # extinction prob at nu=10^-5(->0),10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3(->infty) 
    mono_ext=np.zeros([np.size(d_array), np.size(nu_array)])  # extinction of both species
    
    for i in range(np.size(d_array)):
        death=d_array[i]
        # In absence of species 2
        data=np.loadtxt(str('OneConsumer_Extinction_model%d_death%d.csv'%(model, 10*death)),
                        delimiter=',', skiprows=1)
        mono_ext[i, :]=data
       
    #plot Fig.1A: extinction probability
    rate=np.linspace(-5, 3, 9)
    col_list=['#8dd3c7', '#fb8072', '#80b1d3']
    
    example=[1, 4, 9] # in thes cases, distributions of R(sigma_end) is also available
    for i in range(len(example)):
        plt.plot(rate, mono_ext[example[i], :],  
                 label=str('sensitivity %.1f' %(d_array[example[i]])),
                 marker='D', linewidth=2, color=col_list[i])
        for j in range(np.size(rate)):
            ext1=mono_ext[example[i],j] # probs of extinction and persistence 
            obs1=np.array([ext1, 1-ext1])*10**5 # total observation is 10**5 times
            #estimate 95%HDI 
            hdp_l, hdp_h=SmaplingDIr(obs1, 0)
            #print([hdp_l, hdp_h])
            plt.vlines(rate[j],hdp_l, hdp_h, color=col_list[i])

    plt.xlabel('$\log_{10}$'+'switching rate',fontsize=20)
    plt.ylabel('prob. of extinction',fontsize=20)
    plt.xticks(fontsize=16, ticks=[-4, -2, 0, 2])
    #plt.yticks(fontsize=16, ticks=[0, -0.05, -0.1, -0.15, -0.2])
    plt.legend(loc='lower center', bbox_to_anchor=(.5,1.05),
               ncol=2, fontsize=16)
    plt.savefig('MonoExtinction.pdf',bbox_inches='tight',pad_inches=0.05)
    plt.show()
    
    # plotting abundance distributions (Fig 1B--D)
    for i in range(len(example)):
        path=str('./death%df' %(10*d_array[example[i]]))
        os.chdir(path)
        df=pd.DataFrame({})
        for j in range(np.size(nu_array)):
            fname=str('OneConsumerRTC_switching10^%d_last.csv' %(nu_array[j]))
            data=np.loadtxt(fname, delimiter=',', skiprows=1, dtype=int)[:,1]
            df[str('%d' %(nu_array[j]))]=data
        os.chdir('../')
        
        # plot
        ax=sns.violinplot(data=df, color=col_list[i])
        #ax=sns.stripplot(data=df, color='.2')
        ax.set_xlabel('$\log_{10}$'+'switching rate',fontsize=20)
        ax.set_ylabel('species 1 abundance',fontsize=20)
        plt.ylim(-1, 200)
        plt.xticks(fontsize=16, ticks=[1,3,5,7])
        tname=str('sensitivity %.1f' %(d_array[example[i]]))
        plt.text(x=5, y=175,s=tname, fontsize=16)
        plt.savefig(str('abundance_sensitivity%d.pdf' %(10*d_array[example[i]])),
                 bbox_inches='tight',pad_inches=0.05)
        plt.show()
