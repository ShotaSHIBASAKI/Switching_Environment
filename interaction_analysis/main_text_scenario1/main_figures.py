#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 11:54:38 2020
Plot Figs.2 and 3 in  the main text
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import stats
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




#--------------------------------------------------------------
#Fig.3A: different in extinction probabilities from mono- to co-cultures
#Fig.3B: some exaples  from Fig2A
def Fig3AB (model, val=0):
    # model: scenario of environmental switching:  1,2,or 3
    # val; which species 2 to comp@ete with: in the main text, choose 0. 
    # You can also chose val =1,..., 4 for more sloer grower
    typ=1
    """
    plotting the probabilities of
    diff in extinction
    """
    d_array=np.linspace(0.1, 1.0, 10) # toxin sensitivity
    nu_array=np.linspace(-5, 3, 9)  # log scale in 3d is problematic in matplot
    # extinction prob at nu=10^-5(->0),10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3(->infty) 
    one_ext=np.zeros([np.size(d_array), np.size(nu_array)])  # extinction prob in absence of sepcies 2
    two_ext=np.zeros([np.size(d_array), np.size(nu_array)])  # extinction prob in presence of species 2
    diff=np.zeros([np.size(d_array), np.size(nu_array)])
    for i in range(np.size(d_array)):
        death=d_array[i]
        # In absence of species 2
        path=str('./death%.1f/OneConsumer' %(death))
        os.chdir(path)
        data=np.loadtxt(str('OneConsumer_Extinction_model%d.csv'%(model)),
                        delimiter=',', skiprows=1)
        one_ext[i, :]=data
        os.chdir('../TwoConsumer')
        
        # Extinction prob in presence of species 2
        data=np.loadtxt(str('TwoConsumer_Extinction_model%d_competitor%d.csv'%(model, typ)),
                        delimiter=',', skiprows=1)
        
        two_ext[i, :]=data[val,:] 
        diff[i, :]=one_ext[i,:]- two_ext[i, :]
        os.chdir('../../')
    #plot figures
    #plot Fig.2A
    ax1 = plt.subplot(111)
    heatmap1=ax1.pcolor(diff, cmap='Blues_r')
    #heatmap=ax.pcolor(comp_array, cmap=plt.cm.bwr,vmin=0)
    ax1.xaxis.tick_bottom()
    plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
    plt.ylabel('toxin sensitivity', fontsize=20)
    plt.title('diff. in extinction prob', fontsize=20)
    plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
    plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
    cb1=plt.colorbar(heatmap1, ticks=[-0.15, -0.1, -0.05, 0])
    cb1.ax.tick_params(labelsize=14)
    plt.savefig('DiffExtinction_Heatmap_'+str('_Model%d.pdf'%(model)),
              bbox_inches='tight', pad_inches=0.05)
    plt.show()
    
    #plot Fig.3B
    
    rate=np.linspace(-5, 3, 9)
    col_list=['#8dd3c7', '#fb8072', '#bebada', '#80b1d3']
    lab_list=['sensitivity 0.1', 'sensitivity 0.2', 
              'sensitivity 0.4', 'sensitivity 0.6']
    example=[0,1,3,5]
    for i in range(len(example)):
        plt.plot(rate, diff[example[i], :], color=col_list[i], label=lab_list[i],
                marker='D', linewidth=2)
        #Calculate 95% HDP for each dot using Dirichlet (or beta) distribution
        #As we have 10**5 data, 95% HDP are unobservable in the figure
        for j in range(np.size(rate)):
            ext1=one_ext[example[i],j] # probs of extinction and persistence 
            obs1=np.array([ext1, 1-ext1])*10**5 # total observation is 10**5 times
            ext2=two_ext[example[i],j]
            obs2=np.array([ext2,1-ext2])*10**5
            #estimate 95%HDI 
            sample1=dirichlet.rvs(obs1+1, size=10000)
            sample2=dirichlet.rvs(obs2+1, size=10000)
            sample_diff=sample1[:,0]-sample2[:, 0]
            hdp_l, hdp_h=pm.stats.hpd(sample_diff)
            plt.vlines(rate[j],hdp_l, hdp_h, color=col_list[i])
    plt.xlabel('$\log_{10}$'+'switching rate',fontsize=20)
    plt.ylabel('diff. in extinction prob',fontsize=20)
    plt.xticks(fontsize=16, ticks=[-4, -2, 0, 2])
    plt.yticks(fontsize=16, ticks=[0, -0.05, -0.1, -0.15, -0.2])
    plt.legend(loc='lower center', bbox_to_anchor=(.5,1.05),
               ncol=2, fontsize=16)
    plt.savefig('DiffExtinction_examples.pdf',bbox_inches='tight',pad_inches=0.05)
    plt.show()

#--------------------------------------------------------------

#Fig.3C: probabilities of  exclusion of species
#Fig.3D: some exaples from Fig.2C
def Fig3CD (model, val=0):
    # model: scenario of environmental switching:  1,2,or 3
    # val; which species 2 to comp@ete with: in the main text, choose 0. 
    # You can also chose val =1,..., 4 for more sloer grower
    typ=1
    """
    plotting the probabilities of
    diff in extinction
    """
    d_array=np.linspace(0.1, 1.0, 10) # toxin sensitivity
    nu_array=np.linspace(-5, 3, 9)  # log scale in 3d is problematic in matplot
    # extinction prob at nu=10^-5(->0),10^-4, 10^-3, 10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3(->infty) 
    comp_exc=np.zeros([np.size(d_array),np.size(nu_array)]) # prob of competitive exclusion
    for i in range(np.size(d_array)):
        death=d_array[i]
        # In absence of species 2
        path=str('./death%.1f/TwoConsumer' %(death))
        os.chdir(path)
        #Competitive exclusion
        data=np.loadtxt(str('TwoConsumer_CompetitiveExclusion_model%d_competitor%d.csv'%(model, typ)),
                        delimiter=',', skiprows=1)
        comp_exc[i, :]=data[val,:] 
        os.chdir('../../')
    #plot figures
    
    #plot Fig.3C
    ax1 = plt.subplot(111)
    heatmap1=ax1.pcolor(comp_exc, cmap='Blues',vmin=0)
    #heatmap=ax.pcolor(comp_array, cmap=plt.cm.bwr,vmin=0)
    ax1.xaxis.tick_bottom()
    plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
    plt.ylabel('toxin sensitivity', fontsize=20)
    plt.title('prob. of exclusion of fittest', fontsize=20)
    plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
    plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
    cb1=plt.colorbar(heatmap1, ticks=[0,0.05, 0.1, 0.15])
    cb1.ax.tick_params(labelsize=14)
    plt.savefig('CompetitiveExclusion_Heatmap_'+str('_Model%d.pdf'%(model)),
                bbox_inches='tight', pad_inches=0.05)
    plt.show()
    
    #plot Fig.3D
    toxin=np.linspace(0.1, 1.0, 10)
    Data=np.zeros([3, 10])
    Data[0, :]= comp_exc[:, 0] # nu=10^-5
    Data[1, :]= comp_exc[:, 4] # nu =0.1
    Data[2, :]= comp_exc[:, -1] # nu=10^3
    #plot the data
    #col_list=['#8dd3c7', '#fb8072', '#80b1d3']
    col_list=['grey', 'slategrey', 'black']
    line_list=['dotted', 'dashed', 'solid']
    lab_list=['slow', 'intermediate', 'fast']
    for i in range(np.size(Data, 0)):
        plt.plot(toxin, Data[i, :], color=col_list[i], label=lab_list[i],
                 linestyle=line_list[i], marker='D', linewidth=2)
        """
        Plot 95% HDP. Again, HDPs are too samll to see in the plot
        For simplicity, we use beta distribution to see the probability 
        that competitive exclusion of sp1 by sp2 happens.
        We may also use Dirichlet distribution instead by considering 
         1. both species extinction
         2. both species persistnce
         3. competitive exclusion of sp1 by sp2 and
         4. competitive exclusion of sp2 by sp1
        """
        for j in range(np.size(toxin)):
            comp=Data[i,j]
            a=np.array([comp, 1-comp])*10**5 # total number of simulations is 10**5
            sample=dirichlet.rvs(a+1, size=10000)
            hdp_l, hdp_h=pm.stats.hpd(sample[:,0])
            plt.vlines(toxin[j],hdp_l, hdp_h, color=col_list[i])
    plt.xlabel('toxin sensitivity',fontsize=20)
    plt.ylabel('prob. of exclusion of fittest',fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='best', fontsize=16)
    plt.savefig('CompetitiveExclusion_furthest.pdf',bbox_inches='tight',pad_inches=0.05)
    plt.show()
#-------------------------------------------
def Fig4():
    # plot the sample dynamics where exclusion of fittest happens
    fname_list=['./TimeSeries_switching10^-3_delta0.1_trial2.csv', './TimeSeries_switching10^-1_delta0.2_trial0.csv', 
                './TimeSeries_switching10^1_delta0.4_trial18.csv']
    for i in range(len(fname_list)):
        df=pd.read_csv(fname_list[i])
        df=df.rename(columns={' resource': 'resource',' consumer1': 'species1',
                 ' consumer2': 'species2',' toxin': 'toxin',
                 ' environment': 'environment'})
        df_good=df[df.environment==1]
        df_bad=df[df.environment==-1]
        plt.plot(df.time, df.resource, label='resource', color='#fdb462')
        plt.plot(df.time, df.toxin, label='toxin', color='#80b1d3')
        plt.plot(df.time, df.species1, label='species 1', color='#8dd3c7')
        plt.plot(df.time, df.species2, label='species 2',color='#b3de69')

        plt.legend(bbox_to_anchor=(0.5, 1.15), loc='center', ncol=2, fontsize=18)
        plt.scatter(df_bad.time, 0*df_bad.time, color='k', marker='s', s=30)
        plt.scatter(df_good.time, 0*df_good.time, color='w', marker='s', s=30)
        plt.xlabel('time', fontsize=20)
        plt.xticks(fontsize=16)
        plt.ylabel('amount or abundances', fontsize=20)
        plt.xlim(0, 205)
        plt.ylim(0, 160)
        plt.yticks(fontsize=16)
        plt.text( x=25,y=145, s=str("nu = 10^%d, delta=%.1f" %(-1, 0.2)),fontsize=20)
        gname=str('Ex_exclusion%d.pdf' %(i))
        plt.savefig(gname, bbox_inches='tight',pad_inches=0.05)
        plt.show()
    


#--------------------------------
#Fig.5A-C: both extinction and exclusion of fittest without EFs
def Fig5AC():
   os.chdir('../appendix/Appendix3_constant_env/scenario1')
   CompExcl=np.loadtxt('CompetitiveExclusion_model1.csv',
                        delimiter=',', skiprows=1)
   BothExtinct=np.loadtxt('BothExtinction_model1.csv',
                        delimiter=',', skiprows=1)
   death=np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
   col_list=['#8dd3c7', '#fb8072', '#80b1d3']
   lab_list=['scarce', 'mean', 'abundant']
   width=0.05
   for i in range(len(lab_list)):
  
       plt.bar(death, CompExcl[:, i], color=col_list[i], width=width, align='center',
               hatch='..', label='exclusion of fittest')
       plt.bar(death, BothExtinct[:, i], color=col_list[i], label='both extinction',
               width=width, align='center', bottom=CompExcl[:, i])
   
       plt.xticks(fontsize=16)
       plt.ylabel('probability',fontsize=20)
       plt.yticks(fontsize=16)
       plt.ylim(0, 1)
       plt.xlabel('toxin sensitivity',fontsize=20)
       plt.legend(loc='lower center', bbox_to_anchor=(.5,1.05),
               ncol=2, fontsize=16)
   
       plt.plot(death, CompExcl[:, i], color='k', 
         label=lab_list[i], linewidth=3)
       plt.text(x=0.05, y=0.9, s=lab_list[i], fontsize=20)
    

       for j in range(np.size(death)):
           comp=CompExcl[j,i]            
           a=np.array([comp, 1-comp])*10**5 # total number of simulations is 10**5
           sample=dirichlet.rvs(a+1, size=10000)
           hdp_l, hdp_h=pm.stats.hpd(sample[:,0])
           plt.vlines(death[j],hdp_l, hdp_h, color=col_list[i])
       plt.savefig('ConstantEnv_ver2' + lab_list[i] + '.pdf' 
        ,bbox_inches='tight',pad_inches=0.05)
       plt.show()


    
#------------------------------------------------    

#Fig.4D: Correlation of non-monotonicity and distances between critical toxin sensitivities
def NonMonotonicity(data):
    sign_prev=0
    for i in range(np.size(data)-1):
        
        if data[i+1]-data[i]<-0.01:
            sign_curr=-1 # increasing
        elif data[i+1]-data[i]>0.01:
            sign_curr=1 #decreasing
        else:
            sign_curr=sign_prev
        if sign_prev*sign_curr<-0.1:
            # non-monotonicity
            return 1
        else:
            sign_prev=sign_curr
    #monotnic
    return 0

def Number_NonMonotonic(start, end, Data):
    counter=0
    for j in range(int(start)-1, int(end)):
        counter+=NonMonotonicity(Data[j, :])
    return counter
        
    
def FigA14(peak_harsh, peak_mean, peak_mild):
    """
    This figure was preivously shown in the main text, but in the latest, they are showin in Appendix.
    Need the following three one-dim array that contains critical toxin 
    sensitivities in each scenario without enviornmental scenario
    ----------
    peak_harsh 
    peak_mean 
    peak_mild 
    ----------
    In the main text, the above arrays are
    np.array([0.1, 0.3, 0.1, 0.1, 0.3, 0.1, 0.1]),
    np.array([0.4, 0.4, 0.4, 0.9, 0.5, 0.2, 0.3]),
    np.array([0.8, 1.1, 1.3, 1.2, 0.8, 0.3, 0.8])
    
    Plot the distance of peaks and the range of death rate with non-monotonic effect of environmental siwthcing rate
    Ex. three switching scenarios + four chaging resource supply
    """
    Non_Monotonic=np.zeros([3, 7])  
    # 1st dim: mean-harsh, mild-mean, mild-harsh
    #2nd dim: each scenario
    
    #calculate non-monotonicity in each switchingscenario
    model_array=np.array([1,2,3])

    os.chdir('./Correlation')
    for i in range(np.size(model_array)):
        # non-monotonicity check
        fname=str('DiffExtinction1_Heatmap_Model%d.csv'%(model_array[i]))
        Data=np.loadtxt(fname, delimiter=',', skiprows=0)
        Non_Monotonic[0, i]=Number_NonMonotonic(peak_harsh[i]*10, peak_mean[i]*10, Data)
        Non_Monotonic[1, i]=Number_NonMonotonic(peak_mean[i]*10, min(peak_mild[i]*10,10), Data)
        Non_Monotonic[2, i]=Number_NonMonotonic(peak_harsh[i]*10, min(peak_mild[i]*10, 10), Data)

    
   
    tail_list=['large_xmax', 'large_xmin', 'small_xmax', 'small_xmin']
    for j in range(np.size(tail_list)):
        fname1='DiffExtinction1_Heatmap_Model1_'+tail_list[j]+'.csv'
        Data=np.loadtxt(fname1, delimiter=',', skiprows=0)
        Non_Monotonic[0, j+3]=Number_NonMonotonic(peak_harsh[j+3]*10, peak_mean[j+3]*10, Data)
        Non_Monotonic[1, j+3]=Number_NonMonotonic(peak_mean[j+3]*10, peak_mild[j+3]*10, Data)
        Non_Monotonic[2, j+3]=Number_NonMonotonic(peak_harsh[j+3]*10, peak_mild[j+3]*10, Data)
       
    
    # scatter plot and correlation
    corr, p_val=sp.stats.spearmanr(peak_mean-peak_harsh, Non_Monotonic[0, :])
    plt.scatter(peak_mean-peak_harsh, Non_Monotonic[0, :], s=200,
               marker='o', color='black', label='harsh - mean')
    print(corr, p_val)
    corr, p_val=sp.stats.spearmanr(peak_mild-peak_mean, Non_Monotonic[1, :])
    plt.scatter(peak_mild-peak_mean, Non_Monotonic[1, :], s=100,
               marker='D', color='dimgray', label='mean - mild')
    print(corr, p_val)
    corr, p_val=sp.stats.spearmanr(peak_mild-peak_harsh, Non_Monotonic[2, :])
    plt.scatter(peak_mild-peak_harsh, Non_Monotonic[2, :], s=100,
               marker='x', color='grey', label='harsh - mild')
    print(corr, p_val)
    # add arrow to show the result in the main text
    plt.annotate(text="", xy=[0.35, 2.5], xytext=[0.45, 3.5], 
                 arrowprops=dict(width=2, facecolor='k', 
                                edgecolor='k'))
    plt.xlabel('distance between critical values', fontsize=20)
    plt.ylabel('non-monotonicity', fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='lower center', bbox_to_anchor=(.5, 1.1), ncol=3, fontsize=20)  
    plt.savefig('DistancePeak_Non_monotonicity_ver2.pdf',bbox_inches='tight',pad_inches=0.05)
    plt.show()
    
    print(peak_mean-peak_harsh, Non_Monotonic[0, :])
    print(peak_mild-peak_mean, Non_Monotonic[1, :])
    print(peak_mild-peak_harsh, Non_Monotonic[2, :])
