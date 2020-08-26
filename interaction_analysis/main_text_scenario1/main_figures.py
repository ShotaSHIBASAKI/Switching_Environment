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


#Fig.2A: different in extinction probabilities from mono- to co-cultures
#Fig.2B: some exaples  from Fig2A
def Fig2AB (model, val):
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
    
    #plot Fig.2B
    rate=np.linspace(-5, 3, 9)
    col_list=['#8dd3c7', '#fb8072', '#bebada', '#80b1d3']
    lab_list=['sensitivity 0.1', 'sensitivity 0.2', 
              'sensitivity 0.4', 'sensitivity 0.6']
    example=[0,1,3,5]
    for i in range(len(example)):
        plt.plot(rate, diff[example[i], :], color=col_list[i], label=lab_list[i],
                marker='D', linewidth=2)
    plt.xlabel('$\log_{10}$'+'switching rate',fontsize=20)
    plt.ylabel('prob. of competitive exclusion',fontsize=20)
    plt.xticks(fontsize=16, ticks=[-4, -2, 0, 2])
    plt.yticks(fontsize=16, ticks=[0, -0.05, -0.1, -0.15, -0.2])
    plt.legend(loc='lower center', bbox_to_anchor=(.5,1.05),
               ncol=2, fontsize=16)
    plt.savefig('CompetitiveExclusion_furthest.pdf',bbox_inches='tight',pad_inches=0.05)
    plt.show()

#--------------------------------------------------------------

#Fig.2C: probabilities of competitve exclusion
#Fig.2D: some exaples from Fig.2C
def Fig2CD (model, val):
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
    
    #plot Fig.2C
    ax1 = plt.subplot(111)
    heatmap1=ax1.pcolor(comp_exc, cmap='Blues',vmin=0)
    #heatmap=ax.pcolor(comp_array, cmap=plt.cm.bwr,vmin=0)
    ax1.xaxis.tick_bottom()
    plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
    plt.ylabel('toxin sensitivity', fontsize=20)
    plt.title('prob. of competitive exclusion', fontsize=20)
    plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
    plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
    cb1=plt.colorbar(heatmap1, ticks=[0,0.05, 0.1, 0.15])
    cb1.ax.tick_params(labelsize=14)
    plt.savefig('CompetitiveExclusion_Heatmap_'+str('_Model%d.pdf'%(model)),
                bbox_inches='tight', pad_inches=0.05)
    plt.show()
    
    #plot Fig.2D
    toxin=np.linspace(0.1, 1.0, 10)
    Data=np.zeros([3, 10])
    Data[0, :]= comp_exc[:, 0] # nu=10^-5
    Data[1, :]= comp_exc[:, 4] # nu =0.1
    Data[2, :]= comp_exc[:, -1] # nu=10^3
    #plot the data
    col_list=['#8dd3c7', '#fb8072', '#80b1d3']
    lab_list=['slow', 'intermediate', 'fast']
    for i in range(np.size(Data, 0)):
        plt.plot(toxin, Data[i, :], color=col_list[i], label=lab_list[i],
                marker='D', linewidth=2)
    plt.xlabel('toxin sensitivity',fontsize=20)
    plt.ylabel('prob. of competitive exclusion',fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='best', fontsize=16)
    plt.savefig('CompetitiveExclusion_furthest.pdf',bbox_inches='tight',pad_inches=0.05)
    plt.show()
#-------------------------------------------




#Fig.3A: competitive exclusion without environmental switching
def Fig3A():
    R_array=np.array([50, 125, 200])
    death=np.linspace(0.1, 1.0, 10)
    """
    #use the following code to calculate competitive exclusion probability
    #when the 
    CompExcl=np.zeros([np.size(R_array), np.size(death)])
    for i in range (np.size(R_array)):
        for j in range (np.size(death)):
            fname=str('TwoConsumerRTC_model1_type1_Rsupply%d_death%.1f_last.csv'
                      %(R_array[i], death[j]))#file name
            data=np.loadtxt(fname, delimiter=',', skiprows=1)
            data=data[:, 1:3] # get only species abundances at the end
            count=0
            for k in range(np.size(data,0)):
                if data[k,0]<=0 and data[k,1]>0:
                    # species 1 goers extinct but species 2 persists
                    count+=1
            count=count/np.size(data,0)
            CompExcl[i,j]=count
    """
    #otherwise use the csv file
    os.chdir('../appendix/Appendix3_constant_env/scenario1')
    CompExcl=np.loadtxt('CompetitiveExclusion_model1.csv',
                        delimiter=',', skiprows=1)
    os.chdir('../../../main_text_scenario1')
    #plot the data
    col_list=['#8dd3c7', '#fb8072', '#80b1d3']
    lab_list=['scarce', 'mean', 'abundant']
    for i in range(np.size(R_array)):
        plt.plot(death, CompExcl[:, i], color=col_list[i], label=lab_list[i],
                marker='D', linewidth=2)
    plt.xlabel('toxin sensitivity',fontsize=20)
    plt.ylabel('prob. of competitive exclusion',fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='best', fontsize=16)
    #plt.savefig('CompetitiveExclusion_constantEnv.pdf',bbox_inches='tight',pad_inches=0.05)
    plt.show()
    
#------------------------------------------------    

#Fig.3B: Correlation of non-monotonicity and distances between critical toxin sensitivities
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
        
    
def Fig4B(peak_harsh, peak_mean, peak_mild):
    """
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
