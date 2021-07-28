#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 16:31:00 2020
To plot Fig A.5, 7-11, modify main_figures.py, divertsity.py, or richness.py
Figures in appendices
"""
import sympy as sym
from IPython.display import display
from IPython.display import Math
import numpy as np
import matplotlib.pyplot as plt
import os

def FigA1():
    """
    Analytical calulation in appendix 1 
    with plotting Fig.A1
    """
    #define symbols
    r=sym.Symbol(r'r^*')  # resource
    s=sym.Symbol(r's^*')  # species
    t=sym.Symbol(r't^*')  # toxin
    alpha=sym.Symbol(r'alpha') # dilution rate
    mu=sym.Symbol(r'mu')  # maximum growth rate
    Kr=sym.Symbol(r'K_r')  # median-effect resoruce amount
    Yr=sym.Symbol(r'Y_r')  # biomass yield of resource
    r_mean=sym.Symbol(r'\left<r\right>')  # mean amount of supplied resource
    delta=sym.Symbol(r'delta')  # maximum death rate
    Kt=sym.Symbol(r'K_t') # median-effect toixn  amount
    Yt=sym.Symbol(r'Y_t') # biomass yield of toxin
    t_mean=sym.Symbol(r'\left<t\right>') # mean amount of supplied toxin
    
    #define quadatic equations
    Qx=-r**2+(r_mean-Kr-mu*s/(alpha*Yr))*r+Kr*r_mean
    Qz=-t**2+(t_mean-Kt-delta*s/(alpha*Yt))*t+Kt*t_mean

    solx1, solx2=sym.solve(Qx, r)  # solx1 (solx2) is a negative (positive) root of Qx
    solz1, solz2=sym.solve(Qz, t)  # solx1 (solx2) is a negative (positive) root of Qz
    #display the positive roots
    #display(sym.simplify(solx2))
    #display(sym.simplify(solz2))
    #latex style
    #sym.latex(sym.simplify(solx2))
    #sym.latex(sym.simplify(solz2))
    
    #obtaining y*
    Eq=s-Yt*t+Yr*r-(Yr*r_mean-Yt*t_mean)  # see equilibrium
    Eq_y=Eq.subs([(r, solx2), (t, solz2)])
    #display(Eq_y)
    
    #here an example of parameter set
    Eq_y_subs=sym.simplify(Eq_y.subs([(alpha, 0.1), (mu, 1), (Kr, 100), (Yr, 1), (r_mean, 200),
                    (delta,1.2), (Kt, 100), (Yt, 1), (t_mean,125)]))
    #display(Eq_y_subs)
    s_subs=np.linspace(0, 100, 1000)
    sol=np.zeros((np.size(s_subs)))
    for i in range(np.size(s_subs)):
        sol[i]=Eq_y_subs.subs(s, s_subs[i])
    plt.plot(s_subs, sol,color='b', linestyle='-', linewidth=3)
    plt.plot(s_subs, np.zeros(np.size(s_subs)), color='k', linestyle='--', linewidth=2)
    plt.xlabel(r'$s$', fontsize=18)
    plt.ylabel(r'$f(s)$', fontsize=18)
    plt.xlim(0, 100)
    plt.savefig('Example_species_equilibrium.pdf',bbox_inches='tight', pad_inches=0.1)
    plt.show()
    print(sol[0])
 #-----------------------------   
def FigA2():
    model=1
    """
    Plotting the distritbution of competitive exclusion/ delta P
    model: determine the environmental switching type: deafult is 1
    """
    os.chdir('../main_text_scenario1')
    os.getcwd()
    gname3=str('CompetitiveExclusion_DiffProb_model%d_difference.pdf'%(model))
    d_array=np.arange(0.1, 0.9, 0.1) # sign of species interaction is non-positive
    delta_ext=np.zeros([1, 9*np.size(d_array)])
    comp_excl=np.zeros([1, 9*np.size(d_array)])
    rel_delta=np.array([])
    sum_prob=np.array([])
    for i in range(np.size(d_array)):
        death=d_array[i]
        path=str('./death%.1f/TwoConsumer' %(death))
        os.chdir(path)
        #competitive exclusion probability
        fname=str('TwoConsumer_CompetitiveExclusion_model%d_competitor1.csv' %(model))
        data=np.loadtxt(fname, delimiter=',', skiprows=1)
        comp_excl[:, i*9:i*9+9]=data[0, :] # compare species mu=0.91
        # diff in extinction probabilities
        fname=str('TwoConsumer_Extinction_model%d_competitor1.csv' %(model))
        data=np.loadtxt(fname, delimiter=',', skiprows=1)
        delta_ext[:, i*9:i*9+9]=-data[0, :]
        
        os.chdir('../OneConsumer')
        fname=str('OneConsumer_Extinction_model%d.csv' %(model))
        data=np.loadtxt(fname, delimiter=',', skiprows=1)
        delta_ext[:, i*9:i*9+9]+=data
        os.chdir('../../')
    for j in range(1):
        for k in range(np.size(delta_ext[j, :])):
            if abs(delta_ext[j,k])<=pow(10, -5):
                """
                Need to avoid zero division.
                Note: number of simulation is 10^5
                if abs(delta_ext)< 10^-5, this should come from rounding error.
                """
                delta_ext[j,k]=pow(10, -6)
                if comp_excl[j, k]<pow(10, -5):
                    """
                    if prob of competitive exclusion is also too small due to rounding error,
                    rerio of cmp_excl to delta_ext should be -1.
                    """
                    comp_excl[j, k]=pow(10, -6) 
        
        rel=-comp_excl[j, :]/delta_ext[j, :]
        
        rel_delta=np.concatenate((rel_delta, rel))
    
    plt.hist(rel_delta, bins=16,  range=(0, 1.6), color='b')
    plt.xlabel('- exclusion of the fittest/'+r'$\Delta P\left(s_1(\sigma_{end})=0\right)$', fontsize=18)
    plt.ylabel('frequency', fontsize=18)
    plt.yticks([0,3,6,9,12,15,18], fontsize=16)
    plt.xticks(fontsize=16)
    plt.savefig(gname3, bbox_inches='tight', pad_inches=0.1)
    plt.show()
#--------------------------------

#FigA3

"""
The csv files used here are to huge to upload on Github. Therefore,
we just show the source code.
Please contast SS if you need the original data:
shota.shibasaki@unil.ch
or 
shota-shibasaki@18.alumni.u-tokyo.ac.jp
Instead, we upload csv files that contains alpha, beta, and gamma diversity
on GitHub
Each row: different community
Each column different switching rate 10^-5 -- 10^3
"""
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
        
def Fig3():
    q=1 #  parameter for calculating diversity
    parm_set=100  # number of parameter sets
    iteration=100 # number of iteration using a parameter set
    delta_array=np.array([0.1,0.2, 0.4, 0.6, 1.0])
    species_array=np.array([2,4,6,8,10])
    nu_exp=np.linspace(-5,3,9)
    
    for i in range(np.size(delta_array)):
        delta=delta_array[i]
        path=str('./death%.1f' %(delta))
        os.chdir(path)
        for j in range(np.size(species_array)):
            n=species_array[j]
            path=str('species_number%d' %(n))
            
            #initialize diversity matrix
            alpha_div=np.zeros([iteration, np.size(nu_exp)]) #alpha diversity
            beta_div=np.zeros([iteration, np.size(nu_exp)]) #alpha diversity
            gamma_div=np.zeros([iteration, np.size(nu_exp)]) #alpha diversity
            os.chdir(path)
            for k in range(parm_set):
                os. chdir(str('./parameter%d' %(k+1)))
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
                os.chdir('../')  # end analyzing parameter k data
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
            
            # violin plot over swiching rate
            fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, 
                                                figsize=(10, 4), sharey=False)
            #alpha diversity
            #ax1.set_title('alpha diversity')
            ax1.violinplot(alpha_div)
            alpha_mean=np.mean(alpha_div, 0)
            ax1.plot(np.linspace(1,9,9), alpha_mean, marker='D', color='k')
            ax1.set_xticks([2,4,6,8])
            ax1.set_xticklabels(['$-4$', '$-2$', '$0$', '$2$'], fontsize=16)
            ax1.tick_params(axis='y',labelsize=16)
            
            #beta diversity
            #ax2.set_title('beta diversity')
            if i == np.size(delta_array)-1:
                ax2.set_xlabel(r'$\log_{10}$ switching rate', fontsize=20)
            ax2.violinplot(beta_div)
            beta_mean=np.mean(beta_div, 0)
            ax2.plot(np.linspace(1,9,9), beta_mean, marker='D', color='k')
            ax2.set_xticks([2,4,6,8])
            ax2.set_xticklabels(['$-4$', '$-2$', '$0$', '$2$'],fontsize=16)
            ax2.tick_params(axis='y',labelsize=16)
            
            #gamma diversity
            #ax3.set_title('gamma diversity')
            ax3.violinplot(gamma_div)
            gamma_mean=np.mean(gamma_div, 0)
            ax3.plot(np.linspace(1,9,9), gamma_mean, marker='D', color='k')
            ax3.set_xticks([2,4,6,8])
            ax3.set_xticklabels(['$-4$', '$-2$', '$0$', '$2$'], fontsize=16)
            ax3.tick_params(axis='y',labelsize=16)
            
            plt.savefig(str("ABGdiversities_q%d_species%d_death%d.pdf" %(q, n, 10*delta)),
                        bbox_inches="tight", pad_inches=0.05)
            plt.show()
            """
            # showing only beta div
            plt.violinplot(beta_div)
            plt.plot(np.linspace(1,9,9), beta_mean, marker='D', color='k')
            plt.xticks([2,4,6,8],['$-4$', '$-2$', '$0$', '$2$'],
                       fontsize=24)
            plt.xlabel('$\log_{10}$'+' switching rate', fontsize=28)  
            plt.yticks(fontsize=24)
            #plt.ylabel('beta diversity', fontsize=18)
            plt.savefig(str("Betadiversities_q%d_species%d_death%d.pdf" %(q, n, 10*delta)),
                        bbox_inches="tight", pad_inches=0.05)
            plt.show()
            
            #save csv for beta diversity
            header='switching rate: 10^-5, -4, -3, -2, -1, 0, 1, 2, 3'
            cname=str("beta_div_species%d_death%d.csv" %(n, 10*delta))
            np.savetxt(cname, beta_div, delimiter=',', header=header, fmt='%.6f')
            """
            # save csv files
            Aname=str('AlphaDiversity_death%1.f_%dspecies.csv' %(delta, n))
            np.savetxt(Aname, alpha_div, delimiter=',')
            Bname=str('ABetaDiversity_death%1.f_%dspecies.csv' %(delta, n))
            np.savetxt(Bname, beta_div, delimiter=',')
            Gname=str('AlphaDiversity_death%1.f_%dspecies.csv' %(delta, n))
            np.savetxt(Gname, gamma_div, delimiter=',')
            
            os.chdir("../")
        os.chdir("../")

#--------------------------------
def FigA4():
    #analyzing scenario 1 without environmental switching
    os.chdir('./Appendix3_constant_env/scenario1')
    #Fig A3
    Diff=np.loadtxt('Diff_extinction_model1.csv', delimiter=',', skiprows=0)
    gname_mag='Diff_extinction_FixedEnv.pdf' 
    ax = plt.subplot(111)
    heatmap=ax.pcolor(Diff, cmap='Blues_r')
    #heatmap=ax.pcolor(comp_array, cmap=plt.cm.bwr,vmin=0)
    ax.xaxis.tick_bottom()
    plt.xlabel('resource supply', fontsize=18)
    plt.ylabel('death rate', fontsize=18)
    plt.title('diff. in extinction prob.', fontsize=18)
    plt.xticks([0.5, 1.5, 2.5],['50', '225', '400'], fontsize=14)
    plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=14)
    cb=plt.colorbar(heatmap)
    cb.ax.tick_params(labelsize=14)
    plt.savefig(gname_mag,bbox_inches='tight', pad_inches=0.1)
    plt.show()
    

  
    
def FigA8():
    os.chdir('./Appendix3_constant_env/scenario2_mild')
    comp2=np.loadtxt('CompetitiveExclusion_model2.csv', delimiter=',', skiprows=0)
    os.chdir('../scenario3_mild')
    comp3=np.loadtxt('CompetitiveExclusion_model3.csv', delimiter=',', skiprows=0)
    os.chdir('../')
    toxin=np.linspace(0.1, 2.0, 20)
    plt.plot(toxin, comp2, color='#8dd3c7', label='scenario 2',
                marker='D', linewidth=2)
    plt.plot(toxin, comp3, color='#80b1d3', label='scenario 3',
                marker='D', linewidth=2)
    plt.xlabel('toxin sensitivity',fontsize=20)
    plt.ylabel('prob. of exclusion of fittest',fontsize=20)
    plt.xticks(fontsize=16, ticks=[0.2, 0.6, 1.0, 1.4, 1.8])
    plt.yticks(fontsize=16, ticks=[0, 0.04, 0.08, 0.12, 0.16])
    plt.legend(loc='lower center', bbox_to_anchor=(.5,1.05),
               ncol=2, fontsize=16)
    plt.savefig('CompetitiveExclusion_FixedEnv.pdf',bbox_inches='tight',pad_inches=0.05)
    plt.show()

#--------------------
#different resource supplies
def FigA9to11():
    os.chdir('./Appendix5_change_supply')
    resource_list=['larger_Resource-', 'larger_Resource+',
                   'smaller_Resource-', 'smaller_Resource+']
    for i in range(len(resource_list)):
        os.chdir('./'+resource_list[i])
        nu_array=np.linspace(-5, 3, 9)  # log scale
        if i==1:
            d_array=np.linspace(0.1, 2.0, 20)
            both_ext=np.zeros([np.size(d_array), np.size(nu_array)])
        else:
            d_array=np.linspace(0.1, 1.0, 10) 
            
       
        one_ext=np.zeros([np.size(d_array), np.size(nu_array)])
        two_ext=np.zeros([np.size(d_array), np.size(nu_array)])
        comp_exc=np.zeros([np.size(d_array), np.size(nu_array)])
        for j in range(np.size(d_array)):
            death=d_array[j]
            path1=str('./death%.1f/Model1/OneConsumer' %(death))
            os.chdir(path1)
            one_ext[j, :]=np.loadtxt('OneConsumer_Extinction_model1.csv',
                                     delimiter=',', skiprows=1)
            path2=str('../TwoConsumer')
            os.chdir(path2)
            two_ext[j,:]=np.loadtxt('TwoConsumer_Extinction_model1_competitor1.csv',
                                     delimiter=',', skiprows=1)[0, :]
            comp_exc[j,:]=np.loadtxt('TwoConsumer_CompetitiveExclusion_model1_competitor1.csv',
                                     delimiter=',', skiprows=1)[0, :]
            
            if i==1:
                both_ext[j,:]=np.loadtxt('TwoConsumer_BothExtinction_model1_competitor1.csv',
                                     delimiter=',', skiprows=1)[0, :]
            os.chdir('../../../')
        diff=one_ext - two_ext
        #Diff in extinction
        ax1 = plt.subplot(111)
        heatmap1=ax1.pcolor(diff, cmap='Blues_r')
        #heatmap=ax.pcolor(comp_array, cmap=plt.cm.bwr,vmin=0)
        ax1.xaxis.tick_bottom()
        plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
        plt.ylabel('toxin sensitivity', fontsize=20)
        plt.title('diff. in extinction prob.', fontsize=20)
        plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
        plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
        if i==1:
             plt.yticks([3.5, 7.5, 11.5, 15.5, 19.5],['0.4', '0.8', '1.2', '1.6', '2.0'], fontsize=16)
        else:
             plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
        cb1=plt.colorbar(heatmap1, ticks=[-0.15, -0.1, -0.05, 0])
        cb1.ax.tick_params(labelsize=14)
        plt.savefig('DiffExtinction_Heatmap.pdf',
                  bbox_inches='tight', pad_inches=0.05)
        plt.show()
        
        
        # exclusion
        ax1 = plt.subplot(111)
        heatmap1=ax1.pcolor(comp_exc, cmap='Blues',vmin=0)
        #heatmap=ax.pcolor(comp_array, cmap=plt.cm.bwr,vmin=0)
        ax1.xaxis.tick_bottom()
        plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
        plt.ylabel('toxin sensitivity', fontsize=20)
        if i==1:
             plt.yticks([3.5, 7.5, 11.5, 15.5, 19.5],['0.4', '0.8', '1.2', '1.6', '2.0'], fontsize=16)
        else:
             plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
        plt.title('prob. of exclusion of fittest', fontsize=20)
        plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
        
        cb1=plt.colorbar(heatmap1, ticks=[0,0.05, 0.1, 0.15])
        cb1.ax.tick_params(labelsize=14)
        plt.savefig('CompetitiveExclusion_Heatmap.pdf',
                    bbox_inches='tight', pad_inches=0.05)
        plt.show()
        
        
        if i==1:
            
            # extract positive interaction
            plt.plot(nu_array, diff[9, :], color='k', linewidth=2,
                     marker='D')
            plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
            plt.ylabel('diff. in extinction prob.', fontsize=20)
            plt.xticks([-4, -2, 0, 2], fontsize=16)
            plt.yticks(fontsize=16)
            plt.savefig('Positive_interaction_diff.pdf',
                        bbox_inches='tight', pad_inches=0.05)
            plt.show()
            
            plt.plot(nu_array, comp_exc[9, :], color='k', linewidth=2,
                     marker='D')
            plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
            plt.ylabel('prob.of exclusion of fittest', fontsize=20)
            plt.xticks([-4, -2, 0, 2], fontsize=16)
            plt.yticks(fontsize=16)
            plt.savefig('Positive_interaction_Excl.pdf',
                        bbox_inches='tight', pad_inches=0.05)
            
            plt.show()
            
            plt.plot(nu_array, one_ext[9, :]- both_ext[9, :], color='k', linewidth=2,
                     marker='D')
            plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
            plt.ylabel('mono - both extinction prob.', fontsize=20)
            plt.xticks([-4, -2, 0, 2], fontsize=16)
            plt.yticks(fontsize=16)
            plt.savefig('Positive_interaction_Mono-Both.pdf',
                        bbox_inches='tight', pad_inches=0.05)
            plt.show()
        os.chdir('../')
    
    
#--------------------
#other forms of enviornmental flucutuations
def FigA12():
    os.chdir('./OtherFluctuations/Asymmetric')
    
    d_array=np.linspace(0.1, 1.0, 10) # toxin sensitivity
    nu_array=np.linspace(-5, 3, 9)  # log scale
    for asym in range(1,3):
        os.chdir('./%d' %(asym))
        one_ext=np.zeros([np.size(d_array), np.size(nu_array)])
        two_ext=np.zeros([np.size(d_array), np.size(nu_array)])
        comp_exc=np.zeros([np.size(d_array), np.size(nu_array)])
        
        for i in range(np.size(d_array)):
            death=d_array[i]
            path1=str('./death%.1f/Model1/OneConsumer' %(death))
            os.chdir(path1)
            one_ext[i, :]=np.loadtxt('OneConsumer_Extinction_model1.csv',
                                     delimiter=',', skiprows=1)
            path2=str('../TwoConsumer/type1')
            os.chdir(path2)
            two_ext[i,:]=np.loadtxt('TwoConsumer_Extinction_model1_competitor1.csv',
                                     delimiter=',', skiprows=1)
            comp_exc[i,:]=np.loadtxt('TwoConsumer_CompetitiveExclusion_model1_competitor1.csv',
                                     delimiter=',', skiprows=1)
            os.chdir('../../../../')
        diff=one_ext - two_ext
        #Diff in extinction
        ax1 = plt.subplot(111)
        heatmap1=ax1.pcolor(diff, cmap='Blues_r')
        #heatmap=ax.pcolor(comp_array, cmap=plt.cm.bwr,vmin=0)
        ax1.xaxis.tick_bottom()
        plt.xlabel('$\log_{10}$ basal switching rate', fontsize=20)
        plt.ylabel('toxin sensitivity', fontsize=20)
        plt.title('diff. in extinction prob', fontsize=20)
        plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
        plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
        cb1=plt.colorbar(heatmap1, ticks=[-0.15, -0.1, -0.05, 0])
        cb1.ax.tick_params(labelsize=14)
        plt.savefig(str('DiffExtinction_Asymmetric%d_Heatmap.pdf' %(asym)),
                  bbox_inches='tight', pad_inches=0.05)
        plt.show()
        
        
        # exclusion
        ax1 = plt.subplot(111)
        heatmap1=ax1.pcolor(comp_exc, cmap='Blues',vmin=0)
        #heatmap=ax.pcolor(comp_array, cmap=plt.cm.bwr,vmin=0)
        ax1.xaxis.tick_bottom()
        plt.xlabel('$\log_{10}$ basal switching rate', fontsize=20)
        plt.ylabel('toxin sensitivity', fontsize=20)
        plt.title('prob. of exclusion of fittest', fontsize=20)
        plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
        plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
        cb1=plt.colorbar(heatmap1, ticks=[0,0.05, 0.1, 0.15])
        cb1.ax.tick_params(labelsize=14)
        plt.savefig(str('CompetitiveExclusion_Asymmetric%d_Heatmap.pdf' %(asym)),
                    bbox_inches='tight', pad_inches=0.05)
        plt.show()
        
        os.chdir('../')
    
    
    
    
def FigA13():
    os.chdir('./OtherFluctuations/FourState')
    d_array=np.linspace(0.1, 1.0, 10) # toxin sensitivity
    nu_array=np.linspace(-5, 3, 9)  # log scale
    
    one_ext=np.zeros([np.size(d_array), np.size(nu_array)])
    two_ext=np.zeros([np.size(d_array), np.size(nu_array)])
    comp_exc=np.zeros([np.size(d_array), np.size(nu_array)])
    
    for i in range(np.size(d_array)):
        death=d_array[i]
        path1=str('./death%.1f/Model1/OneConsumer' %(death))
        os.chdir(path1)
        one_ext[i, :]=np.loadtxt('OneConsumer_Extinction_model1.csv',
                                 delimiter=',', skiprows=1)
        path2=str('../TwoConsumer/type1')
        os.chdir(path2)
        two_ext[i,:]=np.loadtxt('TwoConsumer_Extinction_model1_competitor1.csv',
                                 delimiter=',', skiprows=1)
        comp_exc[i,:]=np.loadtxt('TwoConsumer_CompetitiveExclusion_model1_competitor1.csv',
                                 delimiter=',', skiprows=1)
        os.chdir('../../../../')
    diff=one_ext - two_ext
    #Diff in extinction
    ax1 = plt.subplot(111)
    heatmap1=ax1.pcolor(diff, cmap='Blues_r')
    #heatmap=ax.pcolor(comp_array, cmap=plt.cm.bwr,vmin=0)
    ax1.xaxis.tick_bottom()
    plt.xlabel('$\log_{10}$ fluctuation rate', fontsize=20)
    plt.ylabel('toxin sensitivity', fontsize=20)
    plt.title('diff. in extinction prob', fontsize=20)
    plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
    plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
    cb1=plt.colorbar(heatmap1, ticks=[-0.15, -0.1, -0.05, 0])
    cb1.ax.tick_params(labelsize=14)
    plt.savefig('DiffExtinction_Cyclic_Heatmap.pdf',
              bbox_inches='tight', pad_inches=0.05)
    plt.show()
    
    
    # exclusion
    ax1 = plt.subplot(111)
    heatmap1=ax1.pcolor(comp_exc, cmap='Blues',vmin=0)
    #heatmap=ax.pcolor(comp_array, cmap=plt.cm.bwr,vmin=0)
    ax1.xaxis.tick_bottom()
    plt.xlabel('$\log_{10}$ fluctuation rate', fontsize=20)
    plt.ylabel('toxin sensitivity', fontsize=20)
    plt.title('prob. of exclusion of fittest', fontsize=20)
    plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
    plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
    cb1=plt.colorbar(heatmap1, ticks=[0,0.05, 0.1, 0.15])
    cb1.ax.tick_params(labelsize=14)
    plt.savefig('CompetitiveExclusion_Cyclic_Heatmap.pdf',
                bbox_inches='tight', pad_inches=0.05)
    plt.show()
    
    #------------------------------------------------    

#Fig.A14: Correlation of non-monotonicity and distances between critical toxin sensitivities
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
