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
    plt.xlabel('-competitive exclusion/'+r'$\Delta P\left(s_1(T_{end})=0\right)$', fontsize=18)
    plt.ylabel('frequency', fontsize=18)
    plt.yticks([0,3,6,9,12,15,18], fontsize=16)
    plt.xticks(fontsize=16)
    plt.savefig(gname3, bbox_inches='tight', pad_inches=0.1)
    plt.show()

def FigA34():
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
    
    #FigA3
    data=np.loadtxt('BothExtinction_model1.csv', delimiter=',', skiprows=1)
    toxin=np.linspace(0.1, 1, 10)
    col_list=['#8dd3c7', '#fb8072', '#80b1d3']
    lab_list=['scarce', 'mean', 'abundant']
    for i in range(3):
        plt.plot(toxin, data[:, i], color=col_list[i], label=lab_list[i],
                marker='D', linewidth=2)
    plt.xlabel('toxin sensitivity',fontsize=20)
    plt.ylabel('prob. of both exclusion',fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='best', fontsize=16)
    plt.savefig('BothExtinction_FixedEnv.pdf',bbox_inches='tight',pad_inches=0.05)
    plt.show()
  
    
def FigA6():
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
    plt.ylabel('prob. of competitive exclusion',fontsize=20)
    plt.xticks(fontsize=16, ticks=[0.2, 0.6, 1.0, 1.4, 1.8])
    plt.yticks(fontsize=16, ticks=[0, 0.04, 0.08, 0.12, 0.16])
    plt.legend(loc='lower center', bbox_to_anchor=(.5,1.05),
               ncol=2, fontsize=16)
    plt.savefig('CompetitiveExclusion_FixedEnv.pdf',bbox_inches='tight',pad_inches=0.05)
    plt.show()
