#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import stats
import pymc3 as pm
import seaborn as sns
from scipy.stats import dirichlet


def FigA22():
    """
    Analyzing the cases in the absnece of toxin,
    which is equivalent with toxin sensitivity = 0.0
    """
    os.chdir('./OneConsumer')
    mono=np.loadtxt('OneConsumer_Extinction_model1.csv', delimiter=',', skiprows=1)
    os.chdir('../TwoConsumer/type1')
    both=np.loadtxt('TwoConsumer_BothExtinction_model1_competitor1.csv',
                delimiter=',', skiprows=1)[0, :]
    coexist=np.loadtxt('TwoConsumer_Coexistence_model1_competitor1.csv',
                delimiter=',', skiprows=1)[0, :]
    comp_excl=np.loadtxt('TwoConsumer_CompetitiveExclusion_model1_competitor1.csv',
                delimiter=',', skiprows=1)[0, :]
    extinct_co=both+comp_excl # prob that species 1 goes extinct in co-culture
    
    #mono culture
    plt.plot(np.linspace(-5, 3, 9), mono, linewidth=2, marker='o',
            linestyle='-', color='k', label='sensitivity 0.0')
    plt.xlabel('$\log_{10}$'+'switching rate',fontsize=20)
    plt.ylabel('prob. of extinction',fontsize=20)
    plt.ylim(-0.05, 1.05)
    plt.xticks(fontsize=16, ticks=[-4, -2, 0, 2])
    plt.yticks(fontsize=16, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
    plt.legend(loc='lower center', bbox_to_anchor=(.5,1.05),
                   ncol=1, fontsize=16)
    plt.savefig('MonoCulture_NoToxin.pdf',bbox_inches='tight', pad_inches=0.0)
    plt.show()

    # diff in extinction


    plt.plot(np.linspace(-5, 3, 9), mono-extinct_co, linewidth=2, marker='o',
            linestyle='-', color='k', label='sensitivity 0.0')
    plt.xlabel('$\log_{10}$'+'switching rate',fontsize=20)
    plt.ylabel('diff in extinction prob.',fontsize=20)
    plt.ylim(-0.05, 0.0)
    plt.xticks(fontsize=16, ticks=[-4, -2, 0, 2])
    plt.yticks(fontsize=16, ticks=[-0.05, -0.03, -0.01])
    plt.legend(loc='lower center', bbox_to_anchor=(.5,1.05),
                   ncol=1, fontsize=16)
    plt.savefig('DiffExtinction_NoToxin.pdf',bbox_inches='tight', pad_inches=0.0)
    plt.show()

    # competitive exclusion
    plt.plot(np.linspace(-5, 3, 9), comp_excl, linewidth=2, marker='o',
            linestyle='-', color='k', label='sensitivity 0.0')
    plt.xlabel('$\log_{10}$'+'switching rate',fontsize=20)
    plt.ylabel('prob. exclusion of fittest',fontsize=20)
    plt.ylim(0.0, 0.05)
    plt.xticks(fontsize=16, ticks=[-4, -2, 0, 2])
    plt.yticks(fontsize=16, ticks=[0.01, 0.03, 0.05])
    plt.legend(loc='lower center', bbox_to_anchor=(.5,1.05),
                   ncol=1, fontsize=16)
    plt.savefig('CompetitiveExclusion_NoToxin.pdf',bbox_inches='tight', pad_inches=0.0)
    plt.show()
