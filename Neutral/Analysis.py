import os
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import stats
import pymc3 as pm
import seaborn as sns
from scipy.stats import dirichlet


def FigA20and21():
  # load data
  d_array=np.linspace(0.1, 1.0, 10) # toxin sensitivity
  nu_array=np.linspace(-5, 3, 9)  # log scale in 3d is problematic in matplot


  both_ext=np.zeros([np.size(d_array), np.size(nu_array)]) # probability that both  species go extinct
  comp_excl1=np.zeros([np.size(d_array), np.size(nu_array)]) # only species 1 goes extinct and species 2 survives
  coexist=np.zeros([np.size(d_array), np.size(nu_array)]) # coexistence of species 1 and 2
  for i in range(np.size(d_array)):
      d=d_array[i]
      os.chdir(str('./death%.1f' %(d)))
      both_ext[i, :] = np.loadtxt('TwoConsumer_BothExtinction_neutral.csv', delimiter=',',
                                 skiprows=1)
      comp_excl1[i,:] = np.loadtxt('TwoConsumer_CompetitiveExclusion_neutral.csv', delimiter=',',
                                 skiprows=1)
      coexist[i, :] = np.loadtxt('TwoConsumer_Coexistence_neutral.csv', delimiter=',',
                                 skiprows=1)
      os.chdir('../')

  # prob that species 2 goes extinct but species 1 persists
  # this should be similar to comp_excl1
  comp_excl2=1-(both_ext+comp_excl1+coexist)
  
  # Fig A20: analyzing the effect of the initial abundance on the extinction probability in mono-culture
  # by assuming that two species are labeled as the same species
  plt.plot(np.linspace(-5, 3, 9), mono_ext[1, :], linewidth=2, marker='D',
        color='#8dd3c7', label='sensitivity 0.2 (small)')

  plt.plot(np.linspace(-5, 3, 9), mono_ext[4, :], linewidth=2, marker='D',
          color='#fb8072', label='sensitivity 0.5 (small)')

  plt.plot(np.linspace(-5, 3, 9), mono_ext[9, :], linewidth=2, marker='D',
          color='#80b1d3', label='sensitivity 1.0 (small)')

  # larger init size = both extinction prob in neutral
  plt.plot(np.linspace(-5, 3, 9), both_ext[1, :], linewidth=2, marker='o',
          linestyle='--', color='#8dd3c7', label='sensitivity 0.2 (large)')

  plt.plot(np.linspace(-5, 3, 9), both_ext[4, :], linewidth=2, marker='o',
          linestyle='--', color='#fb8072', label='sensitivity 0.5 (large)')

  plt.plot(np.linspace(-5, 3, 9), both_ext[9, :], linewidth=2, marker='o',
          linestyle='--', color='#80b1d3', label='sensitivity 1.0 (large)')

  plt.xlabel('$\log_{10}$'+'switching rate',fontsize=20)
  plt.ylabel('prob. of extinction',fontsize=20)
  plt.ylim(-0.05, 1.05)
  plt.xticks(fontsize=16, ticks=[-4, -2, 0, 2])
  plt.yticks(fontsize=16, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
  plt.legend(loc='lower center', bbox_to_anchor=(.5,1.05),
                 ncol=2, fontsize=16)
  plt.savefig('InitPopSize.pdf',bbox_inches='tight', pad_inches=0.05)
  plt.show()
  
  # A21: Neutral scenarios where two species are identical except for their labels



