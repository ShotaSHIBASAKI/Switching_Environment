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
  
  #check neutrality: prob that species 1 excludes 2 and vice versa should be similar under the neutral scenario
  ax1 = plt.subplot(111)
  heatmap1=ax1.pcolor(abs(comp_excl1-comp_excl2)/np.maximum(comp_excl1, comp_excl2), 
                      cmap='Blues',vmin=0, vmax=1)
  ax1.xaxis.tick_bottom()
  plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
  plt.ylabel('toxin sensitivity', fontsize=20)
  plt.title('diff. in competitive exclusion', fontsize=20)
  plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
  plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
  cb1=plt.colorbar(heatmap1)
  cb1.ax.tick_params(labelsize=14)
  plt.savefig('Diff_CompExclusion_neutral.pdf',  bbox_inches='tight', pad_inches=0.05)
  plt.show()
  # species 2 excludes species 1
  ax2 = plt.subplot(111)
  heatmap2=ax2.pcolor(comp_excl1, cmap='Blues',vmin=0)
  ax2.xaxis.tick_bottom()
  plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
  plt.ylabel('toxin sensitivity', fontsize=20)
  plt.title('competitive exclusion of 1 by 2', fontsize=20)
  plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
  plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
  cb2=plt.colorbar(heatmap2)
  cb2.ax.tick_params(labelsize=14)
  plt.savefig('Exclusion1by2_Neutral.pdf',bbox_inches='tight', pad_inches=0.05)
  plt.show()
  # species 1 excludes species 2
  ax3 = plt.subplot(111)
  heatmap3=ax3.pcolor(comp_excl2, cmap='Blues',vmin=0)
  ax3.xaxis.tick_bottom()
  plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
  plt.ylabel('toxin sensitivity', fontsize=20)
  plt.title('competitive exclusion of 2 by 1', fontsize=20)
  plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
  plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
  cb3=plt.colorbar(heatmap3)
  cb3.ax.tick_params(labelsize=14)
  plt.savefig('Exclusion2by1_Neutral.pdf',bbox_inches='tight', pad_inches=0.05)
  plt.show()
  #coexistenece of two species
  ax4 = plt.subplot(111)
  heatmap4=ax4.pcolor(coexist, cmap='Blues',vmin=0)
  ax4.xaxis.tick_bottom()
  plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
  plt.ylabel('toxin sensitivity', fontsize=20)
  plt.title('coexistence', fontsize=20)
  plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
  plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
  cb4=plt.colorbar(heatmap3)
  cb4.ax.tick_params(labelsize=14)
  plt.savefig('Coexistence_Neutral.pdf',bbox_inches='tight', pad_inches=0.05)
  plt.show()
  # beta diversity
  q=1 #  parameter for calculating diversity
  #iteration=100000 # number of iteration using a parameter set
  delta_array=np.linspace(0.1, 1.0, 10)
  nu_exp=np.linspace(-5,3,9)
  alpha_div=np.zeros([np.size(delta_array), np.size(nu_exp)]) #alpha diversity
  beta_div=np.zeros([np.size(delta_array), np.size(nu_exp)]) #alpha diversity
  gamma_div=np.zeros([np.size(delta_array), np.size(nu_exp)]) #alpha diversity 
    
  for i in range(np.size(delta_array)):
      delta=delta_array[i]
      os.chdir(str('./death%.1f' %(delta)))
          
      for l in range(np.size(nu_exp)):
          #convert aboundance data into csv data
          fname=str("Dynamics_nu%d_last.csv" %(nu_exp[l]))
          data=np.loadtxt(fname, delimiter=',', skiprows=1)
          weight=np.sum(data,1)
          if np.sum(weight)>0:
              weight=weight/np.sum(weight) # normalize
          data_freq=ConvertFreq(data[:,  1:3])
            
                    
          #calculate alpha, beta, gamma diversity
          alpha_div[i,l]=Alpha(data_freq,q,weight)
          gamma_div[i,l]=Gamma(data_freq,q, weight)
          beta_div[i, l]=gamma_div[i,l]/alpha_div[i,l]
                
          
          
      """      
      # showing only beta div
      #plt.violinplot(beta_div[i, :])
      plt.plot(np.linspace(1,9,9), beta_div[i, :], marker='D', color='k')
      plt.xticks([2,4,6,8],['$-4$', '$-2$', '$0$', '$2$'],
                       fontsize=24)
      plt.xlabel('$\log_{10}$'+' switching rate', fontsize=28)  
      plt.yticks(fontsize=24)
      #plt.ylabel('beta diversity', fontsize=18)
      plt.savefig(str("BetaDiversities_neutral_death%d.pdf" %(10*delta)),
                        bbox_inches="tight", pad_inches=0.05)
      plt.show()
            
      os.chdir("../")
      """
  ax5 = plt.subplot(111)
  heatmap5=ax5.pcolor(beta_div, cmap='Blues',vmin=1)
  ax5.xaxis.tick_bottom()
  plt.xlabel('$\log_{10}$ switching rate', fontsize=20)
  plt.ylabel('toxin sensitivity', fontsize=20)
  plt.title('beta diversity', fontsize=20)
  plt.xticks([0.5, 2.5, 4.5, 6.5, 8.5],['-5', '-3', '-1', '1','3'], fontsize=16)
  plt.yticks([1.5, 3.5, 5.5, 7.5, 9.5],['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=16)
  cb5=plt.colorbar(heatmap5)
  cb5.ax.tick_params(labelsize=14)
  plt.savefig('Betda_div_Neutral.pdf',bbox_inches='tight', pad_inches=0.05)
  plt.show()
#---- Functions used in beta diversity analysis---------------  
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


  



