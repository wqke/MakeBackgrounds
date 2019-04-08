import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from numpy import *
import numpy as np
import math
from math import cos,sin,pi
import root_pandas
import pandas as pd
from root_numpy import root2array, rec2array, tree2array
from ROOT import TFile,TChain,TTree
from uncertainties import *



files=['B2DstD0X.root','B2DstDplusX.root','B2DstDsX.root','prompt.root']
labels=[r'$B\rightarrow D^*D^0X$',r'$B\rightarrow D^*D^+X$',r'$B\rightarrow D^*D_sX$','prompt']
columns=['q2_reco','costheta_L_reco','costheta_D_reco','chi_reco','Tau_life_reco']
titles=[r'$q^2$',r'cos$(\theta_L)$',r'cos$(\theta_D)$',r'$\chi$',r'$\tau$ life']
ranges=[[0.,13.],[-1.,1.],[-1.,1.],[-np.pi,np.pi],[0.,2.]]
for j in range(5):
  for i in range(len(files)):
    df=root_pandas.read_root(files[i],columns=columns,key='DecayTree')
    plt.hist(df[columns[j]][~np.isnan(df[columns[j]])],bins=100,label=labels[i],density=True,histtype='step',range=ranges[j])
  plt.title(titles[j])
  plt.legend()
  plt.savefig(columns[j]+'_total.pdf')
  plt.close()
  
  
