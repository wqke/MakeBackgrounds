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

#Background+signal fractions difined wrt the most populated one (Ds mode)
frac={}

frac['signal']=1296/1296.
frac['feed']=0.11*1296/1296.
frac['B2DstDsX']=6835/1296.
frac['B2DstD0X']=445/1296.
frac['B2DstDplusX']=0.245*6835/1296.
frac['prompt']=424/1296.

#Neglect at this stage prompt and feed-down (small fraction)
#The signal is represented by B2Dsttaunu tau->3pinu purely
names=['B2DstD0X','prompt','B2DstDplusX','B2DstDsX','feed','signal']
#labels={label1: r'B$\rightarrow D^* D^0 X$', label2: 'prompt',label3: r'B$\rightarrow D^* D^+ X$',label4: r'B$\rightarrow D^* D_s X$',label5: 'feed-down', label6: 'signal'}
labels=[ r'B$\rightarrow D^* D^0 X$', 'prompt', r'B$\rightarrow D^* D^+ X$', r'B$\rightarrow D^* D_s X$', 'feed-down',  'signal']
colors=['#ef10c3','#cecacd','blue','orange','lightblue','red']
files={'signal':'/home/ke/graphs/signal.root',
      'B2DstD0X':'/home/ke/graphs/B2DstD0X.root',
      'B2DstDplusX':'/home/ke/graphs/B2DstDplusX.root',
      'B2DstDsX':'/home/ke/graphs/Ds.root',
      'feed':'/home/ke/graphs/feed.root',
      'prompt':'/home/ke/graphs/prompt.root'}

M_B = 5.27963
M_Dst = 2.01026

q2_max = (M_B - M_Dst)**2
q2_min = 0.

histo=[]
histo1=[]
histo2=[]
histo3=[]
qc_bin_vals = {}
for i in range(5):
  df=root_pandas.read_root(files[names[i]],columns=['q2_reco','costheta_L_reco','costheta_D_reco','chi_reco'])
  df = df.query("q2_reco>=0. and costheta_L_reco>=-1. and costheta_L_reco<1. and costheta_D_reco>=-1. and costheta_D_reco<1. and chi_reco>=-np.pi and chi_reco<np.pi")
  df=df.sample(n=int(2000000*frac[names[i]]))
  histo.append(df['q2_reco'][~np.isnan(df['q2_reco'])])
  histo1.append(df['costheta_L_reco'][~np.isnan(df['costheta_L_reco'])])
  histo2.append(df['costheta_D_reco'][~np.isnan(df['costheta_D_reco'])])
  histo3.append(df['chi_reco'][~np.isnan(df['chi_reco'])])
      
bin_sample=root_pandas.read_root(files[names[5]],columns=['q2_reco','costheta_L_reco','costheta_D_reco','chi_reco'])    
bin_sample = bin_sample.query("q2_reco>=0. and costheta_L_reco>=-1. and costheta_L_reco<1. and costheta_D_reco>=-1. and costheta_D_reco<1. and chi_reco>=-np.pi and chi_reco<np.pi")
bin_sample=bin_sample.sample(n=int(2000000*frac[names[5]]))
histo.append(bin_sample['q2_reco'][~np.isnan(bin_sample['q2_reco'])])
histo1.append(bin_sample['costheta_L_reco'][~np.isnan(bin_sample['costheta_L_reco'])])
histo2.append(bin_sample['costheta_D_reco'][~np.isnan(bin_sample['costheta_D_reco'])])
histo3.append(bin_sample['chi_reco'][~np.isnan(bin_sample['chi_reco'])])

for b in ["q2_reco"]:
  print b
  qc = pd.qcut(bin_sample[b], q=8, precision=5)
  qc_bins = qc.unique()
  qc_bin_vals[b] = []
  for i in range(0,8):
    qc_bin_vals[b].append(qc_bins[i].left)
    qc_bin_vals[b].append(qc_bins[i].right)
    #Retain unique values then sort
    qc_bin_vals[b] = list(set(qc_bin_vals[b]))
    qc_bin_vals[b].sort()
    print qc_bin_vals[b]    
      
plt.hist(histo,stacked=True,range=[0.,11.],bins=qc_bin_vals["q2_reco"] ,label=labels,color=colors)
plt.legend()
plt.title(r'$q^2$ reco')
plt.savefig('QTOTAL.pdf')
plt.close()

plt.hist(histo1,stacked=True,range=[-1.,1.],bins=qc_bin_vals["q2_reco"],label=labels,color=colors)
plt.legend()
plt.title(r'$cos(\theta_L)$ reco')
plt.savefig('COSTHETAL.pdf')
plt.close()


plt.hist(histo2,stacked=True,range=[-1.,1.],bins=qc_bin_vals["q2_reco"],label=labels,color=colors)
plt.legend()
plt.title(r'$cos(\theta_D)$ reco')
plt.savefig('COSTHETAD.pdf')
plt.close()


plt.hist(histo3,stacked=True,range=[-np.pi,np.pi],bins=qc_bin_vals["q2_reco"],label=labels,color=colors)
plt.legend()
plt.title(r'$\chi$ reco')
plt.savefig('CHI.pdf')
plt.close()



"""
plt.hist(DF['Tau_life_reco'][~np.isnan(DF['Tau_life_reco'])],histtype='step',range=[0.,2.],bins=100)
plt.title(r'$\tau$ life time')
plt.savefig('TAUTOTAL.pdf')
plt.close()
"""





