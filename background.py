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

frac['signal']=1296/6835.
frac['feed']=0.11*1296/6835.
frac['B2DstDsX']=6835/6835.
frac['B2DstD0X']=445/6835.
frac['B2DstDplusX']=0.245*6835/6835.
frac['prompt']=424/6835.

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



histo=[]
histo1=[]
histo2=[]
histo3=[]

for i in range(6):
  df=root_pandas.read_root(files[names[i]],columns=['q2_reco','costheta_L_reco','costheta_D_reco','chi_reco'])
  df=df.sample(n=int(2000000*frac[names[i]]))
  histo.append(df['q2_reco'][~np.isnan(df['q2_reco'])])
  histo1.append(df['costheta_L_reco'][~np.isnan(df['costheta_L_reco'])])
  histo2.append(df['costheta_D_reco'][~np.isnan(df['costheta_D_reco'])])
  histo3.append(df['chi_reco'][~np.isnan(df['chi_reco'])])
plt.hist(histo,stacked=True,range=[0.,11.],bins=8,label=labels,color=colors)
plt.legend()
plt.title(r'$q^2$ reco')
plt.savefig('QTOTAL.pdf')
plt.close()

plt.hist(histo1,stacked=True,range=[-1.,1.],bins=8,label=labels,color=colors)
plt.legend()
plt.title(r'$cos(\theta_L)$ reco')
plt.savefig('COSTHETAL.pdf')
plt.close()


plt.hist(histo2,stacked=True,range=[-1.,1.],bins=8,label=labels,color=colors)
plt.legend()
plt.title(r'$cos(\theta_D)$ reco')
plt.savefig('COSTHETAD.pdf')
plt.close()


plt.hist(histo3,stacked=True,range=[-np.pi,np.pi],bins=8,label=labels,color=colors)
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





