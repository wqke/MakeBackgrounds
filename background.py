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
names=['signal','B2DstDsX','B2DstD0X','B2DstDplusX']
files={'signal':'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstTauNu/3pi_all_Total/model_vars.root',
      'B2DstD0X':'/home/ke/graphs/B2DstD0X.root',
      'B2DstDplusX':'/home/ke/graphs/B2DstDplusX.root',
      'B2DstDsX':'/home/ke/graphs/B2DstDsX.root'}

DF=root_pandas.read_root(files[names[0]],columns=['q2_reco','Tau_life_reco'])
DF=DF.sample(n=int(2000000*frac[names[0]]))
for i in range(1,4):
  df=root_pandas.read_root(files[names[i]],columns=['q2_reco','Tau_life_reco'])
  df=df.sample(n=int(2000000*frac[names[i]]))
  DF=pd.concat([DF, df], ignore_index=True)
plt.hist(DF['q2_reco'][~np.isnan(DF['q2_reco'])],histtype='step',range=[0.,13.],bins=100)
plt.title(r'$q^2$ reco')
plt.savefig('QTOTAL.pdf')
plt.close()

plt.hist(DF['Tau_life_reco'][~np.isnan(DF['Tau_life_reco'])],histtype='step',range=[0.,2.],bins=100)
plt.title(r'$\tau$ life time')
plt.savefig('TAUTOTAL.pdf')
plt.close()






