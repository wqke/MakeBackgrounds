
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

#D* Ds X background
#Fractions defined relative to D* Ds*(=1)
#Fit values
frac_fit={}
frac_fit['Bd2DstDs']=0.594
frac_fit['Bd2DstDsst']=1.
frac_fit['Bd2DstDs1']=0.365

frac_fit['Bu2DststDs1']=0.416*5e-4/(5e-4+0.0177+8e-3)   #fraction of the sum x relative branching fractions
frac_fit['Bu2DststDs']=0.416*0.0177/(5e-4+0.0177+8e-3)
frac_fit['Bu2DststDsst']=0.416*8e-3/(5e-4+0.0177+8e-3)
"""
*The BF are taken to be proportional to the following decays
B+->D1(2420)0 Ds*+ : B -> D* Ds  =8e-3
B+->D1(2420)0 Ds+ : B -> D* Ds*  =0.0177
B+->D1(2420)0 Ds1(2460)+ : B -> D* Ds1  =5e-4
"""   
#SUB MODES
#branching fractions
BF={}
BF['Ds1']={}
BF['Dsst']={}
BF['D10']={}
BF['Dplus']={}
BF['D0']={}
BF['Dst0']={}
BF['Dstmin']={}
BF['eta']={}
BF['etap']={}
BF['rho0']={}
BF['rhoplus']={}
BF['omega']={}
BF['Dsplus']={}


###Ds1(2460) decays
#Ds1(2460)->Ds+ gamma
BF['Ds1']['dsgamma']=0.18
#Ds1(2460)->Ds*+ pi0
BF['Ds1']['dsstpi0']=0.48

###Ds*+ decays
#Ds* -> Ds gamma
BF['Dsst']['dsgamma']=0.299
#Ds* -> Ds pi0
BF['Dsst']['dspi0']=0.058

###D1(2420)0 decays
#D1(2420) -> D*- pi+
BF['D10']['dstpiplus']=1.

###D+ decays
#D+ -> Ks0 3pi
BF['Dplus']['Ks3pi']=0.0297
#D+ -> pi+pi+pi-pi0
BF['Dplus']['3pipi0']=0.0111

###D0 decays
#D0 ->K-pi+pi+pi-
BF['D0']['K3pi']=0.0811
#D0 ->K-pi+pi+pi-pi0
BF['D0']['K4pi']=0.042

###D*0 decays
#D*0 -> Ks0 3pi
BF['Dst0']['d0pi0']=0.647
#D*0 -> D0 gamma
BF['Dst0']['d0gamma']=0.353
  
###D*- decays
#D*- ->D0 pi-
BF['Dstmin']['D0pi']=0.677

###eta decay : eta-> pi+ pi- pi0
BF['eta']['3pi']=0.2292
###eta' decays : eta'->eta pi+ pi- 
BF['etap']['etapipi']=0.426
#eta'->rho0 gamma 
BF['etap']['rhogamma']=0.289
###rho0 decay : rho0 ->pi+pi-
BF['rho0']['2pi']=1.
###rho+ decay : rho+ ->pi+pi0
BF['rhoplus']['2pi']=1.
###omega decay : omega ->pi+pi-pi0
BF['omega']['3pi']=0.892

###Ds+ decays
#Ds+->eta pi+
BF['Dsplus']['etapi']=0.017
#Ds+->(eta->3pi) pi+
#BF['Dsplus']['etapi_3pi']=BF['Dsplus']['etapi']*BF['eta']['3pi']
BF['Dsplus']['etapi']=BF['Dsplus']['etapi']*BF['eta']['3pi']
#Ds+->omega pi+
BF['Dsplus']['omegapi']=2.4e-3
#Ds+->(omega->3pi) pi+
#BF['Dsplus']['omegapi_3pi']=BF['Dsplus']['omegapi']*BF['omega']['3pi']
BF['Dsplus']['omegapi']=BF['Dsplus']['omegapi']*BF['omega']['3pi']
#Ds+->eta rho
BF['Dsplus']['etarho']=0.089
#Ds+->(eta->3pi) (rho->2pi)
#BF['Dsplus']['etarho_5pi']=BF['Dsplus']['etarho']*BF['eta']['3pi']*BF['rho0']['2pi']
BF['Dsplus']['etarho']=BF['Dsplus']['etarho']*BF['eta']['3pi']*BF['rho0']['2pi']
#Ds+->eta' rho+
BF['Dsplus']['etaprho']=0.058
#Ds+->omega rho
BF['Dsplus']['omegarho']=0.028 
#Ds+->rho0(pi+ pi-) rho0 (pi+ pi-) pi+
BF['Dsplus']['5pi']=8e-3  


#Ds+->omega pi+pi+pi-
BF['Dsplus']['omega3pi']=0.016
#Ds+->eta' pi+
BF['Dsplus']['etappi']=0.0394
#Ds+->eta' rho+
BF['Dsplus']['etaprho']=0.058

BF['Dsplus']['etappi_etapipi']=BF['Dsplus']['etappi'] * BF['etap']['etapipi'] * BF['eta']['3pi']
BF['Dsplus']['etappi_rhogamma']=BF['Dsplus']['etappi'] * BF['etap']['rhogamma'] * BF['rho0']['2pi']
BF['Dsplus']['etaprho_etapipi']=BF['Dsplus']['etaprho'] * BF['etap']['etapipi'] * BF['eta']['3pi'] * BF['rhoplus']['2pi']
BF['Dsplus']['etaprho_rhogamma']=BF['Dsplus']['etaprho'] * BF['rhoplus']['2pi'] *BF['etap']['rhogamma'] * BF['rho0']['2pi'] 






############PLOT THE TOTAL HISTOGRAMS############

#BBd2DstDK0  
files=['/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDK0/3pipi0_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDK0/ks3pi_LHCb_Total/model_vars.root']

weights0=[]
values0=[]
df_weights0=[]
for file in files:
  df=root_pandas.read_root(file,columns=['q2_reco','costheta_L_reco','costheta_D_reco','chi_reco'],key='DecayTree')
  components=(file.split('/')[-2]).split('_')
  components=components[:-2]   #extract the sub mode from the file name, remove 'LHCb_Total'
  weight=BF['Dstmin']['D0pi'] * BF['Ds1'][components[0]]
  if components[0]=='dsgamma':
    if len(components)==2:
      weight=weight*BF['Dsplus'][components[1]]
    if len(components)==3:
      weight=weight*BF['Dsplus'][components[1]+'_'+components[2]]
  elif components[0]=='dsstpi0': 
    if len(components)==3:
      weight=weight*BF['Dsst']['dsgamma']*BF['Dsplus'][components[2]]
    if len(components)==4:
      weight=weight*BF['Dsst']['dsgamma']*BF['Dsplus'][components[2]+'_'+components[3]]
  weights0.append(weight)
  values0.append(df.values)
  df_weights0.append(np.ones_like(df.values) )
sum0=sum(weights0)
for i in range(len(weights0)):
  weights0[i]=weights0[i]/sum0   #define the weight with regard to the sum (the proportion)
for i in range(len(weights0)):
  df_weights0[i]=df_weights0[i]*weights0[i]

#Bu2DstDK
files1=['/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DstDK/3pipi0_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DstDK/ks3pi_LHCb_Total/model_vars.root']

weights1=[]
values1=[]
df_weights1=[]
for file in files1:
  df=root_pandas.read_root(file,columns=['q2_reco','costheta_L_reco','costheta_D_reco','chi_reco'],key='DecayTree')
  components=(file.split('/')[-2]).split('_')
  components=components[:-2]   #extract the sub mode from the file name
  if len(components)==1:
    weight=BF['Dstmin']['D0pi']* BF['Dsplus'][components[0]]
  if len(components)==2:
    weight=BF['Dstmin']['D0pi'] * BF['Dsplus'][components[0]+'_'+components[1]]
  weights1.append(weight)
  values1.append(df.values)
  df_weights1.append(np.ones_like(df.values) )

sum1=sum(weights1)
for i in range(len(weights1)):
  weights1[i]=weights1[i]/sum1   #define the weight with regard to the sum (the proportion)

for i in range(len(weights1)):
  df_weights1[i]=df_weights1[i]*weights1[i]

########
########
########
########


#df_merged = pd.concat([df, df2], ignore_index=True)
#df_merged_weights = np.ones_like(df_merged.values) / len(df_merged)
def addlist(a,b):
  lis=[]
  for i in range(len(a)):
    lis.append(a[i]+b[i])
  return lis

def timeslist(a,fraction):
  lis=[]
  for i in range(len(a)):
    lis.append(a[i]*fraction)
  return lis



ranges=[[0.,13.],[-1.,1.],[-np.pi,np.pi],[-1.,1.]]
filenames=['q2','costhetaL','chi','costhetaD']
titles=[r'$q^2$',r'cos$(\theta_L)$',r'$\chi$',r'cos$(\theta_D)$']
totalfile=[files,files1,files2,files3,files4,files5]
totalweights=[df_weights0,df_weights1,df_weights2,df_weights3,df_weights4,df_weights5]
totalvalues=[values0,values1,values2,values3,values4,values5]

binnumber=100

fractions=[frac_fit['Bd2DstDs1'],frac_fit['Bd2DstDs'],frac_fit['Bd2DstDsst'],frac_fit['Bu2DststDs'],frac_fit['Bu2DststDsst'],frac_fit['Bu2DststDs1']]
labels=[r'$B^0 \rightarrow D^{*-} D_{s1}(2460)^+$',r'$B^0 \rightarrow D^{*-} D_s^+$',r'$B^0 \rightarrow D^{*-} D_s^{*+}$',r'$B^+ \rightarrow D_1(2420)^0D_s^+$',r'$B^+ \rightarrow D_1(2420)^0D_s^{*+}$',r'$B^+ \rightarrow D_1(2420)^0D_{s1}(2460)^+$']

for i in range(4):
  all_bin_heights=[[0]*binnumber]
  TOTAL_HEIGHTS=[0]*binnumber
  for k in range(6):
    fig, ax = plt.subplots()
    bin_heights, bin_borders, _=ax.hist(totalvalues[k][0][:,i], weights=totalweights[k][0][:,i], histtype='step',bins=binnumber, range=ranges[i])
    bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2
    plt.close()
    for j in range(1,len(weights3)):
      bin_heights_new, bin_borders_new, _=ax.hist(totalvalues[k][j][:,i], weights=totalweights[k][j][:,i], histtype='step',bins=binnumber, range=ranges[i])
      bin_heights=addlist(bin_heights,bin_heights_new)
      plt.close()
    bin_heights=timeslist(bin_heights,1/(sum(bin_heights)))
    all_bin_heights.append(bin_heights)
    TOTAL_HEIGHTS=addlist(timeslist(bin_heights,fractions[k+1]),TOTAL_HEIGHTS)
    plt.bar(bin_centers,all_bin_heights[k+1],width=bin_centers[1]-bin_centers[0],label=labels[k],bottom=all_bin_heights[k])
  #plt.ylim(bottom=0)
  plt.legend()
  plt.title(titles[i]+r'  (B $\rightarrow$ $D^*$ $D_s$ X)')
  plt.savefig(filenames[i]+'_bar.pdf')
  plt.close()
  plt.close()
  plt.close()
  plt.close()  
  plt.close()
  plt.close()
