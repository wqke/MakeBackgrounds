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

#D* D0 X background
#Fractions defined with regard to the biggest one. B+ and B0 are assumed to have the same amount (so we directly use the branching fractions)
#Fit values
frac_fit={}
frac_fit['Bu2DstD0K0']=3.8e-3/1.06e-2
frac_fit['Bu2DstDst0K0']=9.2e-3/1.06e-2
frac_fit['Bd2DstD0K']=2.47e-3/1.06e-2
frac_fit['Bd2DstDst0K']= 1.06e-2/1.06e-2
frac_fit['Bd2DstDstK0']=8.1e-3/1.06e-2

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

BFerr={}
BFerr['Ds1']={}
BFerr['Dsst']={}
BFerr['D10']={}
BFerr['Dplus']={}
BFerr['D0']={}
BFerr['Dst0']={}
BFerr['Dstmin']={}
BFerr['eta']={}
BFerr['etap']={}
BFerr['rho0']={}
BFerr['rhoplus']={}
BFerr['omega']={}
BFerr['Dsplus']={}



###Ds1(2460) decays
#Ds1(2460)->Ds+ gamma
BF['Ds1']['dsgamma']=0.18
BFerr['Ds1']['dsgamma']=0.04

#Ds1(2460)->Ds*+ pi0
BF['Ds1']['dsstpi0']=0.48
BFerr['Ds1']['dsstpi0']=0.11




###Ds*+ decays
#Ds* -> Ds gamma
BF['Dsst']['dsgamma']=0.299
BFerr['Dsst']['dsgamma']=0.007

#Ds* -> Ds pi0
BF['Dsst']['dspi0']=0.058
BFerr['Dsst']['dspi0']=0.007


###D1(2420)0 decays
#D1(2420) -> D*- pi+
BF['D10']['dstpiplus']=1.
BFerr['D10']['dstpiplus']=0.




###D+ decays
#D+ -> Ks0 3pi
BF['Dplus']['Ks3pi']=0.0297
BFerr['Dplus']['Ks3pi']=0.0011

#D+ -> pi+pi+pi-pi0
BF['Dplus']['3pipi0']=0.0111
BFerr['Dplus']['3pipi0']=8e-4


###D0 decays
#D0 ->K-pi+pi+pi-
BF['D0']['k3pi']=0.0811
BFerr['D0']['k3pi']=

#D0 ->K-pi+pi+pi-pi0
BF['D0']['k3pipi0']=0.042
BFerr['D0']['k3pipi0']=


###D*0 decays
#D*0 -> Ks0 3pi
BF['Dst0']['d0pi0']=0.647
BFerr['Dst0']['d0pi0']=

#D*0 -> D0 gamma
BF['Dst0']['d0gamma']=0.353
BFerr['Dst0']['d0gamma']=

  
###D*- decays
#D*- ->D0 pi-
BF['Dstmin']['D0pi']=0.677
BFerr['Dstmin']['D0pi']=


###eta decay : eta-> pi+ pi- pi0
BF['eta']['3pi']=0.2292
BFerr['eta']['3pi']=
###eta' decays : eta'->eta pi+ pi- 
BF['etap']['etapipi']=0.426
BFerr['etap']['etapipi']=
#eta'->rho0 gamma 
BF['etap']['rhogamma']=0.289
BFerr['etap']['rhogamma']=
###rho0 decay : rho0 ->pi+pi-
BF['rho0']['2pi']=1.
BFerr['rho0']['2pi']=
###rho+ decay : rho+ ->pi+pi0
BF['rhoplus']['2pi']=1.
BF['rhoplus']['2pi']=BF['rhoplus']['2pi']=BF['rhoplus']['2pi']=
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
def inter(lst1, lst2):
    return list(set(lst1) & set(lst2))


columns=[ 'Tau_FD_z',  'Tau_M', 'Tau_E','Tau_P', '3pi_M', '3pi_PZ', 'Tau_m12', 'Tau_m13','Tau_m23',
         'Tau_FD', 'costheta_D_reco','costheta_L_reco','q2_reco','Tau_PZ_reco','Tau_PT',
        'chi_reco', 'Tau_life_reco']

columns=['q2_reco']




#Bu2DstD0K0
files=['/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bu2DstD0K0/k3pi_LHCb_Total/model_vars.root',
       '/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bu2DstD0K0/k3pipi0_LHCb_Total/model_vars.root']

weights0=[]
for file in files:
  components=(file.split('/')[-2]).split('_')
  components=components[:-2]   #extract the sub mode from the file name, remove 'LHCb_Total'
  weight=BF['D0'][components[0]]
  weights0.append(weight)

sum0=sum(weights0)
for i in range(len(weights0)):
  weights0[i]=weights0[i]/sum0   #define the weight with regard to the sum (the proportion)

DF=root_pandas.read_root(files[0],columns=columns,key='DecayTree')
DF=DF.sample(n=int(4000000*weights0[0]*frac_fit['Bu2DstD0K0']))

print 'OK'

for i in range(1,len(files)):
  df=root_pandas.read_root(files[i],columns=columns,key='DecayTree')
  df=df.sample(n=int(4000000*weights0[i]*frac_fit['Bu2DstD0K0']))
  DF=pd.concat([DF, df], ignore_index=True)
  print 'OK'


  
#Bu2DstDst0K0
files1=['/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bu2DstDst0K0/d0gamma_k3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bu2DstDst0K0/d0gamma_k3pipi0_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bu2DstDst0K0/d0pi0_k3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bu2DstDst0K0/d0pi0_k3pipi0_LHCb_Total/model_vars.root']

weights1=[]
for file in files1:
  components=(file.split('/')[-2]).split('_')
  components=components[:-2]   #extract the sub mode from the file name
  weight=BF['Dst0'][components[0]] * BF['D0'][components[1]]
  weights1.append(weight)

sum1=sum(weights1)
for i in range(len(weights1)):
  weights1[i]=weights1[i]/sum1   #define the weight with regard to the sum (the proportion)

for i in range(len(files1)):
  df=root_pandas.read_root(files1[i],columns=columns,key='DecayTree')
  df=df.sample(n=int(4000000*weights1[i]*frac_fit['Bu2DstDst0K0']),random_state=6000)
  DF=pd.concat([DF, df], ignore_index=True)
  print 'OK'


#Bd2DstD0K
files2=['/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bd2DstD0K/k3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bd2DstD0K/k3pipi0_LHCb_Total/model_vars.root']

weights2=[]
for file in files2:
  components=(file.split('/')[-2]).split('_')
  components=components[:-2]   #extract the sub mode from the file name
  weight=BF['D0'][components[0]]
  weights2.append(weight)
sum2=sum(weights2)
for i in range(len(weights2)):
  weights2[i]=weights2[i]/sum2   #define the weight with regard to the sum (the proportion)

for i in range(len(files2)):
  df=root_pandas.read_root(files2[i],columns=columns,key='DecayTree')
  df=df.sample(n=int(4000000*weights2[i]*frac_fit['Bd2DstD0K']),random_state=6000)
  DF=pd.concat([DF, df], ignore_index=True)
  print 'OK'




#Bd2DstDst0K
files3=['/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bd2DstDst0K/d0gamma_k3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bd2DstDst0K/d0gamma_k3pipi0_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bd2DstDst0K/d0pi0_k3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bd2DstDst0K/d0pi0_k3pipi0_LHCb_Total/model_vars.root']


weights3=[]
for file in files3:
  components=(file.split('/')[-2]).split('_')
  components=components[:-2]   #extract the sub mode from the file name
  weight=BF['Dst0'][components[0]] * BF['D0'][components[1]]
  weights3.append(weight)
  
sum3=sum(weights3)
for i in range(len(weights3)):
  weights3[i]=weights3[i]/sum3   #define the weight with regard to the sum (the proportion)

for i in range(len(files3)):
  df=root_pandas.read_root(files3[i],columns=columns,key='DecayTree')
  df=df.sample(n=int(4000000*weights3[i]*frac_fit['Bd2DstDst0K']),random_state=6000)
  DF=pd.concat([DF, df], ignore_index=True)
  print 'OK'





#Bd2DstDstK0
files4=['/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bd2DstDstK0/k3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bd2DstDstK0/k3pipi0_LHCb_Total/model_vars.root']
  
weights4=[]
for file in files4:
  components=(file.split('/')[-2]).split('_')
  components=components[:-2]   #extract the sub mode from the file name
  weight=BF['D0'][components[0]]
  weights4.append(weight)
sum4=sum(weights4)
for i in range(len(weights4)):
  weights4[i]=weights4[i]/sum4   #define the weight with regard to the sum (the proportion)

for i in range(len(files4)):
  df=root_pandas.read_root(files4[i],columns=columns,key='DecayTree')
  df=df.sample(n=int(4000000*weights4[i]*frac_fit['Bd2DstDstK0']),random_state=6000)
  DF=pd.concat([DF, df], ignore_index=True)
  print 'OK'


print len(DF)

DF["hamweight_SM"]=1.

DF["hamweight_T1"]=1.
DF["hamweight_T2"]=1.


#DF.to_root('D0.root', key='DecayTree')





  
"""
                                              
                                              ### HISTOGRAMS ###
ranges=[[-1.,1.],[-np.pi,np.pi],[0.,6.],[0.,13.],[-1.,1.]]
filenames=['costheta_D_reco','chi_reco','Tau_life_reco','q2_reco','costheta_L_reco']
titles=[r'cos$(\theta_D)$',r'$\chi$',r'$\tau$ life',r'$q^2$',r'cos$(\theta_L)$']
binnumber=100
                                              
DF.to_root('D0.root', key='DecayTree')

for i in range(5):
  plt.hist(DF[filenames[i]][~np.isnan(DF[filenames[i]])],histtype='step',bins=binnumber,range=ranges[i])
  plt.ylim(bottom=0)  
  plt.title(titles[i]+r'  (B $\rightarrow$ $D^*$ $D^0$ X)')
  plt.savefig(filenames[i]+'_total.pdf')
  plt.close()
"""
