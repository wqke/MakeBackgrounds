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

#D* Ds X
#Fractions defined relative to D* Ds*(=1)

#Simulation values
frac_sim={}
frac_sim['Bd2DstDs']=0.54
frac_sim['Bd2DstDsst']=1.
frac_sim['Bd2DstDs1']=
frac_sim['Bu2DststDs']=
frac_sim['Bu2DststDsst']=
frac_sim['Bu2DststDs1']=0.39

#Fit values
frac_fit={}
frac_fit['Bd2DstDs']=0.594
frac_fit['Bd2DstDsst']=1.
frac_fit['Bd2DstDs1']=
frac_fit['Bu2DststDs']=
frac_fit['Bu2DststDsst']=
frac_fit['Bu2DststDs1']=0.365


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
BF['Ds1']['dsstpi0']=0.48

###Ds* decays
#Ds* -> Ds gamma
BF['Dsst']['dsgamma']=0.299
#Ds* -> Ds pi0
BF['Dsst']['dspi0']=0.058

###D1(2420)0 decays
#D1(2420) -> D*- pi+
BF['D10']['dstpiplus']='seen'

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
BF['etap']['eta2pi']=0.426
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
BF['Dsplus']['etapi_3pi']=BF['Dsplus']['etapi']*BF['eta']['3pi']
#Ds+->omega pi+
BF['Dsplus']['omegapi']=2.4e-3
#Ds+->(omega->3pi) pi+
BF['Dsplus']['omegapi_3pi']=BF['Dsplus']['omegapi']*BF['omega']['3pi']
#Ds+->eta rho
BF['Dsplus']['etarho']=0.089
#Ds+->(eta->3pi) (rho->2pi)
BF['Dsplus']['etarho_5pi']=BF['Dsplus']['etarho']*BF['eta']['3pi']*BF['rho0']['2pi']
#Ds+->omega rho
BF['Dsplus']['omegarho']=
#Ds+->rho0(pi+ pi-) rho0 (pi+ pi-) pi+
BF['Dsplus']['5pi']=


#Ds+->omega pi+pi+pi-
BF['Dsplus']['omega3pi']=0.016
#Ds+->eta' pi+
BF['Dsplus']['etappi']=0.0394
#Ds+->eta' rho+
BF['Dsplus']['etaprho']=0.058


############PLOT THE TOTAL HISTOGRAMS############

#Bd2DstDs1
files=['/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsgamma_5pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsgamma_etapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsgamma_etappi_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsgamma_etappi_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsgamma_etaprho_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsgamma_etaprho_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsgamma_etarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsgamma_omega3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsgamma_omegapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsgamma_omegarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsstpi0_dsgamma_5pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsstpi0_dsgamma_etappi_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsstpi0_dsgamma_etappi_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsstpi0_dsgamma_etaprho_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsstpi0_dsgamma_etaprho_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsstpi0_dsgamma_etarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsstpi0_dsgamma_omega3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsstpi0_dsgamma_omegapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs1/dsstpi0_dsgamma_omegarho_LHCb_Total/model_vars.root']

#Bd2DstDs
files1=['/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs/5pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs/etapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs/etappi_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs/etappi_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs/etaprho_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs/etaprho_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs/etarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs/omega3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs/omegapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDs/omegarho_LHCb_Total/model_vars.root']

#Bd2DstDsst
files2=['/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dsgamma_5pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dsgamma_etapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dsgamma_etappi_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dsgamma_etappi_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dsgamma_etaprho_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dsgamma_etaprho_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dsgamma_etarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dsgamma_omega3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dsgamma_omegapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dsgamma_omegarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dspi0_5pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dspi0_etapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dspi0_etappi_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dspi0_etappi_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dspi0_etaprho_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dspi0_etaprho_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dspi0_etarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dspi0_omega3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dspi0_omegapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstDsst/dspi0_omegarho_LHCb_Total/model_vars.root']

#Bu2DststDs
files3=['/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs/5pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs/etapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs/etappi_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs/etappi_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs/etaprho_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs/etaprho_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs/etarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs/omega3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs/omegapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs/omegarho_LHCb_Total/model_vars.root']


#Bu2DststDsst
files4=['/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dsgamma_5pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dsgamma_etapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dsgamma_etappi_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dsgamma_etappi_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dsgamma_etaprho_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dsgamma_etaprho_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dsgamma_etarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dsgamma_omega3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dsgamma_omegapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dsgamma_omegarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dspi0_5pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dspi0_etapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dspi0_etappi_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dspi0_etappi_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dspi0_etaprho_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dspi0_etaprho_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dspi0_etarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dspi0_omega3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dspi0_omegapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDsst/dspi0_omegarho_LHCb_Total/model_vars.root']


#Bu2DststDs1
files5=['/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsgamma_5pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsgamma_etapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsgamma_etappi_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsgamma_etappi_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsgamma_etaprho_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsgamma_etaprho_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsgamma_etarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsgamma_omega3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsgamma_omegapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsgamma_omegarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsstpi0_dsgamma_5pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsstpi0_dsgamma_etapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsstpi0_dsgamma_etappi_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsstpi0_dsgamma_etappi_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsstpi0_dsgamma_etaprho_etapipi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsstpi0_dsgamma_etaprho_rhogamma_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsstpi0_dsgamma_etarho_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsstpi0_dsgamma_omega3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsstpi0_dsgamma_omegapi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bu2DststDs1/dsstpi0_dsgamma_omegarho_LHCb_Total/model_vars.root']


#27:-70

df=root_pandas.read_root(files1,columns=['q2_reco','costheta_L_reco','costheta_D_reco','chi_reco'],key='DecayTree')

df_merged = pd.concat([df, df2], ignore_index=True)

# weights
df_weights = np.ones_like(df.values) / len(df_merged)
df2_weights = np.ones_like(df2.values) / len(df_merged)
df_merged_weights = np.ones_like(df_merged.values) / len(df_merged)

plt_range = (df_merged.values.min(), df_merged.values.max())
fig, ax = plt.subplots()
ax.hist(df.values, bins=100, weights=df_weights, color='black', histtype='step', label='df', range=plt_range)
ax.hist(df2.values, bins=100, weights=df2_weights, color='green', histtype='step', label='df2', range=plt_range)
ax.hist(df_merged.values, bins=100, weights=df_merged_weights, color='red', histtype='step', label='Combined', range=plt_range)

ax.margins(0.05)
ax.set_ylim(bottom=0)
#ax.set_xlim([0, 1000])
plt.legend(loc='upper right')
# plt.savefig('output.png')




