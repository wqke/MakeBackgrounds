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


#Fractions defined with regard to the biggest one. B+ and B0 are assumed to have the same amount (so we directly use the branching fractions)
#Fit values
frac_fit={}
frac_fit['2420']=3.03e-3/3.03e-3
frac_fit['2460']=1.01e-3/3.03e-3

frac_fiterr={}
frac_fiterr['2420']=0.2e-3/3.03e-3
frac_fiterr['2460']=0.24e-3/3.03e-3


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
BF['tau']={}

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
BF['Dplus']['ks3pi']=0.0297
#D+ -> pi+pi+pi-pi0
BF['Dplus']['3pipi0']=0.0111

###D0 decays
#D0 ->K-pi+pi+pi-
BF['D0']['k3pi']=0.0811
#D0 ->K-pi+pi+pi-pi0
BF['D0']['kpipi0']=0.042

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
###tau decays 
BF['tau']['3pi']=9.31e-2
BF['tau']['3pipi0']=4.62e-2




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

def inter(lst1, lst2):
    return list(set(lst1) & set(lst2))


columns=['PiP_2_Tau_nu_M2', 'B0_Ori_x_TRUE', 'B_nu_E_TRUE', 'PiP_2_PiM_1_Tau_nu_M2', 'D0_PX_TRUE', 'Dst_PY_TRUE', 'Dst_Ori_z_TRUE', 'B0_End_x', 'B0_End_z', 'Tau_nu_E_TRUE', 'D0_K_E_TRUE', 'B0_End_x_TRUE', 'Dstst_End_z_TRUE', 'Tau_FD_TRUE', 'B0_PX_reco', 'nEvent', 'Tau_P_TRUE', 'PiP_2_PX_TRUE', 'D0_E', 'B0_PT_TRUE', 'D0_M', 'D0_P', 'D0_FD_TRUE', 'Dstst_Tau_M2', 'Dstst_PT_TRUE', 'Dst_FD_TRUE', 'B0_eta_TRUE', 'Dstst_FD_y_TRUE', 'B0_nu_PZ_reco', 'PiP_1_PZ_TRUE', 'B0_PY_TRUE', 'B0_Ori_y_TRUE', 'B0_nu_PX_reco', 'Dst_PZ_TRUE', 'B0_M_TRUE', 'D0_FD', 'D0_Pi_P', 'Dstst_End_x', 'Dstst_End_y', 'Dstst_End_z', 'PiP_1_PiP_2_PiM_1_M2_TRUE', 'D0_Pi_E', 'Dstst_B_nu_M2_TRUE', 'Tau_Ori_y_TRUE', 'Tau_PX_reco', 'Dstst_Ori_z_TRUE', 'Dst_Pi_PZ_TRUE', 'Dst_Pi_E_TRUE', 'Tau_End_z_TRUE', 'Tau_B_nu_M2_TRUE', 'Dstst_Pi_E_TRUE', 'B0_nu_P_reco', 'costheta_L_reco', 'Dstst_FD_z_TRUE', 'D0_eta_TRUE', 'B0_End_y', 'D0_Pi_P_TRUE', 'Tau_gamma', 'Tau_FD_z_TRUE', 'Dstst_Tau_M2_TRUE', 'PiM_1_E_TRUE', 'Tau_beta', 'Dstst_End_x_TRUE', 'PiP_1_PiP_2_M2_TRUE', 'D0_PY_TRUE', 'Tau_Pi1_PZ', 'Tau_Pi1_PX', 'Tau_Pi1_PY', 'Dstst_Pi_PY_TRUE', 'Dst_M', 'Tau_FD_x_TRUE', 'Dst_E', 'Dst_P', 'D0_gamma', 'Tau_Pi3_PX', 'Tau_Pi3_PY', 'Tau_Pi3_PZ', 'Tau_nu_PY_reco', 'PiP_2_PiM_1_M2', 'Tau_FD_z', 'Tau_FD_x', 'B_nu_PY_TRUE', 'Tau_PY_reco', 'Dstst_P_TRUE', 'Dstst_FD_TRUE', 'Tau_E_reco', 'PiP_1_PiP_2_Tau_nu_M2', 'Tau_nu_P_reco', 'Dstst_Pi_PX_TRUE', 'Tau_End_x_TRUE', 'B0_End_z_TRUE', 'PiM_1_Tau_nu_M2', 'Dstst_P', 'Tau_nu_PT_TRUE', 'B_nu_PX_TRUE', 'B_nu_PZ_TRUE', 'PiP_1_P_TRUE', 'D0_Ori_y_TRUE', 'Dst_gamma_TRUE', 'Tau_nu_PX_reco', 'Dstst_FD_z', 'Dst_gamma', 'Tau_Pi3_E', 'Dstst_FD_y', 'B_nu_PT_TRUE', 'Tau_Pi1_E', 'D0_K_P_TRUE', 'PiP_2_PT', 'Tau_Pi2_PZ', 'Tau_Pi2_PY', 'Tau_Pi2_PX', 'PiM_1_PT_TRUE', 'Tau_M', 'Tau_FD_y', 'Tau_E', 'PiP_1_PT', 'Dstst_B_nu_M2', 'D0_Pi_PT', 'D0_Pi_PY', 'D0_Pi_PX', 'D0_Pi_PZ', 'PiP_1_PZ', 'PiP_1_Tau_nu_M2', 'Tau_P', 'PiP_1_PY', 'Dstst_gamma', 'PiP_2_PZ', 'Tau_PY', 'PiP_2_PX', 'D0_K_PT_TRUE', 'DstTau_P_reco', 'Tau_Ori_z', 'Tau_Ori_y', 'Tau_Ori_x', 'Tau_eta_TRUE', 'D0_Ori_x', 'D0_Ori_y', 'D0_Ori_z', 'DstTau_PX_reco', 'Tau_nu_PY_TRUE', '3pi_PZ', '3pi_PY', '3pi_PX', 'D0_Ori_z_TRUE', 'PiP_2_P_TRUE', 'D0_FD_z_TRUE', 'PiM_1_PY', 'PiM_1_PX', 'DstTau_PZ_reco', 'PiM_1_PZ', 'PiM_1_PT', 'Tau_nu_P', 'Dstst_Ori_z', 'Dstst_Ori_y', 'Dstst_Ori_x', 'Tau_nu_E', 'Tau_PT_TRUE', 'PiP_2_PZ_TRUE', 'B0_PZ_reco', 'PiM_1_PY_TRUE', 'Dst_beta_TRUE', 'D0_K_E', 'D0_K_PZ_TRUE', 'Dstst_E_TRUE', 'PiP_1_PT_TRUE', 'B0_PZ_TRUE', 'D0_K_P', 'D0_beta', 'Dst_Pi_P', 'B0_FD_x_TRUE', 'Dstst_PZ_TRUE', 'D0_E_TRUE', 'Dst_Pi_E', 'Dst_beta', 'B0_FD', '3pi_P', 'Dst_FD', 'D0_M_TRUE', 'Dstst_E', 'Dstst_PZ', 'Dstst_PY', 'Dstst_PX', 'Dst_Ori_x', 'Dst_Ori_y', 'Dst_Ori_z', 'Dstst_gamma_TRUE', 'Dstst_PT', 'Tau_angle_max', 'D0_FD_y_TRUE', 'Tau_m12', 'Tau_m13', 'Dst_FD_x_TRUE', 'Tau_gamma_TRUE', '3pi_E', 'Tau_FD_y_TRUE', 'B0_E_TRUE', 'B0_Ori_z_TRUE', 'Dstst_Pi_P_TRUE', 'theta_L_reco', 'Tau_nu_PX_TRUE', 'DstTau_M_reco', 'B0_PX_TRUE', 'D0_End_y_TRUE', 'Dst_PT_TRUE', 'D0_Pi_E_TRUE', 'Dst_eta_TRUE', 'MCorr', 'PiP_2_E', 'B0_nu_PY_reco', 'PiP_2_P', 'Tau_PT', 'PiP_1_PX_TRUE', 'Dst_Pi_PT', 'Dst_Pi_PY', 'Dst_Pi_PX', 'Dst_Pi_PZ', 'Tau_PZ', 'PiP_2_PY', 'Tau_PX', 'Dstst_eta', 'B0_beta_TRUE', 'Dst_End_y_TRUE', 'Dst_Pi_P_TRUE', 'Tau_nu_PZ_reco', 'B0_FD_y', 'Dstst_PX_TRUE', 'Dstst_Ori_y_TRUE', 'Tau_beta_TRUE', 'Dst_PX', 'Tau_FD', 'Dst_PT', 'Dst_PZ', 'PiP_1_PiP_2_Tau_nu_M2_TRUE', 'm2_miss_reco', 'Dst_PY', 'B0_gamma', 'D0_beta_TRUE', 'Dst_E_TRUE', 'costheta_D_reco', 'D0_K_PZ', 'D0_K_PX', 'D0_K_PY', 'Dst_PX_TRUE', 'D0_K_PT', 'D0_End_z', 'D0_End_y', 'D0_End_x', 'Dst_FD_z', 'Dst_FD_y', 'Dst_FD_x', 'B0_P_reco', 'B0_E_reco', 'Tau_m23', 'Dst_P_TRUE', 'B_nu_P', 'Dst_M_TRUE', 'PiP_1_E', 'PiP_1_Tau_nu_M2_TRUE', 'B_nu_E', 'Dst_End_z_TRUE', 'PiP_2_PY_TRUE', 'q2_reco', 'PiP_1_P', 'B0_P_TRUE', 'B0_nu_E_reco', 'Dstst_beta', 'PiP_1_PiP_2_PiM_1_M2', 'Tau_nu_E_reco', 'B0_PT', 'B0_PX', 'B0_PY', 'B0_PZ', 'Dstst_FD_x_TRUE', 'Dst_End_x_TRUE', 'Dstst_Pi_PT_TRUE', 'B_nu_P_TRUE', 'Dst_Pi_PY_TRUE', 'Dst_Ori_x_TRUE', 'B0_eta', 'PiP_1_PiM_1_M2_TRUE', 'Dstst_FD_x', 'Tau_nu_P_TRUE', 'D0_eta', 'D0_End_z_TRUE', 'D0_End_x_TRUE', 'Tau_eta', 'Dst_Pi_PT_TRUE', 'Tau_Ori_z_TRUE', 'Dstst_Pi_PT', 'Dstst_Pi_PZ', 'Dstst_FD', 'Dstst_Pi_PX', 'Dstst_Pi_PY', 'PiP_1_PiP_2_M2', 'D0_Pi_PY_TRUE', '3pi_M', 'Tau_PZ_TRUE', 'Dstst_M', 'PiP_2_PT_TRUE', 'PiP_2_Tau_nu_M2_TRUE', 'PiP_1_PiM_1_Tau_nu_M2_TRUE', 'Dst_Pi_PX_TRUE', 'Dst_eta', 'B0_Ori_z', 'B0_Ori_x', 'B0_Ori_y', 'PiP_1_E_TRUE', 'PiM_1_P_TRUE', 'D0_Pi_PX_TRUE', 'D0_Pi_PZ_TRUE', 'Tau_PX_TRUE', 'Dstst_Pi_E', 'D0_Pi_PT_TRUE', 'Dstst_Pi_P', 'Dstst_Pi_PZ_TRUE', 'DstTau_PY_reco', 'PiP_2_PiM_1_M2_TRUE', 'B0_FD_y_TRUE', 'B_nu_PZ', 'B_nu_PX', 'B_nu_PY', 'Tau_PZ_reco', 'B_nu_PT', 'PiP_2_PiM_1_Tau_nu_M2_TRUE', 'PiM_1_PX_TRUE', 'PiM_1_P', 'D0_K_PY_TRUE', 'chi_reco', 'Tau_life_reco', 'Dstst_PY_TRUE', 'Tau_nu_PZ_TRUE', 'PiM_1_E', 'Dst_FD_z_TRUE', 'PiP_1_PiM_1_M2', 'PiM_1_PZ_TRUE', 'B0_End_y_TRUE', 'D0_PT_TRUE', 'B0_beta', 'Tau_End_y_TRUE', 'D0_Ori_x_TRUE', 'theta_D_reco', 'PiP_1_PX', 'D0_FD_z', 'D0_FD_y', 'D0_FD_x', 'D0_FD_x_TRUE', 'D0_K_PX_TRUE', 'B0_angle_max', 'Dstst_End_y_TRUE', 'B0_P', 'Dstst_eta_TRUE', 'Tau_M_TRUE', 'B0_E', 'Tau_Pi2_E', 'B0_M', 'D0_P_TRUE', 'Tau_End_x', 'Tau_End_y', 'Tau_End_z', 'PiP_1_PY_TRUE', 'B0_PY_reco', 'B0_FD_x', 'B0_FD_z', 'D0_PZ_TRUE', 'coschi_reco', 'Dst_Ori_y_TRUE', 'Dstst_Ori_x_TRUE', 'Tau_nu_PT', 'Tau_nu_PX', 'Tau_nu_PY', 'Tau_nu_PZ', 'Tau_PY_TRUE', 'PiP_1_PiM_1_Tau_nu_M2', 'D0_PT', 'Dstst_beta_TRUE', 'Tau_B_nu_M2', 'PiM_1_Tau_nu_M2_TRUE', 'Tau_P_reco', 'D0_PZ', 'D0_PX', 'D0_PY', 'B0_gamma_TRUE', 'Dst_FD_y_TRUE', 'Dst_End_z', 'Dst_End_y', 'Dst_End_x', 'Tau_Ori_x_TRUE', 'B0_FD_z_TRUE', 'Tau_E_TRUE', 'B0_FD_TRUE', 'D0_gamma_TRUE', 'PiP_2_E_TRUE', 'DstTau_E_reco', 'Dstst_M_TRUE']

############PLOT THE TOTAL HISTOGRAMS############

#Bu2DststTauNu
files=['/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bu2DststTauNu/2420_3pi_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bu2DststTauNu/2420_3pipi0_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bu2DststTauNu/2460_3pipi0_LHCb_Total/model_vars.root',
'/data/lhcb/users/hill/Bd2DstTauNu_Angular/RapidSim_tuples/Bu2DststTauNu/2460_3pi_LHCb_Total/model_vars.root']

weights0=[]
for file in files:
  components=(file.split('/')[-2]).split('_')
  components=components[:-2]   #extract the sub mode from the file name, remove 'LHCb_Total'
  weight=frac_fit[components[0]] * BF['tau'][components[1]]
  weights0.append(weight)

sum0=sum(weights0)
for i in range(len(weights0)):
  weights0[i]=weights0[i]/sum0   #define the weight with regard to the sum (the proportion)



DF=root_pandas.read_root(files[0],columns=columns,key='DecayTree')
Z=len(DF)
#columns=DF.keys()
DF=DF.sample(n=int(Z*weights0[0]),random_state=int(Z/1000))

for i in range(1,len(files)):
  df=root_pandas.read_root(files[i],columns=columns,key='DecayTree')
  df=df.sample(n=int(Z*weights0[i]),random_state=int(Z/1000))
  DF=pd.concat([DF, df], ignore_index=True)




DF["hamweight_SM"]=1.
DF["hamweight_T1"]=1.
DF["hamweight_T2"]=1.


DF.to_root('feed.root', key='DecayTree')


"""                         
                                              ### HISTOGRAMS ###
ranges=[[-1.,1.],[-np.pi,np.pi],[0.,6.],[0.,13.],[-1.,1.]]
filenames=['costheta_D_reco','chi_reco','Tau_life_reco','q2_reco','costheta_L_reco']
titles=[r'cos$(\theta_D)$',r'$\chi$',r'$\tau$ life',r'$q^2$',r'cos$(\theta_L)$']
binnumber=100
                                              
DF.to_root('feed.root', key='DecayTree')

for i in range(5):
  plt.hist(DF[filenames[i]][~np.isnan(DF[filenames[i]])],histtype='step',bins=binnumber,range=ranges[i])
  plt.ylim(bottom=0)  
  plt.title(titles[i]+'  (feed-down)')
  plt.savefig(filenames[i]+'_feed.pdf')
  plt.close()
"""	
