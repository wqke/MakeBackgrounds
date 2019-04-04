import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from numpy import *
import tensorflow as tf
import sys, os
import numpy as np
import math
from math import cos,sin,pi
import root_pandas
from pandas import * 
from tensorflow.python.client import timeline
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

###Ds1(2460) decays
#Ds1(2460)->Ds+ gamma
BF['Ds1_dsgamma']=0.18
BF['Ds1_dsstpi0']=0.48

###Ds* decays
#Ds* -> Ds gamma
BF['Dsst_dsgamma']=0.299
#Ds* -> Ds pi0
BF['Dsst_dspi0']=0.058


###D1(2420)0 decays
#D1(2420) -> D*- pi+
BF['D10_dstpiplus']='seen'


###D0 decays
#D0 -> K- & anything
BF['D0_KminANY']=0.547
#

###D+ decays
#D+ -> Ks0 3pi
BF['Dplus_Ks3pi']=0.0297
#D+ -> pi+pi+pi-pi0
BF['Dplus_3pipi0']=0.0111

###D*0 decays
#D*0 -> Ks0 3pi
BF['Dst0_d0pi0']=0.647
#D*0 -> pi+pi+pi-pi0
BF['Dplus_d0gamma']=0.353


###eta decay : eta-> pi+ pi- pi0
BF['eta_3pi']=0.2292
###eta' decays : eta'->eta pi+ pi- 
BF['etap_eta2pi']=0.426
#eta'->rho0 gamma 
BF['etap_rhogamma']=0.289
###rho0 decay : rho0 ->pi+pi-
BF['rho0_2pi']=1.
###rho+ decay : rho+ ->pi+pi0
BF['rhoplus_2pi']=1.
###omega decay : omega ->pi+pi-pi0
BF['omega_3pi']=0.892

###Ds+ decays
#Ds+->eta pi+
BF['Dsplus_etapi']=0.017
#Ds+->(eta->3pi) pi+
BF['Dsplus_etapi_3pi']=BF['Dsplus_etapi']*BF['eta_3pi']
#Ds+->omega pi+
BF['Dsplus_omegapi']=2.4e-3
#Ds+->(omega->3pi) pi+
BF['Dsplus_omegapi_3pi']=BF['Dsplus_omegapi']*BF['omega_3pi']
#Ds+->eta rho
BF['Dsplus_etarho']=0.089
#Ds+->(eta->3pi) (rho->2pi)
BF['Dsplus_etarho_5pi']=BF['Dsplus_etarho']*BF['eta_3pi']*BF['rho0_2pi']
#Ds+->omega rho
BF['Dsplus_omegarho']=
#Ds+->rho0(pi+ pi-) rho0 (pi+ pi-) pi+
BF['Dsplus_5pi']=


#Ds+->omega pi+pi+pi-
BF['Dsplus_omega3pi']=0.016
#Ds+->eta' pi+
BF['Dsplus_etappi']=0.0394
#Ds+->eta' rho+
BF['Dsplus_etaprho']=0.058







