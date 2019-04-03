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



###Ds+ decays
#Ds+->eta & anything
BF['etaANY']=0.299
#Ds+->eta' & anything
BF['etapANY']=0.103
#Ds+->pi+ & anything
BF['piplusANY']=1.193
#Ds+->pi0 & anything
BF['pi0ANY']=1.234
#Ds+->pi- & anything
BF['piminANY']=0.432
#Ds+->omega & anything
BF['omegaANY']=0.061


#Ds+ -> rho rho 5pi   (pi+ &anything)
BF['5pi']=1.193
#Ds+->eta(->3pi)pi+  (eta &anything)
BF['etapi']=0.017  0.299
#Ds+->eta(something)pi+  (eta' &anything)
BF['etappi_etapipi']=0.103

BF['dsgamma_5pi']=BF['dsgamma']*BF['5pi']
BF['dsgamma_etapi']=BF['dsgamma']*BF['etapi']
BF['dsgamma_etappi_etapipi']=BF['dsgamma']*BF['etappi_etapipi']
BF['dsgamma_etappi_rhogamma']=BF['dsgamma']*BF['etappi_rhogamma']
BF['dsgamma_etaprho_etapipi']=BF['dsgamma']*
BF['dsgamma_etaprho_rhogamma']=BF['dsgamma']*
BF['dsgamma_etarho']=BF['dsgamma']*
BF['dsgamma_omega3pi']=BF['dsgamma']*
BF['dsgamma_omegapi']=BF['dsgamma']*
BF['dsgamma_omegarho']=BF['dsgamma']*
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=





