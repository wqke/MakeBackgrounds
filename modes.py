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
BF['dsgamma_5pi']=
BF['dsgamma_etapi']=
BF['dsgamma_etappi_etapipi']=
BF['dsgamma_etappi_rhogamma']=
BF['dsgamma_etaprho_etapipi']=
BF['dsgamma_etaprho_rhogamma']=
BF['dsgamma_etarho']=
BF['dsgamma_omega3pi']=
BF['dsgamma_omegapi']=
BF['dsgamma_omegarho']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=
BF['dsstpi0_dsgamma']=





