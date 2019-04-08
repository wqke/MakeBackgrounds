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

frac{'signal'}=1296/6835.
frac{'feed'}=0.11*1296/6835.
frac{'B2DstDsX'}=6835/6835.
frac{'B2DstD0X'}=445/6835.
frac{'B2DstDplusX'}=0.245*6835/6835.
frac{'prompt'}=424/6835.












