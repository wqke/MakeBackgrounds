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

frac['signal']=1296/1296.
frac['feed']=0.11*1296/1296.
frac['B2DstDsX']=6835/1296.
frac['B2DstD0X']=445/1296.
frac['B2DstDplusX']=0.245*6835/1296.
frac['prompt']=424/1296.

#Neglect at this stage prompt and feed-down (small fraction)
#The signal is represented by B2Dsttaunu tau->3pinu purely
names=['B2DstD0X','prompt','B2DstDplusX','B2DstDsX','feed','signal']
#labels={label1: r'B$\rightarrow D^* D^0 X$', label2: 'prompt',label3: r'B$\rightarrow D^* D^+ X$',label4: r'B$\rightarrow D^* D_s X$',label5: 'feed-down', label6: 'signal'}
labels=[ r'B$\rightarrow D^* D^0 X$', 'prompt', r'B$\rightarrow D^* D^+ X$', r'B$\rightarrow D^* D_s X$', 'feed-down',  'signal']
colors=['#ef10c3','#cecacd','blue','orange','lightblue','red']
files={'signal':'/home/ke/graphs/signal.root',
      'B2DstD0X':'/home/ke/graphs/D0.root',
      'B2DstDplusX':'/home/ke/graphs/Dplus.root',
      'B2DstDsX':'/home/ke/graphs/Ds.root',
      'feed':'/home/ke/graphs/feed.root',
      'prompt':'/home/ke/graphs/prompt.root'}

M_B = 5.27963
M_Dst = 2.01026
q2_max = (M_B - M_Dst)**2
q2_min = 0.

q2_bins=3
angle_bins = 10  #int((float(200000)*(1.0/3.)/600)**(1.0/3.0))

histo=[]
histo1={}
histo2={}
histo3={}
for i in range(q2_bins):
  histo1['bin'+str(i)]=[]
  histo2['bin'+str(i)]=[]
  histo3['bin'+str(i)]=[]









##

var_bins = {"costheta_D_reco" : angle_bins,
        "costheta_L_reco"  : angle_bins,
        "chi_reco" : angle_bins,
        "q2_reco"  : q2_bins
        }

sample=root_pandas.read_root(files[names[5]],columns=['q2_reco','costheta_L_reco','costheta_D_reco','chi_reco'])    
sample = sample.query("q2_reco>=0. and costheta_L_reco>=-1. and costheta_L_reco<1. and costheta_D_reco>=-1. and costheta_D_reco<1. and chi_reco>=%s and chi_reco<%s" %(-np.pi,np.pi))
sample=sample.sample(n=int(200000*frac[names[5]]))
qc_bin_vals={}
for b in ["q2_reco"]:
  print b
  qc = pd.qcut(sample[b], q=q2_bins, precision=5)
  qc_bins = qc.unique()
  qc_bin_vals[b] = []
  for i in range(0,q2_bins):
    qc_bin_vals[b].append(qc_bins[i].left)
    qc_bin_vals[b].append(qc_bins[i].right)
    #Retain unique values then sort
    qc_bin_vals[b] = list(set(qc_bin_vals[b]))
    qc_bin_vals[b].sort()
    print qc_bin_vals[b]    

      
binning = []
var_type= 'reco'
for b in ["costheta_D_%s" % var_type,"costheta_L_%s" % var_type,"chi_%s" % var_type]:
  for i in range(0,var_bins["q2_%s" % var_type]):
    print "%s %s" % (b,i)
    bin_sample_temp = sample.query("q2_%s > %s and q2_%s <= %s" % (var_type,qc_bin_vals["q2_%s" % var_type][i],var_type,qc_bin_vals["q2_%s" % var_type][i+1]))
    qc = pd.qcut(bin_sample_temp[b], q=var_bins[b], precision=5)
    qc_bins = qc.unique()
    qc_bin_vals["%s_%s" % (b,i)] = []
    for j in range(0,var_bins[b]):
      qc_bin_vals["%s_%s" % (b,i)].append(qc_bins[j].left)
      qc_bin_vals["%s_%s" % (b,i)].append(qc_bins[j].right)
    #Retain unique values then sort
    qc_bin_vals["%s_%s" % (b,i)] = list(set(qc_bin_vals["%s_%s" % (b,i)]))
    qc_bin_vals["%s_%s" % (b,i)].sort()
    print qc_bin_vals["%s_%s" % (b,i)]
 
for i in range(0,var_bins["q2_%s" % var_type]):
  binning.append((qc_bin_vals["costheta_D_%s_%s" % (var_type,i)],qc_bin_vals["costheta_L_%s_%s" % (var_type,i)],qc_bin_vals["chi_%s_%s" % (var_type,i)]))

print 'BINNING =',binning 


### select sample
for i in range(5):
  df=root_pandas.read_root(files[names[i]],columns=['q2_reco','costheta_L_reco','costheta_D_reco','chi_reco'])
  df = df.query("q2_reco>=0. and costheta_L_reco>=-1. and costheta_L_reco<1. and costheta_D_reco>=-1. and costheta_D_reco<1. and chi_reco>=%s and chi_reco<%s" %(-np.pi,np.pi))
  df=df.sample(n=int(200000*frac[names[i]]))
  histo.append(df['q2_reco'][~np.isnan(df['q2_reco'])])
  for j in range(0,var_bins["q2_%s" % var_type]):
    dg=df.query("q2_reco>=%s and q2_reco<%s" %(qc_bin_vals["q2_reco"][j],qc_bin_vals["q2_reco"][j+1]))
    print names[i]+' in bin '+str(j)+' : ',len(dg)
    costhetal=dg['costheta_L_reco'][~np.isnan(dg['costheta_L_reco'])]
    costhetad=dg['costheta_D_reco'][~np.isnan(dg['costheta_D_reco'])]
    chi=dg['chi_reco'][~np.isnan(dg['chi_reco'])]
    histo1['bin'+str(j)].append(costhetal)
    histo2['bin'+str(j)].append(costhetad)
    histo3['bin'+str(j)].append(chi)

bin_sample=root_pandas.read_root(files[names[5]],columns=['q2_reco','costheta_L_reco','costheta_D_reco','chi_reco'])    
bin_sample = bin_sample.query("q2_reco>=0. and costheta_L_reco>=-1. and costheta_L_reco<1. and costheta_D_reco>=-1. and costheta_D_reco<1. and chi_reco>=%s and chi_reco<%s" %(-np.pi,np.pi))
bin_sample=bin_sample.sample(n=int(200000*frac[names[5]]))
histo.append(bin_sample['q2_reco'][~np.isnan(bin_sample['q2_reco'])])
for j in range(0,var_bins["q2_%s" % var_type]):
  dg=bin_sample.query("q2_reco>=%s and q2_reco<%s" %(qc_bin_vals["q2_reco"][j],qc_bin_vals["q2_reco"][j+1]))
  print names[5]+' in bin '+str(j)+' : ',len(dg)
  costhetal=dg['costheta_L_reco'][~np.isnan(dg['costheta_L_reco'])]
  costhetad=dg['costheta_D_reco'][~np.isnan(dg['costheta_D_reco'])]
  chi=dg['chi_reco'][~np.isnan(dg['chi_reco'])]
  histo1['bin'+str(j)].append(costhetal)
  histo2['bin'+str(j)].append(costhetad)
  histo3['bin'+str(j)].append(chi)



####MAKE HISTOGRAMS

plt.hist(histo,stacked=True,bins=qc_bin_vals["q2_reco"] ,label=labels,color=colors)
plt.legend()
plt.title(r'$q^2$ reco')
plt.savefig('QTOTAL.pdf')
plt.close()

for i in range(q2_bins):
  plt.hist(histo1['bin'+str(i)],stacked=True,bins=qc_bin_vals["costheta_L_reco_"+str(i)],label=labels,color=colors)
  #print 'COSTHETAL bins events: ',heights[5][0]-heights[4][0],heights[5][1]-heights[4][1],heights[5][2]-heights[4][2]
  plt.legend()
  plt.title(r'$cos(\theta_L)$ reco - bin '+str(i))
  plt.savefig('COSTHETAL_bin'+str(i)+'.pdf')
  plt.close()

  plt.hist(histo2['bin'+str(i)],stacked=True,bins=qc_bin_vals["costheta_D_reco_"+str(i)],label=labels,color=colors)
  plt.legend()
  plt.title(r'$cos(\theta_D)$ reco - bin '+str(i))
  plt.savefig('COSTHETAD_bin'+str(i)+'.pdf')
  plt.close()

  plt.hist(histo3['bin'+str(i)],stacked=True,bins=qc_bin_vals["chi_reco_"+str(i)],label=labels,color=colors)
  plt.legend()
  plt.title(r'$\chi$ reco - bin '+str(i))
  plt.savefig('CHI_bin'+str(i)+'.pdf')
  plt.close()



