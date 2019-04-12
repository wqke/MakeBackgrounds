# Copyright 2017 CERN for the benefit of the LHCb collaboration
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ===============================================================================

import matplotlib
matplotlib.use('Agg')  #fix python2 tkinter problem

import tensorflow as tf
import numpy as np

import sys, os, math
sys.path.append("../../TensorFlowAnalysis")
os.environ["CUDA_VISIBLE_DEVICES"] = ""   # Do not use GPU

import TensorFlowAnalysis as tfa

from ROOT import TFile, TChain, TH3F
from root_numpy import root2array, rec2array, tree2array

import matplotlib.pyplot as plt
import rootplot.root2matplotlib as r2m
from scipy.stats import norm as sci_norm
from scipy.stats import sem as sem
import matplotlib.mlab as mlab
from root_pandas import to_root, read_root
from uncertainties import *
import pandas as pd
import random
import math


def MakeHistogram(phsp, sample, bins, weights = None, normed = False) : 
  hist = np.histogramdd(sample, bins = bins, range = phsp.Bounds(), weights = weights, normed = normed )
  return hist[0]  # Only return the histogram itself, not the bin boundaries

def MakeHistogram_1D(sample, bins, weights = None, normed = False, density = None) : 
  hist = np.histogram(sample, bins = bins, normed = normed, weights = weights, density = density)
  return hist[0]  # Only return the histogram itself, not the bin boundaries
  
def HistogramNorm(hist) : 
  return np.sum( hist )

def BinnedChi2(hist1, hist2, err) :
	return tf.reduce_sum( ((hist1 - hist2)/err)**2 )


if __name__ == "__main__" :

  #Read RapidSim signal sample for either 3pi mode or 3pipi0 mode
  mode = "Bd2DstTauNu"
  #3pi or 3pipi0
  sub_mode = sys.argv[1]
  #Geometry (all or LHCb)
  geom = sys.argv[2]
  #True or reco angles
  var_type = sys.argv[3]
  #Number of events to run on (in k) - 5, 10, 20, 40, 80
  num_sig = sys.argv[4]
  #Run a toy (Y/N)
  toy = sys.argv[5]
  
 #background fractions
  frac={}

  frac['signal']=1296/1296.
  frac['feed']=0.11*1296/1296.
  frac['B2DstDsX']=6835/1296.
  frac['B2DstD0X']=445/1296.
  frac['B2DstDplusX']=0.245*6835/1296.
  frac['prompt']=424/1296.
	
  M_B = 5.27963
  M_Dst = 2.01026
  
  q2_max = (M_B - M_Dst)**2
  q2_min = M_B - M_Dst
  
  if(var_type=="reco"):
    q2_min = 0.0
 
  #Number of angle bins depending on the signal yield, requiring roughly 50 events per bin
  q2_bins = 3
  angle_bins = int((float(num_sig)*1000.0*(1.0/q2_bins)/50)**(1.0/3.0))
 
  print "NUMBER OF BINS IN EACH ANGLE : %s" % angle_bins
  print "NUMBER OF BINS IN q2 : %s" % q2_bins

  #Binning scheme
  var_bins = {"costheta_D_%s" % var_type: angle_bins,
              "costheta_L_%s" % var_type: angle_bins,
              "chi_%s" % var_type: angle_bins,
              "q2_%s" % var_type: q2_bins
              }
	
  var_range = {"costheta_D_%s" % var_type: (-1.,1.),
               "costheta_L_%s" % var_type: (-1.,1.),
               "chi_%s" % var_type: (-math.pi,math.pi),
               "q2_%s" % var_type: (q2_min,q2_max)
               }
  
  var_titles = {"costheta_D_%s" % var_type: "$\\cos(\\theta_D)$",
                "costheta_L_%s" % var_type: "$\\cos(\\theta_L)$",
                "chi_%s" % var_type: "$\\chi$ [rad]",
                "q2_%s" % var_type: "$q^2$ [GeV$^2/c^4$]"
                }
  	
  # Four body angular phase space is described by 3 angles + q2 for background rejection
  phsp = tfa.RectangularPhaseSpace( ( var_range["costheta_D_%s" % var_type], var_range["costheta_L_%s" % var_type], var_range["chi_%s" % var_type]))

  # TF initialiser
  init = tf.global_variables_initializer()
  sess = tf.Session()
  sess.run(init)
  
  branch_names = ["costheta_D_%s" % var_type, "costheta_L_%s" % var_type, "chi_%s" % var_type, "q2_%s" % var_type]
  

  #Read RapidSim sample used to determine bins
  print "Loading tree"
  bin_file = "/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstTauNu/%s_%s_Total/model_vars_weights.root" % (sub_mode,geom)
   
  #Total rate factor multiplying the PDF
  Rate = tfa.FitParameter("Rate" , 1.0,-100,100)
  
  # Fit parameters
  #I1c = tfa.FitParameter("I1c" , 0.57, -1, 1, 0.01)
  #I1c = (4 - 6 I1s + I2c + 2I2s)/3
  I1s = tfa.FitParameter("I1s" , 0.393099, -100, 100)
  I2c = tfa.FitParameter("I2c" , -0.170286, -100, 100)
  I2s = tfa.FitParameter("I2s" , 0.046821, -100, 100)
  I6c = tfa.FitParameter("I6c" , 0.329175, -100, 100)
  I6s = tfa.FitParameter("I6s" , -0.248686, -100, 100)
  I3  = tfa.FitParameter("I3"  , -0.107168, -100, 100)
  I4  = tfa.FitParameter("I4"  , -0.111871, -100, 100)
  I5  = tfa.FitParameter("I5"  , 0.259955, -100, 100)
  I7  = tfa.FitParameter("I7"  , 0., -100, 100)
  I8  = tfa.FitParameter("I8"  , 0., -100, 100)
  I9  = tfa.FitParameter("I9"  , 0., -100, 100)

  #The floating fraction of Ds
  frac_Ds=tfa.FitParameter("frac_Ds", frac['B2DstDsX'] , -100, 100)

  coeffs = ["I1c","I1s","I2c","I2s","I6c","I6s","I3","I4","I5","I7","I8","I9"]
  	
  #File used to create templates (flat sample with unbinned weights)
  template_file = "/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstTauNu/3pi_%s_Total/model_vars_weights.root" % geom
  
  template_sample = read_root(template_file,"DecayTree",columns=branch_names)
  #Keep 1M events
  #template_sample = template_sample.sample(n=1000000,random_state=9289)
  template_sample = template_sample.query("costheta_D_%s>=-1 and costheta_D_%s<=1 and costheta_L_%s>=-1 and costheta_L_%s<=1 and chi_%s>=-%s and chi_%s<=%s and q2_%s > %s and q2_%s <= %s" % (var_type,var_type,var_type,var_type,var_type,math.pi,var_type,math.pi,var_type,q2_min,var_type,q2_max))
  #Reorder the columns to required order
  template_sample = template_sample[branch_names]

  """
  template_sample_q2 = []
	
  for i in range(0,var_bins["q2_%s" % var_type]):
    template_sample_q2.append(template_sample.query("q2_%s > %s and q2_%s <= %s" % (var_type,qc_bin_vals["q2_%s" % var_type][i],var_type,qc_bin_vals["q2_%s" % var_type][i+1])))
    template_sample_q2[i] = template_sample_q2[i].drop(columns=['q2_%s' % var_type])
    template_sample_q2[i] = template_sample_q2[i].values
  """
	

  ###Ds BACKGROUND
  bkg_file= "/home/ke/graphs/Ds.root"		
  bkg_sample = read_root(bkg_file,"DecayTree",columns=branch_names)
  bkg_sample = bkg_sample.query("costheta_D_%s>=-1 and costheta_D_%s<=1 and costheta_L_%s>=-1 and costheta_L_%s<=1 and chi_%s>=-%s and chi_%s<=%s and q2_%s > %s and q2_%s <= %s" % (var_type,var_type,var_type,var_type,var_type,math.pi,var_type,math.pi,var_type,q2_min,var_type,q2_max))
  #Reorder the columns to required order
  bkg_sample = bkg_sample[branch_names]		
  bkg_sample_q2 = []	
  """
  for i in range(0,var_bins["q2_%s" % var_type]):
    bkg_sample_q2.append(bkg_sample.query("q2_%s > %s and q2_%s <= %s" % (var_type,qc_bin_vals["q2_%s" % var_type][i],var_type,qc_bin_vals["q2_%s" % var_type][i+1])))
    bkg_sample_q2[i] = bkg_sample_q2[i].drop(columns=['q2_%s' % var_type])
    bkg_sample_q2[i] = bkg_sample_q2[i].values	
  """
  #Arrays containing each of the angular weights
  w = {}
  	
  print "Creating weight arrays for each angular term"
  for c in coeffs:
    weight = "w_%s" % c
    branch_names.append(weight)
    w[c] = read_root(template_file,"DecayTree",columns=branch_names)
    #w[c] = w[c].sample(n=1000000,random_state=9289)
    #Weights for each q2 bin
    branch_names.remove(weight)
    """
    for i in range(0,var_bins["q2_%s" % var_type]):
      w["%s_%s" % (c,i)] = w[c].query("costheta_D_%s>=-1 and costheta_D_%s<=1 and costheta_L_%s>=-1 and costheta_L_%s<=1 and chi_%s>=-%s and chi_%s<=%s and q2_%s > %s and q2_%s <= %s" % (var_type,var_type,var_type,var_type,var_type,math.pi,var_type,math.pi,var_type,qc_bin_vals["q2_%s" % var_type][i],var_type,qc_bin_vals["q2_%s" % var_type][i+1]))
      w["%s_%s" % (c,i)] = w["%s_%s" % (c,i)][[weight]]
      w["%s_%s" % (c,i)] = w["%s_%s" % (c,i)].values
      w["%s_%s" % (c,i)] = np.reshape(w["%s_%s" % (c,i)], len(w["%s_%s" % (c,i)]))
    """   

  binning = (4, 4, 8)  	 		
  # List to keep template histograms
  histos = {}
  #Make histogram templates for each angular term
  hist_norm = None
  for c in coeffs:
    print "Creating template for term %s " % (c)
    weight_sample = w["%s" % (c)]
    hist = MakeHistogram(phsp, template_sample, binning,weights = weight_sample)  #bins=?
    if not hist_norm:
      hist_norm = HistogramNorm( hist )
      histos["%s" % (c)] = hist/hist_norm
  ######
  histos_bkg = {}
  #Make histogram templates for each angular term
  hist_bkg_norm = None
  for c in coeffs:
    print "Creating template for term %s " % (c)
    weight_bkg_sample = w["%s" % (c)]
    hist_bkg = MakeHistogram(phsp, bkg_sample, binning,weights = weight_bkg_sample)  #bins=?
    if not hist_bkg_norm:
      hist_bkg_norm = HistogramNorm( hist_bkg )
      histos_bkg["%s" % (c)] = hist_bkg/hist_bkg_norm
  
  #Fit model
  def fit_model(histos):
    pdf = (1.0/3.0)*(4.0 - 6.0*I1s + I2c + 2.0*I2s)*histos["I1c"]
    pdf += I1s*histos["I1s"]
    pdf += I2c*histos["I2c"]
    pdf += I2s*histos["I2s"]
    pdf += I3*histos["I3"]
    pdf += I4*histos["I4"]
    pdf += I5*histos["I5"]
    pdf += I6c*histos["I6c"]
    pdf += I6s*histos["I6s"]
    pdf += I7*histos["I7"]
    pdf += I8*histos["I8"]
    pdf += I9*histos["I9"]
    pdf = Rate*pdf
    pdf += frac_Ds*histos_bkg
    return pdf


  if(toy=="N"):
  	data_file_fit = "/data/lhcb/users/hill/bd2dsttaunu_angular/RapidSim_tuples/Bd2DstTauNu/%s_%s_Total/model_vars_weights.root" % (sub_mode,geom)
  	data_sample_fit = read_root(data_file_fit,"DecayTree",columns=branch_names)
  	data_sample_fit = data_sample_fit.query("costheta_D_%s>=-1 and costheta_D_%s<=1 and costheta_L_%s>=-1 and costheta_L_%s<=1 and chi_%s>=-%s and chi_%s<=%s and q2_%s>=%s and q2_%s<=%s" % (var_type,var_type,var_type,var_type,var_type,math.pi,var_type,math.pi,var_type,q2_min,var_type,q2_max))
  	data_sample_fit = data_sample_fit[branch_names]
  
  	#Randomly sample down to required size
  
  	data_sample_fit = data_sample_fit.sample(n=int(num_sig)*1000,random_state=int(num_sig))
  	
  	#Create datasets for each q2 bin
  	data_sample_fit_q2 = []
  	data_sample_fit_q2_a = []
  	fit_hist = []
  	err_hist = []
  	norm = []
  	
  	for i in range(0,var_bins["q2_%s" % var_type]):
          data_sample_fit_q2.append(data_sample_fit.query("q2_%s > %s and q2_%s <= %s" % (var_type,qc_bin_vals["q2_%s" % var_type][i],var_type,qc_bin_vals["q2_%s" % var_type][i+1])))
          data_sample_fit_q2[i] = data_sample_fit_q2[i].drop(columns=['q2_%s' % var_type])
          data_sample_fit_q2_a.append(data_sample_fit_q2[i].values)
          fit_hist.append(MakeHistogram(phsp, data_sample_fit_q2_a[i], binning[i]))
          err_hist.append(np.sqrt(fit_hist[i] + 0.001))
          norm.append(HistogramNorm(fit_hist[i]))
  	
  else:
  	
  	init_op = tf.initialize_all_variables()
  	sess.run(init_op)
  
  	#Create an instance of the fit PDF, then Poisson vary the values in each bin
  	for i in range(0,var_bins["q2_%s" % var_type]):
          fit_hist.append(sess.run(fit_model(histos,i)))
          #Convert density to number of events
          fit_hist[i] = fit_hist[i]*int(num_sig)*1000*(1.0/var_bins["q2_%s" % var_type])

          fit_hist[i] = np.random.poisson(fit_hist[i])
          
          err_hist.append(np.sqrt(fit_hist[i] + 0.001))
          norm.append(HistogramNorm(fit_hist[i]))

  
  chi2 = []
  result = []
  covmat = []
  
  for i in range(0,var_bins["q2_%s" % var_type]):
  	
    # Define binned Chi2 to be minimised
    chi2.append(BinnedChi2( fit_model(histos,i), fit_hist[i].astype(float)/norm[i], err_hist[i].astype(float)/norm[i] ))
  
    # Run Minuit minimisation
    r, c = tfa.RunMinuit(sess, chi2[i], runHesse=True)
    result.append(r)
    covmat.append(c)
    print result[i]
  	

    #Save covariance matrix
    results_dir = ""
    toy_suf = ""
    if(toy=="N"):
      results_dir = "/home/ke/TensorFlowAnalysis/BinnedResult"
    else:
      results_dir = "/home/ke/TensorFlowAnalysis/BinnedToys"
      toy_rand = random.randint(1,1e10)
      toy_suf = "_%s" % toy_rand
  	
    np.save("%s/cov_%s_%s_%s_%s_q2_%s%s" % (results_dir,sub_mode,geom,var_type,num_sig,i,toy_suf),covmat[i])
  
    #Derived results
    i9=result[i]['I9'][0]
    i8=result[i]['I8'][0]
    i7=result[i]['I7'][0]
    i6s=result[i]['I6s'][0]
    i6c=result[i]['I6c'][0]
    i4=result[i]['I4'][0]
    i5=result[i]['I5'][0]
    i3=result[i]['I3'][0]
    i2s=result[i]['I2s'][0]
    i2c=result[i]['I2c'][0]
    i1s=result[i]['I1s'][0]
    rate=result[i]['Rate'][0]
    (rate,i1s,i2c,i2s,i6c,i6s,i3,i4,i5,i7,i8,i9) = correlated_values([rate,i1s,i2c,i2s,i6c,i6s,i3,i4,i5,i7,i8,i9],covmat[i])
    
    i1c=(4 - 6*i1s + i2c + 2*i2s)/3
    rab=(i1c+2*i1s-3*i2c-6*i2s)/(2*(i1c+2*i1s+i2c+2*i2s))
    rlt= (3*i1c-i2c)/(2*(3*i1s-i2s))
    Gammaq=(3*i1c+6*i1s-i2c-2*i1s)/4.
    afb1=i6c+2*i6s
    afb=(3/8.)*(afb1/Gammaq)
    a3=(1/(np.pi*2))*i3/Gammaq
    a9=(1/(2*np.pi))*i9/Gammaq
    a6s=(-27/8.)*(i6s/Gammaq)
    a4=(-2/np.pi)*i4/Gammaq
    a8=(2/np.pi)*i8/Gammaq
    a5=(-3/4.)*(1-i8-i7-i9-i4-i3-i2s-i1s-i1c-i2c-i6s-i6c)/Gammaq
    a7=(-3/4.)*i7/Gammaq
    para={'RAB':(rab.n,rab.s),'RLT':(rlt.n,rlt.s),'AFB':(afb.n,afb.s),'A6s':(a6s.n,a6s.s),'A3':(a3.n,a3.s),'A9':(a9.n,a9.s),'A4':(a4.n,a4.s),'A8':(a8.n,a8.s),'A5':(a5.n,a5.s),'A7':(a7.n,a7.s), 'I1c': (i1c.n,i1c.s)}
    p = open( "%s/param_%s_%s_%s_%s_q2_%s%s.txt" % (results_dir,sub_mode,geom,var_type,num_sig,i,toy_suf), "w")
    slist=['RAB','RLT','AFB','A6s','A3','A9','A4','A8','A5','A7','I1c']
    for s in slist:
      a=s+" "
      a += str(para[s][0])
      a += " "
      a += str(para[s][1])
      p.write(a + "\n")
    p.close()
    print para
  
    tfa.WriteFitResults(result[i],"%s/result_%s_%s_%s_%s_q2_%s%s.txt" % (results_dir,sub_mode,geom,var_type,num_sig,i,toy_suf))
      #Get final fit PDF
    fit_result = sess.run(fit_model(histos,i))
    
    #1D projections
    fit_hist_proj = {}
    err_hist_proj  = {}
    norm_proj = {}
    fit_result_proj = {}
    data_vals = {}
    branch_names.remove("q2_%s" % var_type)
    for b in branch_names:
  	
      axis = [0,1,2]
      if(b=="costheta_D_%s" % var_type):
        axis.remove(0)
      elif(b=="costheta_L_%s" % var_type):
        axis.remove(1)
      elif(b=="chi_%s" % var_type):
        axis.remove(2)
    	
      if(toy=="N"):
        data_vals["%s_%s" % (b,i)] = data_sample_fit_q2[i][b].values
        #For equi-populated bins
        fit_hist_proj["%s_%s" % (b,i)] = MakeHistogram_1D(data_vals["%s_%s" % (b,i)], (qc_bin_vals["%s_%s" % (b,i)]))
       	#For equal sized bins
        #fit_hist_proj["%s_%s" % (b,i)] = MakeHistogram_1D(data_vals["%s_%s" % (b,i)], var_bins[b])
      else:
        fit_hist_proj["%s_%s" % (b,i)] = np.sum(fit_hist[i], axis=tuple(axis), keepdims=False)
   
      err_hist_proj["%s_%s" % (b,i)] = np.sqrt(fit_hist_proj["%s_%s" % (b,i)])
      norm_proj["%s_%s" % (b,i)] = HistogramNorm(fit_hist_proj["%s_%s" % (b,i)])
      fit_hist_proj["%s_%s" % (b,i)] = fit_hist_proj["%s_%s" % (b,i)].astype(float)/norm_proj["%s_%s" % (b,i)]
      err_hist_proj["%s_%s" % (b,i)] = err_hist_proj["%s_%s" % (b,i)].astype(float)/norm_proj["%s_%s" % (b,i)]
      
      #Binning for equi-populated bins
      bin_centres = []
      bin_width = []
      for j in range(0,len(qc_bin_vals["%s_%s" % (b,i)])-1):
        bin_centres.append(0.5*(qc_bin_vals["%s_%s" % (b,i)][j]+qc_bin_vals["%s_%s" % (b,i)][j+1]))
        bin_width.append(0.5*(qc_bin_vals["%s_%s" % (b,i)][j+1]-qc_bin_vals["%s_%s" % (b,i)][j]))
    
      #Binning for equal sized bins
      #bin_width = 0.5*float(var_range[b][1] - var_range[b][0])/var_bins[b]    
      #bin_centres = []
      #for j in range(0,var_bins[b]):
      #	bin_centres.append(var_range[b][0]+bin_width + j*2*bin_width)
      
      fit_result_proj["%s_%s" % (b,i)] = np.sum(fit_result, axis=tuple(axis), keepdims=False)
  		
      fig,ax = plt.subplots(figsize=(7,7))
      plt.errorbar(bin_centres,fit_result_proj["%s_%s" % (b,i)]/np.sum(fit_result_proj["%s_%s" % (b,i)]),xerr=bin_width,edgecolor=None,ls='none',alpha=0.5,color='b',label="Fit")
      plt.errorbar(bin_centres,fit_hist_proj["%s_%s" % (b,i)],yerr=err_hist_proj["%s_%s" % (b,i)],ls='none',color='k',markersize='3',fmt='o',label="Data")
      
      plt.ylabel("Density")
      plt.xlabel(var_titles[b])
      plt.legend(loc='lower right')
      
      y_min,y_max = ax.get_ylim()
      plt.ylim(0.0,y_max*1.05)
      plt.show()

      if(toy=="N"):
        fig.savefig('/home/ke/TensorFlowAnalysis/BinnedFigs/%s_%s_%s_%s_%s_q2_%s.pdf' % (b,sub_mode,geom,var_type,num_sig,i))
    branch_names.append("q2_%s" % var_type)
    

   

    #Unrolled 1D plot of all bins
    fit_result_1d = fit_result.ravel()
    
    data_norm = fit_hist[i].astype(float)/norm[i]
    data_norm_1d = data_norm.ravel()
    
    err_norm = err_hist[i].astype(float)/norm[i]
    err_norm_1d = err_norm.ravel()
    
    x_max = angle_bins**3
    x = np.linspace(0,x_max-1,x_max)
    
    fig,ax = plt.subplots(figsize=(15,5))
    
    plt.bar(x,fit_result_1d,edgecolor=None,color='r',alpha=0.5,label="Fit")
    plt.errorbar(x,data_norm_1d,yerr=err_norm_1d,ls='none',color='k',markersize='3',fmt='o',alpha=0.8,label="Data")
    
    plt.ylabel("Density")
    plt.xlabel("Bin number")
    plt.xlim(-1,x_max)
    
    plt.legend()
    
    plt.tight_layout()
    plt.show()
    if(toy=="N"):
      fig.savefig('/home/ke/TensorFlowAnalysis/BinnedFigs/Fit_%s_%s_%s_%s_q2_%s.pdf' % (sub_mode,geom,var_type,num_sig,i))
      
  		
    #Pull plot	
    pull = (data_norm_1d - fit_result_1d)/err_norm_1d
    
    fig,ax = plt.subplots(figsize=(15,5))
    
    plt.bar(x,pull,edgecolor='navy',color='royalblue',fill=True)
    
    plt.ylabel("Pull ($\sigma$)")
    plt.xlabel("Bin number")
    plt.xlim(-1,x_max)
    plt.ylim(-5,5)
    
    plt.tight_layout()
    plt.show()
    if(toy=="N"):
      fig.savefig('/home/ke/TensorFlowAnalysis/BinnedFigs/Pull_%s_%s_%s_%s_q2_%s.pdf' % (sub_mode,geom,var_type,num_sig,i))
  	
    #Histogram of the pull values with a fit
    fig,ax = plt.subplots(figsize=(7,7))
      
    pull_bins = int(np.sqrt(x_max))
    n, hist_bins, patches = plt.hist(pull,bins=pull_bins,range=(-5,5),histtype='step',color='navy',normed=True)
    
    plt.xlabel("Pull value ($\\sigma$)")
    plt.ylabel("Fraction of bins")
    
    mu = pull.mean()
    mu_err = sem(pull)
    sigma = pull.std()
    
    plt.title("$\\mu_{Pull} = %.3f \\pm %.3f$, $\\sigma_{Pull} = %.3f$" % (mu,mu_err,sigma))
    
    plt.show()
    if(toy=="N"):
      fig.savefig('/home/ke/TensorFlowAnalysis/BinnedFigs/Pull_Hist_%s_%s_%s_%s_q2_%s.pdf' % (sub_mode,geom,var_type,num_sig,i))
