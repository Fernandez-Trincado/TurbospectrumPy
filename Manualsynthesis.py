#!/usr/bin/python
############################################################################################################################################################
# Author: J. G. Fernandez-Trincado
# Date: October 03, 2016
# This program compute the chemical abundances element by element
# Example: python Manualsynthesis.py 2M16011638-1201525.fits s4500_g+1.5_m1.0_t02_x2_z-1.50_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod -1.30 
#              APOGEE star: 2M16011638-1201525.fits      
#         Atmosphere model: s4500_g+1.5_m1.0_t02_x2_z-1.50_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod
#         Star metallicity: -1.30
############################################################################################################################################################
import numpy as np
import scipy as sc
import pylab as plt
import sys
import os
import pyfits 
from scipy.interpolate import interp1d
from scipy import interpolate
from pylab import *
import matplotlib.lines as lines
############################################################################################################################################################
# Input data in .fit format 
input_fit_  = sys.argv[1]   # Input spectrum = 2M16011638-1201525.fits (example)
model_      = sys.argv[2]   # Input atmosphere model = s4500_g+1.5_m1.0_t02_x2_z-1.50_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod (example) 
metall      = sys.argv[3]
#seq         = sys.argv[4]
init_elem   = 'Al'
solar       = {'Al':6.37, 'Mg':7.53, 'N':7.78}                          # Asplund 2005
elem        = {'C':6, 'N':7, 'O':8, 'Na':11, 'Mg': 12, 'Al':13, 'Mg': 12} # Atomic weight
############################################################################################################################################################
#Teff = input_fit_[4].data['FPARAM'][0][0] # Teff 
#Logg = input_fit_[4].data['FPARAM'][0][1] # LOGG 
#MH   = input_fit_[4].data['FPARAM'][0][3] # [M/H]
#GridTeff = [3500, 3600, 3700, 3800, 3900, 4000]
#GridLogg = [0,0.5,1.0,1.5,2.0,2.5,3.0]
#GridMH      = [0.00,0.50]
#GridMHstr   = ['0.00','0.50']
#for i in np.arange(len(GridTeff)-1):
#
#	if (Teff >= GridTeff[i] ) & (Teff < GridTeff[i+1]):
#		diff1 = np.abs(Teff - GridTeff[i])
#		diff2 = np.abs(Teff - GridTeff[i+1])
#		if diff1 < diff2: Teffend = GridTeff[i]	
#		else: Teffend = GridTeff[i+1]
#	else: pass
#for j in np.arange(len(GridLogg)-1):
#	
#	if (Logg >= GridLogg[j]) & (Logg < GridLogg[j+1]):
#		diff1 = np.abs( Logg - GridLogg[j])
#		diff2 = np.abs( Logg - GridLogg[j+1])
#		if diff1 < diff2: Loggend = GridLogg[j]
#		else: Loggend = GridLogg[j+1]
#	else: pass
#for k in np.arange(len(GridMH)-1):
#
#
#	if (MH >= GridMH[k]) & ( MH < GridMH[k+1]):
#		diff1 = np.abs(MH - GridMH[k])	
#		diff2 = np.abs(MH - GridMH[k+1])
#		if diff1 < diff2: MHend = str(GridMHstr[k])
#		else: MHend = str(GridMHstr[k+1])
#	else: pass
############################################################################################################################################################
# Creating the Script ...
def file_(inputabu, xelem):
	defelem = elem[xelem]
	f = open('Manual.com','w')
	f.write("#!/bin/csh -f \n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
	f.write("\n")
	f.write("date\n")
	f.write("set mpath=models\n")
	f.write("\n")
	#f.write("foreach MODEL (s"+str(Teffend)+"_g+"+str(Loggend)+"_m1.0_t02_x2_z+"+str(MHend)+"_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod)\n")
	f.write("foreach MODEL ("+model_+")")
	f.write("\n")
	f.write("\n")
	f.write("set Xabu = "+str(inputabu)+"\n")
	f.write("\n")
	f.write("set lam_min    = '15000'\n")
	f.write("set lam_max    = '17000'\n")
	f.write("\n")
	f.write("set deltalam   = '0.05'\n")
	f.write("set METALLIC   = '"+str(float(metall))+"'\n")
	f.write("set TURBVEL    = '2.3'\n")
	f.write("set SUFFIX     = _${lam_min}-${lam_max}_xit${TURBVEL}_${Xabu}.spec\n")
	f.write("set result     = ${MODEL}${SUFFIX}\n")
	f.write("\n")
	f.write("# \n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
	f.write("# ABUNDANCES FROM THE MODEL ARE NOT USED !!!\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                
	f.write("\n")
	f.write("../exec-gf-v15.1/babsma_lu << EOF\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                           
	f.write("'LAMBDA_MIN:'  '${lam_min}'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
	f.write("'LAMBDA_MAX:'  '${lam_max}'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
	f.write("'LAMBDA_STEP:' '${deltalam}'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                
	f.write("'MODELINPUT:' '$mpath/${MODEL}'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                             
	f.write("'MARCS-FILE:' '.true.'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
	f.write("'MODELOPAC:' 'contopac/${MODEL}opac'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                        
	f.write("'METALLICITY:'    '${METALLIC}'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                             
	f.write("'ALPHA/Fe   :'    '0.00'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("'HELIUM     :'    '0.00'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("'R-PROCESS  :'    '0.00'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("'S-PROCESS  :'    '0.00'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("'INDIVIDUAL ABUNDANCES:'   '1'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                              
	f.write(str(defelem)+"  ${Xabu}\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
	f.write("'XIFIX:' 'T'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
	f.write("$TURBVEL\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("EOF\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
	f.write("\n")
	f.write("########################################################################\n")                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("\n")
	f.write("../exec-gf-v15.1/bsyn_lu <<EOF\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                              
	f.write("'LAMBDA_MIN:'     '${lam_min}'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                              
	f.write("'LAMBDA_MAX:'     '${lam_max}'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                              
	f.write("'LAMBDA_STEP:'    '${deltalam}'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                             
	f.write("'INTENSITY/FLUX:' 'Flux'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("'COS(THETA)    :' '1.00'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("'ABFIND        :' '.false.'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
	f.write("'MODELOPAC:' 'contopac/${MODEL}opac'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                        
	f.write("'RESULTFILE :' 'syntspec/${result}'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                         
	f.write("'METALLICITY:'    '${METALLIC}'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                             
	f.write("'ALPHA/Fe   :'    '0.00'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("'HELIUM     :'    '0.00'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("'R-PROCESS  :'    '0.00'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("'S-PROCESS  :'    '0.00'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("'INDIVIDUAL ABUNDANCES:'   '1'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                              
	f.write(str(defelem)+"  ${Xabu}\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
	f.write("'ISOTOPES : ' '0'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
	f.write("'NFILES   :' '3'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
	f.write("DATA/Hlinedata\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
	f.write("linelists/turboatoms.20150714\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                               
	f.write("linelists/turbomolec.20150714\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                               
	f.write("'SPHERICAL:'  'T'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
	f.write("  30\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
	f.write("  300.00\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("  15    \n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("  1.30  \n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("EOF\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("########################################################################\n")                                                                                                                                                                                                                                                                                                                                                                                                                    
	f.write("date\n")
	f.write("end\n")
	f.close()
	os.system('chmod +x Manual.com')
	os.system('./Manual.com')
	os.system("rm Manual.com")
############################################################################################################################################################
# Step 2:
# Convertion of the spectrum from vac to air ...

def vac2air(wave,sdssweb = False):
	if sdssweb: return wave/(1.+2.735182*10.**-4.+131.4182/wave**2.+2.76249*10.**8./wave**4.)
	else:       return wave/(1.+0.05792105/(238.0185-(10000./wave)**2.)+0.00167917/(57.362-(10000./wave)**2.))

images2 = pyfits.open(sys.argv[1])

wl = np.array([])
for i in np.arange(8575):  wl = np.append(wl, 4.179 + 6.E-6*i)

wl    = 10.**(wl)
flux  = images2[1].data
wlnew = vac2air(wl,sdssweb = False)

input_obs   = flux 

#spec = open(sys.argv[1]+'_vac.dat','w') # file with wavelenght in vac ...
#for lam in np.arange(len(flux)): spec.write(str(wlnew[lam])+"\t"+str(flux[lam])+"\n")
#spec.close()
############################################################################################################################################################
# Step 3:


init_        = np.linspace(float(solar[init_elem])-2.,float(solar[init_elem]+2.), 10.) # change this for the correct elment
windows      = sc.genfromtxt('windows_apogee_manual_vac.dat', dtype=str)
mask         = (windows[:,0] == init_elem)  
windows_min  = windows[mask,2]
windows_max  = windows[mask,3] 
windows_N    = int(len(windows_max))


for seq in np.arange(windows_N):

	Abu, Chi2abu = np.array([]), np.array([])
	
	Nwin        = seq #float(seq) - 1
	
	fig = plt.figure(int(seq)+1,(25,15))
	for i in np.arange(len(init_)):
		file_(init_[i], init_elem) # change this for the correct element 
		os.system("ls syntspec/"+str(model_)+"_15000-17000_xit2.3_"+str(init_[i])+".spec > tempsys.temp")
		os.system("./syntspec/manualfaltbo < tempsys.temp")
		os.system("rm tempsys.temp")
		os.system("rm syntspec/"+str(model_)+"_15000-17000_xit2.3_"+str(init_[i])+".spec")
		os.system("mv out.conv out."+str(init_[i])+"_.dat")
		abunfit = sc.genfromtxt("out."+str(init_[i])+"_.dat") # Synthetic spectrum ... 

		input_model = abunfit[:,1] 
		mask  = (input_obs   > np.median(input_obs) )
		mask2 = (input_model > np.median(input_model) )
		obs   = np.mean(input_obs[mask])
		model = np.mean(input_model[mask2])
		diff1 = (obs - model)
		if diff1 < 0:
			factor  = np.abs(diff1)
			newflux = input_model - factor
		else:   
			factor  = np.abs(diff1)
			newflux = input_model + factor

		maskabu = (abunfit[:,0]>= float(windows_min[Nwin])) & (abunfit[:,0]<= float(windows_max[Nwin]))
		syntheX, syntheY = abunfit[maskabu,0], newflux[maskabu] # abunfit[maskabu,1] 
		maskobs = (wlnew >= float(windows_min[Nwin])) & (wlnew <= float(windows_max[Nwin]))
		obsX,  obsY  = wlnew[maskobs], images2[1].data[maskobs]
		maskobs = (obsX > np.min(syntheX)) &  (obsX < np.max(syntheX))
		obsX, obsY = obsX[maskobs], obsY[maskobs]
		f2 = interp1d(syntheX, syntheY, kind='cubic')
		Chi2 = np.sum(((obsY - f2(obsX))*(obsY - f2(obsX))) / f2(obsX)) 
		ax = fig.add_subplot(1,2,1)
		ax.plot(syntheX/1E4, syntheY, '--')
		ax.plot(obsX/1E4, obsY, 's', color='black')
		min_ab, max_ab = np.min(syntheX), np.max(syntheX)
		x              = np.linspace(min_ab/1E4, max_ab/1E4, len(wl[mask]))
		locs,labels    = xticks()
		xticks(locs, map(lambda x: "%g" % x, locs))
		Abu     = np.append(Abu, init_[i])
		Chi2abu = np.append(Chi2abu, Chi2)
		ax.set_xlabel('Wavelength [Angstrom]',fontsize=25)
		ax.set_ylabel('Flux', fontsize=25)
		os.system("rm out."+str(init_[i])+"_.dat")
		ax.tick_params(labelsize = 20)
	
	
	ax2 = fig.add_subplot(1,2,2)
	ax2.plot(Abu, Chi2abu, 's', color='red')
	x_new = np.linspace(np.min(init_),np.max(init_),50000)
	f = interpolate.interp1d(Abu, Chi2abu, 'quadratic')
	ax2.plot(x_new, f(x_new),'-',color='black')
	mask2nd = (f(x_new) == float(np.min(f(x_new))))
	ax2.axvline(x=float(x_new[mask2nd][0]), color='black')
	ax2.set_ylabel(r'$\chi^2$',fontsize=25)
	ax2.set_xlabel("A("+init_elem+")",fontsize = 25)
	ax2.set_title("A("+init_elem+") = "+str(float(x_new[mask2nd][0])),fontsize=25)
	#ax.add_subplot(1,2,1)
	file_(float(x_new[mask2nd][0]), init_elem)
	os.system("ls syntspec/"+str(model_)+"_15000-17000_xit2.3_"+str(float(x_new[mask2nd][0]))+".spec > tempsys.temp")
	os.system("./syntspec/manualfaltbo < tempsys.temp")
	os.system("rm tempsys.temp")
	os.system("rm syntspec/"+str(model_)+"_15000-17000_xit2.3_"+str(float(x_new[mask2nd][0]))+".spec")
	os.system("mv out.conv out."+str(float(x_new[mask2nd][0]))+"_.dat")
	abunfit = sc.genfromtxt("out."+str(float(x_new[mask2nd][0]))+"_.dat") # Synthetic spectrum ... 

	input_model = abunfit[:,1]
	mask  = (input_obs   > np.median(input_obs) )
	mask2 = (input_model > np.median(input_model) )
	obs   = np.mean(input_obs[mask])
	model = np.mean(input_model[mask2])
	diff1 = (obs - model)
	if diff1 < 0:
		factor  = np.abs(diff1)
		newflux = input_model - factor
	else:   
		factor  = np.abs(diff1)
		newflux = input_model + factor

	maskabu = (abunfit[:,0]>= float(windows_min[Nwin])) & (abunfit[:,0]<= float(windows_max[Nwin]))
	syntheX, syntheY = abunfit[maskabu,0], newflux[maskabu] # abunfit[maskabu,1]
	ax.plot(syntheX/1E4, syntheY, '-',lw=3, color='grey')
	os.system("rm out."+str(float(x_new[mask2nd][0]))+"_.dat")
	ax2.tick_params(labelsize = 20)
	
	plt.savefig(sys.argv[1]+"_elem_"+str(init_elem)+"_win_"+str(int(Nwin)+1)+"_.png")



