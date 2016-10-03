#!/usr/bin/python

# Author: J. G. Fernandez-Trincado
# Date: October 03, 2016
# This program compute the chemical abundances element by element

import numpy as np
import scipy as sc
import pylab as plt
import sys
import os
import pyfits 
from scipy.interpolate import interp1d
from scipy import interpolate

# Input data in .fit format 
input_fit_  = sys.argv[1]   # Input spectrum = 2M16011638-1201525.fits (example)
model_      = sys.argv[2]   # Input atmosphere model = s4500_g+1.5_m1.0_t02_x2_z-1.50_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod (example) 

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

	elem = {'C':6, 'N':7, 'O':8, 'Na':11, 'Mg': 12, 'Al':13, 'Mg': 12}

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
	f.write("set METALLIC   = '-1.300'\n")
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

########################################################################################################################
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

#spec = open(sys.argv[1]+'_vac.dat','w') # file with wavelenght in vac ...
#for lam in np.arange(len(flux)): spec.write(str(wlnew[lam])+"\t"+str(flux[lam])+"\n")
#spec.close()
########################################################################################################################
# Step 3:

init_elem = 'Al'

solar     = {'Al':6.37, 'Mg':7.53, 'N':'7.78'} # Asplund 2005 
init_     = np.linspace(solar[init_elem]-2.,solar[init_elem]+2., 10.) # change this for the correct elment
windows   = sc.genfromtxt('windows_apogee_manual_vac.dat', dtype=str)
mask      = (windows[:,0] == init_elem)  

windows_min = windows[mask,2]
windows_max = windows[mask,3] 
windows_N   = int(len(windows_max))
Nwin        = 1

ab = open('abu','a')

Abu, Chi2abu = np.array([]), np.array([])

for i in np.arange(len(init_)):

	file_(init_[i], init_elem) # change this for the correct element 

	os.system("ls syntspec/"+str(model_)+"_15000-17000_xit2.3_"+str(init_[i])+".spec > tempsys.temp")
	os.system("./syntspec/manualfaltbo < tempsys.temp")
	os.system("rm tempsys.temp")
	os.system("rm syntspec/"+str(model_)+"_15000-17000_xit2.3_"+str(init_[i])+".spec")
	os.system("mv out.conv out."+str(init_[i])+"_.dat")

	abunfit = sc.genfromtxt("out."+str(init_[i])+"_.dat") # Synthetic spectrum ... 

	maskabu = (abunfit[:,0]>= float(windows_min[Nwin])) & (abunfit[:,0]<= float(windows_max[Nwin]))  
	syntheX, syntheY = abunfit[maskabu,0], abunfit[maskabu,1] 

	maskobs = (wlnew >= float(windows_min[Nwin])) & (wlnew <= float(windows_max[Nwin]))
	obsX,  obsY  = wlnew[maskobs], images2[1].data[maskobs]

	f2 = interp1d(syntheX, syntheY, kind='cubic')

	Chi2 = np.sum(((obsY - f2(obsX))*(obsY - f2(obsX))) / f2(obsX)) 

	plt.subplot(1,2,1)

	plt.plot(syntheX, syntheY, '--')
#	plt.plot(obsX, f2(obsX), 'x', color='black')
	plt.plot(obsX, obsY, 's', color='black')

#	plt.subplot(1,2,2)
#	plt.plot(init_[i], Chi2, 's', color='red')

	Abu     = np.append(Abu, init_[i])
	Chi2abu = np.append(Chi2abu, Chi2)

	ab.write(str(init_[i])+'\t'+str(Chi2)+'\n')

	plt.xlabel('Wavelength [Angstrom]',fontsize=25)
	plt.ylabel('Flux', fontsize=25)

	os.system("rm out."+str(init_[i])+"_.dat")


ab.close()

plt.subplot(1,2,2)
plt.plot(Abu, Chi2abu, 's', color='red')

x_new = np.linspace(np.min(init_),np.max(init_),50000)

f = interpolate.interp1d(Abu, Chi2abu, 'quadratic')
plt.plot(x_new, f(x_new),'-',color='black')
mask2nd = (f(x_new) == float(np.min(f(x_new))))
plt.axvline(x=float(x_new[mask2nd][0]), color='black')
plt.ylabel(r'Chi$^2$',fontsize=25)
plt.xlabel("A("+init_elem+")",fontsize = 25)

plt.subplot(1,2,1)
file_(float(x_new[mask2nd][0]), init_elem)

os.system("ls syntspec/"+str(model_)+"_15000-17000_xit2.3_"+str(float(x_new[mask2nd][0]))+".spec > tempsys.temp")
os.system("./syntspec/manualfaltbo < tempsys.temp")
os.system("rm tempsys.temp")
os.system("rm syntspec/"+str(model_)+"_15000-17000_xit2.3_"+str(float(x_new[mask2nd][0]))+".spec")
os.system("mv out.conv out."+str(float(x_new[mask2nd][0]))+"_.dat")

abunfit = sc.genfromtxt("out."+str(float(x_new[mask2nd][0]))+"_.dat") # Synthetic spectrum ... 

maskabu = (abunfit[:,0]>= float(windows_min[Nwin])) & (abunfit[:,0]<= float(windows_max[Nwin]))
syntheX, syntheY = abunfit[maskabu,0], abunfit[maskabu,1]
plt.title("A("+init_elem+") = "+str(float(x_new[mask2nd][0])),fontsize=25)
plt.plot(syntheX, syntheY, '-',lw=3, color='grey')


plt.show()











