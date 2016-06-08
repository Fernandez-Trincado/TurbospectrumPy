#!/usr/bin/python
# Athor: J. G. Fernandez-Trincado
# Any comments, suggestions, and your complaints, are welcome.
# Contact me: jfernandez@obs-besancon.fr and/or jfernandezt87@gmail.com

import numpy as np
import scipy as sc
import pylab as plt
import sys
import pyfits 


# Input data in .fit format 

input_fit_ = pyfits.open('2M17352602-1800213.fits')

Teff = input_fit_[4].data['FPARAM'][0][0] # Teff 
Logg = input_fit_[4].data['FPARAM'][0][1] # LOGG 
MH   = input_fit_[4].data['FPARAM'][0][3] # [M/H]

GridTeff = [3500, 3600, 3700, 3800, 3900, 4000]
GridLogg = [0,0.5,1.0,1.5,2.0,2.5,3.0]
GridMH      = [0.00,0.50]
GridMHstr   = ['0.00','0.50']

for i in np.arange(len(GridTeff)-1):

	if (Teff >= GridTeff[i] ) & (Teff < GridTeff[i+1]):
		diff1 = np.abs(Teff - GridTeff[i])
		diff2 = np.abs(Teff - GridTeff[i+1])
		if diff1 < diff2: Teffend = GridTeff[i]	
		else: Teffend = GridTeff[i+1]
	else: pass

for j in np.arange(len(GridLogg)-1):
	
	if (Logg >= GridLogg[j]) & (Logg < GridLogg[j+1]):
		diff1 = np.abs( Logg - GridLogg[j])
		diff2 = np.abs( Logg - GridLogg[j+1])
		if diff1 < diff2: Loggend = GridLogg[j]
		else: Loggend = GridLogg[j+1]
	else: pass

for k in np.arange(len(GridMH)-1):


	if (MH >= GridMH[k]) & ( MH < GridMH[k+1]):
		diff1 = np.abs(MH - GridMH[k])	
		diff2 = np.abs(MH - GridMH[k+1])
		if diff1 < diff2: MHend = str(GridMHstr[k])
		else: MHend = str(GridMHstr[k+1])
	else: pass


# Creating the Script ...
f = open('test.com','w')
f.write("#!/bin/csh -f \n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
f.write("\n")
f.write("date\n")
f.write("set mpath=models\n")
f.write("\n")
f.write("foreach MODEL (s"+str(Teffend)+"_g+"+str(Loggend)+"_m1.0_t02_x2_z+"+str(MHend)+"_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod)\n")
f.write("\n")
f.write("set Mgabu = 8.71185\n")
f.write("\n")
f.write("set lam_min    = '15000'\n")
f.write("set lam_max    = '17000'\n")
f.write("\n")
f.write("set deltalam   = '0.05'\n")
f.write("set METALLIC   = '+0.327'\n")
f.write("set TURBVEL    = '2.3'\n")
f.write("set SUFFIX     = _${lam_min}-${lam_max}_xit${TURBVEL}_${Mgabu}.spec\n")
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
f.write("6  ${Mgabu}\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
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
f.write("6  ${Mgabu}\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
f.write("'ISOTOPES : ' '0'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
f.write("'NFILES   :' '3'\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
f.write("DATA/Hlinedata\n")                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
f.write("linelists/privative_line_list\n") You need authorization for this line list, sorry I can not help you with this                                                                                                                                                                                                                                                                                                                                                                                                                                                             
f.write("linelists/privative_line_list\n") You need authorization for this line list, sorry I can not help you with this                                                                                                                                                                                                                                                                                                                                                                                                                                                         
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

# Step 2:
# Convertion of the spectrum from vac to air ...

def vac2air(wave,sdssweb = False):
	if sdssweb: return wave/(1.+2.735182*10.**-4.+131.4182/wave**2.+2.76249*10.**8./wave**4.)
	else:       return wave/(1.+0.05792105/(238.0185-(10000./wave)**2.)+0.00167917/(57.362-(10000./wave)**2.))

# input_ = sc.genfromtxt(sys.argv[1])












