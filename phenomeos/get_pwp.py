#!/usr/bin/python

#__doc__ = 'Constructs a piecewise polytrope EoS from its parameters'
#__usage__ = 'get_pwp.py log10p1 [Gamma1,Gamma2,Gamma3]'
#__author__ = 'philippe.landry@ligo.org'
#__date__ = '02-2019'

import numpy as np
from .constants import *

def rho(i,p,K,Gamma): # rest mass density in each segment
	
	return (p/K[i])**(1./Gamma[i])
		
def mu(i,p,a,rho,K,Gamma): # total energy density in each segment
	
	return (1.+a[i])*rho(i,p,K,Gamma) + p/(Gamma[i]-1.)
	
def findseg(p,seg,divp): # identify to which segment given pressure belongs
	
	test = np.full(seg,p) - divp
		
	i = 0 
	while (test[i] > 0. and i < seg-1):
		
		i = i+1
				
	return i

def pwpoly(log10p1, Gamma, div=0, stp=1e2, pmax=1e17, pmin=1e0, plist=False):

# DEFINE SEGMENTS AND DIVIDING DENSITIES

	seg = len(Gamma) # number of segments in piecewise polytrope
	divrho = np.zeros(seg)

	if div == 0:
	
		for i in range(seg): # default dividing rest-mass densities [g/cm^3]
		
			divrho[i] = 10.**(14.7+i*0.3)
		
	else:
	
		divrho = div # use chosen dividing densities instead
		divrho = np.append(divrho,1e30)
		
# DEFINE CORE PARAMETERS

	K = np.zeros(seg)
	divp = np.zeros(seg)
	a = np.zeros(seg) # below taken from Read+ PRD 79 (2009)
	
	divp[0] = 10.**log10p1/c**2 # pressure at first dividing density in [g/cm^3]
	K[0] = divp[0]/(divrho[0]**Gamma[0]) # constant K in first segment
	rhomin = (pmin/K[0])**(1./Gamma[0]) # lowest rest mass density
	a[0] = - pmin/(rhomin*(Gamma[0]-1.)) # constant a in first segment
	
	for i in range(1,seg): # constants in each segment, plus dividing pressures
	
		K[i] = divp[i-1]/(divrho[i-1]**Gamma[i])
		a[i] = mu(i-1,divp[i-1],a,rho,K,Gamma)/divrho[i-1] - 1. - K[i]*divrho[i-1]**(Gamma[i]-1.)/(Gamma[i]-1.)	
		divp[i] = K[i]*divrho[i]**Gamma[i]

# EXPORT EOS DATA
		
	rhodat = []
	mudat = []
	pdat_out = []
	
	if len(plist) > 1:
	
		pdat = plist
	
	else:
	
		pdat = np.logspace(np.log10(pmin),np.log10(pmax),stp) # pressure in g/cm^3
	
	for count,p in enumerate(pdat):
	
		mupt = mu(findseg(p,seg,divp),p,a,rho,K,Gamma)
		rhopt = rho(findseg(p,seg,divp),p,K,Gamma)
		if count > 1 and (p - pdat[count-1])/(mupt - mudat[count-1]) > c**2: break
		
		pdat_out.append(p)
		mudat.append(mupt) # total energy density in g/cm^3
		rhodat.append(rhopt) # baryon density in g/cm^3

	return rhodat, mudat, pdat_out	
