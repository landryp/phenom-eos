#!/usr/bin/python

#__doc__ = 'Constructs a piecewise polytrope EoS from its parameters'
#__usage__ = 'get_pwp.py log10p1 [Gamma1,Gamma2,Gamma3]'
#__author__ = 'philippe.landry@ligo.org'
#__date__ = '02-2019'

import numpy as np
from constants import c

def pwpoly(log10p1, Gamma, div=0, stp=1e2, pmax=1e17, pmin=1e1):

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
	a = np.zeros(seg)
	
	divp[0] = 10.**log10p1/c**2 # pressure at first dividing density in [g/cm^3]
	K[0] = divp[0]/(divrho[0]**Gamma[0])
	rhomin = (pmin/K[0])**(1./Gamma[0])
	a[0] = - pmin/(rhomin*(Gamma[0]-1.))
	
	def rho(i,p):
	
		return (p/K[i])**(1./Gamma[i])
		
	def mu(i,p):
	
		return (1.+a[i])*rho(i,p) + p/(Gamma[i]-1.)
	
	for i in range(1,seg):
	
		K[i] = divp[i-1]/(divrho[i-1]**Gamma[i])
		a[i] = mu(i-1,divp[i-1])/divrho[i-1] - 1. - K[i]*divrho[i-1]**(Gamma[i]-1.)/(Gamma[i]-1.)	
		divp[i] = K[i]*divrho[i]**Gamma[i]

	print K[0], K[1], K[2], K[3]
# EXPORT EOS DATA

	def findseg(p):
	
		test = np.full(seg,p) - divp
		
		i = 0 
		while (test[i] > 0. and i < seg-1):
		
			i = i+1
				
		return i
		
	rhodat = []
	mudat = []
	
	pdat = np.logspace(np.log10(pmin),np.log10(pmax),stp)
	
	for p in pdat:
	
		mupt = mu(findseg(p),p)
		rhopt = rho(findseg(p),p)
		
		mudat.append(mupt) # total energy density in g/cm^3
		rhodat.append(rhopt) # baryon density in g/cm^3

	return rhodat, mudat, pdat	
