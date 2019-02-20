#!/usr/bin/python

#__doc__ = 'Construct a linearized strange quark matter EoS from its bag-model parameters'
#__usage__ = 'get_sqm.py B cqm2'
#__author__ = 'philippe.landry@ligo.org'
#__date__ = '02-2019'

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interp

def sqm(B, cqm2, stp=1e2, pmax=1e17, eps=1e-4): # bag constant B in eV, SQM sound speed squared in c^2

# DEFINE CORE PARAMETERS

	def mu(p): # total energy density in g/cm^3 as a function of pressure in g/cm^3
	
		return 4.*B**4 + p/cqm2
		
	def rho(p): # rest-mass energy density in g/cm^3, matching rho and mu at surface p=0
	
		num = 4.**cqm2*(mu(p)+p)*cqm2
		denom = cqm2*B**4
		power = 1./(1.+cqm2)
	
		return B**4*(num/denom)**power

# EXPORT EOS DATA

	rhodat = []
	mudat = []

	pdat = np.logspace(np.log10(eps),np.log10(pmax),stp)
	
	for p in pdat:
	
		mupt = mu(p)
		rhopt = rho(p)
		
		mudat.append(mupt) # total energy density in g/cm^3
		rhodat.append(rhopt) # baryon density in g/cm^3

	return rhodat, mudat, pdat
