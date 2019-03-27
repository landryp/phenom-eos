#!/usr/bin/python

#__doc__ = 'Construct a linearized strange quark matter EoS from its bag-model parameters'
#__usage__ = 'get_sqm.py B cqm2'
#__author__ = 'philippe.landry@ligo.org'
#__date__ = '02-2019'

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interp

def sqm(B, cqm2, stp=1e2, pmax=1e17, eps=1e-4): # bag constant B in MeV, SQM sound speed squared in c^2

	MeVc2tog = 1.78266191e-27
	fmtocm = 1e-13
	fac = MeVc2tog/fmtocm**3

# DEFINE CORE PARAMETERS

	def mu(p): # total energy density in g/cm^3 as a function of pressure in g/cm^3
	
		return 4.*B*fac + p/cqm2 # from Han+Steiner 1810.10967
		
	def rho(p): # rest-mass energy density in g/cm^3, matching rho and mu at surface p=0
	
		num = 4.**cqm2*(mu(p)+p)*cqm2
		denom = cqm2*B*fac
		power = 1./(1.+cqm2)
	
		return B*fac*(num/denom)**power # follows from 1st law of thermodynamics

# EXPORT EOS DATA

	rhodat = []
	mudat = []

	pdat = np.logspace(np.log10(eps),np.log10(pmax),stp) # pressure in g/cm^3
	
	for p in pdat:
	
		mupt = mu(p)
		rhopt = rho(p)
		
		mudat.append(mupt) # total energy density in g/cm^3
		rhodat.append(rhopt) # baryon density in g/cm^3

	return rhodat, mudat, pdat
