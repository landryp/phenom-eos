#!/usr/bin/python

#__doc__ = 'Constructs a linearized strange quark matter EoS from its bag-model parameters'
#__usage__ = 'get_sqm.py B cqm2'
#__author__ = 'philippe.landry@ligo.org'
#__date__ = '02-2019'

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interp

def sqm(B, cqm2, stp=1e2, pmax=1e1, eps=1e-4, pts=10):

# DEFINE CONSTANTS

	c = 2.99792458e10 # speed of light in cgs
	G = 6.67408e-8 # Newton's constant in cgs
	Msun = 1.3271244e26/G # solar mass in cgs
	rhonuc = 2.7e14 # nuclear density in cgs

# DEFINE CORE PARAMETERS

	def mu(p): # total energy density
	
		return 4.*B**4 + p/cqm2
		
	def rho(p): # rest-mass energy density in units of central value
	
		integral, err = integrate.quad(lambda p1: 1./(mu(p1)+p1), pc, p)
	
		return np.exp(integral/cqm2)

# EXPORT EOS DATA

	rhodat = []
	pdat = []
	mudat = []
	xpts = []

	xdat = np.linspace(eps,xmax,pts)
	
	for x in xdat:
	
		ppt = p0*np.exp(x)
				
		if ppt > 0.:
		
			mupt = mu(x)
		
			if mupt > 0.:
			
				rhopt = rho(x)
			
				if rhopt > 0.:
		
					rhodat.append(rhopt*rhonuc) # baryon density in g/cm^3
					pdat.append(ppt*rhonuc) # pressure in g/cm^3
					mudat.append(mupt*rhonuc) # total energy density in g/cm^3
					xpts.append(x)	

	rhofunc = interp.interp1d(xpts,rhodat,kind='linear',bounds_error=False,fill_value=None)
	pfunc = interp.interp1d(xpts,pdat,kind='linear',bounds_error=False,fill_value=None)
	mufunc = interp.interp1d(xpts,mudat,kind='linear',bounds_error=False,fill_value=None)
	
	xdat = np.linspace(eps,xmax,stp)

	return rhofunc(xdat), mufunc(xdat), pfunc(xdat)

