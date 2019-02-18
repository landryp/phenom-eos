#!/usr/bin/python

#__doc__ = 'Constructs a spectral EoS from its parameters'
#__usage__ = 'get_spec.py g0 g1 g2 g3'
#__author__ = 'philippe.landry@ligo.org'
#__date__ = '02-2019'

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interp

def spec(g0, g1, g2, g3, stp=1e2, xmax=1e1, eps=1e-4, pts=10):

# DEFINE CONSTANTS

	c = 2.99792458e10 # speed of light in cgs
	G = 6.67408e-8 # Newton's constant in cgs
	Msun = 1.3271244e26/G # solar mass in cgs
	rhonuc = 2.7e14 # nuclear density in cgs

# FIX CRUST-CORE INTERFACE AND CALCULATE INTERFACE VALUES

	rho0 = 0.5 # interface density in units of rhonuc
	
	p0 = pSLy(rho0)[0]
	mu0 = muSLy(rho0)[0]

# DEFINE CORE PARAMETERS   # Shouldn't need to specify scales here!

	def G(x):
	
		term0 = g0
		term1 = g1*x
		term2 = g2*x**2
		term3 = g3*x**3
	
		return np.exp(term0+term1+term2+term3)
		
	def rho(x): # rest-mass energy density in units of arbitrary scale rho0
	
		integral, err = integrate.quad(lambda x1: 1./G(x1), 0., x)
	
		return np.exp(integral)
	
	def mu(x): # total energy density in units of corresponding scale mu0
	
		integral, err = integrate.quad(lambda x1: rho0*np.exp(x1)/(rho(x1)*G(x1)), 0., x)
	
		return rho(x) + (p0/mu0)*rho(x)*integral

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

