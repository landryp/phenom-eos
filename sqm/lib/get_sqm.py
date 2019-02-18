#!/usr/bin/python

#__doc__ = 'Constructs a linearized strange quark matter EoS from its bag-model parameters'
#__usage__ = 'get_sqm.py B cqm2'
#__author__ = 'philippe.landry@ligo.org'
#__date__ = '02-2019'

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interp

def sqm(B, cqm2, stp=1e2, xmax=1e1, eps=1e-4, xi=-1.):

# DEFINE CONSTANTS

	c = 2.99792458e10 # speed of light in cgs
	G = 6.67408e-8 # Newton's constant in cgs
	Msun = 1.3271244e26/G # solar mass in cgs
	rhonuc = 2.7e14 # nuclear density in cgs

# DEFINE CORE PARAMETERS

	def p(x): # pressure
	
		return rhonuc*np.exp(x)

	def mu(x): # total energy density
	
		return 4.*B**4 + p(x)/cqm2
		
	def rho(x): # rest-mass energy density
	
		integral, err = integrate.quad(lambda x1: p(x1)/(mu(x1)+p(x1)), xi, x)
	
		return rhonuc*np.exp(integral/cqm2)  # note: rhonuc scale is wrong here

# EXPORT EOS DATA

	rhodat = []
	pdat = []
	mudat = []

	xdat = np.linspace(eps,xmax,stp)
	
	for x in xdat:
	
		ppt = p(x)
		mupt = mu(x)
		rhopt = rho(x)
		
		rhodat.append(rhopt) # baryon density in g/cm^3
		pdat.append(ppt) # pressure in g/cm^3
		mudat.append(mupt) # total energy density in g/cm^3

	return rhodat, mudat, pdat
