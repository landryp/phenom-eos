#!/usr/bin/python

#__doc__ = 'Constructs a spectral EoS from its parameters'
#__usage__ = 'get_spec.py g0 g1 g2 g3'
#__author__ = 'philippe.landry@ligo.org'
#__date__ = '02-2019'

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interp
from constants import *

def spec(g0, g1, g2, g3, pref=False, muref=False, stp=1e2, xmax=6.57, xmin=-4., xlist=False):

# DEFINE CORE PARAMETERS
	
	def G(x): # adiabatic index
	
		term0 = g0
		term1 = g1*x
		term2 = g2*x**2
		term3 = g3*x**3
	
		return np.exp(term0+term1+term2+term3) # here and below taken from Lindblom PRD 82 (2010)
		
	def z(x): # auxiliary variable (related to rho)
	
		integral, err = integrate.quad(lambda x1: 1./G(x1), 0., x)
	
		return np.exp(-integral)

	if pref == False: # reference pressure (where spectral EoS starts)
	
		pref = 0.5*rhonuc # if no reference pressure given, use half nuclear density
		
	if muref == False: # if no reference density given, calculate by extrapolating spectral fit to low pressures
			
		mupref, err = integrate.quad(lambda x1: np.exp(x1)*z(x1)/G(x1), xmin, 0.) 
		muref = mupref*pref # total energy density at reference pressure in g/cm^3

	def mu(x): # total energy density in g/cm^3
	
		integral, err = integrate.quad(lambda x1: np.exp(x1)*z(x1)/G(x1), 0., x)
	
		return (muref+pref*integral)/z(x)
		
	rhoref = muref # match rho and mu at reference pressure
		
	def rho(x): # rest-mass energy density in g/cm^3
	
		return rhoref/z(x)

# EXPORT EOS DATA

	rhodat = []
	pdat = []
	mudat = []
	
	if len(xlist) > 1:
	
		xdat = xlist
	
	else:
	
		xdat = np.linspace(xmin,xmax,stp) # x = log10(p/pref)
	
	for x in xdat:
	
		ppt = pref*np.exp(x)
		mupt = mu(x)
		rhopt = rho(x)
				
		pdat.append(ppt) # pressure in g/cm^3
		mudat.append(mupt) # total energy density in g/cm^3
		rhodat.append(rhopt) # baryon density in g/cm^3

	return rhodat, mudat, pdat
