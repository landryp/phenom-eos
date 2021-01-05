#!/usr/bin/python

#__doc__ = 'Constructs a spectral EoS from its parameters'
#__usage__ = 'get_spec.py g0 g1 g2 g3'
#__author__ = 'philippe.landry@ligo.org'
#__date__ = '02-2019'

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interp
from .constants import *

# DEFINE CORE PARAMETERS

def G(x, g0, g1, g2, g3): # adiabatic index
	
	term0 = g0
	term1 = g1*x
	term2 = g2*x**2
	term3 = g3*x**3
	
	return np.exp(term0+term1+term2+term3) # here and below taken from Lindblom PRD 82 (2010)
	
def z(x,g0,g1,g2,g3): # auxiliary variable (related to rho)
	
	integral, err = integrate.quad(lambda x1: 1./G(x1,g0,g1,g2,g3), 0., x)
	
	return np.exp(-integral)
	
def mu(x,g0,g1,g2,g3,muref,pref): # total energy density in g/cm^3
	
	integral, err = integrate.quad(lambda x1: np.exp(x1)*z(x1,g0,g1,g2,g3)/G(x1,g0,g1,g2,g3), 0., x)
	
	return (muref+pref*integral)/z(x,g0,g1,g2,g3)
	
def rho(x,g0,g1,g2,g3,rhoref): # rest-mass energy density in g/cm^3
	
	return rhoref/z(x,g0,g1,g2,g3)

def spec(g0, g1, g2, g3, pref=False, muref=False, stp=1e2, xmax=18., xmin=-4., xlist=False):

	if pref == False: # reference pressure (where spectral EoS starts)
	
		pref = 6e11 # if no reference pressure given, use value from Carney+ 2018
		
	if muref == False: # if no reference density given, get from sly eos's total energy density at pref

		muref = 1.3e14 # from sly eos
	
	rhoref = muref # match rho and mu at reference pressure

# EXPORT EOS DATA

	rhodat = []
	pdat = []
	mudat = []
	
	if len(xlist) > 1: xdat = xlist
	else: xdat = np.linspace(xmin,xmax,stp) # x = ln(p/pref)
	
	for count,x in enumerate(xdat):
	
		ppt = pref*np.exp(x)
		mupt = mu(x,g0,g1,g2,g3,muref,pref)
		rhopt = rho(x,g0,g1,g2,g3,rhoref)
		if count > 1 and (ppt - pdat[count-1])/(mupt - mudat[count-1]) > c**2: break
		if rhopt > 20*rhonuc: break
				
		pdat.append(ppt) # pressure in g/cm^3
		mudat.append(mupt) # total energy density in g/cm^3
		rhodat.append(rhopt) # baryon density in g/cm^3

	return rhodat, mudat, pdat
	return rhodat, mudat, pdat
