#!/usr/bin/python

#__doc__ = 'Constructs a constant sound-speed EoS from its parameters'
#__usage__ = 'get_css.py ptr Deltae cqm2'
#__author__ = 'philippe.landry@ligo.org'
#__date__ = '02-2019'

import numpy as np
from constants import *

def css(params, stp=1e2, pmax=1e17, plist=False):

# DEFINE CORE PARAMETERS

	[ptr, mutr, Dmu, cqm2] = params
	
	def mu(p): # total energy density in [g/cm^3]
	
		return mutr + Dmu*mutr + (p-ptr)/cqm2	# taken from Han+Steiner 1810.10967
		
	def rho(p): # rest-mass energy density in [g/cm^3]
	
		factor1 = -ptr+cqm2*(Dmu*mutr+mutr)
		factor2 = (1.+cqm2)*p - ptr + cqm2*(Dmu*mutr+mutr)
		
		power1 = cqm2/(1.+cqm2)
		power2 = 1./(1.+cqm2)
	
		return factor1**power1*factor2**power2/cqm2 # follows from 1st law of thermodynamics

# EXPORT EOS DATA

	rhodat = []
	mudat = []
	
	if len(plist) > 1:
	
		pdat = plist
	
	else:
	
		pdat = np.logspace(np.log10(ptr),np.log10(pmax),stp) # pressure in g/cm^3
	
	for p in pdat:
	
		mupt = mu(p)
		rhopt = rho(p)
		
		mudat.append(mupt) # total energy density in g/cm^3
		rhodat.append(rhopt) # baryon density in g/cm^3

	return rhodat, mudat, pdat	
