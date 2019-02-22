#!/usr/bin/python

#__doc__ = 'Constructs a piecewise polytrope EoS from its parameters'
#__usage__ = 'get_pwp.py log10p1 Gamma1 Gamma2 Gamma3'
#__author__ = 'philippe.landry@ligo.org'
#__date__ = '02-2019'

import numpy as np
from constants import rhonuc, G, c, Msun
rho0 = rhonuc

def pwpoly(logp1, G1, G2, G3, stp=1e2, rhomax=1e1, eps=1e-4):

# DEFINE SLY CRUST PARAMETERS

	rhoSLy = np.array([2.44033979e7,3.78358138e11,2.62780487e12])/rho0 # div densities
	GSLy = np.array([1.58424999,1.28732904,0.62223344,1.35692395]) # adiabatic indices
	KSLy = np.array([6.80109613e-9,1.06186086e-6,5.32696797e1,3.99873692e-8])*rho0**(GSLy-1.) # prop consts
	
# DEFINE CORE PARAMETERS

	rho1 = 10.**14.7/rho0 # core div densities
	rho2 = 10.**15./rho0

# LOCATE CRUST-CORE INTERFACE

	p1 = 10.**logp1

	p1a = p1/(rho0*c**2.) # dimensionless p1
	rhocr = (p1a/(KSLy[3]*rho1**G1))**(1./(GSLy[3]-G1)) # crust-core interface density
	
	rhoSLy = np.append(rhoSLy,rhocr)

# CONSTRUCT PRESSURE-DENSITY PROFILE

	rhoUni = np.array([rhoSLy[0],rhoSLy[1],rhoSLy[2],rhocr,rho1,rho2,1e18]) # unified eos div densities
	GUni = np.array([GSLy[0],GSLy[1],GSLy[2],GSLy[3],G1,G2,G3]) # adiabatic indices of unified eos
	KUni = np.array([KSLy[0],KSLy[1],KSLy[2],KSLy[3],p1a/rho1**G1,p1a/rho1**G2,p1a*((rho2/rho1)**G2)/rho2**G3]) # prop consts of unified eos

	def hi(i,rho): # enthalpies of each segment

		if i==0:

			return rho + KUni[0]/(GUni[0]-1.)*rho**GUni[0]

		else:

			return (1. + hi(i-1,rhoUni[i-1])/rhoUni[i-1]-1.-KUni[i]/(GUni[i]-1.)*rhoUni[i-1]**(GUni[i]-1.))*rho + KUni[i]/(GUni[i]-1.)*rho**GUni[i]

	aUni = np.zeros(7) # integration constant related to internal energy in each segment
	HUni = np.zeros(7) # rescaled enthalpy in each segment

	for j in range(1,7):

		aUni[j] = hi(j-1,rhoUni[j-1])/rhoUni[j-1]-1.-KUni[j]/(GUni[j]-1.)*rhoUni[j-1]**(GUni[j]-1.)
		HUni[j] = aUni[j] + GUni[j]/(GUni[j]-1.)*KUni[j]*rhoUni[j]**(GUni[j]-1.)

	HUni[0] = aUni[0] + GUni[0]/(GUni[0]-1.)*KUni[0]*rhoUni[0]**(GUni[0]-1.)

	if rhomax > rho2: # choose core parameters based on location of rhomax rel to div densities

		aCor = aUni[-1]
		GCor = GUni[-1]
		KCor = KUni[-1]

	elif (rho1 < rhomax) and (rhomax <= rho2):

		aCor = aUni[-2]
		GCor = GUni[-2]
		KCor = KUni[-2]

	elif rhomax <= rho1 :

		aCor = aUni[-3]
		GCor = GUni[-3]
		KCor = KUni[-3]

	hc = 1. + aCor + (GCor/(GCor-1.))*KCor*rhomax**(GCor-1.) # central values of fluid vars
	muc = (1. + aCor)*rhomax + (1./(GCor-1.))*KCor*rhomax**GCor
	pc = KCor*(((hc-1.)-aCor)/(KCor*(GCor/(GCor-1.))))**(GCor/(GCor-1.))

	thetaUni = np.zeros(7) # dividing rescaled enthalpies

	for j in range(0,7):

		thetaUni[j] = HUni[j]/(hc-1.)

	def pUni(i,theta): # pressure, mass density, tot energy density in each segment

		return KUni[i]*((theta*(hc-1.)-aUni[i])/(KUni[i]*(GUni[i]/(GUni[i]-1.))))**(GUni[i]/(GUni[i]-1.))

	def rho1Uni(i,theta):

		return ((theta*(hc-1.)-aUni[i])/(KUni[i]*(GUni[i]/(GUni[i]-1.))))**(1./(GUni[i]-1.))

	def muUni(i,theta):

		return rho1Uni(i,theta)*(1.+(aUni[i]+theta*(hc-1.)/(GUni[i]-1.))/(GUni[i]/(GUni[i]-1.)))

	def p(theta): # unified pressure function (of theta)

		return np.piecewise(theta,[theta <= thetaUni[0], (theta > thetaUni[0]) & (theta <= thetaUni[1]), (theta > thetaUni[1]) & (theta <= thetaUni[2]), (theta > thetaUni[2]) & (theta <= thetaUni[3]), (theta > thetaUni[3]) & (theta <= thetaUni[4]), (theta > thetaUni[4]) & (theta <= thetaUni[5]), (theta > thetaUni[5]) & (theta <= thetaUni[6]) ],[lambda theta=theta, i=j: pUni(i,theta) for j in np.arange(0,7)])

	def mu(theta): # unified tot energy density function (of theta)

		return np.piecewise(theta,[theta <= thetaUni[0], (theta > thetaUni[0]) & (theta <= thetaUni[1]), (theta > thetaUni[1]) & (theta <= thetaUni[2]), (theta > thetaUni[2]) & (theta <= thetaUni[3]), (theta > thetaUni[3]) & (theta <= thetaUni[4]), (theta > thetaUni[4]) & (theta <= thetaUni[5]), (theta > thetaUni[5]) & (theta <= thetaUni[6]) ],[lambda theta=theta, i=j: muUni(i,theta) for j in np.arange(0,7)])
		
		
	def rho(theta): # unified tot energy density function (of theta)

		return np.piecewise(theta,[theta <= thetaUni[0], (theta > thetaUni[0]) & (theta <= thetaUni[1]), (theta > thetaUni[1]) & (theta <= thetaUni[2]), (theta > thetaUni[2]) & (theta <= thetaUni[3]), (theta > thetaUni[3]) & (theta <= thetaUni[4]), (theta > thetaUni[4]) & (theta <= thetaUni[5]), (theta > thetaUni[5]) & (theta <= thetaUni[6]) ],[lambda theta=theta, i=j: rho1Uni(i,theta) for j in np.arange(0,7)])

# EXPORT EOS DATA

	thetadat = np.logspace(np.log10(0.+eps),np.log10(1.-eps),stp)
	rhodat = rho(thetadat)*rho0 # baryon density in g/cm^3
	pdat = p(thetadat)*rho0 # pressure in g/cm^3
	mudat = mu(thetadat)*rho0 # total energy density in g/cm^3

	return rhodat, mudat, pdat
	
