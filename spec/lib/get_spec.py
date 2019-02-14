#!/usr/bin/python

#__doc__ = 'Constructs a spectral EoS from its parameters'
#__usage__ = 'get_pwp.py g1 g2 g3 g4'
#__author__ = 'philippe.landry@ligo.org'
#__date__ = '10-2018'

import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interp

def spec(g0, g1, g2, g3, stp=1e2, xmax=1e1, eps=1e-4, pts=10):

# DEFINE CONSTANTS

	c = 2.99792458e10 # speed of light in cgs
	G = 6.67408e-8 # Newton's constant in cgs
	Msun = 1.3271244e26/G # solar mass in cgs
	rhonuc = 2.7e14 # nuclear density in cgs

# BUILD SLY CRUST
	
	K = [0,6.80110e-9,1.06186e-6,5.32697e1,3.99874e-8]
	Gamma = [1,1.58425,1.28733,0.62223,1.35692]
	rhodiv = [0,2.44034e7,3.78358e11,2.62780e12,1e20]
	a = np.zeros(5)
	a[1] = 0
	a[2] = (rhodiv[1]+K[1]*rhodiv[1]**Gamma[1]/(Gamma[1]-1.)+a[1]*rhodiv[1])/rhodiv[1]-1.-K[1]*rhodiv[1]**(Gamma[1]-1.)/(Gamma[1]-1)
	a[3] = (rhodiv[2]+K[2]*rhodiv[2]**Gamma[2]/(Gamma[2]-1.)+a[2]*rhodiv[2])/rhodiv[2]-1.-K[2]*rhodiv[2]**(Gamma[2]-1.)/(Gamma[2]-1)
	a[4] = (rhodiv[3]+K[3]*rhodiv[3]**Gamma[3]/(Gamma[3]-1.)+a[3]*rhodiv[3])/rhodiv[3]-1.-K[3]*rhodiv[3]**(Gamma[3]-1.)/(Gamma[3]-1)

	def pSLy(rho):
	
		if isinstance(rho,list) or isinstance(rho,np.ndarray):	
	
			rho = np.asarray(rho)
		
		else:
	
			rho = np.asarray([rho])
	
		pwpoly = [K[1]*(rho*rhonuc)**Gamma[1]/rhonuc,K[2]*(rho*rhonuc)**Gamma[2]/rhonuc,K[3]*(rho*rhonuc)**Gamma[3]/rhonuc,K[4]*(rho*rhonuc)**Gamma[4]/rhonuc]
		
		conds = [(rhodiv[0]/rhonuc <= rho) & (rho < rhodiv[1]/rhonuc),(rhodiv[1]/rhonuc <= rho) & (rho < rhodiv[2]/rhonuc),(rhodiv[2]/rhonuc <= rho) & (rho < rhodiv[3]/rhonuc),(rhodiv[3]/rhonuc <= rho) & (rho < rhodiv[4]/rhonuc)]
				
		return np.piecewise(rho,conds,[lambda rho=rho, i=j: pwpoly[i] for j in np.arange(0,4)])
	
	def muSLy(rho):
	
		if isinstance(rho,list) or isinstance(rho,np.ndarray):	
	
			rho = np.asarray(rho)
		
		else:
	
			rho = np.asarray([rho])
	
		pwpolymu = [(K[1]*(rho*rhonuc)**Gamma[1]/(Gamma[1]-1.)+(1.+a[1])*rho*rhonuc)/rhonuc,(K[2]*(rho*rhonuc)**Gamma[2]/(Gamma[2]-1.)+(1.+a[2])*rho*rhonuc)/rhonuc,(K[3]*(rho*rhonuc)**Gamma[3]/(Gamma[3]-1.)+(1.+a[3])*rho*rhonuc)/rhonuc,(K[4]*(rho*rhonuc)**Gamma[4]/(Gamma[4]-1.)+(1.+a[4])*rho*rhonuc)/rhonuc]
	
		conds = [(rhodiv[0]/rhonuc <= rho) & (rho < rhodiv[1]/rhonuc),(rhodiv[1]/rhonuc <= rho) & (rho < rhodiv[2]/rhonuc),(rhodiv[2]/rhonuc <= rho) & (rho < rhodiv[3]/rhonuc),(rhodiv[3]/rhonuc <= rho) & (rho < rhodiv[4]/rhonuc)]
	
		return np.piecewise(rho,conds,[lambda rho=rho, i=j: pwpolymu[i] for j in np.arange(0,4)])
		
# FIX CRUST-CORE INTERFACE AND CALCULATE INTERFACE VALUES

	rho0 = 0.5 # interface density in units of rhonuc
	
	p0 = pSLy(rho0)[0]
	mu0 = muSLy(rho0)[0]

# DEFINE CORE PARAMETERS

	def G(x):
	
		term0 = g0
		term1 = g1*x
		term2 = g2*x**2
		term3 = g3*x**3
	
		return np.exp(term0+term1+term2+term3)
		
	def rho(x):
	
		integral, err = integrate.quad(lambda x1: 1./G(x1), 0., x)
	
		return rho0*np.exp(integral)
	
	def mu(x):
	
		const = mu0*rho(x)/rho0
		integral, err = integrate.quad(lambda x1: rho0*np.exp(x1)/(rho(x1)*G(x1)), 0., x)
	
		return const + p0*rho(x)*integral/rho0

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

