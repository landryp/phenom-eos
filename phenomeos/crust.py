#!/usr/bin/python

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
from .constants import rhonuc

def crust(eos,pts=1e2,srange=[3.8e11,3.8e14],rhoi=0.28,rholist='0'):

	# LOAD EOS DATA

	eosdatx = eos[0]/rhonuc # rho in units of rhonuc
	eosdaty = eos[1]/rhonuc # p in units of rhonuc
	eosintp = interp1d(eosdatx,eosdaty,kind='linear',bounds_error=False,fill_value=0.)
	
	def p(rho):
	
		return eosintp(rho)
		
	eosdaty2 = eos[2]/rhonuc # mu in units of rhonuc
	eosintp2 = interp1d(eosdatx,eosdaty2,kind='linear',bounds_error=False,fill_value=0.)
	
	def mu(rho):
	
		return eosintp2(rho)
		
	# BUILD SLY CRUST
	
	K = [0,6.80110e-9,1.06186e-6,5.32697e1,3.99874e-8]
	Gamma = [1,1.58425,1.28733,0.62223,1.35692]
	rhodiv = [0,2.44034e7,3.78358e11,2.62780e12,1e20]
	a = np.zeros(5)
	a[1] = 0
	a[2] = (rhodiv[1]+K[1]*rhodiv[1]**Gamma[1]/(Gamma[1]-1.)+a[1]*rhodiv[1])/rhodiv[1]-1.-K[1]*rhodiv[1]**(Gamma[1]-1.)/(Gamma[1]-1)
	a[3] = (rhodiv[2]+K[2]*rhodiv[2]**Gamma[2]/(Gamma[2]-1.)+a[2]*rhodiv[2])/rhodiv[2]-1.-K[2]*rhodiv[2]**(Gamma[2]-1.)/(Gamma[2]-1)
	a[4] = (rhodiv[3]+K[3]*rhodiv[3]**Gamma[3]/(Gamma[3]-1.)+a[3]*rhodiv[3])/rhodiv[3]-1.-K[3]*rhodiv[3]**(Gamma[3]-1.)/(Gamma[3]-1)
	
	for i in range(2,5): # K values in Read+ table are rounded, so crust EoS slightly discontinuous unless corrected
	
		K[i] = K[i-1]*rhodiv[i-1]**Gamma[i-1]/(rhodiv[i-1]**Gamma[i])

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

	# FIND CORE-CRUST INTERFACE
	
	def corecrust(p,pSLy):
	
		def intersect(rho):
		
			return abs(p(rho)-pSLy(rho))
	
		res = minimize_scalar(intersect,bounds=(srange[0]/rhonuc,srange[1]/rhonuc),method='bounded')
		rhocr = res.x
		residual = res.fun
	
		return rhocr, residual

	rhocr, residual = corecrust(p,pSLy)

#	print 'Crust-core interface found at {0} rho_nuc, with absolute pressure residual of {1} rho_nuc'.format(rhocr,residual)
	
	# BUILD UNIFIED EOS
	
	def puni(rho):
	
		if isinstance(rho,list) or isinstance(rho,np.ndarray):	
	
			rho = np.asarray(rho)
		
		else:
	
			rho = np.asarray([rho])
	
		plist = [pSLy(rho),p(rho)]
	
		return np.piecewise(rho,[rho < rhocr, rho >= rhocr],[lambda rho=rho, i=j: plist[i] for j in np.arange(0,2)])
		
	def muuni(rho):
	
		if isinstance(rho,list) or isinstance(rho,np.ndarray):	
	
			rho = np.asarray(rho)
		
		else:
	
			rho = np.asarray([rho])
	
		mulist = [muSLy(rho),mu(rho)]
	
		return np.piecewise(rho,[rho < rhocr, rho >= rhocr],[lambda rho=rho, i=j: mulist[i] for j in np.arange(0,2)])

	if len(rholist) > 1:
		rholist = [rhopt/rhonuc for rhopt in rholist if rhopt/rhonuc <= rhocr]
		pts = len(rholist)
	else:
		rholist = np.logspace(np.log10(rhoi/rhonuc),np.log10(rhocr),pts)

	mudat = np.zeros(pts)
	pdat = np.zeros(pts)
	
	for i in np.arange(pts):
	
		rho = rholist[i]
		
		mudat[i] = muuni(rho)
		pdat[i] = puni(rho)
	
	return [rholist,mudat,pdat]
