	# -*- coding: iso-8859-1 -*

"""
This file defines the functions used to implement the formalism from Section 3.5 and Appendix 1 of Ranjan & Sasselov (2015)
"""

##************************************************************************************************************************####Global definitions

#Import globally useful libraries
import numpy as np
import scipy as sp
from scipy import interpolate as interp
import matplotlib.pyplot as plt
import scipy.integrate
import math as math
import pdb
import cookbook as cb #This is a user-defined library which contains Python cookbook codes useful to binning and rebinning.

#Define globally used constants
g=9.81e2 #cm/s2
k=1.3806e-16 #ergs/K
amu2g=1.6605e-24 #amu to g
bar2bayres=1.0e6 #bar to bayres; 1 bayres= 1 dyne/cm2
bayres2bars=1.0e-6 #bayres to bar
R_e=6.378e8 #radius of Earth in cm
hc=1.98645e-9 # value of h*c in erg-nm
erg2eV=6.24150934e11 #1 erg in eV
eV2erg=1.60217657e-12 #1 eV in erg

#Set atmospheric parameters assumed in our calculations
T=290. #Kelvin, atmosphere temperature
mu_amu=0.1*(12.+2*15.9994)+0.9*(2*14.0067)#mean molecular weight of atmosphere in amu. Corresponds to the 0.9 bar N2/0.1 bar CO2 atmosphere ad-hoc model.

#Derived atmospheric parameters
mu=mu_amu*amu2g #mean molecule weight of atmosphere in g
H=k*T/(mu*g) #scale height of atm in cm







##************************************************************************************************************************##
def ComputeHCNAbundance(coldensity_co2, z0, plotdiagnostics):
	"""	
	This function estimates the partial pressure of HCN at the planet surface (bottom of the atmosphere).
	
	Inputs:
	coldensity_co2: column density of CO2 shielding the HCN in cm**-2.
	z0: the height of the HCN layer, in cm. For z>z0 the HCN abundance is 0. 
	plotdiagnostics: if 1, shows diagnostic plots showing intermediate steps of calculation. If 0, does not.
	
	Output:
	Surface partial pressure of HCN, in bar.
	
	Assumptions:
	-1 bar Isothermal atmosphere (T=290) in hydrostatic equilibrium, with mean molecular weight corresponding to 0.9 bar N2 and 0.1 bar CO2. This reflects the ad-hoc model of Rugheimer et al (2015).
	-Assumes sole source of HCN is geochemical. Following Zahnle et al (1986), assumes 0.1 of the methane flux gets converted to HCN. We take the methane flux to be given by the Emmanuel & Ague (2007) estimate for abiotic methane flux to the atmosphere for the modern Earth.
	-Assumes sole sink of HCN is photolysis. Assumes all absorptions lead to photolysis based on Lee et al (1980)
	-Solar signal for photolysis calculation is attenuated by the CO2 column density taken as input. A layer-by-layer approach is not used.
	"""
	###########################################################################################################################################
	### Read in young Sun spectrum.
	importeddata=np.genfromtxt('./claire_3d9Ga.dat', skip_header=1, skip_footer=0)
	stellar_wav=importeddata[:,0] #nm
	stellar_flux=importeddata[:,1] #erg/s/cm**2/nm
	stellar_flux=stellar_flux*stellar_wav/hc #Convert stellar flux to units of photons/nm/s/cm^2.

	stellar_lefts, stellar_rights=cb.get_bin_edges(stellar_wav) #Extract wavelength bins of stellar spectrum.

	
	###Read in HCN absorption of Nuth et al (1982)
	importeddata=np.genfromtxt('./HCN_cxs.dat', skip_header=1, skip_footer=22)
	hcn_wav=np.reshape(importeddata, -1)/10. #convert from A to nm
	del importeddata
	importeddata=np.genfromtxt('./HCN_cxs.dat', skip_header=23, skip_footer=0)
	hcn_xc=np.reshape(importeddata, -1) #cm**2
	
	hcn_func=interp.interp1d(hcn_wav, hcn_xc, kind='linear') #This function provides the HCN cross-section in cm**2 as a function of wavelength in nm. It linearly interpolates the measurements of Nuth et al (1982) to do so. 
	
	
	###Read in CO2 absorption spectrum compiled by Huestis & Berkowitz (2010)
	importeddata=np.genfromtxt('./CO2_HuestisBerkowitz(2010)_300K_0.1254-201.6nm(evaluation).txt')
	co2_wav=importeddata[:,0] #nm
	co2_xc=importeddata[:,1] #cm2
	
	co2_func=interp.interp1d(co2_wav, co2_xc, kind='linear') #This function provides the CO2 cross-section in cm**2 as a function of wavelength in nm. It linearly interpolates the dataset of Huestis & Berkowitz (2010) to do so.


	########################################################################################################################################################
	###Select down the stellar input data we have to the wavelengths we have coverage for.
	min_lambda=np.min(hcn_wav)
	max_lambda=np.max(co2_wav)

	stellarwavinds=(stellar_lefts >= min_lambda) & (stellar_rights <= max_lambda)
	stellar_lefts=stellar_lefts[stellarwavinds]
	stellar_rights=stellar_rights[stellarwavinds]
	stellar_wav=stellar_wav[stellarwavinds]
	stellar_flux=stellar_flux[stellarwavinds]
	
	stellar_binwidths=stellar_rights-stellar_lefts #widths of each of the stellar spectrum bins.
	
	###Find mean absorption cross-section of CO2 and HCN in each of the stellar spectrum bins.
	co2_xc_binned=np.zeros(np.shape(stellar_wav))
	hcn_xc_binned=np.zeros(np.shape(stellar_wav))
	
	for ind in range(0, len(stellar_wav)):
		left=stellar_lefts[ind]
		right=stellar_rights[ind]
		
		co2_xc_binned[ind]=scipy.integrate.quad(co2_func, left, right)[0]/(right-left) 
		hcn_xc_binned[ind]=scipy.integrate.quad(hcn_func, left, right)[0]/(right-left)
		#returns cm2/molecule*nm/nm, average cross-section across that band
	
	
	#######################################################################################################################################################
	###Compute atmospheric attenuation due to CO2 & resultant photolysis rate
	#Compute optical depth of CO2
	tau_co2=co2_xc_binned*coldensity_co2

	#Compute emergent spectrum after CO2 shielding
	emergent_flux=stellar_flux*np.exp(-tau_co2)

	#Compute absorptions/s/nm/particle
	absrateperparticleperwav=emergent_flux*hcn_xc_binned #units: s**-1 nm**-1

	#Compute absorptions/s/particle
	B=np.sum(absrateperparticleperwav*stellar_binwidths)

	###Compare photolysis sink to geochemical source.
	#Methane input from geochemical sources, based on work of Emmanuel & Ague (2007)
	S0_ch4=2.7e27 #s**-1

	#Following Zahnle+1986, assume 0.1 of CH4 particles are converted to HCN. 
	S0_hcn=0.1*S0_ch4 #s**-1

	#Compute HCN surface pressure 
	hcn_n0=(S0_hcn)/(np.pi*R_e**2.0*(B)*H*(1.0-np.exp(-z0/H))) #number density
	hcn_P0=hcn_n0*k*T*bayres2bars #surface pressure of Convert to bars

	#Compute optical depth of the HCN.
	tau_hcn=hcn_xc_binned*hcn_n0*H*(1.0-np.exp(-z0/H))


	##########################################################################################################################################################
	###Diagnostic plots for intermediates of computation
	if plotdiagnostics==1:
		###Plot absorption cross-sections and check that finding their average values has been done properly
		fig1, (ax1,ax2)=plt.subplots(2, figsize=(8.0, 6.0))

		ax1.plot(hcn_wav,hcn_xc, color='red',marker='.', label='Data')
		ax1.bar(stellar_lefts, hcn_xc_binned, width=stellar_binwidths, color='maroon', label='Binned')
		ax1.set_yscale('log')
		ax1.legend(loc=0)
		ax1.set_xlim([np.min(stellar_lefts),np.max(stellar_rights)])
		ax1.set_ylabel('HCN Absorption (cm$^2$)', fontsize=14)
		ax1.set_xlabel('Wavelength (nm)')
		ax1.yaxis.grid(True)
		ax1.xaxis.grid(True)

		ax2.plot(co2_wav, co2_xc, color='green',marker='.', label='Data')
		ax2.bar(stellar_lefts, co2_xc_binned, width=stellar_binwidths, color='blue', label='Binned')
		ax2.set_yscale('log')
		#ax2.set_ylim([1.e-21, 1.e-15])
		ax2.set_xlim(ax1.get_xlim())
		ax2.legend(loc=0)
		ax2.set_ylabel('CO2 Absorption (cm$^2$)', fontsize=14)
		ax2.set_xlabel('Wavelength (nm)')
		ax2.yaxis.grid(True)
		ax2.xaxis.grid(True)
				
		###Plot attenuation of stellar flux.
		fig2, (ax1,ax2, ax3)=plt.subplots(3, figsize=(8.0, 11.0))

		ax1.plot(stellar_wav,stellar_flux, color='red',marker='.', label='Stellar Input')
		ax1.bar(stellar_lefts, stellar_flux, width=stellar_binwidths, color='red',label='Binned', alpha=0.5)
		ax1.set_yscale('log')
		ax1.legend(loc=0) #4, prop={'size':10.5}
		ax1.set_ylim([1.e-4*np.max(stellar_flux), 1.e2*np.max(stellar_flux)])
		ax1.set_xlim([np.min(stellar_lefts), np.max(stellar_rights)])
		ax1.set_title('Incident Stellar Flux', fontsize=14)
		ax1.set_xlabel('Wavelength (nm)', fontsize=14)
		ax1.set_ylabel('Flux (photons/s/nm/cm2)', fontsize=14)
		ax1.yaxis.grid(True)
		ax1.xaxis.grid(True)
		
		ax2.bar(stellar_lefts, emergent_flux, width=stellar_binwidths, color='blue',label='After CO2 Shielding', alpha=0.5)
		ax2.set_yscale('log')
		ax2.legend(loc=0) #4, prop={'size':10.5}
		ax2.set_ylim([1.e-68*np.max(emergent_flux), 1.e2*np.max(emergent_flux)])
		ax2.set_title('Emergent Flux (Post CO2 Shielding)', fontsize=14)
		ax2.set_xlabel('Wavelength (nm)', fontsize=14)
		ax2.set_ylabel('Flux (photons/s/nm/cm2)', fontsize=14)
		ax2.yaxis.grid(True)
		ax2.xaxis.grid(True)
		ax2.set_xlim(ax1.get_xlim())

		
		ax3.bar(stellar_lefts, tau_co2, width=stellar_binwidths, color='red', label='Optical Depth of CO2', alpha=0.5)
		ax3.bar(stellar_lefts, tau_hcn, width=stellar_binwidths,color='blue', label='Optical Depth of CH4', alpha=0.5)
		ax3.set_yscale('log')
		ax3.legend(loc=0) #4, prop={'size':10.5}
		ax3.set_title('Optical Depths', fontsize=14)
		ax3.set_xlabel('Wavelength (nm)', fontsize=14)
		ax3.set_ylabel('Optical Depth', fontsize=14)
		ax3.yaxis.grid(True)
		ax3.xaxis.grid(True)
		ax3.set_xlim(ax1.get_xlim())

		plt.show()

	return hcn_P0






###************************************************************************************************************************##
def ComputeCH4Abundance(coldensity_co2, z0, plotdiagnostics):
	"""	
	This function estimates the partial pressure of CH4 at the planet surface (bottom of the atmosphere).
	
	Inputs:
	coldensity_co2: column density of CO2 shielding the CH4 in cm**-2.
	z0: the height of the methane layer, in cm. For z>z0 the CH4 abundance is 0. 
	plotdiagnostics: if 1, shows diagnostic plots showing intermediate steps of calculation. If 0, does not.
	
	Output:
	Surface partial pressure of CH4, in bar.
	
	Assumptions:
	-1 bar Isothermal atmosphere (T=290) in hydrostatic equilibrium, with mean molecular weight corresponding to 0.9 bar N2 and 0.1 bar CO2. This reflects the ad-hoc model of Rugheimer et al (2015).
	-Assumes sole source of CH4 is geochemical. We take the methane flux to be given by the Emmanuel & Ague (2007) estimate for abiotic methane flux to the atmosphere for the modern Earth.
	-Assumes sole sink of CH4 is photolysis. Take all absorptions to lead to photolysis at these wavelengths
	-Solar signal for photolysis calculation is attenuated by the CO2 column density taken as input. A layer-by-layer approach is not used.
	"""
	
	###############################################################################################################
	### Read in Ranjan/Claire Young Earth model
	importeddata=np.genfromtxt('./claire_3d9Ga.dat', skip_header=1, skip_footer=0)
	stellar_wav=importeddata[:,0] #nm
	stellar_flux=importeddata[:,1] #erg/s/cm**2/nm
	stellar_flux=stellar_flux*stellar_wav/hc #Convert stellar flux to units of photons/nm/s/cm^2.

	stellar_lefts, stellar_rights=cb.get_bin_edges(stellar_wav) #Extract wavelength bins of stellar spectrum.
	
	###Read in Methane Absorption from Au et al (1993)
	ch4file2='./CH4_Au(1993)_298K_5.6-165nm(e,e).txt' #Some key differences with more modern dataset
	del importeddata
	importeddata=np.genfromtxt(ch4file2)
	ch4_wav=importeddata[:,0] #nm
	ch4_xc=importeddata[:,1] #cm**2
		
	ch4_func=interp.interp1d(ch4_wav, ch4_xc, kind='linear') #returns ch4 xc in cm2 as a function of nm via linear interpolation
	
	###Read in CO2 Absorption
	co2file2='./CO2_HuestisBerkowitz(2010)_300K_0.1254-201.6nm(evaluation).txt' #Consistent with earlier dataset, larger
	del importeddata
	importeddata=np.genfromtxt(co2file2)
	co2_wav=importeddata[:,0] #nm
	co2_xc=importeddata[:,1] #cm2
	
	co2_func=interp.interp1d(co2_wav, co2_xc, kind='linear') #returns co2 xc in cm2 as a function of nm via linear interpolation
	
	
	
	###############################################################################################################
	###Select down the stellar input data we have to the wavelengths we have coverage for.
	min_lambda=np.min(ch4_wav)
	max_lambda=np.max(ch4_wav)

	stellarwavinds=(stellar_lefts >= min_lambda) & (stellar_rights <= max_lambda)
	stellar_lefts=stellar_lefts[stellarwavinds]
	stellar_rights=stellar_rights[stellarwavinds]
	stellar_wav=stellar_wav[stellarwavinds]
	stellar_flux=stellar_flux[stellarwavinds]
	
	stellar_binwidths=stellar_rights-stellar_lefts #widths of the stellar spectrum bins
	
	###Find mean absorption cross-section of CO2 and CH4 in each of the stellar spectrum bins.
	co2_xc_binned=np.zeros(np.shape(stellar_wav))
	ch4_xc_binned=np.zeros(np.shape(stellar_wav))
	
	for ind in range(0, len(stellar_wav)):
		left=stellar_lefts[ind]
		right=stellar_rights[ind]
		
		co2_xc_binned[ind]=scipy.integrate.quad(co2_func, left, right)[0]/(right-left) 
		ch4_xc_binned[ind]=scipy.integrate.quad(ch4_func, left, right)[0]/(right-left)
		#returns cm2/molecule*nm/nm, average cross-section across that band
	
	
	###############################################################################################################
	###Computer atmospheric attenuation due to CO2 & resultant photolysis rate

	###Compute optical depth of CO2
	tau_co2=co2_xc_binned*coldensity_co2

	###Compute emergent spectrum after CO2 shielding
	emergent_flux=stellar_flux*np.exp(-tau_co2)

	###Compute absorptions/s/nm/particle
	absrateperparticleperwav=emergent_flux*ch4_xc_binned #units: s**-1 nm**-1

	###Compute absorptions/s/particle	
	B=np.sum(absrateperparticleperwav*stellar_binwidths) #integrate over nm

	###Compare photolysis sink to geochemical source
	#Methane input from geochemical sources, based on work of Emmanuel & Ague (2007)
	S0_g=7.3e4 #g/s, from literature
	S0=2.7e27 #s^-1, calculated
	
	#Compute methane surface pressure
	ch4_n0=(S0)/(np.pi*R_e**2.0*(B)*H*(1.0-np.exp(-z0/H)))
	ch4_P0=ch4_n0*k*T*bayres2bars #Convert to bars

	#Compute optical depth of methane
	tau_ch4=ch4_xc_binned*ch4_n0*H*(1.0-np.exp(-z0/H))
	
	
	##############################################################################################################
	###Plot absorption cross-sections and check that finding their average values has been done properly

	if plotdiagnostics==1:
		fig1, (ax1,ax2)=plt.subplots(2, figsize=(8.0, 6.0))

		ax1.plot(ch4_wav,ch4_xc, color='red',marker='.', label='Data')
		ax1.bar(stellar_lefts, ch4_xc_binned, width=stellar_binwidths, color='maroon', label='Binned')
		ax1.set_yscale('log')
		ax1.legend(loc=0)
		ax1.set_xlim([np.min(stellar_lefts),np.max(stellar_rights)])
		ax1.set_ylabel('CH4 Absorption (cm$^2$)', fontsize=14)
		#ax1.set_xlabel('Wavelength (nm)')
		ax1.yaxis.grid(True)
		ax1.xaxis.grid(True)

		ax2.plot(co2_wav, co2_xc, color='green',marker='.', label='Data')
		ax2.bar(stellar_lefts, co2_xc_binned, width=stellar_binwidths, color='blue', label='Binned')
		ax2.set_yscale('log')
		ax2.set_ylim([1.e-21, 1.e-15])
		ax2.set_xlim(ax1.get_xlim())
		ax2.legend(loc=0)
		ax2.set_ylabel('CO2 Absorption (cm$^2$)', fontsize=14)
		ax2.set_xlabel('Wavelength (nm)')
		ax2.yaxis.grid(True)
		ax2.xaxis.grid(True)
				
		###Plot attenuation of stellar flux.				
		fig2, (ax1,ax2, ax3)=plt.subplots(3, figsize=(8.0, 11.0))

		ax1.plot(stellar_wav,stellar_flux, color='red',marker='.', label='Stellar Input')
		ax1.bar(stellar_lefts, stellar_flux, width=stellar_binwidths, color='red',label='Binned', alpha=0.5)
		ax1.set_yscale('log')
		ax1.legend(loc=0)
		ax1.set_ylim([1.e-4*np.max(stellar_flux), 1.e2*np.max(stellar_flux)])
		ax1.set_xlim([np.min(stellar_lefts), np.max(stellar_rights)])
		ax1.set_title('Incident Stellar Flux', fontsize=14)
		#ax1.set_xlabel('Wavelength (nm)', fontsize=14)
		ax1.set_ylabel('Flux (photons/s/nm/cm2)', fontsize=14)
		ax1.yaxis.grid(True)
		ax1.xaxis.grid(True)
		
		ax2.bar(stellar_lefts, emergent_flux, width=stellar_binwidths, color='blue',label='After CO2 Shielding', alpha=0.5)
		ax2.set_yscale('log')
		ax2.legend(loc=0)
		ax2.set_ylim([1.e-68*np.max(emergent_flux), 1.e2*np.max(emergent_flux)])
		ax2.set_title('Emergent Flux (Post CO2 Shielding)', fontsize=14)
		#ax2.set_xlabel('Wavelength (nm)', fontsize=14)
		ax2.set_ylabel('Flux (photons/s/nm/cm2)', fontsize=14)
		ax2.yaxis.grid(True)
		ax2.xaxis.grid(True)
		ax2.set_xlim(ax1.get_xlim())

		
		ax3.bar(stellar_lefts, tau_co2, width=stellar_binwidths, color='red', label='Optical Depth of CO2', alpha=0.5)
		ax3.bar(stellar_lefts, tau_ch4, width=stellar_binwidths,color='blue', label='Optical Depth of CH4', alpha=0.5)
		ax3.set_yscale('log')
		ax3.legend(loc=0)
		ax3.set_title('Optical Depths', fontsize=14)
		ax3.set_xlabel('Wavelength (nm)', fontsize=14)
		ax3.set_ylabel('Optical Depth', fontsize=14)
		ax3.yaxis.grid(True)
		ax3.xaxis.grid(True)
		ax3.set_xlim(ax1.get_xlim())
		plt.show()

	return ch4_P0






##************************************************************************************************************************##
def getCO2coldensity_CH4(ch4_surfpressure_target, z0, guess):
	"""
	This function computes the CO2 column density required to permit the buildup of CH4 to a given level
	
	Inputs:
	-Desired CH4 surface pressure, in bar.
	-Height z0 of the CH4 layer, in cm. For z>z0, the CH4 mixing ratio is 0.
	-guess: an initial guess for the CO2 column density, in cm**-2
	
	Output:
	-The CO2 column density, in cm**-2.
	
	Note:
	-The Python solver is very sensitive to the initial guess. The initial guess needs to be close to the true value. 
	"""
	import scipy.optimize
	def func(colden):
		return ComputeCH4Abundance(colden, z0, 0)-ch4_surfpressure_target
	
	return scipy.optimize.newton(func, guess)


##************************************************************************************************************************##


def getCO2coldensity_HCN(hcn_surfpressure_target, z0, guess):
	"""
	This function computes the CO2 column density required to permit the buildup of HCN to a given level
	
	Inputs:
	-Desired HCN surface pressure, in bar.
	-Height z0 of the HCN layer, in cm. For z>z0, the HCN mixing ratio is 0.
	-guess: an initial guess for the CO2 column density, in cm**-2
	
	Output:
	-The CO2 column density, in cm**-2.
	
	Note:
	-The Python solver is very sensitive to the initial guess. The initial guess needs to be close to the true value. 
	"""
	import scipy.optimize
	def func(colden):
		return ComputeHCNAbundance(colden, z0, 0)-hcn_surfpressure_target
	
	return scipy.optimize.newton(func, guess)



##************************************************************************************************************************##
def coldensity2pressure(coldensity):
	"""
	This function computes the gas partial pressure corresponding to a given gas column density by integrating up the atmospheric column to infinity from that gas partial pressure. It assumes, as before, an isothermal atmosphere in hydrostatic equilibrium.
	
	Note: This function can be used in conjunction with getCO2coldensity_XXX to estimate the partial pressure of CO2 required to provide the column density of CO2 required to permit a given feedstock gas to build up to a desired abundance for z<z0. However, this implicitly assumes all the CO2 attenuation of the incoming solar flux is occuring at z>z0; CO2 shielding for z<z0 is neglected. Hence, the CO2 partial pressures extracted using this method should be regarded as upper bounds. 
	
	Input: 
	-column density, in cm**-2
	
	Output:
	-Pressure in bar.
	"""
	pressure=(k*T*coldensity/H)*bayres2bars
	return pressure

##************************************************************************************************************************##
def getsurfacepressure(Pz, z):
	"""
	Convert the pressure of a given gas (e.g. CO2) at an altitude z to the pressure at the surface. Assumes isothermal atmosphere in hydrostatic equilibrium.
	
	Inputs:
	-Partial pressure at altitude
	-altitude in cm. 
	
	Outputs:
	-Surface partial pressure
	"""
	surfacepressure=Pz*np.exp(z/H)
	return surfacepressure




###************************************************************************************************************************##
##Test cases, to ensure code is running properly. Uncomment to run.

#print 'CH4'
#print ComputeCH4Abundance(2.59e20, 50.e5, 1) #should be ~1.e-6
#y=getCO2coldensity_CH4(1., 17.e5, 7.40e20) #should be equal to 7.41e20
#print y
#z=coldensity2pressure(y) #should be qual to 3.57e-5
#print z
#print getsurfacepressure(z, 17.e5) #should be equal to 2.77e-4

#print 'HCN'
#print ComputeHCNAbundance(4.78e22, 17.e5, 0) #should be ~1.e-6
#y=getCO2coldensity_HCN(1., 1.e5, 2.0e23) #should be equal to 2.05e23
#print y
#z=coldensity2pressure(y) #should be equal to 9.90e-3
#print z
#print getsurfacepressure(z, 1.e5) #should be equal to 1.12e-2
