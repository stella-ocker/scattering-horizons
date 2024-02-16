from numpy import *
from matplotlib.pyplot import *
#from astropy.cosmology import Planck15 as cosmo
from astropy import constants as const
from astropy import units as u
from astropy.coordinates import Distance
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from colossus.halo import mass_so

cosmo1 = cosmology.setCosmology('planck18')

def dH(z):
    c = const.c 
    H0 = cosmo1.Hz(0.)
    numer = c.to('km/s')/H0
    denom = sqrt((cosmo1.Om0*(1.+z)**3.)+cosmo1.Ode0)
    return (numer/denom).value

def nhalo(logmass1,logmass2,z):
	''' Number density of halos for given redshift, mass bin '''
	mass_bin = 10.**linspace(logmass1,logmass2,10)
	hmf = mass_function.massFunction(mass_bin,z,q_in='M', q_out='dndlnM', mdef = '200c', model = 'tinker08')
	# units of hmf are (Mpc/h)^{-3}, but seems like integrating over mass in M_sun/h would change this??
	density = trapz(hmf,log(mass_bin))
	return density

def intersection_number(zlow,zhigh,logmass1,logmass2):

	#Hz = cosmo1.Hz(zhigh)
	H0 = cosmo1.Hz(0.)
	h = H0/100. # H0 in units of 100 km/s/Mpc
	c = const.c # speed of light in m/s

	mean_mass = (logmass1+logmass2)/2.

	z_arr = linspace(zlow,zhigh,200)
	nhalo_arr = []
	r200_arr = []
    
	for j in range(len(z_arr)):
		z_dens = nhalo(logmass1,logmass2,z_arr[j]) # not sure about units
		r200 = mass_so.M_to_R(10.**mean_mass,z_arr[j],'200c') # kpc/h units
		nhalo_arr.append(z_dens)
		r200_arr.append(r200*10.**(-3.)) # Mpc/h units
	nhalo_arr = array(nhalo_arr)
	r200_arr = array(r200_arr)

    #rdisk = r200_arr*(17./0.236)*10.**(-3.)

	int_halo = (pi*(2.*r200_arr)**2.) * nhalo_arr * dH(z_arr) / (1.+z_arr)
	Nval_halo = trapz(int_halo,z_arr)

	return Nval_halo


