from numpy import *
from matplotlib.pyplot import *
from astropy.cosmology import Planck15 as cosmo
import math

# random constants and conversion factors

kpc2cm = 3.086*(10.**(21.)) # kiloparsecs to centimeters
Msun = 1.99*(10.**(33.)) # mass of Sun in grams

# all of the functions we need

def Gscatt(dso,dlo,dls,leff): # d's in Gpc, leff in Mpc
	return (dls*dlo)/(leff*dso)

def virial_radius(Mhalo,z1): 
	kpc2cm = 3.086*(10.**(21.)) # kiloparsecs to centimeters
	# virial radius for a galaxy with total mass (Mhalo) at a redshift (z1)
	rho_c = cosmo.critical_density(z1) # critical density in grams cm^-3
	rho_c = rho_c.value
	rhoc = rho_c*(kpc2cm**3.) # need to convert cm^-3 to kpc^-3
	return (Mhalo/((4*pi/3)*200.*rhoc))**(1./3.)

def c200(Mhalo,z1):
	# concentration parameter (K_c in my paper)
	# depends on the halo mass and redshift
	Hz = cosmo.H(z1) # Hubble constant at redshift z1
	h = Hz.value/100 # Hubble constant in units of 100 km/s/Mpc
	return 4.67*(Mhalo/(10.**(14.) * h**(-1.) * Msun))**(-0.11)

def ne_halo(rint,Mtot,rvcm,conc,fb,mue,muh,mp,Om,Ob,r200): 

	# electron density as a function of distance rint from halo center

	a = 2
	y0 = 2

	y = linspace(0.1,conc,100) #linspace(0.1,7.7,100) # y = c(r/r200)
	massfac = (rvcm/conc)**3. * (4.*pi*y**2.) / (y**(1-a) * (y0+y)**(2+a))
	mint = trapz(massfac,x=y)
	rhoh = Mtot/mint

	r = rint
	ynew = conc*r/r200
	mnfwr = rhoh*fb*(Ob/Om)/((ynew**(1.-a)) * (y0+ynew)**(2.+a))
	ne_mnfw = mue*mnfwr/(mp*muh)

	necut = ne_mnfw*0.5*(1-tanh((r-2.*r200)/20.)) # rolls off mNFW profile at twice virial radius
    
	return necut

def Leff(r200,Rperp): # kpc
	# effective path length through halo at a given distance from the center (Rperp) and a given virial radius r200
	return 2.*sqrt((r200)**2. - (Rperp**2.)) # cuts the halo off at r200

def dm_halo(z1,zsource,Mhalo,Rperp):

	# input Mtot in solar mass units
	Mtot = Msun*Mhalo # total mass of halo

	Om = cosmo.Om(z1)
	Ob = cosmo.Ob(z1)

	r200 = virial_radius(Mtot,z1) 
	conc = c200(Mtot,z1)
	rvcm = r200*kpc2cm # r200 in cm
	fb = 0.75 # fraction of baryonic matter in halo that's not part of galactic disk
	mue = 1.167
	muh = 1.3
	mp = 1.67*10.**(-24.) # grams

	dso = cosmo.angular_diameter_distance(zsource).value * 10.**(-3.) # Gpc
	dlo = cosmo.angular_diameter_distance(z1).value * 10.**(-3.)
	dls = cosmo.angular_diameter_distance_z1z2(z1,zsource).value * (10.**-3.)

	leff = Leff(2.*r200,Rperp) # gives leff in kpc
	Gsc = Gscatt(dso,dlo,dls,leff*10.**(-3.)) # distances in Gpc, leff in Mpc -- only valid for distant galaxies

	y = linspace(0.,leff/2)
	rtr = sqrt(Rperp**2. + y**2.) # gives the actual distance between the FRB path thru halo & the halo center (this is actually just half of the path thru the sphere)
	ne_r = ne_halo(rtr,Mtot,rvcm,conc,fb,mue,muh,mp,Om,Ob,r200) # gives density profile along half of the FRB path thru the sphere
	dm_ih = 2.*trapz(ne_r,rtr*10.**3.) # total DM contribution along the line-of-sight thru halo (takes DM from half the halo and multiplies by two b/c symmetry)

	return dm_ih,Gsc,leff # return leff for some code and not others

# ALL YOU NEED TO DO IS CALL THIS FUNCTION BELOW
# z1 is the lens galaxy redshift, zsource is the FRB source redshift, Mhalo is the halo mass in solar mass units, 
# and Rperp is the distance between the FRB line-of-sight and the galaxy center in kpc

#dm_halo(z1,zsource,Mhalo,Rperp)




