
from numpy import *
from matplotlib.pyplot import *
from os import getpid
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy import stats
from scipy import integrate
from scipy.interpolate import interp1d
from datetime import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import subprocess
from subprocess import Popen, call, PIPE
import string
from astropy import constants as const
from astropy import units as u
from astropy.coordinates import Distance
#from colossus.cosmology import cosmology
#cosmo = cosmology.setCosmology('planck18')
from astropy.cosmology import Planck15 as cosmo
from DMhalo import virial_radius

### evaluates tau and DM for one LOS through disk

Msun = 1.99*(10.**(33.)) # mass of Sun in grams

def density(nh,zh,z,rh,r):
    z_comp = (cosh(abs(z)/zh))**(-2.)
    r_comp = (cosh(abs(r)/rh))**(-2.)
    return nh*z_comp*r_comp

def density_rolloff(rh,r):
	return cosh(abs(r)/rh)**(-2.)

def elliptical_cutoff(ri,zi,rcut,zcut):
    return (ri/rcut)**2. + (zi/zcut)**2.

def lognorm(mun,sigman): # get mu, sigma of lognormal distributions for Ftilde
    mu = log(mun**2./sqrt(mun**2. + sigman**2.))
    sigma = sqrt(log(1.+(sigman/mun)**2.))
    #pdf = (exp(-(log(x) - mu)**2 / (2 * sigma**2)) / (x * sigma * sqrt(2 * pi)))
    return mu,sigma

def Ftilde_z(Ftilde_0,z):
	# Ftilde z dependence re-defined based on Madau and Dickinson (2014)
    num = (1.+z)**(2.7)
    denom = 1. + ((1.+z)/2.9)**(5.6)
    return Ftilde_0*num/denom

def tauDM(DMl,ftilde,dls,dlo,dso,nu,zl,L):
	# Gscatt has factor of 2 included
	geom = 2.*dls*dlo/(dso*L) # all distances in same units
	nu4 = nu**4.
	zfac = (1.+zl)**(-3.)
	prefac = 48.03*10.**(-6.) # convert nanoseconds to milliseconds
	return prefac*nu4*zfac*ftilde*geom*(DMl**2.) 

def tauDM_host(DMl,ftilde,nu,zl): # in galaxy frame (not redshift corrected)
	geom = 1.
	nu4 = nu**4.
	prefac = 48.03*10.**(-6.) # convert nanoseconds to milliseconds
	return prefac*nu4*ftilde*geom*(DMl**2.) 

def Gscatt_full(dsl,dso,L):
    x = L/dso
    y = dsl/dso
    num = 1. - (2.*x/3.) + (2.*y/x)*(1-y-x)
    denom = 1. - (2.*x/3.)
    G = num/denom
    if G<1.:
        G = 1.
    return G

def tauDM_near(DMl,ftilde,dls,dlo,dso,nu,L):
	# Gscatt has factor of 2 included
	geom = Gscatt_full(dls,dso,L)
	nu4 = nu**4.
	prefac = 48.03*10.**(-6.) # milliseconds
	return prefac*nu4*ftilde*geom*(DMl**2.) 

def check_gal_type(gal_type,r200,lz):

	# returns density model parameters for given galaxy type and halo virial radius
	# fiducial halo virial radius is scaled with redshift appropriately

	if gal_type==0: # spiral galaxy

		# draw mid-plane density and Ftilde from distributions
		nh_thin = abs(random.normal(0.2,0.07)) # just in case it returns a negative number (unlikely)
		nh_thick = abs(random.normal(0.015,0.005))
		Ftilde_thin = abs(random.lognormal(lognorm(1.,0.5)[0],lognorm(1.,0.5)[1]))
		Ftilde_thick = abs(random.lognormal(lognorm(0.003,0.001)[0],lognorm(0.003,0.001)[1])) 

		# fiducial length parameters
		zh_thick = 1600. # pc
		rh_thick = 1*10**4.
		zh_thin = 200.
		rh_thin = 0.5*10**4.
		rcut_thin = 2.*10**4. # pc, 20 kpc
		rcut_thick = 2.*(10.**4.) # 20 kpc
		zcut_thin = 5000. # pc
		zcut_thick = 5000. # pc
		r200_ref = virial_radius(Msun*1.5*10.**(12.),lz)*1000. # pc

		# re-scale length parameters by halo virial radius
		rh_thin2 = r200*(rh_thin/r200_ref)
		rh_thick2 = r200*(rh_thick/r200_ref)
		zh_thin2 = r200*(zh_thin/r200_ref)
		zh_thick2 = r200*(zh_thick/r200_ref)

		rcut_thin2 = r200*(rcut_thin/r200_ref)
		rcut_thick2 = r200*(rcut_thick/r200_ref)
		zcut_thin2 = r200*(zcut_thin/r200_ref)
		zcut_thick2 = r200*(zcut_thick/r200_ref)

		return nh_thin,nh_thick,Ftilde_thin,Ftilde_thick,rh_thin2,rh_thick2,zh_thin2,zh_thick2,rcut_thin2,rcut_thick2,zcut_thin2,zcut_thick2

	if gal_type==1: # elliptical galaxy

		# draw mid-plane density and Ftilde from distributions
		nh = abs(random.normal(0.015,0.005)) # just in case it returns a negative number (unlikely)
		Ftilde = abs(random.lognormal(lognorm(0.003,0.001)[0],lognorm(0.003,0.001)[1]))

		# fiducial length parameters
		zh = 5000.
		rh = 5000.
		rcut = 20.*1000. # pc, 20 kpc
		zcut = 20.*1000.
		r200_ref = virial_radius(Msun*1.5*10.**(12.),lz)*1000. # pc

		# re-scale length parameters by halo virial radius
		rh2 = r200*(rh/r200_ref) 
		zh2 = r200*(zh/r200_ref)

		rcut2 = r200*(rcut/r200_ref)
		zcut2 = r200*(zcut/r200_ref)

		return nh,Ftilde,rh2,zh2,rcut2,zcut2

	if gal_type==2: # dwarf galaxy

		# draw mid-plane density and Ftilde from distributions
		nh = abs(random.normal(0.05,0.017)) # just in case it returns a negative number (unlikely)
		Ftilde = abs(random.lognormal(lognorm(0.2,0.5)[0],lognorm(0.2,0.5)[1]))

		# fiducial length parameters
		zh = 3000.
		rh = 3000.
		rcut = 10.*1000. # pc, 10 kpc 
		zcut = 10.*1000.
		r200_ref = virial_radius(Msun*10.**(9.8),lz)*1000

		# re-scale length parameters by halo virial radius
		rh2 = r200*(rh/r200_ref) # virial radius for SMC-mass halo
		zh2 = r200*(zh/r200_ref)

		rcut2 = r200*(rcut/r200_ref)
		zcut2 = r200*(zcut/r200_ref)

		return nh,Ftilde,rh2,zh2,rcut2,zcut2


def DM_disk(rs,incl,lz,sz,nu,thetas,r200,gal_type=0,Ftildez=False,host_gal=False,host_path=0.): 

	# input a single inclination angle, single impact parameter, lens redshift, source redshift, 
	# ref. freq., azimuthal angle, virial radius, galaxy type, whether to do Ftilde(z), whether
	# to treat as host galaxy, if host galaxy how far FRB source is within host
	
	# two-component density model for spiral galaxy disk

	# grab density model parameters 
	nh_thin2,nh_thick2,Ftilde_thin,Ftilde_thick,rh_thin2,rh_thick2,zh_thin2,zh_thick2,rcut_thin2,rcut_thick2,zcut_thin2,zcut_thick2 = check_gal_type(gal_type,r200,lz)

	if Ftildez==True:
		Ftilde_thin = Ftilde_z(Ftilde_thin,lz) # (pc^2 km)^{-1/3}
		Ftilde_thick = Ftilde_z(Ftilde_thick,lz) 

	#thetas = 0. # 0 = x-axis, 90 = y-axis
	dso = cosmo.angular_diameter_distance(sz).value*(10.**6.) # pc
	dlo = cosmo.angular_diameter_distance(lz).value*(10.**6.)
	dls = cosmo.angular_diameter_distance_z1z2(lz,sz).value*(10.**6.)
	#dls = dso - dlo

	#if host_gal==False:
	# calculate path length through galaxy at impact parameter of zero, to use for Gscatt
	xs0 = 0.
	ys0 = 0.
	#zs = linspace(-20.,20.,10000)*1000.  # pc -- observer's axis to integrate along, edge-on z0 = 0, face-on zs = z0
	zs = linspace(-rcut_thick2,rcut_thick2,10000)

	if host_gal==True:
		zs = linspace(host_path,rcut_thick2,5000)

	x0 = xs0
	y0 = ys0*cos(deg2rad(incl)) + zs*sin(deg2rad(incl)) # edge-on y0 --> zs
	z0 = -ys0*sin(deg2rad(incl)) + zs*cos(deg2rad(incl)) # edge-on z0 --> -ys --> 0.

	r0 = sqrt((x0**2.) + (y0**2.))
	cut_thin = elliptical_cutoff(r0,z0,rcut_thin2,zcut_thin2)
	cut_thick = elliptical_cutoff(r0,z0,rcut_thick2,zcut_thick2)

	# path length at center of galaxy
	path_thin = max(zs[(cut_thin<1.)]) - min(zs[(cut_thin<1.)])
	path_thick = max(zs[(cut_thick<1.)]) - min(zs[(cut_thick<1.)])

	xs = rs*cos(deg2rad(thetas))
	ys = rs*sin(deg2rad(thetas))
		
	x0 = xs
	y0 = ys*cos(deg2rad(incl)) + zs*sin(deg2rad(incl)) # edge-on y0 --> zs
	z0 = -ys*sin(deg2rad(incl)) + zs*cos(deg2rad(incl)) # edge-on z0 --> -ys --> 0.

	r0 = sqrt((x0**2.) + (y0**2.))

	cut_thin = elliptical_cutoff(r0,z0,rcut_thin2,zcut_thin2)
	if len(cut_thin[cut_thin<1.])>0:
		ne_thin = density(nh_thin2,zh_thin2,z0,rh_thin2,r0)
		dm_thin = trapz(ne_thin[(cut_thin<1.)],zs[(cut_thin<1.)])
		tau_thin = tauDM(dm_thin,Ftilde_thin,dls,dlo,dso,nu,lz,path_thin)
		if host_gal==True:
			tau_thin = tauDM_host(dm_thin,Ftilde_thin,nu,lz)
	elif len(cut_thin[cut_thin<1.])==0:
		dm_thin = 0.
		tau_thin = 0.

	cut_thick = elliptical_cutoff(r0,z0,rcut_thick2,zcut_thick2)
	if len(cut_thick[cut_thick<1.])>0:
		ne_thick = density(nh_thick2,zh_thick2,z0,rh_thick2,r0)
		dm_thick = trapz(ne_thick[(cut_thick<1.)],zs[(cut_thick<1.)])
		tau_thick = tauDM(dm_thick,Ftilde_thick,dls,dlo,dso,nu,lz,path_thick)
		if host_gal==True:
			tau_thick = tauDM_host(dm_thick,Ftilde_thick,nu,lz)
	elif len(cut_thick[cut_thick<1.])==0:
		dm_thick = 0.
		tau_thick = 0.

	dmint = dm_thin + dm_thick
	tau_simp = tau_thick + tau_thin

	#except:

	#	dmint = 0.
	#	tau_simp = 0.

	return dmint,tau_simp


def DM_ellipdwarf(rs,incl,lz,sz,nu,thetas,r200,gal_type,Ftildez=False,host_gal=False,host_path=0.): # input a single inclination angle, single impact parameter, lens redshift, source redshift, ref. freq.
	
	# single component density model for elliptical and dwarf galaxies

	# grab density model parameters depending on gal_type
	nh2,Ftilde,rh2,zh2,rcut2,zcut2 = check_gal_type(gal_type,r200,lz)

	if Ftildez==True:
		Ftilde = Ftilde_z(Ftilde,lz) # (pc^2 km)^{-1/3}

	dso = cosmo.angular_diameter_distance(sz).value*(10.**6.) # pc
	dlo = cosmo.angular_diameter_distance(lz).value*(10.**6.)
	dls = cosmo.angular_diameter_distance_z1z2(lz,sz).value*(10.**6.)

	#if host_gal==False:
	# calculate path length through galaxy at impact parameter of zero, to use for Gscatt
	xs0 = 0.
	ys0 = 0.
	zs = linspace(-rcut2,rcut2,10000) # pc -- observer's axis to integrate along, edge-on z0 = 0, face-on zs = z0

	if host_gal==True:
		zs = linspace(host_path,rcut2,5000)

	x0 = xs0
	y0 = ys0*cos(deg2rad(incl)) + zs*sin(deg2rad(incl)) # edge-on y0 --> zs
	z0 = -ys0*sin(deg2rad(incl)) + zs*cos(deg2rad(incl)) # edge-on z0 --> -ys --> 0.

	r0 = sqrt((x0**2.) + (y0**2.))
	cut = elliptical_cutoff(r0,z0,rcut2,zcut2)

	# path length at center of galaxy
	path = max(zs[(cut<1.)]) - min(zs[(cut<1.)])

	xs = rs*cos(deg2rad(thetas))
	ys = rs*sin(deg2rad(thetas))
		
	x0 = xs
	y0 = ys*cos(deg2rad(incl)) + zs*sin(deg2rad(incl)) # edge-on y0 --> zs
	z0 = -ys*sin(deg2rad(incl)) + zs*cos(deg2rad(incl)) # edge-on z0 --> -ys --> 0.

	r0 = sqrt((x0**2.) + (y0**2.))

	cut = elliptical_cutoff(r0,z0,rcut2,zcut2)
	if len(cut[cut<1.])>0:
		ne = density(nh2,zh2,z0,rh2,r0)
		dm = trapz(ne[(cut<1.)],zs[(cut<1.)])
		tau = tauDM(dm,Ftilde,dls,dlo,dso,nu,lz,path)
		if host_gal==True:
			tau = tauDM_host(dm,Ftilde,nu,lz)
	elif len(cut[cut<1.])==0:
		dm = 0.
		tau = 0.

	return dm,tau

r'''

def DM_disk_near(rs,incl,dlo,sz,nu,thetas,rdisk): # input a single inclination angle, single impact parameter, lens distance, source redshift, ref. freq.
	# very nearby objects (e.g. LMC, SMC) -- input a lens distance rather than a redshift, no Ftilde(z), set the disk diameter explicitly rather than scaling to virial radius

	Ftilde_thin = 1.
	Ftilde_thick = 0.003

	rh_thin2 = rdisk*(rh_thin/20000.)
	rh_thick2 = rdisk*(rh_thick/20000.)
	zh_thin2 = rdisk*(zh_thin/20000.)
	zh_thick2 = rdisk*(zh_thick/20000.)

	rcut_thin2 = rdisk*(rcut_thin/20000.)
	rcut_thick2 = rdisk*(rcut_thick/20000.)
	zcut_thin2 = rdisk*(zcut_thin/20000.)
	zcut_thick2 = rdisk*(zcut_thick/20000.)

	nh_thin2 = rh_thin2*(nh_thin/rh_thin)
	nh_thick2 = rh_thick2*(nh_thick/rh_thick)

	dso = cosmo.angular_diameter_distance(sz).value*(10.**6.) # pc
	dls = dso - dlo

	dm_tot = []
	tau_tot = []
	path_tot = []

	# calculate path length through galaxy at impact parameter of zero, to use for Gscatt
	xs0 = 0.
	ys0 = 0.
	#zs = linspace(-20.,20.,10000)*1000.  # pc -- observer's axis to integrate along, edge-on z0 = 0, face-on zs = z0
	zs = linspace(-rcut_thick2,rcut_thick2,10000)

	x0 = xs0
	y0 = ys0*cos(deg2rad(incl)) + zs*sin(deg2rad(incl)) # edge-on y0 --> zs
	z0 = -ys0*sin(deg2rad(incl)) + zs*cos(deg2rad(incl)) # edge-on z0 --> -ys --> 0.

	r0 = sqrt((x0**2.) + (y0**2.))
	cut_thin = elliptical_cutoff(r0,z0,rcut_thin2,zcut_thin2)
	cut_thick = elliptical_cutoff(r0,z0,rcut_thick2,zcut_thick2)

	# path length at center of galaxy
	path_thin = max(zs[(cut_thin<1.)]) - min(zs[(cut_thin<1.)])
	path_thick = max(zs[(cut_thick<1.)]) - min(zs[(cut_thick<1.)])

	xs = rs*cos(deg2rad(thetas))
	ys = rs*sin(deg2rad(thetas))
		
	x0 = xs
	y0 = ys*cos(deg2rad(incl)) + zs*sin(deg2rad(incl)) # edge-on y0 --> zs
	z0 = -ys*sin(deg2rad(incl)) + zs*cos(deg2rad(incl)) # edge-on z0 --> -ys --> 0.

	r0 = sqrt((x0**2.) + (y0**2.))

	cut_thin = elliptical_cutoff(r0,z0,rcut_thin2,zcut_thin2)
	if len(cut_thin[cut_thin<1.])>0:
		ne_thin = density(nh_thin2,zh_thin2,z0,rh_thin2,r0)
		dm_thin = trapz(ne_thin[(cut_thin<1.)],zs[(cut_thin<1.)])
		tau_thin = tauDM_near(dm_thin,Ftilde_thin,dls,dlo,dso,nu,path_thin) 
	elif len(cut_thin[cut_thin<1.])==0:
		dm_thin = 0.
		tau_thin = 0.

	cut_thick = elliptical_cutoff(r0,z0,rcut_thick2,zcut_thick2)
	if len(cut_thick[cut_thick<1.])>0:
		ne_thick = density(nh_thick2,zh_thick2,z0,rh_thick2,r0)
		dm_thick = trapz(ne_thick[(cut_thick<1.)],zs[(cut_thick<1.)])
		tau_thick = tauDM_near(dm_thick,Ftilde_thick,dls,dlo,dso,nu,path_thick)
	elif len(cut_thick[cut_thick<1.])==0:
		dm_thick = 0.
		tau_thick = 0.

	dmint = dm_thin + dm_thick
	tau_simp = tau_thick + tau_thin

	return dmint,tau_simp
r'''



