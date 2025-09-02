from numpy import *
import random
from matplotlib.pyplot import *
import matplotlib.cm as cmx
import matplotlib.colors as colors
import os
import sys
#import calcticks
import datetime
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord,SkyOffsetFrame,Distance
from astropy.visualization import quantity_support
from astropy.time import Time
from matplotlib.patches import Ellipse
#from astropy.cosmology import Planck15 as cosmo
from astropy import constants as const
from astropy.coordinates import Distance
from astropy.table import Table
import matplotlib as mpl
from intersection_probabilities import *
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from colossus.halo import mass_so
from tauDMdisk import *

## HALOS ##

def tau_halo(Ftilde,DM_l,obs_freq): 
    Gsc = 1.
    return (48.03*10.**(-6.))*Ftilde*(DM_l**2.)*(obs_freq**(-4.))*Gsc # milliseconds

# random constants and conversion factors

kpc2cm = 3.086*(10.**(21.)) # kiloparsecs to centimeters
Msun = 1.99*(10.**(33.)) # mass of Sun in grams

def virial_radius(Mhalo,z1): 
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

def tauDM_halo(z1,Mhalo,Rperp,nu):

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

    rtr = linspace(Rperp,2.*r200) # LOS for inside host galaxy starts at impact parameter and goes to edge of halo
    ne_r = ne_halo(rtr,Mtot,rvcm,conc,fb,mue,muh,mp,Om,Ob,r200) # density profile along FRB LOS
    dm_ih = trapz(ne_r,rtr*10.**3.) 

    # calculate tau 
    Ftilde_ldraw = random.lognormal(-9.556913962256155,0.8325546111576977) # corresponds to mean and sigma of 10**(-4) for normal distribution
    Ftilde_l = Ftilde_z(Ftilde_ldraw,z1) 
    tau_ih = tau_halo(Ftilde_l,dm_ih,nu)

    return dm_ih, tau_ih # output dm and tau in galaxy frame, apply z correction later 

def nhalo(logmass1,logmass2,z):
    # calculates number density of halos in given mass bin as a function of redshift
    mass_bin = 10.**linspace(logmass1,logmass2,10)
    hmf = mass_function.massFunction(mass_bin,z,q_in='M', q_out='dndlnM', mdef = '200c', model = 'tinker08')
    # units of hmf are (Mpc/h)^{-3}
    density = trapz(hmf,log(mass_bin))
    return density

# need to parametrize the (intrinsic) GSMF for star-forming & quiescent galaxies, using McLeod et al. (2021)

def gsmf(logM,z): # input mass in units of log10(M/Msun), redshift 

    a1 = 10.55 
    a2 = 0.
    a3 = -0.16
    a4 = 0.12
    a5 = -1.45
    a6 = -0.08
    a7 = -2.43
    a8 = -0.17
    a9 = -0.08
    a10 = -2.94
    a11 = -0.22

    Mstar = a1 + a2*z
    alpha1 = a3 + a4*z
    alpha2 = a5 + a6*z
    phistar1 = 10.**(a7 + a8*z + a9*(z**2.))
    phistar2 = 10.**(a10 + a11*z)

    phistar_param = phistar1*10.**(alpha1*(logM-Mstar)) + phistar2*10.**(alpha2*(logM-Mstar))
    phi_logM = log(10.) * exp(-10.**(logM-Mstar)) * 10.**(logM-Mstar) * phistar_param # number density per log(M/Msun) per Mpc^3
    
    return phi_logM

# define a function to draw stellar mass and galaxy type

def stellar_mass_type(Nsamples,logM_min,logM_max,z_ref):
    # input number of samples, min and max stellar mass bins, reference redshift
    
    # calculate maximum likelihood -- corresponds to minimum mass bin 
    max_like = gsmf(logM_min,z_ref)
    
    masses = []
    gal_types = []
    while len(masses) < Nsamples:
        
        logM = random.uniform(logM_min,logM_max) # draw candidate stellar mass

        # calculate GSMF
        phi_logM = gsmf(logM,z_ref)
        if (z_ref < 2.) and (logM<=10.7):
            phi_active = phi_logM # star-forming galaxies
            phi_passive = 0.1*phi_logM # quiescent galaxies
            likelihood_list = array([phi_active,phi_passive])
            gal_type = argmax(likelihood_list) # 0 = star-forming, 1 = quiescent
        
        if (z_ref > 2.) and (logM<=10.7):
            phi_active = phi_logM
            phi_passive = 10.**(-5.)*phi_logM 
            likelihood_list = array([phi_active,phi_passive])
            gal_type = argmax(likelihood_list) 
    
        if logM>10.7:
            # equal probability of passive or active
            likelihood_list = array([phi_logM])
            gal_type = random.binomial(1.,0.5) # equal probability of 0 or 1

        if (gal_type==0) and (logM<8.):
            gal_type = 2 # set dwarf galaxies 
        
        # Accept randomly.
        u = random.uniform(0.,max_like)
        if u < max(likelihood_list):
            masses.append(logM)
            gal_types.append(gal_type)
        
    return array(masses),array(gal_types)

def stell_to_halomass(Mstar,z): # input Mstar in units of Msun (not log!)
    
    B = 11.79
    mu = 0.20
    C = 0.046
    nu = -0.38
    D = 0.709
    eta = -0.18
    F = 0.043
    E = 0.96
    
    MA_z = 10.**(B + z*mu)
    Az = C*(1+z)**nu
    gammaz = D*(1+z)**eta
    betaz = F*z + E
    
    # difficult to invert the function, so calculate array of possible halo masses and find the one that matches
    
    Mh_tests = 10.**arange(10.,15.1,0.05)
    Mstar_tests = []
    for i in range(len(Mh_tests)):
        Mstar_Mh_ratio = 2.*Az / ((Mh_tests[i]/MA_z)**(-betaz) + (Mh_tests[i]/MA_z)**(gammaz))
        Mstar_tests.append(Mstar_Mh_ratio*Mh_tests[i])
    Mstar_tests = array(Mstar_tests)
    Mstar_resids = abs(log10(Mstar_tests) - log10(Mstar)) # residual of log mass
    Mstar_true = argmin(Mstar_resids)
    Mhalo_true = Mh_tests[Mstar_true]
    
    return Mhalo_true

def host_galaxies_eval(frb_zs):

    nu = 1.
    tautot_arr = zeros(361*181)
    dmtot_arr = zeros(361*181)
    galtype_arr = zeros(361*181)

    # old procedure:
    # assign probability of each galaxy type
    # draw halo mass from HMF, again keeping dwarfs < 10^10 Msun
    # evalute DM & scattering within each host 

    #galtype_choices = array([0,1,2])
    #galtype_pdf = array([0.85,0.05,0.1])
    #galtype_pdf = array([0.60,0.01,0.39])


    #dwarf_massbins = arange(9.,10.,0.1)
    #spirellip_massbins = arange(10.,15.,0.5)

    #hmf_dwarf = mass_function.massFunction(10.**dwarf_massbins,frb_zs,q_in='M', q_out='dndlnM', mdef = '200c', model = 'tinker08')
    #hmf_norm = sum(hmf_dwarf)
    #dwarf_masspdf = hmf_dwarf/hmf_norm

    #hmf_massive = mass_function.massFunction(10.**spirellip_massbins,frb_zs,q_in='M', q_out='dndlnM', mdef = '200c', model = 'tinker08')
    #hmf_norm = sum(hmf_massive)
    #spirellip_masspdf = hmf_massive/hmf_norm

    # new procedure:
    # draw stellar mass and galaxy type using GSMF
    # convert stellar mass to halo mass
    # evaluate DM & scattering, limiting impact parameter and phi to assumed FRB locations

    gal_masses, gal_types = stellar_mass_type(len(tautot_arr),7.,12.5,frb_zs)

    for i in range(len(tautot_arr)):

        # draw galaxy type, inclination angle, azimuthal angle
        #gal_type = random.choice(galtype_choices,p=galtype_pdf)

        galmass = gal_masses[i]
        gal_type = gal_types[i]

        halo_mass = stell_to_halomass(10.**galmass,frb_zs)
        r200_halo = virial_radius((halo_mass)*Msun,frb_zs)

        i_arr = linspace(0.,90.,500)
        i_prob = sin(deg2rad(i_arr))
        norm = sum(i_prob)
        i_prob = i_prob/norm
        i_samp = random.choice(i_arr,p=i_prob)
        phi_samp = random.uniform(0.,90.)
        #phi_samp = 0. # forces LOS to thin disk? not exactly but a test case

        if gal_type==0:
            #galmass = random.choice(spirellip_massbins,p=spirellip_masspdf)
            #r200_halo = virial_radius((10.**galmass)*Msun,frb_zs)
            r200_fid = virial_radius((1.5*10.**12.)*Msun,frb_zs)
            rdisk = 0.2*(r200_halo/r200_fid) # restricted locations
            #rdisk = 19.*(r200_halo/r200_fid) # anywhere in the disk
            r_samp = random.uniform(0.,rdisk)
            host_path = random.uniform(-rdisk,rdisk)
            dmhalo,tauhalo = tauDM_halo(frb_zs,halo_mass,r_samp,nu)
            dmdisk,taudisk = DM_disk(r_samp*1000.,i_samp,frb_zs,frb_zs,nu,phi_samp,r200_halo*1000.,
                    gal_type=gal_type,Ftildez=True,host_gal=True,host_path=host_path)

        if gal_type==1:
            #galmass = random.choice(spirellip_massbins,p=spirellip_masspdf)
            #r200_halo = virial_radius((10.**galmass)*Msun,frb_zs)
            r200_fid = virial_radius((1.5*10.**12.)*Msun,frb_zs)
            rdisk = 0.2*(r200_halo/r200_fid) # restricted locations
            #rdisk = 19.*(r200_halo/r200_fid) # anywhere in the disk
            r_samp = random.uniform(0.,rdisk)
            host_path = random.uniform(-rdisk,rdisk)
            dmhalo,tauhalo = tauDM_halo(frb_zs,halo_mass,r_samp,nu)
            dmdisk,taudisk = DM_ellipdwarf(r_samp*1000.,i_samp,frb_zs,frb_zs,nu,phi_samp,r200_halo*1000.,
                    gal_type=gal_type,Ftildez=True,host_gal=True,host_path=host_path)

        if gal_type==2:
            #galmass = random.choice(dwarf_massbins,p=dwarf_masspdf)
            #r200_halo = virial_radius((10.**galmass)*Msun,frb_zs)
            r200_fid = virial_radius((10.**9.8)*Msun,frb_zs)
            rdisk = 0.2*(r200_halo/r200_fid) # restricted locations
            #rdisk = 19.*(r200_halo/r200_fid) # anywhere in the disk
            r_samp = random.uniform(0.,rdisk)
            host_path = random.uniform(-rdisk,rdisk)
            dmhalo,tauhalo = tauDM_halo(frb_zs,halo_mass,r_samp,nu)
            dmdisk,taudisk = DM_ellipdwarf(r_samp*1000.,i_samp,frb_zs,frb_zs,nu,phi_samp,r200_halo*1000.,
                    gal_type=gal_type,Ftildez=True,host_gal=True,host_path=host_path)

        tautot_arr[i],dmtot_arr[i],galtype_arr[i] = (tauhalo+taudisk) * (1.+frb_zs)**(-3.), (dmhalo+dmdisk)/(1.+frb_zs), gal_type

    return tautot_arr,dmtot_arr,galtype_arr


print('zs = 5')
tautot_arr,dmtot_arr,galtype_arr = host_galaxies_eval(frb_zs=5.)
print('saving')
savez('hostgal_allsky_zs5_3gal_stellmass_restricted',tautot_arr=tautot_arr,dmtot_arr=dmtot_arr,galtype_arr=galtype_arr)

print('zs = 1')
tautot_arr,dmtot_arr,galtype_arr = host_galaxies_eval(frb_zs=1.)
print('saving')
savez('hostgal_allsky_zs1_3gal_stellmass_restricted',tautot_arr=tautot_arr,dmtot_arr=dmtot_arr,galtype_arr=galtype_arr)

print('zs = 0.5')
tautot_arr,dmtot_arr,galtype_arr = host_galaxies_eval(frb_zs=0.5)
print('saving')
savez('hostgal_allsky_zs0point5_3gal_stellmass_restricted',tautot_arr=tautot_arr,dmtot_arr=dmtot_arr,galtype_arr=galtype_arr)


