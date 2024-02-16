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
from DMhalo import *
from astropy.table import Table
import matplotlib as mpl
from tauDMdisk import *
from intersection_probabilities import *
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from colossus.halo import mass_so

def tau_halo(Ftilde,DM_l,z_l,obs_freq,Gsc): # d's in Gpc, leff in Mpc
    # multiple Gscatt by two to correct
    return 48.03*Ftilde*(DM_l**2.)*((1.+z_l)**(-3.))*(obs_freq**(-4.))*2.*Gsc # microseconds

def gal_redshift(mean_gal_z,zlow,zhigh):
    # generate intervening galaxy redshifts from Poisson distribution, within a given redshift bin ending at zhigh
    
    n = [random.random() for j in range(500)]
    _lambda = mean_gal_z
    inter_event_z = array([-log(1.0 - n[j]) / _lambda for j in range(len(n))])
    absolute_event_z = zlow + cumsum(inter_event_z)
    bad_events = where(absolute_event_z>=zhigh)[0]
    if len(bad_events)>0:
        absolute_event_z = delete(absolute_event_z,bad_events)
    
    return absolute_event_z

def halo_to_stellmass(Mh,z): # input Mh in units of Msun (not log!)
    # stellar mass to halo relation
    
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
    
    Mstar_Mh_ratio = 2.*Az / ((Mh/MA_z)**(-betaz) + (Mh/MA_z)**(gammaz))
    
    return Mstar_Mh_ratio*Mh

def draw_gal_type(Mh,z): # input Mh in units of Msun (not log!)

    # determine galaxy type based on stellar mass 

    Mstar = log10(halo_to_stellmass(Mh,z)) # log10 Mstar/Msun
    if Mstar < 8.: # dwarf galaxy stellar mass range
        gal_type = 2
    if 8. <= Mstar <= 10.5: # spiral galaxy stellar mass range
        gal_type = 0
    if Mstar > 10.5: # spiral or elliptical equally likely
        gal_type = random.binomial(1.,0.5)

    return gal_type

def intervening_galaxies(frb_zs,massbin_arr):
    
    # five mass bins, five redshift bins
    # calculate mean number of galaxies intersected in each mass bin
    # calculate redshifts of those galaxies
    mean_gals_arr = zeros(len(massbin_arr)-1)
    z_arr = arange(0.,frb_zs,1.)

    int_galz_arr = []
    for i in range(len(mean_gals_arr)):
    
        # mean number of galaxies in each redshift bin
        # CONSIDER WHETHER TO USE SMALLER REDSHIFT BINS GIVEEN STEEP Z-DEPENDENCE FOR Z<1
        mass1 = massbin_arr[i]
        mass2 = massbin_arr[i+1]
        mean_galsz1 = intersection_number(z_arr[0],z_arr[1],mass1,mass2)
        mean_galsz2 = intersection_number(z_arr[1],z_arr[2],mass1,mass2)
        mean_galsz3 = intersection_number(z_arr[2],z_arr[3],mass1,mass2)
        mean_galsz4 = intersection_number(z_arr[3],z_arr[4],mass1,mass2)
    
        gal_z1 = gal_redshift(mean_galsz1,z_arr[0],z_arr[1])
        gal_z2 = gal_redshift(mean_galsz2,z_arr[1],z_arr[2])
        gal_z3 = gal_redshift(mean_galsz3,z_arr[2],z_arr[3])
        gal_z4 = gal_redshift(mean_galsz4,z_arr[3],z_arr[4])
    
        allgal_z = concatenate((gal_z1,gal_z2,gal_z3,gal_z4)) # all galaxy redshifts in one array
        # exclude redshifts greater than source or less than 0.001 
        bad_events1 = where(allgal_z>=z_arr[-1])[0]
        bad_events2 = where(allgal_z<0.001)[0]
        bad_events = concatenate((bad_events1,bad_events2))
        if len(bad_events)>0:
            allgal_z = delete(allgal_z,bad_events)
    
        int_galz_arr.append(allgal_z) # redshifts for all galaxies in each mass bin
    
    int_galz_arr = array(int_galz_arr,dtype='object')
    
    tautot = 0.
    dmtot = 0.
    #dmhalo_arr = []
    #leff_arr = []

    for i in range(len(int_galz_arr)):
    
        int_galz = int_galz_arr[i]
        Mhalo = mean(massbin_arr[i])
    
        for j in range(len(int_galz)):
        
            zl = int_galz[j] # lens redshift
            #print('zl/zs',zl/frb_zs)
            r200_halo = virial_radius((10.**Mhalo)*Msun,zl)
            Ftilde_ldraw = random.lognormal(-9.556913962256155,0.8325546111576977) # corresponds to mean and sigma of 10**(-4) for normal distribution
            Ftilde_l = Ftilde_z(Ftilde_ldraw,zl) 
            #Ftilde_l = Ftilde_ldraw # no z dependence 
            #Ftilde_l = 10.**(-4.) # fix Ftilde at a single constant value
        
            # draw impact parameter, inclination angle, azimuthal angle
            phi_samp = random.uniform(0.,90.)

            rtest = linspace(0.,2.*r200_halo-5.,500) # impact parameter in kpc
            r_prob = rtest**2.
            norm = sum(r_prob)
            r_prob = r_prob/norm
            r_samp = random.choice(rtest,p=r_prob) # only draws impact parameters within 2*r200

            i_arr = linspace(0.,90.,500)
            i_prob = sin(deg2rad(i_arr))
            norm = sum(i_prob)
            i_prob = i_prob/norm
            i_samp = random.choice(i_arr,p=i_prob)
        
            # calculate tau from halo
            dmhi,ghi,leff = dm_halo(zl,frb_zs,10.**Mhalo,r_samp) # impact parameter in kpc
            tauhi = tau_halo(Ftilde_l,dmhi,zl,1.,ghi) # nu = 1 GHz, tau in microseconds

            #dmhalo_arr.append(dmhi) # to check whether DM values make sense
            #leff_arr.append(leff)
        
            #print('Mhalo',Mhalo)
            #print('r/r200',r_samp/r200_halo)
            #print('DMhalo',dmhi)
            #print('Gscatt_halo',ghi)
            #print('tau_halo',tauhi)
            #print('')

            r'''
            # old code with galaxy type tied to halo mass
            # calculate tau from galaxy ISM, all lengths in pc, for different galaxy types
            if Mhalo < 10.:
                gal_type = 2
                dmdi,taudi = DM_ellipdwarf(r_samp*1000.,i_samp,zl,frb_zs,1.,phi_samp,r200_halo*1000.,gal_type=gal_type,Ftildez=False) # tau in milliseconds
                #dmdi,taudi = 0.,0. # ignore dwarf  ISM

            if Mhalo >= 10.:
                gal_type = random.binomial(1.,0.5) # equal probability of 0 or 1
                if gal_type == 0:
                    dmdi,taudi = DM_disk(r_samp*1000.,i_samp,zl,frb_zs,1.,phi_samp,r200_halo*1000.,gal_type=gal_type,Ftildez=False) # tau in milliseconds
                if gal_type==1:
                    dmdi,taudi = DM_ellipdwarf(r_samp*1000.,i_samp,zl,frb_zs,1.,phi_samp,r200_halo*1000.,gal_type=gal_type,Ftildez=False) # tau in milliseconds
            r'''

            gal_type = draw_gal_type(10.**Mhalo,zl)
            if gal_type == 0:
                dmdi,taudi = DM_disk(r_samp*1000.,i_samp,zl,frb_zs,1.,phi_samp,r200_halo*1000.,gal_type=gal_type,Ftildez=True) # tau in milliseconds
            if gal_type == 1:
                dmdi,taudi = DM_ellipdwarf(r_samp*1000.,i_samp,zl,frb_zs,1.,phi_samp,r200_halo*1000.,gal_type=gal_type,Ftildez=True) # tau in milliseconds
            if gal_type == 2:
                dmdi,taudi = DM_ellipdwarf(r_samp*1000.,i_samp,zl,frb_zs,1.,phi_samp,r200_halo*1000.,gal_type=gal_type,Ftildez=True) # tau in milliseconds

            #taudi = 0. # no tau from disks

            # all tau are in observer frame
            tau_l = (tauhi*10.**(-3.)) + taudi # milliseconds

            tautot = tautot + tau_l
            dmtot = dmtot + ((dmhi+dmdi)/(1.+zl)) # observer frame
        
    return tautot,dmtot #,array(leff_arr)

frb_zs = 5.
massbin_arr = arange(9.,15.,1.) # to include dwarf galaxies
#massbin_arr = arange(10.,16.,1) # to exclude dwarf galaxies


# full grid run 
tautot_arr = zeros(361*181)
dmtot_arr = zeros(361*181)
for k in range(len(tautot_arr)):
	print(k,' out of ',len(tautot_arr))
	tautot_arr[k],dmtot_arr[k] = intervening_galaxies(frb_zs,massbin_arr)


# testing code
#tautot_arr = zeros(10*10)
#dmtot_arr = zeros(361*181)
#for k in range(len(tautot_arr)):
#    print(k,' out of ',len(tautot_arr))
#    tautot_arr[k],dmtot_arr[k] = intervening_galaxies(frb_zs,massbin_arr)

#tautot_arr,leff_arr = intervening_galaxies(frb_zs,massbin_arr)
#hist(log10(tautot_arr.flatten()))
#xlabel('Halo Path Length (kpc)')
#ylabel('Counts')
#yscale('log')
#show()

print('saving data')
savez('/Users/stellaocker/Desktop/intervening_galaxies/intgal_allsky_zs5_3galtypes_dm_stellmass_noFz',tautot_arr=tautot_arr,dmtot_arr=dmtot_arr)


