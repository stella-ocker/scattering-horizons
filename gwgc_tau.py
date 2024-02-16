from numpy import *
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
from astropy.cosmology import Planck15 as cosmo
from astropy import constants as const
from astropy.coordinates import Distance
from DMhalo import *
from astropy.table import Table
import matplotlib as mpl
from tauDMdisk import *

def tau(Ftilde,DM_l,z_l,obs_freq,Gsc): # Gscatt dimensionless
    return (48.03*10.**(-3.))*Ftilde*(DM_l**2.)*((1.+z_l)**(-3.))*(obs_freq**(-4.))*Gsc # microseconds

infile = 'gwgc_cat.fit' # distance > 0.5 Mpc
t = Table.read(infile)
glon = t['_Glon']
glat = t['_Glat']
ra = t['_RAJ2000']
dec = t['_DEJ2000']
BMAG = t['bmag_lc']
e_BMAG = t['e_BMAG']
a = t['a']
e_a = t['e_a']
b = t['b']
e_b = t['e_b']
ba = t['b_a']
pa = t['PA']
dist = t['Dist']
e_dist = t['e_Dist']
TT = t['TT']

incl_deg = rad2deg(arccos(ba))

if glon.min() >=  0:
    glon[where(glon>180)] -= 360

ellipse_centers = []
for i in range(len(glon)):
    ellipse_centers.append((glon[i],glat[i]))
ellipse_centers = array(ellipse_centers)

# calculate halo mass, scaling by B-band magnitude and mass of Milky Way
MW_ratio = (1.3*10.**(12.))/(-20.8) # mass over absolute B-band magnitude
gal_mass = BMAG*MW_ratio
gal_mass_err = e_BMAG*MW_ratio # only accounting for error in B-band mag of other galaxies
# re-set mass of M31, M33
gal_mass[677] = 1.5*10.**(12.)
gal_mass[1703] = 10.**(11.72)

#Ftilde_halo = 10.**(-4.) 
#Ftilde_halo = abs(random.lognormal(-0.9556,0.83255)) # corresponds to mean and sigma of 10**(-4) for normal distribution

zls = dist*cosmo.H0/(const.c.to('km/s'))

# go through each galaxy, calculate 2*virial radius and re-scale to an angular size given its distance
r200 = []
diams = []
for i in range(len(gal_mass)):
    mvir = virial_radius(Msun*gal_mass[i],zls[i])*2. # twice virial radius, kpc
    r200.append(mvir)
    theta_diam = mvir*10**(-3.)/dist[i] # distance in Mpc, mvir in kpc, gives angle in radians
    diam_deg = rad2deg(theta_diam)
    diams.append(diam_deg)
r200 = array(r200)
diams = array(diams)

# grab all galaxies > 1 deg in diameter (1/2 deg for plotting higher res)
big_gals = where((2*diams)>0.5)[0]
print('max distance ',max(dist[big_gals]))
print('number of galaxies ',len(big_gals))

frb_zs = 5.
frb_dist = cosmo.angular_diameter_distance(frb_zs).value

# create coordinate grid on which to calculate tau
coord_gridx = arange(-180.,181.,0.5) 
coord_gridy = arange(-90.,91.,0.5)
dm_grid = zeros((len(coord_gridx),len(coord_gridy)))
tau_grid = zeros((len(coord_gridx),len(coord_gridy)))

r'''
for i in range(len(coord_gridx)):
    print('l = ',coord_gridx[i])
    for j in range(len(coord_gridy)):
        for k in range(len(big_gals)):
            ind = big_gals[k]
            test_coords = ellipse_centers[ind]
            impact_deg = sqrt((test_coords[0]-coord_gridx[i])**2. + (test_coords[1]-coord_gridy[j])**2.)
            impact_dist = deg2rad(impact_deg)*dist[ind]*10.**(3.) # kpc
            vir_kpc = virial_radius(Msun*gal_mass[ind],zls[ind]) 
            # calculate halo contribution
            if impact_dist < 2.*vir_kpc:
                DMi,Gold,leff = dm_halo(zls[ind],frb_zs,gal_mass[ind],impact_dist) 
                Gi = Gscatt_full(dist[ind],frb_dist,leff*10.**(-3.)) # distances in Mpc, leff output in kpc
                taui = tau(Ftilde_z(zls[ind],Ftilde_halo),DMi,zls[ind],1.,Gi) # nu = 1 GHz, tau in us
                dm_grid[i,j] = dm_grid[i,j] + DMi # redshift correction is negligible for DM
                tau_grid[i,j] = tau_grid[i,j] + taui
            # calculate disk contribution -- check which galaxy type
            if gal_mass[ind] < 10.**(10.):
                disk_cutoff = vir_kpc*(10./39.) # dwarf galaxy cutoff
                gal_type = 2
            if gal_mass[ind] >= 10.**(10.):
                disk_cutoff = vir_kpc*(20./236.) # spiral/elliptical cutoff
                if ind==677 or ind==1703: # force M31,M33 to be disks
                    gal_type = 0
                else:
                    gal_type = random.binomial(1.,0.5)
            #disk_cutoff = vir_kpc*(20./236.)
            if impact_dist < disk_cutoff:
                phi_samp = random.uniform(0.,90.)
                if gal_type==0:
                    DMd,taud = DM_disk(impact_dist*1000.,incl_deg[ind],zls[ind],frb_zs,1.,phi_samp,vir_kpc*(10**3.),gal_type=gal_type,Ftildez=True) 
                if gal_type==1:
                    DMd,taud = DM_ellipdwarf(impact_dist*1000.,incl_deg[ind],zls[ind],frb_zs,1.,phi_samp,vir_kpc*(10**3.),gal_type=gal_type,Ftildez=True) 
                if gal_type==2:
                    DMd,taud = DM_ellipdwarf(impact_dist*1000.,incl_deg[ind],zls[ind],frb_zs,1.,phi_samp,vir_kpc*(10**3.),gal_type=gal_type,Ftildez=True) 
                # taud in ms, need to convert to mus
                taud = taud*1000.
                dm_grid[i,j] = dm_grid[i,j] + DMd
                tau_grid[i,j] = tau_grid[i,j] + taud
r'''

# swap looping order so that every galaxy halo has a single n0, Ftilde
for k in range(len(big_gals)):

    print('loop ',k,' out of ',len(big_gals))
    ind = big_gals[k]
    test_coords = ellipse_centers[ind]
    vir_kpc = virial_radius(Msun*gal_mass[ind],zls[ind]) 
    Ftilde_halo = random.lognormal(-9.556913962256155,0.8325546111576977) # draw value of Ftilde for the halo

    # check what galaxy type -- base on TT rather than mass 
    if TT[ind]>=9: # dwarf
        disk_cutoff = vir_kpc*(10./39.) # dwarf galaxy cutoff
        gal_type = 2
    if TT[ind]<=0: # elliptical
        disk_cutoff = vir_kpc*(20./236.) # spiral/elliptical cutoff
        gal_type = 1
    if TT[ind]>0 and TT[ind]<9: # spiral
        disk_cutoff = vir_kpc*(20./236.) # spiral/elliptical cutoff
        gal_type = 0
    if ind==677 or ind==1703:
        Ftilde_halo = 10.**(-4.) # fix Ftilde for M31,M33 

    #if gal_mass[ind] >= 10.**(10.):
    #    disk_cutoff = vir_kpc*(20./236.) # spiral/elliptical cutoff
    #    if ind==677 or ind==1703: # force M31,M33 to be disks
    #        gal_type = 0
    #        Ftilde_halo = 10.**(-4.)
    #    else:
    #        gal_type = random.binomial(1.,0.5)

    # loop through coordinate grid to calculate DM, tau
    for i in range(len(coord_gridx)):
        for j in range(len(coord_gridy)):
            impact_deg = sqrt((test_coords[0]-coord_gridx[i])**2. + (test_coords[1]-coord_gridy[j])**2.)
            impact_dist = deg2rad(impact_deg)*dist[ind]*10.**(3.) # kpc
            # calculate halo contribution
            if impact_dist < 2.*vir_kpc:
                DMi,Gold,leff = dm_halo(zls[ind],frb_zs,gal_mass[ind],impact_dist) 
                Gi = Gscatt_full(dist[ind],frb_dist,leff*10.**(-3.)) # distances in Mpc, leff output in kpc
                taui = tau(Ftilde_z(zls[ind],Ftilde_halo),DMi,zls[ind],1.,Gi) # nu = 1 GHz, tau in us
                dm_grid[i,j] = dm_grid[i,j] + DMi # redshift correction is negligible for DM
                tau_grid[i,j] = tau_grid[i,j] + taui
            # calculate ISM contribution -- still allows n0, Ftilde to be different within each galaxy! maybe can get away with it because disk angular sizes so small? 
            if impact_dist < disk_cutoff:
                phi_samp = random.uniform(0.,90.)
                if gal_type==0:
                    DMd,taud = DM_disk(impact_dist*1000.,incl_deg[ind],zls[ind],frb_zs,1.,phi_samp,vir_kpc*(10**3.),gal_type=gal_type,Ftildez=True) 
                if gal_type==1:
                    DMd,taud = DM_ellipdwarf(impact_dist*1000.,incl_deg[ind],zls[ind],frb_zs,1.,phi_samp,vir_kpc*(10**3.),gal_type=gal_type,Ftildez=True) 
                if gal_type==2:
                    DMd,taud = DM_ellipdwarf(impact_dist*1000.,incl_deg[ind],zls[ind],frb_zs,1.,phi_samp,vir_kpc*(10**3.),gal_type=gal_type,Ftildez=True) 
                # taud in ms, need to convert to mus
                taud = taud*1000.
                dm_grid[i,j] = dm_grid[i,j] + DMd
                tau_grid[i,j] = tau_grid[i,j] + taud


dm_grid = swapaxes(dm_grid,0,1)
tau_grid = swapaxes(tau_grid,0,1)

savez('/Users/stellaocker/Desktop/gwgc_allsky_tauDM_zs5_3galtypes_halfdegres_TT',dm_grid=dm_grid,tau_grid=tau_grid)


