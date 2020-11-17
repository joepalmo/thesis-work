''''
    
    AUTHOR:
        Joe Palmo

'''





################### IMPORT STATEMENTS #########################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
import seaborn as sb
from astropy import constants as const
import random
import astropy.constants as const
import math
from tqdm import tqdm
import extinction as ex
import pdb
import glob
import scipy.optimize as optimization
import scipy.interpolate as sinterp

################### Mdot Calculations #########################

def unlog(x):
    return 10**x    
    
def flux_to_luminosity(flux, distance):
    '''
    This function will turn a flux estimate (erg/s/cm^2) into a luminosity using distance
    '''
    lum = flux * (4*np.pi*(dist*const.pc.to('cm').value)**2)
    return lum
    

def luminosity_to_flux(lum, dist):
    '''
    This function will turn a flux estimate (erg/s/cm^2) into a luminosity using distance
    '''
    flux = lum / (4*np.pi*(dist*const.pc.to('cm').value)**2)
    return flux
    
    
def Lacc_to_Mdot(Lacc, mass, radius, Rin=5):
    '''
    This function will turn an accretion luminosity into a mass accretion rate estimate using the widely
    accepted relationship. The inputs of the function are:
    
    Lacc (float) : Accretion luminosity [Lsun]
    mass (float) : mass of object [Msun]
    radius (float) : radius of object [Rsun]
    
    Optional:
    Rin (float) : magnetospheric radius (default 5 [Rsun])
    
    Outputs:
    Mdot (float) : mass accretion rate [Msun/yr]
    '''
    Mdot = (1-(radius/Rin))**-1*(radius*const.R_sun.value*Lacc*const.L_sun.value)/(const.G.value*mass*const.M_sun.value) * (365*24*3600) / (const.M_sun.value)


    return Mdot
    
    
def Mdot_to_Lacc(Mdot, mass, radius, Rin=5):
    '''
    This function will turn a mass accretion rate estimate into an accretion luminosity using the widely
    accepted relationship.
    
    Inputs:
    Mdot (float) : mass accretion rate [Msun/yr]
    mass (float) : mass of object [Msun]
    radius (float) : radius of object [Rsun]
    
    Optional:
    Rin (float) : magnetospheric radius (default 5 [Rsun])
    
    Outputs:
    Lacc (float) : Accretion luminosity [Lsun]
    '''
    Lacc = Mdot*(1-(radius/Rin))*(const.G.value*mass*const.M_sun.value)*(const.M_sun.value)/(radius*const.R_sun.value)/(365*24*3600)/const.L_sun.value


    return Lacc
    
    
def UVExcess_to_Mdot(UVexcess, bc, dist, mass, radius, Av, Rin=5, ):
    '''
    This function will transform a UV Excess flux value into a mass accretion rate estimate by following 
    the process described in Herczeg 2008. 
    
    Inputs:
    UVexcess - UV continuum excess flux [erg/(s*cm^2)]
    bc - bolometric correction
    dist - distance to object [pc]
    mass - mass of object [Msun]
    radius - radius of object [Rsun]
    
    Optional:
    Rin - magnetospheric radius (default 5 [Rsun])
    
    Outputs:
    Mdot - mass accretion rate [Msun/yr]
    '''
    #Extinction correction of flux
    deredUV = ex.remove(Av, UVexcess)
    
    #use the bolometric correction factor to scale the UV excess flux to total accretion flux
    total_excess = UderedUV*bc
    
    #accretion flux to accretion luminosity
    Lacc = flux_to_luminosity(total_excess, dist)
    
    #convert accretion luminosity to solar luminosity
    Lacc = Lacc / const.L_sun.to('erg/s').value
    
    #accretion luminosity to Mdot
    Mdot = Lacc_to_Mdot(Lacc, mass, radius, Rin)
    
    return Mdot
    
    
def Mdot_to_UVExcess(Mdot, bc, dist, mass, radius, Rin=5):
    '''
    This function will transform a mass accretion rate estimate to a UV Excess flux value by following 
    the process described in Herczeg 2008. 
    
    Inputs:
    Mdot - mass accretion rate [Msun/yr]
    bc - bolometric correction
    dist - distance to object [pc]
    mass - mass of object [Msun]
    radius - radius of object [Rsun]
    
    Optional:
    Rin - magnetospheric radius (default 5 [Rsun])
    
    Outputs:
    UVexcess - UV continuum excess flux [erg/(s*cm^2)]
    '''
    #Mdot to accretion luminosity
    Lacc = Mdot_to_Lacc(Mdot, mass, radius, Rin)
    
    #convert accretion luminosity to erg/s
    Lacc = Lacc * const.L_sun.to('erg/s').value
    
    #accretion luminosity to accretion flux
    total_excess = luminosity_to_flux(Lacc, dist)
    
    #use the bolometric correction factor to scale the total accretion flux to UV excess flux
    UVexcess = total_excess / bc
    
    #Extinction correction of flux
    redUV = ex.apply(Av, UVexcess)
    
    return redUV
    
    
def UbandExcess_to_Mdot(Uexcess, dist, mass, radius, Av, Rin=5, unc=False):
    '''
    This function will transform a U band flux value into a mass accretion rate estimate by following 
    the process described in Robinson 2019. 
    
    Inputs:
    Uexcess - U band continuum excess flux [erg/(s*cm^2)]
    dist - distance to object [pc]
    mass - mass of object [Msun]
    radius - radius of object [Rsun]
    
    Optional:
    Rin - magnetospheric radius (default 5 [Rsun])
    
    Outputs:
    Mdot - mass accretion rate [Msun/yr]
    '''
    #Extinction correction of flux
    deredU = ex.remove(Av, Uexcess) 
    
    #U-band flux to Lu
    Lu = flux_to_luminosity(deredU, dist)
    
    #convert Lu to solar luminosity
    Lu = Lu / const.L_sun.to('erg/s').value
    
    #Lu to Lacc using Robinson paper -- natural logarithms
    #uncertainties 0.03 for each constant
    if unc == False:
        logLacc = 0.93 *np.log(Lu) + 0.5 
    else:
        logLacc = (0.93+np.random.normal(0.03)) *np.log(Lu) + (0.5+np.random.normal(0.03)) 
    Lacc = np.exp(logLacc)
    
    #accretion luminosity to Mdot
    Mdot = Lacc_to_Mdot(Lacc, mass, radius, Rin)
    
    return Mdot
    


def lineflux_to_Mdot(flux, dist, mass, radius, Rin=5, line=None, A=None, B=None):
    '''
    This function will turn a line flux into a mass accretion rate estimate using the Lacc-Lline fits derived
    by Alcala et al 2017. 
    
    Inputs:
    flux (float) : line flux [erg/(s*cm^2)]?
    dist (float) : distance to object [pc]
    mass (float) : mass of object [Msun]
    radius (float) : radius of object [Rsun]
    
    Optional:
    Rin (float) : magnetospheric radius (default 5 [Rsun])
    line (str) : type of line. for now acceptable inputs are H-alpha: 'Ha', Pa-beta: 'Pab', and Br-gamma: 'Brg'
    a (float) : If you want to input the parameters for your own line flux vs Lacc relationship
    b (float) : If you want to input the parameters for your own line flux vs Lacc relationship
    
    Outputs:
    Mdot (float) : mass accretion rate [Msun/yr]
    '''
    
    #a & b values pulled directly from the paper
    
    if line == None:
        a = A
        b = B
    elif line == 'Ha':
        a=1.13 #+/- 0.05
        b=1.74  #+/- 0.19
    elif line == 'Pab':
        a=1.06 #+/- 0.07
        b=2.76  #+/- 0.34
    elif line == 'Brg':
        a=1.19 #+/- 0.10
        b=4.02  #+/- 0.51
    else:
        print('Line not found.')
        return
    
    #find Lline in erg/s
    Lline = flux * (4*np.pi*(dist*const.pc.to('cm').value)**2)
    #convert to solar luminosity
    Lline = Lline / const.L_sun.to('erg/s').value
    
    #Find Lacc using Alcala relationships
    logLacc = a*np.log10(Lline)+b
    #solar luminosity
    Lacc = unlog(logLacc)
    
    Mdot = Lacc_to_Mdot(Lacc, mass, radius)
    
    return Mdot


################### Object Parameter Estimation #########################


#Need to have these models handy
models = glob.glob('StellarParams/Baraffe*txt')
mesamodels = glob.glob('StellarParams/MESA_*.txt')



def to_SpTyNum(spectype):
    '''
    This function takes a spectral type in the form - 'letternumber', i.e. 'M5'.
    It then translates that spec type into a spectral type number identifier used for interpolation.
    '''
    spty_dict = {'B' : 0,
                 'A' : 1,
                 'F' : 2,
                 'G' : 3,
                 'K' : 4,
                 'M' : 5}
    
    letter = spectype[0]
    number = spectype[1:]
    
    sptynum = spty_dict[letter]*10 + int(number)
    
    return sptynum
    
    
def SpTy_to_Teff(spectype):
    '''
    This function will take a numerical spectral type identifier, and using interpolation from the tables in
    Herczeg and Hillenbrand 2014 it will calculate an effective temperature.
    '''
    SpTyNum, Teff, SpTy = np.genfromtxt('StellarParams/HerczegHillenbrand_SpTyTeff_Numericized.txt', skip_header=1, dtype='str').T
    SpTyNum, Teff = [float(x) for x in SpTyNum], [float(y) for y in Teff]
    
    spl = sinterp.UnivariateSpline(SpTyNum, Teff)
    
    teff = spl(spectype)
    
    return teff
    
    
def Teff_to_params(Teff, age):
    '''
    This function will take an effective temperature and an age, and by interpolating the Baraffe 2015 models 
    it will return a mass and radius estimate for the object.
    
    Inputs:
    Teff (float) - effective temperature in Kelvin
    age (float) - age of object in Myr
    
    Outputs:
    m (float) - interpolated mass of object in solar masses
    r (floar) - interpolated radius of object in solar radii
    
    '''
    #Find the most accurate model given age
    ages = [1, 2, 3, 8, 10]
    age_diff = []
    for a in ages:
        age_diff.append(np.abs(age-a))
    closest_age = ages[age_diff.index(min(age_diff))]
    
    model = None
    for m in models:
        if (str(closest_age)+'Myr') in m:
            model = m
            
    mass, teff, loglum, logg, radius = np.loadtxt(model, skiprows=1).T
    
    p_mass = np.polyfit(teff, mass, 10)
    f_mass = np.poly1d(p_mass)
    
    p_radius = np.polyfit(teff, radius, 10)
    f_radius = np.poly1d(p_radius)
    
    m = f_mass(Teff)
    r = f_radius(Teff)
    
    return m, r
    
    
def mass_to_Teff(massobj, age):
    '''
    This function will take a mass and an age, and by interpolating the Baraffe 2015 models 
    it will return a temperature estimate for the object.
    '''
    #Find the most accurate model given age
    ages = [1, 2, 3, 8, 10]
    age_diff = []
    for a in ages:
        age_diff.append(np.abs(age-a))
    closest_age = ages[age_diff.index(min(age_diff))]
    
    model = None
    for m in models:
        if (str(closest_age)+'Myr') in m:
            model = m
            
    mass, teff, loglum, logg, radius = np.loadtxt(model, skiprows=1).T
    
    p_teff = np.polyfit(mass, teff, 10)
    f_teff = np.poly1d(p_teff)
    
    t = f_teff(massobj)
    
    return t
    
    
def Teff_to_SpTy(teff):
    '''
    This function will take an effective temperature, and using interpolation from the tables in
    Herczeg and Hillenbrand 2014 it will calculate a numerical spectral type identifier.
    '''
    SpTyNum, Teff, SpTy = np.genfromtxt('StellarParams/HerczegHillenbrand_SpTyTeff_Numericized.txt', skip_header=1, dtype='str').T
    SpTyNum, Teff = [float(x) for x in SpTyNum], [float(y) for y in Teff]
    SpTyNum.reverse()
    Teff.reverse()
    
    spl = sinterp.UnivariateSpline(Teff, SpTyNum)
    
    spty = spl(teff)
    
    return spty
    
################### Uncertainty Distributions #########################


def gaussian(x, mu, sigma):
    g = (1/(sigma*np.sqrt(2*np.pi)))*np.exp((-1/2)*(x-mu)**2/sigma)
    return g
    
def cdf(g, dx):
    cdf = np.cumsum(g) * dx
    return cdf
    
    
def inverse_transform_sampling(x, f, nsamples):
    '''
    This function performs the inverse transform sampling method to sample from a probability distribution 
    function (PDF) using uniform samples between 0 and 1. It does this by calculating the cumulative distribution 
    function (CDF) of the normalized PDF, inverting it, then sampling along the probability axis from 0 to 1.
    
    Inputs:
    x (array-like) - the x-values used to generate f
    f (array-like) - the PDF
    nsamples (int) - the number of inverse transform sampled values needed, i.e. the size of the output
    
    Outputs:
    samples (array-like) - an inverse transform sampled distribution which resembles the PDF
    '''
    
    normalized = f / np.trapz(f, x=x)
    dx = x[1]-x[0]
    c = cdf(normalized, dx)
    inv = sinterp.interp1d(c, x)
    r = np.random.rand(nsamples)
    samples = inv(r)
    
    return np.array(samples)