# -*- coding: utf-8 -*-
# Copyright (C) 2022 by Jason Koglin, Los Alamos National Laboratory
# All rights reserved. BSD 3-clause License.
# This file is part of the catemis package. Details of the copyright and
# user license can be found in the 'LICENSE' file distributed with the
# package.

"""
Blackbody Methods
"""

def stefan_boltzmann_radiance(temp):
    """
    Radiance from Stefan-Boltzmann T^4 power relationship
    """
    from scipy import constants
    sigma = constants.Stefan_Boltzmann
        
    return sigma*temp**4


def stefan_boltzmann_tempearture(power=2400,
                                area=None,
                                emissivity=0.5,
                                radiated_fraction=1):
    """
    Temperature for a simplified blackbody that follows
    the Stefan-Boltzmann T^4 power relationship

    Parameters
    ----------
    power : float
        Power in W. The default is 2400.
    area : float
        Area in cm2. The default is 6.5" circle.
    emissivity : float
        Nominal emissivity of object. Default = 0.25
        This is typically wavelength dependent.
    radiated_fraction : float
        Fraction of power radiated (some may be conductive). 
        Default = 0.9
        
    Returns
    -------
    float
        Temperature in K

    """
    import numpy as np
    from scipy import constants
    sigma = constants.Stefan_Boltzmann
    if area is None:
        area = np.pi*(6.5/2*2.54)**2 # cm2
        
    if isinstance(power, list):
        power = np.array(power)
        
    return (power*radiated_fraction
            /(emissivity*sigma*area/1e4))**(1./4.)


def planck_pdf(wavelength, T):
    """
    Planck spectral radiance pdf distribution for fitting
    Units of W/srad/m3
    Note that an extra pi is needed to give the spectral radiant exitance
    
    Parameters
    ----------
    wavelength : array
        Wavelength [m]
    T : scalar
        Temperature [k]
    
    """
    from scipy.constants import h,k,c
    import numpy as np
    return 2*h*c**2 / (wavelength**5 
                 * (np.exp(h*c / (wavelength*k*T)) - 1))


def planck_fit_pdf(wavelength, ydata, sigma, **kwargs):
    """
    Planck fitting method using scipy.optimize.curve_fit
    """
    from scipy.optimize import curve_fit
    
    popt, pcov = curve_fit(planck_pdf, wavelength, ydata, 
                           sigma=sigma, **kwargs)

    return popt, pcov



def planck_distribution(wavelength, T, norm):
    """
    Planck distribution for fitting
    
    Parameters
    ----------
    wavelength : array
        Wavelength [m]

    T : scalar
        Temperature [k]

    norm : scalar
        Overall Normalization factor
    
    """
    from scipy.constants import h,k,c
    import numpy as np
    return norm*(2*h*c**2 / (wavelength**5 
                 * (np.exp(h*c / (wavelength*k*T)) - 1)))


def planck_fit(wavelength, ydata, sigma, **kwargs):
    """
    Planck fitting method using scipy.optimize.curve_fit
    """
    from scipy.optimize import curve_fit
    
    popt, pcov = curve_fit(planck_distribution, wavelength, ydata, 
                           sigma=sigma, **kwargs)

    return popt, pcov


def planck_distribution_offset(wavelength, T, norm, offset):
    """
    Planck distribution for fitting
    
    Parameters
    ----------
    wavelength : array
        Wavelength [m]

    T : scalar
        Temperature [k]

    norm : scalar
        Overall Normalization factor

    offset : scalar
        Offset after normalization
    
    """
    from scipy.constants import h,k,c
    import numpy as np
    return norm*(2*h*c**2 / (wavelength**5 
                 * (np.exp(h*c / (wavelength*k*T)) - 1))) + offset


def planck_fit_offset(wavelength, ydata, sigma, **kwargs):
    """
    Planck fitting method using scipy.optimize.curve_fit
    """
    from scipy.optimize import curve_fit
        
    popt, pcov = curve_fit(planck_distribution_offset, wavelength, ydata, 
                           sigma=sigma, **kwargs)

    return popt, pcov


def planck_two_temperature(wavelength, T1, T2, norm1, norm2):
    """
    Two tempearture planck distribution
    """
    return  planck_distribution(wavelength, T1, norm1) \
           +planck_distribution(wavelength, T2, norm2) 


def planck_fit2(wavelength, ydata, sigma, p0=(1500., 5000., 0.002, 1.e-9), **kwargs):
    """
    Planck fitting method using scipy.optimize.curve_fit
    """
    from scipy.optimize import curve_fit
        
    popt, pcov = curve_fit(planck_two_temperature, wavelength, ydata, 
                           p0=p0, sigma=sigma, **kwargs)

    return popt, pcov


def planck_two_temperature0(wavelength, T1, T2, norm1, norm2, offset):
    """
    Two tempearture planck distribution
    """
    return  planck_distribution(wavelength, T1, norm1) \
           +planck_distribution(wavelength, T2, norm2) + offset


def planck_fit20(wavelength, ydata, sigma, p0=(1500., 5000., 0.002, 1.e-9, 0), **kwargs):
    """
    Planck fitting method using scipy.optimize.curve_fit
    """
    from scipy.optimize import curve_fit
        
    popt, pcov = curve_fit(planck_two_temperature0, wavelength, ydata, 
                           p0=p0, sigma=sigma, **kwargs)

    return popt, pcov


def planck(wavelength, temperature, units='m'):
    """
    Planck distribution in units of W/sr/m3
    
    planck = (2hc^2/λ^5)/(exp(hc/λkT)−1)
    
    see also:  astropy.modeling.blackbody.BlackBody1D
    
    Parameters
    ----------
    wavelength : scalar, list or array
        Wavelength [m]
    temperature : scalar, list or array
        Temperature [k]
    
    """
    import numpy as np
    import xarray as xr
    from scipy import constants
    from scipy.constants import physical_constants

    eV2K_pc = physical_constants.get('electron volt-kelvin relationship')
    eV2K = eV2K_pc[0]

    if units == 'm':
        wave_scale = 1
    elif units == 'cm':
        wave_scale = 1e-2
    elif units == 'mm':
        wave_scale = 1e-3
    elif units == 'um':
        wave_scale = 1e-6
    elif units == 'nm':
        wave_scale = 1e-9
    elif units == 'A':
        wave_scale = 1e-10
    else:
        print('Invalid units {:} -- setting units = "m"'.format(units))
        units = 'm'
        wave_scale = 1
            
    awavelength = np.array(wavelength)
    if awavelength.size == 1:
        awavelength = np.array([float(awavelength)])
    
    alambda = awavelength*wave_scale
    
    atemperature = np.array(temperature)
    if atemperature.size == 1:
        atemperature = np.array([float(atemperature)])
    
    adata = np.zeros((atemperature.size,awavelength.size))
    for itemp,temp in enumerate(atemperature):
        b = constants.h*constants.c/(alambda*constants.k*temp)
        adata[itemp] = 2.0*constants.h*constants.c**2/alambda**5 \
                  / (np.exp(b) - 1.0) 

    da = xr.DataArray(adata, dims=['kelvin','wavelength'],
                      coords={'wavelength': awavelength,
                              'kelvin': atemperature})

    da.attrs['long_name'] = 'Planck Distribution'
    da.attrs['units'] = 'W/sr/m3'
    da.coords['kelvin'].attrs['long_name'] = 'Absolute Temperature'
    da.coords['kelvin'].attrs['units'] = 'K'
    da.coords['wavelength'].attrs['long_name'] = 'Wavelength'
    da.coords['wavelength'].attrs['units'] = units

    da.coords['energy'] = (('wavelength'),
                           alambda*constants.eV/(constants.h*constants.c))
    da.coords['energy'].attrs['long_name'] = 'Photon Energy'
    da.coords['energy'].attrs['units'] = 'eV'

    da.coords['celcius'] = (('kelvin'), atemperature-273.)
    da.coords['celcius'].attrs['long_name'] = 'Temperature'
    da.coords['celcius'].attrs['units'] = 'C'

    da.coords['eV'] = (('kelvin'), atemperature/eV2K)
    da.coords['eV'].attrs['long_name'] = 'Temperature'
    da.coords['eV'].attrs['units'] = 'eV'    
    
    return da


def powerlaw_distribution_offset(xdata, alpha, norm=1, offset=0):
    """
    Powerlaw distribution with offset for fitting
    
    Parameters
    ----------
    xdata : array
        Array of X data
        
    alpha : scalar
        Power 

    norm : scalar
        Overall Normalization factor

    offset : scalar
        Offset after normalization
    
    """
    import numpy as np
    return norm * np.power(np.abs(xdata), alpha) + offset


def powerlaw_fit_offset(xdata, ydata, **kwargs):
    """
    Powerlaw fitting method using scipy.optimize.curve_fit
    """
    from scipy.optimize import curve_fit
        
    popt, pcov = curve_fit(powerlaw_distribution_offset,
                           xdata, ydata, 
                           **kwargs)

    return popt, pcov


def powerlaw_distribution(xdata, alpha, norm=1):
    """
    Powerlaw distribution for fitting
    
    Parameters
    ----------
    xdata : array
        Array of X data
        
    alpha : scalar
        Power 

    norm : scalar
        Overall Normalization factor
    """
    import numpy as np
    return norm * np.power(np.abs(xdata), alpha)


def powerlaw_fit(xdata, ydata, **kwargs):
    """
    Powerlaw fitting method using scipy.optimize.curve_fit
    """
    from scipy.optimize import curve_fit
        
    popt, pcov = curve_fit(powerlaw_distribution,
                           xdata, ydata, 
                           **kwargs)

    return popt, pcov
