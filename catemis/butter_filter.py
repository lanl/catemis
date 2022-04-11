# -*- coding: utf-8 -*-
# Copyright (C) 2022 by Jason Koglin, Los Alamos National Laboratory
# All rights reserved. BSD 3-clause License.
# This file is part of the catemis package. Details of the copyright and
# user license can be found in the 'LICENSE' file distributed with the
# package.

"""
Butterworth Filter Methods
"""

def butter_bandpass(fs, lowcut=None, highcut=None, order=10):
    """
    Butterworth Bandpass Filter Method
    
    Parameters
    ----------
    fs : float
        Sample Rate
    lowcut : float
        Low pass frequency cutoff
    highcut : float
        High pass frequency cutoff
    order : int
        Butterworth filter order [Default = 10, i.e., 10th order Butterworth filter]

    Notes
    -----
    http://scipy-cookbook.readthedocs.io/items/ButterworthBandpass.html
    """
    from scipy.signal import butter
    nyq = 0.5 * fs
    if lowcut and highcut:
        high = highcut / nyq
        low = lowcut / nyq
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
    elif highcut:
        high = highcut / nyq
        sos = butter(order, high, btype='low', output='sos')
    elif lowcut:
        low = lowcut / nyq
        sos = butter(order, low, btype='high', output='sos')
    else:
        print('Error -- must supply lowcut, highcut or both')
        sos = None

    return sos

def butter_bandpass_filter(da, lowcut=None, highcut=None, 
                           order=10, dim=None, safe=None):
    """
    Buterworth high, low or band pass filter.
    
    Uses scipy.signal.sosfiltfilt -- A forward-backward digital filter using cascaded 
    second-order sections.

    Parameters
    ----------
    da : xr.DataArray
        Data array
        
    fs : float
        Sample Rate

    lowcut : float
        Low pass frequency cutoff

    highcut : float
        High pass frequency cutoff

    order : int
        Butterworth filter order [Default = 10, i.e., 10th order Butterworth filter]

    safe : bool
        Safe filter:  checks that highcut and lowcut are valid

    Notes
    -----
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sosfiltfilt.html
    http://scipy-cookbook.readthedocs.io/items/ButterworthBandpass.html
    """
    if not dim:
        if len(da.dims) == 1:
            dim = da.dims[0]
        elif 'time' in da.dims:
            dim = 'time'
        else:
            raise Exception('Error: must specify dim for DataArray')

    xaxis = da[dim].data    
    fs = 1/(xaxis[1] - xaxis[0])

    if safe:
        if highcut and fs < highcut*2.:
            return da
    
        if lowcut and fs < lowcut*2.:
            return da
        
    from scipy.signal import sosfiltfilt
    sos = butter_bandpass(fs, lowcut=lowcut, highcut=highcut, order=order)
    dfilt = sosfiltfilt(sos, da.data)
    dout = da.copy()
    dout.data = dfilt
    
    return dout



