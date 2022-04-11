# -*- coding: utf-8 -*-
# Copyright (C) 2022 by Jason Koglin, Los Alamos National Laboratory
# All rights reserved. BSD 3-clause License.
# This file is part of the catemis package. Details of the copyright and
# user license can be found in the 'LICENSE' file distributed with the
# package.

"""
Controls Information
"""

def read_controls_csv(filename=None, 
                      daystr=None,
                      data_path=None, 
                      time_offset=-0.0195, 
                      nroll=7,
                      **kwargs):
    """
    Read controls data file
    
    Currently this has only been used for 2020-07-01 cathode camera
    and pyrometer calibration.
    
    Default file is _{daystr}.csv in the catemis/data folder
    
    Returns
    -------
    xarray.Dataset
    
    Parameters
    ----------
    filename : str or Path
        Name of file. 
    data_path : str or Path
        File path.  Default is catemis/data folder
    time_offset : float
        Offset time of Omega pyrometer compared to other data.
        Default = -0.0195 
    nroll : int
        Number of points for rolling mean of controls data
    """
    import numpy as np
    import pandas as pd
    from pathlib import Path
    import datetime
    
    if not data_path:
        data_path = Path(__file__).parent/'data'
    
    if filename is None:
        filename = data_path/'controls_{:}.csv'.format(daystr)
        
    df = pd.read_csv(filename, index_col=0, skiprows=1)
    dhist = df.to_xarray()
    
    dheader = pd.read_csv(filename, index_col=0, skiprows=0, nrows=1)
    dnames = dheader.T.to_dict()['Timestamp']
    
    dunits = {
     'Cathode_Current': 'A',
     'Cathode_Voltage': 'V',
     'Cathode_Power': 'kW',
     'Omega': 'C',
     }
    
    dhist.attrs['daystr'] = str(daystr)
    dhist['hour'] = (dhist.Timestamp-dhist.Timestamp.astype(int))*24 \
        +time_offset
   
    year = int(daystr[:4])
    month = int(daystr[4:6])
    day = int(daystr[6:])
    dt0 = datetime.datetime(year, month, day)
    datetimes = [np.datetime64(dt0 + datetime.timedelta(seconds=hour*60*60))
                 for hour in dhist.hour.values]

    dhist.coords['datetime'] = (('hour'), datetimes)
   
    for name, key in dnames.items():
        dhist[key].attrs['long_name'] = name
        dhist[key].attrs['units'] = dunits.get(key,'')

    dhist = dhist.swap_dims({'Timestamp': 'hour'}).sortby('hour')
        
    if nroll:
        dhist = dhist.rolling(hour=nroll).mean()

    dhist['Cathode_Power'] *= 1000
    dhist['Cathode_Power'].attrs['units'] = 'W'
    
    dhist['Cathode_VA'] = dhist.Cathode_Voltage*dhist.Cathode_Current
    dhist['Cathode_VA'].attrs['long_name'] = 'Cathode Power (V*A)'
    dhist['Cathode_VA'].attrs['units'] = 'W'

    return dhist


def make_pyro_correction(da,
                  fTc=20000,
                  fTo=32,
                  **kwargs):
    """
    Make Pyrometer temperature correction with emissivity correction 
    based on Eq 8 and fTc and fTo parameters.

    Parameters
    ----------
    da : array
        Temperature array
    fTc : float
        Omega IR2C emissivity Tc correction factor
    fTo : float
        Omega IR2C emissivity To correction factor
        
    """
        
    return 1/(1/(da+273+fTo)+1/fTc)-273
