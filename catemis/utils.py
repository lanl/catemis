# -*- coding: utf-8 -*-
# Copyright (C) 2022 by Jason Koglin, Los Alamos National Laboratory
# All rights reserved. BSD 3-clause License.
# This file is part of the catemis package. Details of the copyright and
# user license can be found in the 'LICENSE' file distributed with the
# package.

"""
Utilities
"""

def add_metadata(ds, cat=None,
                 **kwargs):
    """
    Add version, date_created and time_created metadata to xarray.Dataset.
    If cat is set then add metadate to that DataSet data_var or coord.
    Any additional kwargs will be added as key, value pairs.
    
    Parameters
    ----------
    cat : str
        Catagory within DataSet
        
    """
    import numpy as np
    from . import __version__
    
    try:
        if cat:
            ds[cat].attrs['version'] = __version__
            day_created, time_created = \
                np.datetime64('now').astype(str).split('T')
            ds[cat].attrs['date_created'] = day_created
            ds[cat].attrs['time_created'] = time_created
            for attr, val in kwargs.items():
                ds[cat].attrs[attr] = val
            
        else:
            ds.attrs['version'] = __version__
            day_created, time_created = \
                np.datetime64('now').astype(str).split('T')
            ds.attrs['date_created'] = day_created
            ds.attrs['time_created'] = time_created
            for attr, val in kwargs.items():
                ds.attrs[attr] = val

    except:
        print('... Failed adding metadata', ds)

    return ds


def get_local_datetime(
        datetime_utc,
        zone='us/mountain'):
    """
    Get the local datetime from 

    Parameters
    ----------
    datetime : np.datetime64
        UTC datetime.
    zone : str, optional
        Time Zone. The default is 'us/mountain'.

    Returns
    -------
    np.datetime64
        Local datetime
    """
    import pandas as pd
    import pytz

    tz = pytz.timezone(zone)
    tstamp = pd.Timestamp(datetime_utc)
    try:
        datetime = (pd.Timedelta(tz.utcoffset(tstamp), units='ns')
                    +datetime_utc).to_datetime64()
    except:
        print('WARNTING:  Could not convert {:} to local datetime'
              .format(datetime_utc))
        datetime = tstamp

    return datetime


def get_daystr_path(daystr, 
                    quiet=True, 
                    **kwargs):
    """
    Get thet path for the input day.
    This points to the files for the paper in the catemis/data folder,
    but this function will generally be facility specific.

    Returns
    -------
    str or Path
        Path name or object if as_path is set

    """
    from pathlib import Path
    if not isinstance(daystr, str):
        daystr = str(daystr)

    else:
        # make sure no - or / between YYYYMMDD format
        daystr = daystr.replace('-', '').replace('/','')

    year = daystr[0:4]
    month = daystr[4:6]
    day = daystr[6:8]
    datestr = '-'.join([year,month,day])

    try:
        path = Path(__file__).parent/'data'/daystr

    except:
        path = Path('.')/'data'/daystr

    if not path.is_dir():
        if not quiet:
            print('...  get_daystr_path for {:}...'.format(datestr))
            print('... WARNING: path {:} is not available'.format(str(path)))
            print(kwargs)
            print('...  get_daystr_path ...')
    
        return None

    return path


def get_daystr_files(daystr, path=None, globstr='*', **kwargs):
    """
    Get list of files for input day

    Parameters
    ----------
    daystr : str
        string giving the date in the format yyyymmdd
    path : Path object
        Path for files
    globstr : str
        Format of files to search under path

    """
    from pathlib import Path
    
    if not path:
        path = get_daystr_path(daystr, **kwargs)
    
    if path is None:
        return []
        
    path = Path(path)
        
    pfiles = list(path.glob(globstr))

    return pfiles


def _get_shot_min(complete=False, **kwargs):
    """
    Get minimum valid shot number.  
    Default is 1200 since lower numbered shots are only unique for 
    a specific day and only used for machine warmup or diagnostics tests
    
    If complete keyword is set then use all shots for that day

    Parameters
    ----------
    complete : bool
        Including warmup shots 

    """
    if complete:
        shot_min = 1

    else:
        shot_min = 1200
        
    return shot_min


def _set_dtypes(ds, 
        int_types=['shot'],
        bool_types=['pulsed_shot']):
    """
    Set data types.
    """
    import numpy as np
    for attr in ds.coords:
        try:
            if attr in bool_types:
                ds.coords[attr] = ds[attr].astype(bool)
    
            elif attr in int_types:
                ds.coords[attr] = ds[attr].astype(np.int64)
    
            elif isinstance(ds.coords[attr].values.flatten()[0], bytes):
                ds.coords[attr] = ds.coords[attr].astype(str)
    
            elif str(ds.coords[attr].dtype) == 'object':
                ds.coords[attr] = ds.coords[attr].astype(str)

        except:
            print('... failed setting dtype for {:}'.format(attr))

    for attr in ds.data_vars:
        if isinstance(ds[attr].values.flatten()[0], bytes):
            ds[attr] = ds[attr].astype(str)

    return ds


def set_byte_to_str(ds):
    """
    Set byte data types to unicode.
    This should be just a temporary fix until the netcdf readers
    are fixed to read this correctly.  Seems only to be an issue
    on Windows and not on MacOS.
    
    https://github.com/pydata/xarray/issues/2059
    
    """
    for attr in ds.coords:
        if isinstance(ds.coords[attr].values.flatten()[0], bytes):
            ds.coords[attr] = ds.coords[attr].astype(str)
    for attr in ds.data_vars:
        if isinstance(ds[attr].values.flatten()[0], bytes):
            ds[attr] = ds[attr].astype(str)

    return ds


def list_to_str(runs, portable=False):
    """
    Convert list of numbers (e.g., runs/shots) to string representation.
    e.g.,

    [3,4,5,9,12,13,14,15] becomes '3:5,9,12:15'

    Parameters
    ----------
    portable : bool
        Use '-' and '_' instead of ':' and ',' for use in sequential
        steps and individual numbers
    """
    import numpy as np
    
    if portable:
        strmulti = '-'
        strsingle = '_'

    else:
        strmulti = ':'
        strsingle = ','

    runs = np.array(sorted(set(runs)))
    runsteps = runs[1:] - runs[0:-1]
    runstr = '{:}'.format(runs[0])
    for i in range(len(runsteps)):
        if i == 0:
            if runsteps[i] == 1:
                grouped_runs = True

            else:
                grouped_runs = False

            if i == len(runsteps)-1:
                runstr += '{:}{:}'.format(strsingle, runs[i+1])

        elif i == len(runsteps)-1:
            if grouped_runs:
                if runsteps[i] == 1:
                    runstr += '{:}{:}'.format(strmulti, runs[i+1])

                else:
                    runstr += '{:}{:}{:}{:}'.format(strmulti, runs[i],
                                                    strsingle, runs[i+1])

            else:
                runstr += '{:}{:}:{:}'.format(strsingle, runs[i], runs[i+1])

        elif i > 0:
            if runsteps[i] == 1:
                if not grouped_runs:
                    runstr += '{:}{:}'.format(strsingle, runs[i])

                grouped_runs = True

            else:
                if grouped_runs:
                    runstr += '{:}{:}'.format(strmulti, runs[i])

                else:
                    runstr += '{:}{:}'.format(strsingle, runs[i])

                grouped_runs = False

    return runstr


def str_to_array(runstr, portable=False):
    """
    Convert string representing a sequence of numbers
    (e.g., run/shot) to list.

    '3:5,9,12:15' becomes [3,4,5,9,12,13,14,15]

    Parameters
    ----------
    portable : bool
        Use '-' and '_' instead of ':' and ',' for use in sequential
        steps and individual numbers

    """
    import numpy as np
    
    if portable:
        strmulti = '-'
        strsingle = '_'

    else:
        strmulti = ':'
        strsingle = ','

    runs = []
    for item in runstr.split(strsingle):
        runrange = item.split(strmulti)
        if len(runrange) == 2:
            for run in range(int(runrange[0]), int(runrange[1])+1):
                runs.append(run)

        else:
            runs.append(int(runrange[0]))

    return np.array(sorted(runs))
