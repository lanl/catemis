# -*- coding: utf-8 -*-
# Copyright (C) 2022 by Jason Koglin, Los Alamos National Laboratory
# All rights reserved. BSD 3-clause License.
# This file is part of the catemis package. Details of the copyright and
# user license can be found in the 'LICENSE' file distributed with the
# package.

"""
FAR Pyrometer Analysis Methods
"""

def load_pyro(log_num=None, 
        pyro='FMP2',
        pyro_path=None, 
        latest=None,
        daystr=None,
        log_path=None,
        sample_path=None,
        spec_path=None,
        sequence=None,
        load_waveforms=True,
        status=None,
        **kwargs):
    """
    Load pyro data.
    Calculate two-band temperature using pyro_temperature_twoband method.
    
    Returns
    -------
    xarray.Dataset
        Pyrometer data
    
    Parameters
    ----------
    log_num : int or str
        Log number provides Log file (with .log extension)
        as well as the Sample folder name.
        Default is automatically found using daystr

    pyro : str
        Name of pyrometer.  Default = 'FMP2'
    
    daystr : str
        string that gives the date in the format yyyymmdd.

    pyro_path : Path
        Base name of pyrometer data path
       
    log_path : Path
        Path of log files.  Default is pyro_path/pyro/'Log'
       
    sample_path : Path
        Path of sample files.  Default is pryo_path/pryo/'Sample'

    spec_path : Path
        Path of spectrum files.  Default is sample_path/log_num

    load_waveforms : bool
        Load waveforms along with log file.
        Default = True

    latest : bool or int
        Get only latest values
        Default = False
        
    sequence : list
        List of sequence numbers to load.
        Default is to load all sequences in log file for day.
        
    ispec_max : int
        First ispec_max bins used for calculating background level.
        The mean over this range is used as dspec0 in far_temperature_fit
        
    status : pn.pane.Markdown.object, optional
        Used to provide status in panel apps.          

    """
    import numpy as np
    import pandas as pd
    import xarray as xr
    from tqdm import tqdm
    from pathlib import Path
    
    from . import utils

    kwargs.get('folder', None)
    kwargs.get('as_path', None)

    if not pyro_path:
        pyro_path = Path(__file__).parent/'far'
        
    else:
        pyro_path = Path(pyro_path)
        
    log_path = pyro_path/pyro/'Log'

    if daystr:
        if isinstance(daystr, int):
            daystr = str(daystr)

        year = daystr[0:4]
        month = daystr[4:6]
        day = daystr[6:8]
        datestr = '{:}-{:}-{:}'.format(year,month,day)        

    if not log_num:
        log_nums = sorted([f.parts[-1].split('.')[0] \
                           for f in list(log_path.glob('*.log'))], \
                          reverse=True)
        
        if daystr:
            df = None
            # Find the log number with the given datestr
            for ilog, log_num in enumerate(log_nums):
                log_name = '{:}/{:}.log'.format(log_path,log_num)
                try:
                    df, pyro_version = _get_log(log_name, pyro=pyro) 
                    if df is not None and datestr in df.Date.values:
                        break
                    
                except:
                    # catches bad log file
                    print('... missed error for log_num', log_num)
                
            if df is None or datestr not in df.Date.values:
                raise Exception('ERROR: {:} not found in {:} logs'.format(
                                datestr, pyro))

        else:
            log_num = log_nums[-1]
            latest = True

    if sample_path is None:
        sample_path = log_path.parent/'Sample'

    if spec_path is None:
        spec_path = sample_path/str(log_num)

    if isinstance(spec_path, str):
        spec_path = Path(spec_path)

    log_name = '{:}/{:}.log'.format(log_path,log_num)
    df, pyro_version = _get_log(log_name, pyro=pyro)

    if latest:
        if latest is True:
            datestr = str(df[-1:].Date.values[0])
            df = df.where(df.Date == datestr).dropna()
            df['Sequence'] = df.Sequence.astype(int)

        else:
            df = df[-latest:]

    if daystr:
        df = df.where(df.Date == datestr).dropna()
        df['Sequence'] = df.Sequence.astype(int)

    # Make into xarray DataArray
    dseq = df.to_xarray() \
             .set_coords(['Sequence']) \
             .swap_dims({'index':'Sequence'})
    dftime = pd.to_datetime(dseq.Date.values+'T'+dseq.Time.values)
    dseq.coords['datetime'] = (('Sequence'),dftime)
    dseq.coords['year'] = (('Sequence'),dftime.year)
    dseq.coords['year'].attrs['long_name'] = 'Year'
    dseq.coords['year'].attrs['units'] = 'Y'
    dseq.coords['month'] = (('Sequence'),dftime.month)
    dseq.coords['month'].attrs['long_name'] = 'Month'
    dseq.coords['month'].attrs['units'] = 'M'
    dseq.coords['day'] = (('Sequence'),dftime.day)
    dseq.coords['day'].attrs['long_name'] = 'Day'
    dseq.coords['day'].attrs['units'] = 'D'
    dseq['Temp'].attrs['long_name'] = 'Temperature'
    dseq['Temp'].attrs['units'] = 'C'

    dseq = dseq.set_coords(['Date','Time','Signal','Tol'])
    dseq.coords['Tol'].attrs['long_name'] = 'Tolerance'
    dseq.coords['Tol'].attrs['units'] = 'C'
    dseq.coords['Signal'].attrs['long_name'] = 'Signal'
    dseq.coords['Signal'].attrs['units'] = ''
    
    if pyro == 'FMP2':
        dseq['Average'].attrs['long_name'] = 'Average Temp'
        dseq['Average'].attrs['units'] = 'C'
        
    if 'BB' in dseq:
        dseq = dseq.set_coords(['BB'])
        dseq.coords['BB'].attrs['long_name'] = 'Blackbody Type'    
    
    import time
    ahour = []
    for Time in dseq.Time.data:
        x = time.strptime(Time.split('.')[0],'%H:%M:%S')  
        ahour.append(x.tm_hour+x.tm_min/60+x.tm_sec/60**2)

    dseq.coords['hours'] = (('Sequence'), ahour)
    dseq.coords['hours'].attrs['long_name'] = 'Hour of day'
    dseq.coords['hours'].attrs['units'] = 'h'
    dseq.coords['log_num'] = (('Sequence'),
                              np.repeat(int(log_num), dseq.Sequence.size))
    dseq.attrs['pyro_path'] = str(pyro_path)
    dseq.attrs['sample_path'] = str(sample_path)
    dseq.attrs['log_path'] = str(log_path)
    dseq.attrs['log_num'] = log_num
    dseq.attrs['pyro'] = pyro

    spec_cols = load_pyro_spectrum(pyro=pyro, 
                                   pyro_version=pyro_version, 
                                   get_spec_cols=True)
    usecols = load_pyro_spectrum(pyro=pyro, 
                                 pyro_version=pyro_version,
                                 get_usecols=True)
    dseq.attrs['spec_cols'] = spec_cols
    dseq.attrs['usecols'] = usecols
    
    if sequence:
        sel_seq = [seq for seq in sequence if seq in dseq.Sequence]
        if sel_seq:
            dseq = dseq.sel(Sequence=sel_seq)
    
    if not load_waveforms:
        return dseq
    
    if load_waveforms:
        ada = []
        aseq = []
        desc = pyro+' spectra'
        if daystr:
            desc += ' '+daystr
            
        iseq = 0
        nseq = dseq.Sequence.size
        for seq in tqdm(dseq.Sequence.values, desc=desc):
            spec_name = spec_path/'Temp{:03}.dat'.format(seq)
            iseq += 1
            if status is not None and np.mod(iseq, 5) == 0:
                status = '... Loading Sequence {:} of {:}'.format(iseq, nseq)
            
            if spec_name.is_file():
                da = load_pyro_spectrum(spec_name, pyro_version=pyro_version, 
                                        pyro=pyro)
                # check to make sure not empty
                if da.shape[0] > 0:
                    ada.append(da.to_xarray())
                    aseq.append(seq)
    
        if aseq == []:
            print('Failed to find and load {:} spectra data...'
                  .format(pyro))
            print('Returning only {:} temperature measurements'
                  .format(pyro))

            return dseq.swap_dims({'Sequence': 'datetime'})
    
        try:        
            if not daystr or latest:
                print('merging...')

            ds = xr.concat(ada, dim='Sequence')
            ds.coords['Sequence'] = aseq
            ds.coords['wavelength'].attrs['long_name'] = 'Wavelength'
            ds.coords['wavelength'].attrs['units'] = 'nm'
            if 'raw' in ds:
                ds['raw'].attrs['long_name'] = 'Raw Spectrum'
            
            if 'eraw' in ds:
                ds = ds.set_coords('eraw')
                ds.coords['eraw'].attrs['long_name'] = 'Error in Raw Spectrum'
                
            ds['intensity'].attrs['long_name'] = 'Blackbody Intensity'
            ds['t_exposure'].attrs['long_name'] = 'Exposure Time'
            ds['t_exposure'].attrs['units'] = 's'

            if dseq.Sequence.size != ds.Sequence.size:
                print('... merging {:} of {:} with waveforms'.format(
                        ds.Sequence.size, dseq.Sequence.size ))
                dseq = dseq.sel(Sequence=ds.Sequence).merge(ds)

            else:
                dseq = dseq.merge(ds)

            dseq.attrs['pyro'] = pyro
            dseq.attrs['pyro_version'] = pyro_version
            
        except:
            print('Failed merging waveforms...')
            raise Exception(log_num)

        ds = dseq
        t0str = '2010-01-01'
        t0 = np.datetime64(t0str)
        ds.coords['time'] = ((ds.datetime-t0)/1.e9).astype(float)
        ds.coords['time'].attrs['long_name'] = 'Time since {:}'.format(t0str)
        ds.coords['time'].attrs['units'] = 's'
    
        # delete 'Date_cal' info to avoid merge error
        attrs = list(ds.data_vars)
        for attr in attrs:
            try:
                if ds[attr].dims[0] == 'Date_cal':
                    del ds[attr] 
            except:
                pass

        if 'Date_cal' in ds.coords:
            del ds.coords['Date_cal']
    
        ds = ds.swap_dims({'Sequence': 'datetime'})
        ds.coords['Time'].attrs['long_name'] = 'Time'
        ds.coords['Date'].attrs['long_name'] = 'Date'
        ds.coords['datetime'].attrs['long_name'] = 'Date'
        ds.coords['Sequence'].attrs['long_name'] = 'Sequence Number'
        ds.coords['hours'].attrs['long_name'] = 'Hour of day'
        ds.coords['hours'].attrs['units'] = 'h'

        adaystr = [''.join(date.split('-')) for date in ds.Date.values]
        ds.coords['daystr'] = (('datetime'), adaystr)
        ds.coords['Year'] = (ds.year+ds.month/100).astype(str)
        ds.coords['Day'] = ds.day.astype(str)

        ds.attrs['pyro_version'] = pyro_version
        ds = utils.add_metadata(ds, build_method='load_pyro')
            
        return ds
    
    
def _get_log(log_name, pyro):
    """
    Get pyrometer log file information.
    
    Find header rows that start with 'Seq'
    in some files there may be more than one if pyrometer was 
    restarted before beginning to take data.
    
    Parameters
    ----------
    log_name : str or Path
        Full file name of log file
        
    pyro : str
        Pyrometer type (either 'FMP1' of 'FMP2')
        
    Returns
    -------
    df : pandas.DataFrame
        Pyrometer log data
        
    pyro_version : str
        Pyrometer version from log data
    
    """
    import pandas as pd
    
    dlines = []
    pyro_version = None
    with open(log_name) as f:
        for iline,line in enumerate(f):
            if line.startswith('Seq'):
                dlines.append(iline)
                
            elif line.startswith('Version'):
                pyro_version = line.split()[1]

    if len(dlines) == 0:
        return None, None

    skiprows = dlines[-1]+2

    if pyro == 'FMP1':
        names=['Sequence','Date','Time','Temp','Tol','Signal']

    elif pyro == 'FMP2':
        if pyro_version == 'FMP4_45L':
            # 2016 FMP2 pyrometer installed at DARTH
            # on 2021-05-05 
            names=['Sequence','Date','Time','Temp',
                   'Tol','Signal','Average',
                   'V0','V1','V2','V3','V4']

        else:
            names=['Sequence','Date','Time','Temp',
                   'Tol','Signal','BB','Average']

    else:
        raise('pyro must be either "FMP1" or "FMP2"')

    # Load in table of sequence info
    df = pd.read_csv(log_name, names=names, skiprows=skiprows, \
                     delimiter="\s+")    

    # Skip last line if starts with 'Session Ended'
    dcut = df[['Sequence']].apply(lambda x: x.astype(str).str.isnumeric())
    df = df.where(dcut.Sequence).dropna()
    df.Sequence = df.Sequence.astype(int)
    df.Temp = df.Temp.astype(float)

    if pyro_version is None:
        pyro_version = ''

    return df, pyro_version


def pyro_temperature_twoband(ds, 
                             twoband_pars=None, 
                             low_slice=None,
                             high_slice=None,
                             ispec_max=20,
                             **kwargs):
    """
    Two band pyrometer approximation.
    e.g., FMP1 twoband_pars = [ 108.1, -245.2,  987.3,  256.1]
    e.g., FMP2 twoband_pars = [367.8, 321.0]
    
    Parameters
    ----------
    twoband_pars : list
        polynomial's coefficients (in decreasing powers) used
        to transform low/high spectrum ratio to temperature
        using np.poly1d one-dimensional polynomial class
        
    low_slice : slice
        wavelength slice used to calculate low band sum
        Defaults:  FMP1 = slice(1275,1300),  FMP2 = slice(750,800)

    high_slice : slice
        wavelength slice used to calculate high band sum
        Defaults:  FMP1 = slice(1610,1630), FMP2 = slice(880,900)
        
    """
    import numpy as np
    
    from . import utils

    pyro = ds.attrs['pyro']
   
    if twoband_pars is None:
        print('... See docs for starting twoband_pars values for FMP1/FMP2')
        raise Exception('ERROR:  twoband_pars must be supplied.')

    if pyro == 'FMP1':
        if not low_slice:
            low_slice = slice(1275,1300)

        if not high_slice:
            high_slice = slice(1610,1630)

    elif pyro == 'FMP2':
        if not low_slice:
            low_slice = slice(750,800)

        if not high_slice:
            high_slice = slice(880,900)

    else:
        if not low_slice and not high_slice:
            print('... unknown pyro = {:}'.format(pyro))
            raise Exception('ERROR: low_slice and high_slice required')

    spec_back = ds.raw.isel(wavelength=slice(0,ispec_max)) \
        .mean(dim='wavelength')
    spec_corr = (ds.raw-spec_back)
    
    ds['spec_sum_low'] = spec_corr.sel(wavelength=low_slice) \
                                  .sum(dim='wavelength')
    ds['spec_sum_low'].attrs['long_name'] = 'Corrected low band'
    ds['spec_sum_low'].attrs['start_wavelength'] = low_slice.start
    ds['spec_sum_low'].attrs['stop_wavelength'] = low_slice.stop
        
    ds['spec_sum_high'] = spec_corr.sel(wavelength=high_slice) \
                                  .sum(dim='wavelength')
    ds['spec_sum_high'].attrs['long_name'] = 'Corrected high band'
    ds['spec_sum_high'].attrs['start_wavelength'] = high_slice.start
    ds['spec_sum_high'].attrs['stop_wavelength'] = high_slice.stop

    ds['spec_lh'] = ds.spec_sum_low/(ds.spec_sum_high)
    ds['spec_lh'].attrs['long_name'] = 'Corrected low/high band ratio'

    ds['twoband_temp'] = (ds.spec_lh.dims, 
                          np.poly1d(twoband_pars)(ds['spec_lh']))
    ds['twoband_temp'].attrs['long_name'] = 'Temperature (2-bands)'
    ds['twoband_temp'].attrs['units'] = 'C'
    ds['twoband_temp'].attrs['low_slice_start'] = low_slice.start 
    ds['twoband_temp'].attrs['low_slice_stop'] = low_slice.stop
    ds['twoband_temp'].attrs['high_slice_start'] = high_slice.start 
    ds['twoband_temp'].attrs['high_slice_stop'] = high_slice.stop
    ds['twoband_temp'].attrs['twoband_pars'] = twoband_pars

    ds = utils.add_metadata(ds, 'twoband_temp', 
                            build_method='pyro_temperature_twoband')

    return ds


def load_pyro_spectrum(spec_name=None,
                       usecols=None,
                       spec_cols=None,
                       pyro=None,
                       get_spec_cols=False,
                       get_usecols=False,
                       pyro_version=None,
                       **kwargs):
    """
    Load in a single pyro spectrum file into a pandas DataFrame
    
    Parameters
    ----------
    spec_name : str
        Complete name of spectrum file
        
    usecols : list
        List of column numbers to use in loading spectrum tab delimited file.
        
    spec_cols : list
        List of names of usecols columns for creating pandas.DataFrame.
        
    pyro : str
        FAR pyrometer name (e.g., FMP2)
    
    get_spec_cols : bool
        Returns spec_cols for specified pyro
        
    get_usecols : bool
        Returns usecols for specified pyro
        
    pyro_version : str
        Pyrometer version used to select default spec_cols and usecols 
        for version of FMP2 that start with FMP4.
    """
    import numpy as np
    import pandas as pd
    
    if not spec_cols and not pyro:
        pyro='FMP2'
    
    if not pyro_version:
        pyro_version = ''
    
    if pyro == 'FMP1':
        if spec_cols is None:
            spec_cols = ['wavelength','raw_intensity','background',
                         'intensity','t_exposure']

        if usecols is None:
            usecols=[0,1,2,3,4]


    if pyro == 'FMP2' and pyro_version.startswith('FMP4'):
        if spec_cols is None:
            spec_cols = ['wavelength','raw_intensity','background',
                         'intensity','far_emissivity']

        if usecols is None:
            usecols=[0,1,2,3,4]

    elif pyro == 'FMP2':
        if spec_cols is None:
            spec_cols = ['wavelength','raw_intensity',
                         'background','intensity',
                         'col5','col6','col7',
                         't_exposure','col9']
                
        if usecols is None:
            usecols=[0,1,2,3,7]

    else:
        if spec_cols is None:
            raise('provide spec_cols or set pyro = "FMP1" or "FMP2"')
    
    if get_spec_cols:
        return spec_cols
    
    if get_usecols:
        return usecols
    
    if spec_name is None:
        raise('ERROR, must provide spec_name')
    
    da = pd.read_csv(spec_name, names=spec_cols, skiprows=1, 
             sep='\t', usecols=usecols, index_col=0)

    da['raw'] = da.raw_intensity - da.background
    da['eraw'] = da.background.pipe(np.sqrt)

    if pyro == 'FMP2' and pyro_version.startswith('FMP4'):
        try:
            with open(spec_name, 'r') as f:
                first_line = f.readline()

            vals0 = first_line.rstrip('\n').split('\t')

            da['t_exposure'] = np.repeat(float(vals0[5]), da.shape[0])
            da['internal_tempA'] = np.repeat(float(vals0[6]), da.shape[0])
            da['internal_tempB'] = np.repeat(float(vals0[7]), da.shape[0])
            da['internal_tempC'] = np.repeat(float(vals0[8]), da.shape[0])

        except:
            print('WARNING:  Failed adding internal temps for FMP2 version {:}'
                  .format(pyro_version))
            
    return da


def far_temperature(ds, pyro=None, 
                    dcalib=None,
                    npar=2, 
                    ispec_max=None,
                    temp_min=300, 
                    emissivity=None,
                    wave_fit=None,
                    pyro_version='',
                    status=None,
                    **kwargs):
    """
    Calculate FAR temperature from pyrometer Dataset output by 
    load_pyro method.
    Defaults for dcalib and emissivity arrays are provided
    for DARHT pyro = FMP1 or FMP2.
    Uses far_temperature_fit to fit temperature.
    
    Notes
    -----
    If wave_fit data_var is not in ds then the wavelength range 
    in nm will be used:
    FMP2:  600 - 900 (500 - 900 for new FMP4 pyrometer version):  
    FMP1:  1000 - 1330 and 1430 - 1636
        
    Total ranges [nm]:
    FMP2:  454 - 966
    FMP1:  875 - 1636

    Parameters
    ----------    
    pyro : str
        Pyrometer name.  'FMP1' or 'FMP2'

    dcalib : xarray.DataArray
        calibration vs wavelength.  Default from far_pyro.load_fmp_calibration
        
    emissivity : xarray.DataArray
        emissivity vs wavelength.  Default from cathode.load_emissivity
        if emissivity is False then no emissivity will be used, 
        otherwise the 
        
    wave_fit : xarray.DataArray
        wavelength fit mask.
        
    pyro_version : str
        pyrometer version is automaticaly provided in load_pyro method
        for newer FAR pyrometers.
        
    ispec_max : int
        First ispec_max bins used for calculating background level.
        The mean over this range is used as dspec0 in far_temperature_fit

    Returns
    -------
    xarray.Dataset
        

    """
    import traceback
    import numpy as np
    import pandas as pd
    import xarray as xr
    from tqdm import tqdm
    
    from . import utils

    
    for attr in ['planck','planck_popt','planck_pcov','par0','par1']:
        if attr in ds:
            del ds[attr]
            
    if not pyro:
        pyro = ds.attrs.get('pyro')
        if not pyro:
            raise Exception('pyro keyword must be supplied or in attrs')
    
    if not pyro_version:
        pyro_version = ds.attrs.get('pyro_version')
        if pyro_version is None:
            pyro_version = ''
            
    # Set ispec_max for additional background subtraction
    # for the old FMP2 pyrometer, which has no pyrometer version
    # number in the data files.
    # The background measurement in FMP1 and the new FMP2 
    # with pyro_version = 'FMP4_45L' have more accurate background
    # levels in the readout spectra and do not need additional
    # background subtraction
    if ispec_max is None and pyro == 'FMP2' and pyro_version == '':
        ispec_max = 20
        
    if dcalib is None:
        estr = 'ERROR: dcalib must now be supplied for planck far_temperature'
        raise Exception(estr)
        
    else:
        if isinstance(dcalib, pd.DataFrame) or isinstance(dcalib, pd.Series):
            dcalib = dcalib.to_xarray()
            
        if isinstance(dcalib, xr.Dataset) and 'dcalib' in dcalib:
            dcalib = dcalib.dcalib
            
        ds.coords['dcalib'] = (('wavelength'), dcalib.interp_like(ds).values)

    if emissivity is None:
        # no emissivity will be used
        demiss = np.zeros(ds.wavelength.size)
        
    else:
        if isinstance(emissivity, pd.DataFrame) \
                or isinstance(emissivity, pd.Series):
            emissivity = emissivity.to_xarray()
            
        if isinstance(emissivity, xr.Dataset) and 'emissivity' in emissivity:
            emissivity = emissivity.emissivity
            
        demiss = emissivity.sel(wavelength=ds.wavelength, method='nearest')
        # if emissivity set then change dcalib and set dcalib0 to original
        # calibration without emissivity correction
        ds.coords['emissivity'] = (('wavelength'), demiss.values)
        ds.coords['dcalib0'] = (('wavelength'), ds.dcalib.values)
        ds.coords['dcalib'] = ds.emissivity*ds.dcalib


    if 'decalib' not in ds.coords:
        ds.coords['decalib'] = ds.dcalib * ds.eraw / ds.raw
    
    if wave_fit is None:
        ds.coords['wave_fit'] = (('wavelength'), np.ones(ds.wavelength.size))
        
    else:
        if isinstance(wave_fit, pd.DataFrame) \
                or isinstance(wave_fit, pd.Series):
            wave_fit = wave_fit.to_xarray()
            
        if isinstance(wave_fit, xr.Dataset) and 'wave_fit' in wave_fit:
            wave_fit = wave_fit.wave_fit
            
        ds.coords['wave_fit'] = wave_fit.sel(wavelength=ds.wavelength, 
                                             method='nearest')
        ds.coords['wave_fit'] = ds.coords['wave_fit'] \
            .astype(bool).fillna(False)

    nwv = ds.wavelength.size
    ndt = ds.datetime.size
    tdim = ds.raw.dims[0]

    if ispec_max:
        dspec0 = ds.raw.isel(wavelength=slice(0,ispec_max)) \
            .mean(dim='wavelength')
        ds.coords['dspec0'] = dspec0
    else:
        ds.coords['dspec0'] = ((tdim), np.zeros((ndt)))

    ds.coords['dspec0'].attrs['long_name'] = pyro+' spectrum offset'

    ds['planck'] = ((tdim,'wavelength'), np.zeros((ndt,nwv)))
    ds['planck_popt'] = ((tdim,'par0'), np.zeros((ndt,npar)))
    ds['planck_pcov'] = ((tdim,'par0','par1'), np.zeros((ndt,npar,npar)))
    ds = ds.set_coords(['Temp'])
    
    nraw = ds[tdim].size
    for i in tqdm(range(nraw), desc='Fit {:} Planck'.format(pyro)):
        if status is not None and np.mod(i+1, 5) == 0:
            status = '... Calculating Planck {:} of {:}'.format(i, nraw)

        da = ds.raw[i]
        if da.Temp > temp_min:
            try:
                daf = far_temperature_fit(da, npar=npar)
                ds['planck'][i] = daf.planck_fit
                ds['planck_popt'][i] = daf.popt
                ds['planck_pcov'][i] = daf.pcov
                
            except:
                try:
                    iseq = int(da.Sequence[i])
                except:
                    iseq = -1
                print('... failed calibrating Sequence {:}'.format(iseq))
                traceback.print_exc()                

    ds.coords['planck_norm'] = ds['planck_popt'].sel(par0=1)
    ds.coords['planck_norm'].attrs['long_name'] = 'Fit Normalization' 
    ds.coords['planck_enorm'] = ds['planck_pcov'].sel(par0=1, par1=1) \
        .pipe(np.sqrt)
    ds.coords['planck_enorm'].attrs['long_name'] = 'Fit Normalization Error'  
        
    ds.coords['planck_temp'] = ds['planck_popt'].sel(par0=0)-273
    ds.coords['planck_temp'].attrs['units'] = 'C' 
    ds.coords['planck_temp'].attrs['long_name'] = 'Fit Temperature' 

    ds.coords['planck_etemp'] = ds['planck_pcov'].sel(par0=0, par1=0) \
        .pipe(np.sqrt)
    ds.coords['planck_etemp'].attrs['units'] = 'C' 
    ds.coords['planck_etemp'].attrs['long_name'] = 'Fit Temperature Error' 

    ds = utils.add_metadata(ds, 'planck_temp', 
                            build_method='far_temperature')    
    
    return ds


def far_temperature_fit(da, 
        norm=1,
        npar=2,
        dspec0=None,
        temp0=None,
        **kwargs):
    """
    Fit the temperature for a pyrometer DataArray used in 
    far_temperature method
    
    Coordinates and attributes of input da DataArray include
    - wavelength: in nm
    - eraw: error in raw measurement (from FAR output)
    - dcalib: calibration vs wavelength
    - decalib: estimate of error in calibration
    - wave_fit: bool array as cut for fitted wavelength regions
    defined by default in far_temperature
    - Temp: temperature measurement used as starting point C
    
    Fit algorithms based on npar fit parameter number
    - 1: planck_fit_pdf, planck_pdf
    - 2: planck_fit, planck_distribution
    - 3: planck_fit_offset, planck_distribution_offset
    
    Parameters
    ----------
    da : xarray.DataArray
        Pyrometer spectral array
    
    dspec0 : float or xarray.DataArray
        Offset applied to initail da spectrum
    
    norm : float
        Normalization factor
        
    npar : int
        Number of fit parameters.  Default=2
        Option 3 is used if not 1 or 2
    
    """
    import numpy as np
    from . import planck
    
    vwave0 = (da.wavelength*1e-9)
    if dspec0 is None:
        dspec0 = da.dspec0.astype(float).values

    if 'decalib' not in da.coords:
        da.coords['decalib'] = da.dcalib * da.eraw / da.raw
    
    dacorr0 = da-dspec0
    ydata0 = (1./norm)*(dacorr0/da.dcalib)
    sigma0 = (1./norm)*(np.sqrt(da.eraw**2/da**2 \
                     +da.decalib**2/da.dcalib**2) \
                     *da/da.dcalib)
    
    dcut = (ydata0.isnull() == False)
    
    if 'wave_fit' in da.coords and int(da.wave_fit.sum()) < da.wave_fit.size:
        # add mask cut if available
        dcut = (dcut) & (da.wave_fit>0)
        
    vwave0 = vwave0.where(dcut, drop=True)
    ydata0 = ydata0.where(dcut, drop=True)
    sigma0 = sigma0.where(dcut, drop=True)
    
    if temp0 is None:
        if 'Temp' in da.coords:
            temp0 = float(da.Temp.values)+273
        else:
            temp0 = 1400
        
    if npar == 1:
        p0 = [temp0]
        bounds = ([400], [5000])

        popt0, pcov0 = planck.planck_fit_pdf(vwave0,ydata0,sigma0, 
                             p0=p0,bounds=bounds)

        pfit0 = planck.planck_pdf(da.wavelength.values*1e-9, *popt0)

    elif npar == 2:
        p0 = [temp0, 1]
        bounds = ([400,0], [5000,10000])

        popt0, pcov0 = planck.planck_fit(vwave0,ydata0,sigma0, 
                                         p0=p0,bounds=bounds)

        pfit0 = planck.planck_distribution(da.wavelength.values*1e-9, 
                                                  *popt0)

    else:
        p0 = [temp0,1,0]
        bounds = ([400,0,-1e9], [5000,10000,1e9])

        popt0, pcov0 = planck.planck_fit_offset(vwave0,ydata0,sigma0, 
                                         p0=p0[0:npar],bounds=bounds[0:npar])

        pfit0 = planck.planck_distribution_offset(da.wavelength.values*1e-9, 
                                                  *popt0)

    da.coords['planck_fit'] = (('wavelength'), 
                               norm*pfit0*da.dcalib.values+dspec0)
    da.attrs['popt'] = popt0
    da.attrs['pcov'] = pcov0

    return da


def load_wave_fit(file=None, path=None):
    """
    Load FAR mask file

    Parameters
    ----------
    file : str
        file name
    
    path : str or Path
        file path.  Default is 'catemis/data'
    
    Examples
    --------
    
    >>> from catemis import far_pyro
    >>> da = far_pyro.load_wave_fit('FMP1_masc.csv')
    """
    import pandas as pd
    from pathlib import Path

    if path:
        file = Path(path)/file
        
    elif not Path(file).is_file():
        try:
            file = Path(__file__).parent/'data'/file
            
        except:
            file = Path('.')/'data'/file
            
    da = pd.read_csv(file, index_col=0)
            
    return da


def load_far_crx(crx_file=None, crx_path=None):
    """
    Load FAR pyrometer crx file

    Parameters
    ----------
    crx_file : str
        crx file name
    
    crx_path : str or Path
        crx file path.  Default is 'catemis/data'
    
    Examples
    --------
    
    >>> from catemis import far_pyro
    >>> da = far_pyro.load_far_crx('FMP2.crx')
    """
    import pandas as pd
    from pathlib import Path

    if crx_path:
        file = Path(crx_path)/crx_file
        
    else:
        try:
            file = Path(__file__).parent/'data'/crx_file
            
        except:
            file = Path('.')/'data'/crx_file
            
    spec_cols = ['wavelength','crx']
    usecols = [0,1]
    df = pd.read_csv(file, names=spec_cols, skiprows=1, 
             sep='\t', usecols=usecols, index_col=0)
    
    da = df.to_xarray()
    da['dcalib'] = 1/da.crx
    return da


def _correct_fmp2_time_with_fmp1(ds2, ds1,
                                 attr='twoband_temp',
                                 temp_thresh=900,
                                 quiet=True,
                                 **kwargs):
    """
    The FMP2 computer absolute time is not reliable.
    This method corrects the FMP2 time using FMP1.

    Parameters
    ----------
    ds2 : xarray.Dataset
        Pyrometer data that needs time correction
    
    ds1 : xarray.Dataset
        Pyrometer data to align to
    
    attr : str
        Temperature variable to use for aligning the two pyrometer times
        
    temp_thresh : number
        Temperature threshold to use for selecting temperatures points
        to align
    """
    from tqdm import tqdm
    import numpy as np
    import pandas as pd
    
    if 'hours0' not in ds2.coords: 
        # Save original hours time info for FMP2
        ds2.coords['hours0'] = ds2.hours.copy() 
        ds2.coords['hours0'].attrs['long_name'] = 'Hours (corrected)'
        ds2.coords['hours0'].attrs['units'] = 'h'    

    else:
        # Reset to original hours time info for FMP2
        ds2.coords['hours'] = ds2.hours0

    tdim = 'hours'
    tdim0 = ds2.raw.dims[0]
    nt = ds2[tdim].size
    ds2.coords['inds'] = ((tdim0), range(nt))
    ds2.coords['hshift'] = ((tdim0), np.zeros(nt))
    ds2.coords['hshift'].attrs['long_name'] = 'Time offset'
    ds2.coords['hshift'].attrs['units'] = 'h'    
    
    dates = sorted(list(set(ds2.Date.values) & set(ds1.Date.values)))
    ndates = len(dates)
    vshift0 = np.arange(-0.2,0.2,.001)
    atshift = np.zeros((ndates,vshift0.size))
    
    for idate, date in tqdm(enumerate(dates), desc='FMP2 time correct'):
        if idate >= 8888:
            continue
        ds1sel = ds1.swap_dims({tdim0:'Date'}).sel(Date=date)
        ds1sel = ds1sel.swap_dims({'Date':'hours'})
        ds2sel = ds2.where(ds2.Date==date, drop=True).copy()
        ds2sel = ds2sel.swap_dims({tdim0:'hours'})
        
        da2o = ds2sel[attr].where(ds2sel[attr]>temp_thresh, drop=True)
        da1o = ds1sel[attr].where(ds1sel[attr]>temp_thresh, drop=True)
        t2 = da2o.hours.max()-da2o.hours.min()
        t1 = da1o.hours.max()-da1o.hours.min()
        da1 = ds1sel[attr]
        
        if np.abs(t2-t1) < 0.3:
            d2scale = float(da1o.mean()/da2o.mean())
            hshift0 = float(da2o.hours[1]-da1o.hours[1])
        
        elif da1o.hours.min() > (da2o.hours.min()+1):
            ds2sel = ds2sel.where(ds2sel.hours>da1o.hours.min(), drop=True)
            da2o = ds2sel[attr].where(ds2sel[attr]>900, drop=True)
            d2scale = 1.
            hshift0 = float(da2o.hours[1]-da1o.hours[1])
            if hshift0 > 2:
                hshift0 = 0
            
        else:
            d2scale = 1.
            hshift0 = 0
    
        vtshift = vshift0+hshift0
            
        try:
            itry = 0
            try_again = True
            while try_again:
                ashift = {}
    
                for tshift in vtshift:
                    da2 = ds2sel[attr].sel(hours=(da1.hours+tshift),
                                           method='nearest')*d2scale
                    ashift[tshift] = float(np.sqrt(((da1-da2.values)**2)
                                                   .mean()))
    
                hshift = pd.Series(ashift).idxmin()
                if not quiet:
                    print(idate, date, itry, hshift, 
                          vtshift.min()+0.02, vtshift.max()-0.02)
    
                itry +=1
                if itry > 8:
                    hshift = hshift0
                    try_again = False
    
                elif hshift <= vtshift.min()+0.02:
                    vtshift = vshift0+hshift
    
                elif hshift >=  vtshift.max()-0.02:
                    vtshift = vshift0+hshift
    
                else:
                    try_again = False
    
            atshift[idate] = list(ashift.values())
    
        except:
            print('failed', idate, date)
            
        inds = ds2.where(ds2.Date == date, drop=True).inds
        ds2.coords['hshift'][inds] = hshift
        
    ds2['hours'] = ds2['hours0']-ds2['hshift']
    
    # correct time offset 
    ds2['datetime0'] = ds2.datetime.copy()
    ds2['datetime0'].attrs['long_name'] = 'Date (uncorrected)'
    
    doffset = [pd.DateOffset(hours=hshift) for hshift in ds2.hshift.values]
    ds2['datetime'] = [ds2['datetime0'].to_pandas()[i] - val \
                       for i,val in enumerate(doffset)]
    
    for attr in ['datetime','hours']:
        ds2[attr].attrs = ds1[attr].attrs
    
    return ds2


def add_ramp_flags(ds, 
        temp_thresh=900.,
        temp_min=650.,
        time_up=1.0,
        time_down=0.25,
        attr='twoband_temp', 
        quiet=False,
        **kwargs):
    """
    Add flags to data to roughly define ramped state.
    Calculates relative time in hours
    
    ramp_state = -1 for below thesh_min
    ramp_state = 0 for ramping up
    ramp_state = 1 for ramped to operational temperature
    ramp_state = 2 for ramping down
    ramp_state > 2 for unstable

    Parameters
    ----------
    attr : str
        Attribute that is used for temperature.  Default = 'twoband_temp'
    
    temp_thresh : float
        Temperature ramping threshold.  Default = 900 C
        
    temp_min : float
        Temperature minimum threshold.  Default = 650 C
        
    time_up : float
        Extra time after reaching theshold that is considered ramping up [h].
        Default = 1 h

    time_down : float
        Extra time before dropping back to theshold that is considered 
        ramping down [h].  Default = 0.25 h
        

    """
    from tqdm import tqdm
    import numpy as np

    if attr not in ds:
        estr = 'Unable to use {:} for FAR pyrometer flags'.format(attr)
        print('WARNING: {:}'.format(estr))
        attr = 'Temp'
        print('... trying {:}'.format(attr))
        if attr not in ds:
            estr = 'Unable to use {:} for FAR pyrometer flags.'.format(attr)
            raise Exception('ERROR: {:}'.format(estr))
    
    tdim = 'hours'
    tdim0 = ds.raw.dims[0]
    nt = ds[tdim].size
    ds.coords['inds'] = ((tdim0), range(nt))
    dates = ds['Date'].to_pandas().unique()
    
    ds.coords['ramp_state'] = ((tdim0), np.zeros(nt))
    ds.coords['ramp_state'].attrs['long_name'] = 'Ramp State'    
    
    ds.coords['hmin'] = ((tdim0), np.zeros(nt))
    ds.coords['hmin'].attrs['long_name'] = 'Min time'
    ds.coords['hmin'].attrs['units'] = 'h'
    
    ds.coords['hmax'] = ((tdim0), np.zeros(nt))
    ds.coords['hmax'].attrs['long_name'] = 'Max time'
    ds.coords['hmax'].attrs['units'] = 'h'
    
    for idate, date in tqdm(enumerate(dates), \
                            desc='add flags', disable=quiet):
        dsel = ds.where(ds.Date == date, drop=True).copy() 
        if dsel[attr].max() > temp_thresh:
            dselc = dsel.where(dsel[attr] > temp_thresh, drop=True)
            tdown = (dsel.hours >= dselc.hours.max()-time_down)
            tcut = (dsel.hours > dselc.hours.min()+time_up) \
                 & (dsel.hours < dselc.hours.max()-time_down)
            tunstable = tcut & (dsel[attr] < temp_thresh)
            tval = tcut*1+tdown*2+tunstable*3
            inds = dsel.inds.astype(int)
            ds.coords['ramp_state'][inds] = dsel.coords['ramp_state'] \
                        .where(dsel[attr] > temp_min).fillna(-1)
            ds.coords['ramp_state'][inds] = tval
            ds.coords['hmin'][inds] = dsel.hours.min()
            ds.coords['hmax'][inds] = dsel.hours.max()

    ds.coords['hrel'] = ds['hours']-ds['hmin']
    ds.coords['hrel'].attrs['long_name'] = 'Relative time'
    ds.coords['hrel'].attrs['units'] = 's'
    
    return ds
