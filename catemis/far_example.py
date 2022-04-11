# -*- coding: utf-8 -*-
# Copyright (C) 2022 by Jason Koglin, Los Alamos National Laboratory
# All rights reserved. BSD 3-clause License.
# This file is part of the catemis package. Details of the copyright and
# user license can be found in the 'LICENSE' file distributed with the
# package.

"""
Examples for building calibration files
"""

def build_crx_default_file(pyro='FMP2', pyro_version=None, **kwargs):
    """
    The load_fmp_calibration method is used to build the dcalib0 calibration waveforms
    These are just the inverse of the crx files used by FAR
    but with a different normalization. 
    The default crx files for FMP1 and FMP2 were built using this method
    """
    from pathlib import Path
    default_path = Path(__file__).parent
    da0 = load_fmp_calibration(pyro, pyro_version=pyro_version, **kwargs)
    dcrx0 = (1/da0.reset_coords().dcalib0).to_dataframe()
    crx_file = str(default_path/'data'/(pyro+'.crx'))
    dcrx0.to_csv(crx_file, sep='\t', header=False)
    
    
def load_fmp_calibration(pyro, 
                         icalib0=None,
                         spec_files=None,
                         npar=2,
                         norm0=1,
                         ispec_max=20,
                         offset=0,
                         pyro_version=None,
                         refit_dcalib=True,
                         dcalib0=None,
                         dtrans0=None,
                         dspec0=None,
                         **kwargs):
    """
    Load FAR FMP1 and FMP2 calibration data 
    
    2019-11-16:  Initially calibration of 2018-2019 data was made 
    based on scaling FMP1 calibration temperature performed 
    in July 2016 from 1150 C to 1093 C 
    (i.e., pyro_path/'b890cal/Sample/16070712/Temp100.DAT')
    so that it would match the calibration of FMP2 at 1150 C
    (i.e., fmp2_path/'b309cal/Sample/16070611/Temp025.DAT').
    After processing the 2018-2019 FMP1 and FMP2 data in this way,
    the two pyrometers generally had the same operational temperature 
    of ~1150 C, but showed a clear bias where FMP2 changed faster than
    FMP1 around this operational temperature.
    
    2019-11-18:  Subsequently the July 2016 FMP1 calibration temperature
    of 1150 C was assumed and the temperature for the
    '18060409/Temp93071.DAT' FMP2 spectra was set to be 1200 C 
    and used to infer the FMP2 calibration.
    This temperature of 1200 C is consistent with the FMP1 temperature
    obtained from the FMP1 spectra '18060409/Temp12238060.DAT' using the 
    2016 FMP1 calibration.

    FMP1 Calibration data notes:
    - 2012-11-11 measurements at 800 and 1008 C in far/one folder 
    have consistent calibration
    - 2012-11-13 measurement in far/b890cal folder needs to be adjusted 
    back to 1008 C to be consistent with previous measurments
    - 2016-07-07 measurement in far/b890cal folder needs is 1150 C 
    which is consistent with other calibration measurements
    - 2012-10-29 measurement in far/one folder needs to be 3080 C 
    (2610 C in log, which appears to have uncalibrated numbers)
    
    Relevant files have been copied into darhtio/far/FMP1 
    and darhtio/far/FMP2 subfolders to keep in git repository
    
    See also examples/far_pyro_calibration.py for example of loading
    and plotting data.
    
    Parameters
    ----------
    offset : float
        Temperature offset
        
    Examples
    --------
    >>> from catsim import far_pyro
    >>> ds1cal = far_pyro.load_fmp_calibration(pyro='FMP1')
    
    """
    from pathlib import Path
    import numpy as np
    import xarray as xr
    
    from . import planck, butter_filter, far_pyro

    if pyro_version is None:
        pyro_version = ''
    
    if spec_files is None:
        if pyro == 'FMP2':
            fmp2_path = Path(__file__).parent/'far'/pyro/'Sample'
            spec_files = [
                (1004, fmp2_path/'06022100/Temp010.DAT'), 
                (1488, fmp2_path/'06022814/Temp1648.DAT'), # orig calibration 
                (1970, fmp2_path/'06022814/Temp020.DAT'), # From LA1-far.xls
                (1101, fmp2_path/'06031400/Temp001.DAT'),  
                (1982, fmp2_path/'11031612/Temp001.DAT'), 
                (1150, fmp2_path/'11071100/Temp001.DAT'), 
                (1150, fmp2_path/'12050817/Temp13013.DAT'), 
                (1150, fmp2_path/'16070611/Temp025.DAT'), 
                (1150, fmp2_path/'16070600/Temp001.DAT'), 
                (1150, fmp2_path/'18043009/Temp246688.DAT'), # 2018 cal
                (1150, fmp2_path/'18060409/Temp93071.DAT'),  # 2018 cal
                ]
            if not icalib0:
                icalib0 = 7

        elif pyro == 'FMP1':
            fmp1_path = Path(__file__).parent/'far'/pyro/'Sample'
            spec_files = [
                (1005.6, fmp1_path/'06120000/Temp001.DAT'), # 2006 original 
                (1008.5, fmp1_path/'12111300/Temp001.DAT'),  
                (800, fmp1_path/'12111315/Temp53528.DAT'),
                (1008, fmp1_path/'12111315/Temp240.DAT'),  # 2012 cal 
                (1150, fmp1_path/'16070712/Temp100.DAT'), # 2016 cal
                ]
            if not icalib0:
                icalib0 = 4
            
        else:
            raise Exception('ERROR:  Only pyro = FMP1 or FMP2 valid')

    elif icalib0 is None:
        raise Exception('ERROR:  icalib0 must be supplied with spec_files')
    
    adaf = []
    for items in spec_files:
        temp0 = items[0]
        spec_file = items[1]
        try:
            daf = far_pyro.load_pyro_spectrum(str(spec_file), 
                     pyro=pyro, 
                     pyro_version=pyro_version).to_xarray()
            daf.coords['cal_temp'] = (('calibration'), [temp0+offset])
            adaf.append(daf)

        except:
            import traceback
            traceback.print_exc()
            print('cannot load {:}'.format(spec_file))
    
    dscal = xr.concat(adaf, dim='calibration').dropna(dim='wavelength')
        
    nwv = dscal.wavelength.size
    ncal = dscal.calibration.size
    
    for attr in ['planck','planck_popt','planck_pcov','par0','par1']:
        if attr in dscal:
            del dscal[attr]
    
    dscal['trans'] = (('calibration','wavelength'), np.ones((ncal,nwv), dtype=float))
    dscal['corr'] = (('calibration','wavelength'), np.zeros((ncal,nwv)))
    dscal['ecorr'] = (('calibration','wavelength'), np.zeros((ncal,nwv)))
    dscal['dcalib'] = (('calibration','wavelength'), np.zeros((ncal,nwv)))
    dscal['planck'] = (('calibration','wavelength'), np.zeros((ncal,nwv)))
    dscal['planck0'] = (('calibration','wavelength'), np.zeros((ncal,nwv)))    
    dscal['pfit0'] = (('calibration','wavelength'), np.zeros((ncal,nwv)))    
    dscal['planck_popt'] = (('calibration','par0'), np.zeros((ncal,npar)))
    dscal['planck_pcov'] = (('calibration','par0','par1'),
                            np.zeros((ncal,npar,npar)))
    dscal.coords['calibration_file'] = (('calibration'), \
                [str(items[1]) for items in spec_files])

    for ical, items in enumerate(spec_files):
        temp0 = items[0]
        spec_file = items[1]
        
        if len(items) > 2:
            trans_factor = items[2]

        else:
            trans_factor = 1.
            
        if isinstance(trans_factor, xr.DataArray):
            trans_factor = trans_factor.interp_like(dscal.wavelength) \
                    .ffill('wavelength').bfill('wavelength')
            dscal['trans'][ical] *= trans_factor

        else:
            dscal['trans'][ical] *= trans_factor

            
    # update to be similar to far_termperature
    if 'wave_fit' not in dscal.coords:
        if pyro == 'FMP1':
            dscal.coords['wave_fit'] = (abs(dscal.wavelength-1380)>80) \
                                     & (dscal.wavelength>1000)

        else:
            if pyro_version.startswith('FMP4'):
                dscal.coords['wave_fit'] = (dscal.wavelength>600) \
                                         & (dscal.wavelength<980)

            else:
                dscal.coords['wave_fit'] = (dscal.wavelength>600) \
                                         & (dscal.wavelength<980)

    if dspec0 is not False and ispec_max:
        dscal.coords['dspec0'] = dscal.raw \
            .isel(wavelength=slice(0,ispec_max)).mean(dim='wavelength')

    else:
        dscal.coords['dspec0'] = (('calibration'), np.repeat(0,ncal))
        
    dscal.coords['dspec0'].attrs['long_name'] = pyro+' spectrum offset'

    if dtrans0 is None:
        dtrans0 = dscal['trans'][icalib0]

    else:
        dtrans0 = dtrans0.interp_like(dscal.wavelength)
        
    dscal['dtrans0'] = dtrans0
    dscal['dtrans0'].attrs['long_name'] = pyro+' calibration tranmission'

    for i in range(ncal):
        da = dscal.raw[i]
        dtrans = dscal.trans[i]

        temp0 = float(dscal.cal_temp[i])
        pfit0 = planck.planck_distribution(da.wavelength.values*1e-9, 
                                           temp0+273, norm0)
        dcalib = ((da-da.dspec0)/pfit0)/dtrans*dtrans0
        if pyro == 'FMP2':
            highcut = 0.08
            dacm = butter_filter.butter_bandpass_filter(dcalib, 
                                                        highcut=highcut)
            filt_cut = (dscal.wavelength>500)
            dcalib = dacm.where(filt_cut).bfill('wavelength')
            
        dscal['dcalib'][i] = dcalib    
    
    if dcalib0 is not None:
        dcalib0 = dcalib0.interp_like(dscal.wavelength)
        
    else:
        dcalib0 = dscal['dcalib'][icalib0]

    dscal.coords['dcalib0'] = (('wavelength'), dcalib0.values)
    dscal.coords['dcalib0'].attrs['long_name'] = pyro+' calibration'
    dscal.coords['decalib0'] = dscal.dcalib0 * dscal.eraw / dscal.raw
    
    # Fit calibration spectra using 
    for i in range(ncal):
        da = dscal.set_coords(['eraw','background']).raw[i].squeeze()
        dtrans = dscal.trans[i]
        da.coords['dcalib'] = da.coords['dcalib0']/dtrans*dtrans0
        da.coords['decalib'] = da.coords['decalib0']/dtrans*dtrans0
        daf = far_pyro.far_temperature_fit(da, npar=npar, temp0=da.cal_temp.values)
        dscal['planck'][i] = daf.planck_fit
        dscal['planck_popt'][i] = daf.popt
        dscal['planck_pcov'][i] = daf.pcov
        norm = daf.popt[1]
        pfit0 = planck.planck_distribution(da.wavelength.values*1e-9, 
                                           da.cal_temp.values+273, norm)
        dscal['pfit0'][i] = pfit0
        rawnorm = (da-da.dspec0).where(da.wave_fit).sum().values
        fitnorm = (dscal.pfit0[i]*dcalib0).where(da.wave_fit).sum().values
        dscal['planck0'][i] = float(rawnorm/fitnorm)*pfit0*dcalib0 \
            +da.dspec0.astype(float).values

    # Redo calibration with temperature fits assuming dcalib0
    if refit_dcalib:
        for i in range(ncal):
            da = dscal.raw[i]
            dspec0 = float(da.isel(wavelength=slice(0,ispec_max))
                           .mean(dim='wavelength'))
            tempK = float(dscal.planck_popt[i][0])
            norm = float(dscal.planck_popt[i][1])
            pfit0 = planck.planck_distribution(da.wavelength.values*1e-9, 
                                               tempK, norm)
            dscal['dcalib'][i] = ((da-dspec0)/pfit0).values

    dscal.coords['planck_temp'] = dscal['planck_popt'].sel(par0=0)-273
    dscal.coords['planck_norm'] = dscal['planck_popt'].sel(par0=1)
    dscal.coords['planck_temp'].attrs['units'] = 'C' 
    dscal.coords['planck_temp'].attrs['long_name'] = 'Fit Temperature' 

    dscal.coords['daystr'] = (('calibration'), 
                               [Path(f).parts[-2] for f 
                                in dscal.calibration_file.values])
    dscal.coords['year'] =   (('calibration'),
                               [2000+int(a[0:2]) for a 
                                in dscal.daystr.values])
    dscal.coords['month'] =  (('calibration'),
                               [int(a[2:4]) for a 
                                in dscal.daystr.values])
    dscal.coords['day'] =    (('calibration'),
                               [int(a[4:6]) for a 
                                in dscal.daystr.values])
    dscal.coords['file'] =   (('calibration'), 
                               [Path(f).parts[-1] for f 
                                in dscal.calibration_file.values])

    dscal.attrs['pyro'] = pyro
    dscal.attrs['icalib0'] = icalib0

    return dscal



    
