# -*- coding: utf-8 -*-
# Copyright (C) 2022 by Jason Koglin, Los Alamos National Laboratory
# All rights reserved. BSD 3-clause License.
# This file is part of the catemis package. Details of the copyright and
# user license can be found in the 'LICENSE' file distributed with the
# package.

"""
Cathode Camera Emissivity and Response Model
"""

# This dictionary will be used as the default if cat_defaults
# keyword is not passed to get_cat_defaults in 
# load_cathode or load_catcam 
# Can be updated with with changes to catcam config versus date
cat_defaults0 = {
    20180101: {
        'desc': 'New cathode camera',
        'shot_null': 1,
        'daystr_null': '20190314',
        'ry': 53,
        'rx': -11,
        'y0': 630,
        'x0': 550,
        'dcat': 490,
        'orientation': 0,
        'focal': 100,
        'rcat': 6.5/2*25.4, # 6.5 inch diameter
        'Zcam': (30.7+10)*25.4,
        'psize': 0.00375,
        },
    }

def get_cat_defaults(daystr=None, cat_defaults=None, **kwargs):
    """
    Find the default settings based on date
    """
    idef = 0
    if daystr:    
        if not cat_defaults:
            cat_defaults = cat_defaults0

        while idef < len(cat_defaults0)-1 and \
                int(daystr) > list(cat_defaults.keys())[idef+1]:
            idef += 1
    
    return list(cat_defaults0.values())[idef]
    

def load_cathode(daystr,
        shot_null = None,
        daystr_null = None,
        complete=True,
        raw=False,
        **kwargs):   
    """
    Load a set of Cathode Camera images and related cathode information.
    Passes load_catacam parameters to load and process cathode camera images.

    Returns
    -------
    xarray.Dataset

    Parameters
    ----------
    daystr : str
        string that gives the date in the format yyyymmdd
        shot number may also be passed, which will load all shots
        for the day of that shot

    complete : bool
        Including warmup shots 

    raw : bool
        Return raw images.
        Default = False (return aspect corrected and gauss filtered)
        
    shot_null : int
        Shot number for null imaage subtraction
        
    daystr_null : str
        String that gives the date (yyyymmdd) for the shot_null
        
    ds : xarray.Dataset
        cathode dataset

    """    
    import numpy as np

    da = load_catcams(daystr, raw=raw, complete=complete, **kwargs)

    attrs = da.attrs
    
    params = get_cat_defaults(daystr, **kwargs).copy()
    if not shot_null:
        shot_null = params.get('shot_null')

    if not daystr_null:
        daystr_null = params.get('daystr_null')

    if daystr_null:
        da_null = load_catcam(shot=shot_null, daystr=daystr_null)
        da['CatNull'] = da_null.reset_coords(drop=True)
        da['CatNull'].attrs['long_name'] = 'Cathode Pedestal'
        da['CatNull'].attrs['units'] = 'ADU'
        da['CatNull'].attrs['shot_null'] = shot_null
        da['CatNull'].attrs['daystr_null'] = daystr_null
        
        # subtract null
        da = da-da['CatNull']
        da.attrs = attrs
        da.attrs['desc'] = 'Cathode after pedestal subtraction'
    
    else:
        print('No null subtracted')

    da.coords['CAT_R'] = (da.CAT_X**2+da.CAT_Y**2).pipe(np.sqrt)
    da.coords['CAT_theta'] = np.arctan2(da.CAT_Y, da.CAT_X)
    
    da.attrs['long_name'] = 'Cathode'
    da.attrs['units'] = 'ADU'

    ds = da.rename('CatCam').to_dataset()
    ds.attrs['CathodeArea'] = np.pi*(ds.CatCam.CatDiameter/10/2)**2 # cm2
    ds.attrs['CathodeArea_units'] = 'cm2'

    return ds


def add_pyro_data(ds, dpyro, **kwargs):
    """
    Add pyrometer data to cathode data using far_pyro.far_temperature
    method (passes any keyword argument to this method).
    Closest pyrometer measurement to cathode shot is used.
    
    Parameters
    ----------
    ds : xr.Dataset
        Cathode dataset, e.g., ds = cathode.load_cathode('20200701')
        
    dpyro : xr.Dataset
        Pyrometer dataset, e.g., dpyro = ds1 = far_pyro.load_pyro('20200701') 
    """
    from . import far_pyro
    
    pyro = dpyro.attrs['pyro']
    dpyros = dpyro.sel(datetime=ds.datetime, method='nearest')
    dpyros = far_pyro.far_temperature(dpyros, **kwargs)

    pyro_coords0 = ['Temp', 'eraw', 'emissivity',
                   'background', 'intensity', 'Tol', 'Signal', 'BB',
                   'dcalib', 'dcalib0', 'decalib', 'dspec0', 'wave_fit']
    pyro_coords0 = pyro_coords0 + ['twoband_temp', 'planck', 
                                   'planck_popt', 'planck_pcov',
                                   'planck_temp', 'planck_etemp',
                                   'planck_norm', 'planck_enorm']

    pyro_coords = [a for a in pyro_coords0 if a in dpyros]
    pyro_attrs = ['raw'] + pyro_coords
    pyro_dict = {attr: '{:}_{:}'.format(pyro, attr)
                 for attr in pyro_attrs+['wavelength']}
    dpyrosel = dpyros.reset_coords()[pyro_attrs] \
            .set_coords(pyro_coords) \
            .rename(pyro_dict)

    ds = ds.merge(dpyrosel)
            
    return ds


def calculate_catcam_temperature(ds, 
        calfactor=2.85e-15, 
        catcam_integration=0.080,
        cathode_radius_cm=None,
        vtemperature=None,
        **kwargs):
    """
    Calculate cathode temperature from cathode camera data

    Currently DAAAC camera integration setting is not available in 'flatfile'
    images and is not otherwise read in from any DAAAC settings.
    This is typically never changed and assumed to always be 80 msec.
    
    Returns
    -------
    xarray.Dataset

    Parameters
    ----------
    calfactor : float
        Calibration factor for acceptance and photons to ADU
    catcam_integration : float
        Camera integration time [sec]
    cathode_radius_cm : float
        Cathode radius [cm].  Default is calculated from ds.CatCam.CatDiameter
    vtemperature : array
        Temperature array for model.  Default = np.arange(100,4001,1)
    
    """
    import numpy as np
    from scipy import constants
    from . import planck

    # Load camera response
    df = read_camera_response(**kwargs)
    
    # Load emissivity model
    dsb = load_emissivity(**kwargs)
    
    if vtemperature is None:
        vtemperature = np.arange(800,1201,1)
    
    ds.coords['temperature'] = vtemperature
    ds.coords['temperature'].attrs['long_name'] = 'Temperature'
    ds.coords['temperature'].attrs['units'] = 'C'

    # Make planck model
    ds.coords['CatCam_wavelength'] = df.index.astype(float).values
    ds.coords['CatCam_wavelength'].attrs['long_name'] = 'Wavelength'
    ds.coords['CatCam_wavelength'].attrs['units'] = 'nm'
    ds.coords['CatCam_dwavelength'] = (('CatCam_wavelength'), 
                                       ds.CatCam_wavelength.pipe(np.gradient))
    ds['CatCam_planck_pdf'] = planck.planck_pdf(ds.CatCam_wavelength*1e-9, 
                                          ds.temperature+273)
    ds['CatCam_planck_pdf'].attrs['long_name'] = 'Planck Spectral Radiance'
    ds['CatCam_planck_pdf'].attrs['units'] = '$W/srad/m^3$'
    ds['CatCam_emissivity'] = dsb.emissivity \
                                 .rename({'wavelength': 'CatCam_wavelength'}) \
                                 .interp_like(ds.CatCam_wavelength)
                             
    ds.coords['QE_color'] = (('color', 'CatCam_wavelength'),
                             np.zeros((3, ds.CatCam_wavelength.size)))
    ds.coords['color'] = ['Blue', 'Green', 'Red']
    ds.coords['QE_color'].loc['Blue'] = df.B.values
    ds.coords['QE_color'].loc['Green'] = df.G.values
    ds.coords['QE_color'].loc['Red'] = df.R.values
    ds.coords['QE_color'].attrs['long_name'] = 'QE color CCD'
    ds.coords['QE_mono'] = (('CatCam_wavelength'), df.M.values)
    ds.coords['QE_mono'].attrs['long_name'] = 'QE mono CCD'
    
    ds['CatCam_Mean'] = ds.CatCam.where(ds.CatCut).mean(dim=['CAT_X','CAT_Y'])

    if not cathode_radius_cm:
        cathode_radius_cm = ds.CatCam.CatDiameter/2/10
    
    cathode_area = np.pi*(cathode_radius_cm)**2 # cm2

    ds['CatCam_pixelADU'] = (
            calfactor
            *catcam_integration
            *cathode_area/1e4
            *(ds['CatCam_wavelength']/1e9/constants.h/constants.c)
            *ds['CatCam_dwavelength']*1e-9
            *ds['QE_mono']
            *ds['CatCam_emissivity']
            *ds['CatCam_planck_pdf']
        ).sum(dim='CatCam_wavelength')
    
    ds['CatTemp'] = (('shot', 'CAT_Y', 'CAT_X'),
        np.interp(ds.CatCam, ds.CatCam_pixelADU, ds.temperature))
    ds['CatTemp'].attrs['long_name'] = 'Cathode Camera Temperature'
    ds['CatTemp'].attrs['units'] = 'C'
    ds.coords['CatTemp_mean'] = ds.CatTemp.where(ds.CatCut) \
                                  .mean(dim=['CAT_X','CAT_Y'])
    
    return ds


def load_catcams(daystr=None, 
                 files=None,
                 globstr='CAT_CAM.*.tif',
                 complete=False,
                 parallel=False,
                 num_workers=8,
                 **kwargs):
    """
    Load set of catcam images for a day
    
    Returns
    -------
    xarray.Dataset

    Parameters
    ----------
    daystr : str
        string that gives the date in the format yyyymmdd
    files : list
        list of files
    globstr : str
        Format of files to search under path
    parallel : bool
        Load image files in parallel using xyzpy
    num_workers : int
        Number of workers for parallel processing
        
    """
    import numpy as np
    from . import utils
    
    if daystr is None:
        raise Exception('daystr in yyyymmdd format must be provided')

    if isinstance(daystr, int):
        daystr = str(daystr)

    # make sure the max number of workers is less than 64
    # otherwise multiprocessing error occurs
    # at least in Python 3.8
    num_workers_max = 40
    num_workers = min([num_workers, num_workers_max])

    if not files:
        if not daystr:
            estr = 'If list of files is not supplied then daystr must be given'
            raise Exception('ERROR: {:}'.format(estr))
        
        files = utils.get_daystr_files(daystr, 
                                       globstr=globstr, 
                                       **kwargs)
        
    if not files:
        return

    atimes = {}
    for p in files:
        shot = int(p.parts[-1].split('.')[1])
        datetime_utc = np.datetime64(p.stat().st_mtime_ns, 'ns')
        datetime = utils.get_local_datetime(datetime_utc)
        atimes[shot] = datetime
 
    if parallel:
        print('... loading catcam images in parallel')
        try:
            import xyzpy
            
        except:
            print('WARNING: xyzpy is required for parallel load in load_cathode')
            print('... using standard sequential method to load images')
            parallel = False
                  
        print
        constants = {'daystr': str(daystr)}
        constants.update(**kwargs)
        combos = {'shot': list(atimes.keys())}

        r = xyzpy.Runner(load_catcam, var_names=['CAT_CAM'])
        ds = r.run_combos(combos,
                          parallel=parallel,
                          constants=constants,
                          num_workers=num_workers)

    if not parallel:
        from tqdm import tqdm
        import xarray as xr
        print('... loading catcam images sequentially')
        ads = []
        shots = list(atimes.keys())
        for shot in tqdm(shots):
            da = load_catcam(shot, daystr=daystr, **kwargs)
            ads.append(da)
    
        ds = xr.concat(ads, dim='shot')        
        ds.coords['shot'] = shots

    if not complete:
        shot_min = utils._get_shot_min(complete=complete)
        
        ds = ds.where(ds.shot>=shot_min, drop=True)
        
    ds.coords['daystr'] = (('shot'), [daystr for shot in ds.shot])
    ds.coords['day'] = ds.shot*0 + int(daystr[6:8])
    ds.coords['month'] = ds.shot*0 + int(daystr[4:6])
    ds.coords['year'] = ds.shot*0 + int(daystr[0:4])
    ds['datetime'] = (('shot'), [atimes.get(shot) 
                                 for shot in ds.shot.values])
    
    ds = ds.sortby('shot')

    attrs = load_catcam(daystr=daystr, get_attrs=True)
    ds.attrs.update(**attrs)

    # Make a cut to on the cathode that leaves off edge pixels
    xdim = ds.dims[1]
    ydim = ds.dims[2]
    ds.coords['CatCut'] = ((
                           (((ds/ds.mean(dim=[xdim,ydim])-1)*100)
                            .mean(dim='shot')
                            .rolling({ydim:3}).mean()
                            .rolling({xdim:3}).mean()
                           ) >0 )
                           .rolling({ydim:3}).mean()
                           .rolling({xdim:3}).mean()
                           ) == 1
    ds.coords['CatCut'].attrs['long_name'] = 'Cathode Image Cut'
    for attr in list(ds.attrs.copy()):
        if ds.attrs[attr] is None:
            del ds.attrs[attr]

    ds = utils.add_metadata(ds, build_method='catcam.load_catcams')
        
    return ds


def load_catcam(
        shot=None,
        fname=None,
        daystr=None,
        raw=False,
        params={},
        get_attrs=None,
        **kwargs):
    """
    Load DARHT Axis-II Cathode Camera image and apply aspect correction.
    
    The Cathode is viewed from above at a ~50 degree angle 
    with the camera horizontally offset from the beam axis by ~10 degrees.
    The Cathode camera image is transformed using the 
    skimage.transform.warp method 
    based on Euler angles using a rotation matrix generated with the 
    scipy.spatial.transform.Rotation.from_euler method.
    The camera pixel intensities are corrected inversely with 
    the square of the relative distance from a position on the cathode 
    to the camera focal plane using these same transformations. 
    
    Returns
    -------
    xarray.DataArray
        Load and automatically calibrate cathode image.
    
    Parameters
    ----------
    shot : number
        Shot number
    raw : bool
        Return raw images.
        Default = False (return aspect corrected and gauss filtered)
    rcat : float
        Cathode radius [cm].  Default = 20.32 cm (i.e., 8")
    orientation : double, optional
        Major axis orientation in clockwise direction as radians.
        Default = pi radians
    sigma : number
        Standard deviation for Gaussian kernel 
        used in scipy gaussian filter
    params : dict
        Image processing parameters.  Defaault from get_cat_defaults method
    desc : str
        Camera description:  Default = 'New cathode camera'
    shot_null : int
        Shot number for null imaage subtraction
    daystr_null : str
        String that gives the date (yyyymmdd) for the shot_null
    ry : number
        Rotaion in y axis [degrees].  Default = 53
    rx : number
        Rotation in x axis [degrees].  Default = -11 
    y0: number
        Cathode y center pixel.  Default = 630
    x0' : number
        Cathode x center pixel.  Default = 550
    dcat : number
        Distance to cathode camera [cm].  Default = 490
    orientation : number
        Rolll orientation of camera [degrees].  Default = 0
    focal : number
        Effective focal distance of camera [cm].  Default = 100
    rcat : number
        Radius of camera [cm].  Default = 6.5/2*25.4 (6.5 inch diameter)
    Zcam : number
        Z position of camera [cm].  Default = (30.7+10)*25.4,
    psize : number
        Camera pixel size [mm].  Default = 0.00375
        
    """
    import numpy as np
    import xarray as xr
    from scipy.spatial.transform import Rotation
    from skimage import transform
    
    from .image import open_image

    if not daystr:
        raise Exception('ERROR:  daystr must be supplied')

    if not params and daystr:
        params = get_cat_defaults(daystr, **kwargs).copy()
        
    params.update(**kwargs)

    sigma = params.get('sigma', 0)
    transpose = params.get('transpose', False)
    orientation = params.get('orientation', 0.)

    if not fname:
        fname = 'CAT_CAM.{:}.tif'.format(shot)

    ry = params.get('ry', 53)
        
    rx = params.get('rx', -11)
        
    focal = params.get('focal', 100)
        
    rcat = params.get('rcat', 6.5/2*25.4) # 6.5 inch

    width = rcat*2
    height = rcat*2
    
    Zcam = params.get('Zcam', (30.7+10)*25.4)

    extra_pixels = params.get('extra_pixels', 10)
        
    dbx = extra_pixels
    dby = dbx

    dcat = params.get('dcat', 490)
    
    psize = params.get('psize', 0.00375)
    
    x0 = params.get('x0', 550)
    
    y0 = params.get('y0', 630)

    Ycam = Zcam/np.arctan(np.deg2rad(ry))
    Xcam = Zcam*np.arctan(np.deg2rad(rx))

    attrs = {
        'cameraX': Xcam,
        'cameraY': Ycam,
        'cameraZ': Zcam,
        'focal': focal,
        'CatDiameter': rcat*2,
        'CatPixelDiameter': dcat,
        'cameraPixel': psize,
        'cameraPixelX0': x0,
        'cameraPixelY0': y0,
        'cameraElAngle': ry,
        'cameraAzAngle': rx,
        'orientation': orientation,
        'transpose': transpose,
        'filter_sigma': sigma,
        'shot': shot,
        'daystr': daystr,
        }

    if get_attrs:
        return attrs

    da = open_image(fname, daystr=daystr, **kwargs)
    if raw:
        da = da.rename({'x':'CAT_Xind', 'y':'CAT_Yind'})
        da.attrs['long_name'] = 'Cathode'
        da.attrs['units'] = 'ADU'
        return da

    im = da.values
    if transpose:
        im = im.T
    
    if orientation:
        im = transform.rotate(im, orientation, resize=True)
    
    output_shape = (int(im.shape[0]*2),
                    int(im.shape[1]*2))
    
    smatrix = np.array([[focal, 0, width/2],
                        [0, focal, height/2],
                        [0, 0, 1]])
    tmatrix = Rotation.from_euler('zyx', [0,ry,rx], degrees=True).as_matrix()
    rmatrix = smatrix*tmatrix*np.linalg.inv(smatrix)
    tf_img = transform.warp(im, rmatrix,
                            output_shape=output_shape, 
                            mode='wrap')
    
    dx = dcat+dbx*2
    dy = dcat+dby*2
    
    aimg = tf_img[y0-dby:y0+dcat+dby,
                  x0-dbx:x0+dcat+dbx]

    if sigma:
        from scipy.ndimage import gaussian_filter
        aimg = gaussian_filter(aimg, sigma, mode='nearest')
        
    vx = np.linspace(-width/2*dx/dcat, width/2*dx/dcat, dx)
    vy = np.linspace(-height/2*dy/dcat, height/2*dy/dcat, dy)
    da = xr.DataArray(aimg,
                      dims=['CAT_Y','CAT_X'],
                      coords={'CAT_Y':vy,'CAT_X':vx})
    da.attrs['long_name'] = 'Cathode'
    da.attrs['units'] = 'ADU'
    da.coords['CAT_Y'].attrs['long_name'] = 'Y'
    da.coords['CAT_Y'].attrs['units'] = 'mm'
    da.coords['CAT_X'].attrs['long_name'] = 'X'
    da.coords['CAT_X'].attrs['units'] = 'mm'

    # Effective Z distance relative to camera pointing from offset position
    daZ = da.CAT_Y*np.tan(np.deg2rad(-ry)) \
        + da.CAT_X*np.tan(np.deg2rad(-rx))
    Dcam2 = (Xcam-da.CAT_Y)**2+(Ycam-da.CAT_X)**2+(Zcam-daZ)**2
    da.coords['cameraD'] = np.sqrt(Dcam2)
    da.coords['cameraD'].attrs['long_name'] = 'Pixel to camera distance'
    da.coords['cameraD'].attrs['units'] = 'mm'

    # Normalize image
    da.coords['CatCam_corr'] = (float(Dcam2.sel(CAT_X=0, CAT_Y=0, 
                                method='nearest'))/Dcam2)
    da *= 2**16*da.CatCam_corr
    da.attrs.update(**attrs)

    return da


def read_camera_response(ccd_response_file=None, 
                         navg=3, 
                         **kwargs):
    """
    Read default ccd camera response shown in Figure 12 
    Data was taken from a figure using rough approximation.

    Returns
    -------
    pandas.DataFrame

    Parameters
    ----------
    filename : str
        Sensys filename.
    navg : integer
        Number of point to apply rolling average to response file.

    """
    from pathlib import Path

    if ccd_response_file:
        if not Path(ccd_response_file).is_file():
            raise Exception('ERROR:  Invalid Camera Response File {:}'
                            .format(ccd_response_file))

    else:
        basefile = 'ccd_response_typical.csv'
        try:
            filename = Path(__file__).parent/'data'/basefile

        except:
            filename = Path('.')/'data'/basefile
            
    import pandas as pd
    df = pd.read_csv(filename, index_col=0)

    if navg:
        df = df.rolling(navg).mean().fillna(0)

    return df/100


def load_emissivity(emissivity_file=None, 
                    aparam=0.35, 
                    bparam=0.24,
                    vtemperature_emissivity=None,
                    vwavelength_emissivity=None,
                    **kwargs):
    """
    Load starting point emissivity function and apply parametric correction 
    of the form
    
    emissivity = emissivity0*aparm + bparam
    
    where emissivity0 is the read from the csv file.
    The default emissivity_file provides the Spectroradiometer data at 25C 
    shown in Figure 3, and the default parameters applied to this function 
    provide the emissivity model in Figure 3.
    
    Returns
    -------
    xarray.Dataset

    Parameters
    ----------
    filename : str
        Filename for basic emissivity function in csv format
    aparam : number
        Scale factor applied to the emissivity0 function.  Default=0.35
    bparam : number
        Offset factor applied to the emissivity0 function.  Default=0.24
    vtemperature_emissivity : array
        Temperature array for model.  Default = np.arange(100,4001,1)
    vwavelength_emissivity : array
        Wavelength array for model.  Default = np.arange(10,20000,1)
    
    """
    import numpy as np
    import pandas as pd
    import xarray as xr
    from pathlib import Path
    
    from . import planck

    if emissivity_file:
        if not Path(emissivity_file).is_file():
            raise Exception('ERROR:  Invalid Emissivity File {:}'
                            .format(emissivity_file))
        
    else:
        basefile = 'emissivity.csv'
        try:
            emissivity_file = Path(__file__).parent/'data'/basefile
            
        except:
            emissivity_file = Path('.')/'data'/basefile
    
    dfemissivity = pd.read_csv(emissivity_file, index_col=0)
    
    if not vtemperature_emissivity:
        vtemperature_emissivity = np.arange(100,4001,1)
        
    if not vwavelength_emissivity:
        vwavelength_emissivity = np.arange(10,20000,1)
    
    dsb = xr.Dataset()
    dsb.coords['wavelength'] = vwavelength_emissivity
    dsb.coords['wavelength'].attrs['long_name'] = 'Wavelength'
    dsb.coords['wavelength'].attrs['units'] = 'nm'
    dsb['emissivity0'] = (('wavelength'), 
                         np.interp(dsb.wavelength,
                                   dfemissivity.index,
                                   dfemissivity.emissivity))
    dsb['emissivity0'] = dsb.emissivity0.rolling(wavelength=250).mean() \
            .shift(wavelength=-125).bfill(dim='wavelength') \
            .ffill(dim='wavelength')
    dsb['emissivity'] = dsb['emissivity0']*aparam+bparam 

    dsb.coords['temperature'] = vtemperature_emissivity 
    dsb.coords['temperature'].attrs['long_name'] = 'Temperature'
    dsb.coords['temperature'].attrs['units'] = 'C'
    dsb.coords['dwavelength'] = (('wavelength'), 
                                 dsb.wavelength.pipe(np.gradient))
    dsb['planck_pdf'] = planck.planck_pdf(dsb.wavelength*1e-9, 
                                          dsb.temperature+273)
    dsb['planck_pdf'].attrs['long_name'] = 'Planck Spectral Radiance'
    dsb['planck_pdf'].attrs['units'] = '$W/srad/m^3$'

    return dsb


def calculate_powermodel_tempearture(ds,
        cat_radiated=0.617, 
        **kwargs):
    """
    Apply cathode power model to the cathode camera data
    
    Power model = Cathode_Voltage * Cathode_Current * cat_radiated
    
    Returns
    -------
    xarray.Dataset
    
    Parameters
    ----------
    cat_radiated : float
        Fraction of cathode power radiated
    """
    import numpy as np

    dsb = load_emissivity()
    dsb.coords['model_power'] = np.pi*ds.CathodeArea/1e4*(
            dsb.planck_pdf*dsb.emissivity*dsb.dwavelength*1e-9
            ).sum(dim='wavelength') 
    dsb.coords['Cathode_power'] = dsb.model_power/cat_radiated
    
    dtemp = dsb.temperature \
        .swap_dims({'temperature': 'Cathode_power'}) \
        .interp(Cathode_power=ds.Cathode_VA).astype(float)
    
    ds['Model_temp'] = (('shot'), dtemp.values)
    ds['Model_temp'].attrs['long_name'] = 'Cathode Temperature Model'
    ds['Model_temp'].attrs['units'] = 'C'
    
    return ds


