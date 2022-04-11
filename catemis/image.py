# -*- coding: utf-8 -*-
# Copyright (C) 2022 by Jason Koglin, Los Alamos National Laboratory
# All rights reserved. BSD 3-clause License.
# This file is part of the catemis package. Details of the copyright and
# user license can be found in the 'LICENSE' file distributed with the
# package.

"""
DAAAC Image Access
"""

def open_image(fname,
               ftype=None,
               name=None,
               shot=None,
               daystr=None,
               path=None, 
               img_split=None, flipud=None, fliplr=None,
               xrange=None, yrange=None,
               **kwargs):
    """
    Read camera image into xarray DataArray

    Parameters
    ----------
    fname : str
        file name
    path : str or Path object
        path name
    name : str
        name of returned DataArray.
        Automatically hosen from base filename
    img_split : (int,int)
        Split data into images
        if ftype is auto detected as simacon image then img_split = (4,4)
    flipud : bool
        Flip image up-down (default for simacon images)
    fliplr : bool
        Flip image left-right
    xrange : list
        x axis cut for image [xmin,xmax] 
    yrange : list
        y axis cut for image [ymin,ymax] 
        
    Returns
    -------
    xarray.DataArray
        Image Array

    Examples
    --------
    
    >>> from catemis import image
    >>> da = image.open_image('CAT_CAM.1.tif')

    """
    import imageio
    import xarray as xr
    import numpy as np
    from pathlib import Path

    fpath = Path(fname)
    pparts = fpath.parts[-1].split('.')
    if not name:
        name = pparts[0]

    try:
        if not shot:
            if pparts[1].isnumeric():
                shot = pparts[1]
                shot = int(shot)

    except:
        print('... cannot determine shot from file name {:}'.format(pparts))

    if path is not None or len(fpath.parts) == 1:
        from . import utils
        if path is None:
            if not daystr:
                raise Exception('ERROR:  daystr or path must be supplied')
            path = utils.get_daystr_path(daystr)
            
        if not isinstance(path, Path):
            path = Path(path)

        if not path.is_dir():
            raise Exception('ERROR: Invalid file path {:}'.format(
                str(path)))

        file_name = path/fpath
        
    else:
        file_name = fpath
        if not isinstance(file_name, Path):
            file_name = Path(file_name)
            
        path = file_name.parent
        
    fbase = file_name.parts[-1]
            
    if not file_name.is_file():
        raise Exception('ERROR: No image file {:}'.format(
            str(file_name)))

    fext = fname.split('.')[-1]
    if ftype == None and fext == 'spe':
        ftype = 'spe'
        
    reader = imageio.get_reader(str(file_name))
    header = reader.get_meta_data()
    compression = header.pop('compression', None)

    nimg = reader.get_length()
    if nimg == 1:
        data = reader.get_data(0)
        dims = ['y','x']
        
        if not ftype:
            if header.get('software') == 'SPECIALISED - OS':
                ftype = 'simacon'
            else:
                ftype = ''
    
        # check for Simicon and set split image and flip up-down
        if ftype == 'simacon':
            if img_split is None:
                img_split = (4,4)
            if flipud is None:
                flipud = True

        if xrange:
            data = data[:,xrange[0]:xrange[1]]

        if yrange:
            data = data[yrange[0]:yrange[1]]

        if flipud:
            data = np.flipud(data)

        if fliplr:
            data = np.fliplr(data)            
                
        if img_split:
            shape = data.shape
            ny = img_split[0]
            nx = img_split[1]
            nimg = nx*ny
            ysize = int(shape[0]/ny)
            xsize = int(shape[1]/nx)
            
            # sort images in order left to right, top to bottom
            adata = np.zeros((ny*nx,ysize,xsize))
            for ix in range(nx):
                for iy in range(ny):
                    adata[ix+(ny-iy-1)*nx] = \
                            data[iy*ysize:(iy+1)*ysize,
                                 ix*xsize:(ix+1)*xsize]
                            
            dims = ['image','y','x']
            data = adata

    else:
        print('... need to implement multi image reader')
        data = np.array([reader.get_data(i) for i in range(nimg)])
        dims = ['image','y','x']
        
    da = xr.DataArray(data, dims=dims, name=name)
    da.attrs.update(**header)
    da.attrs['shot'] = shot
    da.attrs['name'] = name
    da.attrs['ext'] = fext
    da.attrs['file_type'] = ftype
    da.attrs['path'] = str(path)
    da.attrs['file'] = str(fbase)

    return da
