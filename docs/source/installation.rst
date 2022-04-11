============
Installation
============

Install release from gitlab repository::

    $ pip install -e git+ssh://git@github.com/lanl/catemis.git@v0.3.0#egg=catemis

See the package README.txt for more details on 
how to install and configure Anaconda and the catemis package
including recommended way to setup a conda environment:
    
    - https://github.com/lanl/catemis

Required dependencies
---------------------

Anaconda Defaults:
    - Python 3.8 or higher
    - imageio: library for reading and writing image data
    - matplotlib: Comprehensive 2-D plotting (and some 3-D plotting)
    - numpy: Base N-dimensional array package
    - pandas: Data structures and plotting
    - scikit-image: image processing algorithms
    - scipy: Fundamental library for scientific computing

Additional primary packages:
    - xarray: N-D labeled arrays and datasets
    - h5netcdf: library for reading and writing netCDF4 files

Optional dependencies
---------------------

Optional dependencies for interactive web based applications.
    - notebook: Jupyter notebook
    - panel: library for creating custom interactive web apps and dashboards
    - param: library used for providing parameters used in panel
    - bokeh: library for creating interactive visualizations for web browsers
    - holoviews: data analysis and visualization library
    - hvplot: high-level plotting API built on HoloViews 
    - xyzpy: methods for parallel processing using xarray

Repository
----------
Git repository located at https://github.com/lanl/catemis

Please submit issues and bug reports through the github project page. 
For additional information and questions contact Jason Koglin at koglin@lanl.gov

