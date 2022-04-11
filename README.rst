=======
catemis
=======

Python package for DARHT Axis-II Dispenser-Cathode Emissivity and Temperature analysis

Analysis and data used in 2021 Weapons Engineering Symposium and Journal paper for
"DARHT Axis-II Dispenser-Cathode Emissivity and Temperature"

* LA-UR-21-25647 (paper)
* LA-UR-21-28231 (presentation)


Python Installation
-------------------

* Windows 10 Anaconda Python 3.8 Installation

     Download Anaconda (https://www.anaconda.com/download/)

     Click on Anaconda3-...-Windows-x86_64.exe file and intall with default settings

* Open 'Anaconda Prompt'  (Use 'Type here to search' or find in Start\Anaconda3)

* For new cond install, create .condarc file

     (base) C:\\Users\\ZZZZZZ> conda config --add channels defaults

* Locate .condarc file with where command (For MacOS, Linux and GitBash use 'which' instead of 'where')

     (base) C:\\Users\\ZZZZZZ> where .condarc

* Update proxy info as needed in .condarc file (WARNING:  spaces matter for this file)::


Code and Applications Installation 
----------------------------------

* Create and activate a conda environment

    (base) C:\\Users\\ZZZZZZ> conda create --name catemis python=3.8

    (base) C:\\Users\\ZZZZZZ> conda activate catemis

* Change directory to the folder where yout python packages are located (recommended to use 'python' in your base user directory; create directory if needed).

    (catemis) C:\\Users\\ZZZZZZ> cd python

* Install release from github repository using pip with the ‘-e’ flag to install a copy of the git repository code locally.

    (catemis) C:\\Users\\ZZZZZZ\\python> pip install -e git+ssh://git@github.com/lanl/catemis.git@v0.3.0#egg=catemis

* Example:  start the Cathode viewer application

     (catemis) C:\\Users\\ZZZZZZ\\python\\catemis\\panel panel serve --show "cathode_viewer.ipynb" --port=5051

* When using --show as in this example, the default browser will automatically open up with the URL

     http://localhost:5051/cathode_viewer


* In each new prompt window, activate the conda environment for the catemis package setup and start 
  other applications with a unique port number, e.g., for the FAR Pyrometer viewer application

     (base) C:\\Users\\ZZZZZZ\\python\\catemis\\panel  conda activate catemis

     (catemis) C:\\Users\\ZZZZZZ\\python\\catemis\\panel panel serve --show "far_viewer.ipynb" --port=5052


API Installation
----------------

The python packages included in the requirements.txt file are needed for the web based application.
A smaller set of packages defined in the requirement-api.txt file are needed to use the 
methods for loading and processing spectra and images. 
The xyzpy package can be optinally installed for parallel image loading in the cathode.load_catcams method.

License
-------
catemis is distributed as open-source software under a BSD 3-Clause License. LANL Copyright No. C21111. 
(see the ``LICENSE`` file for details).

© 2022. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.


