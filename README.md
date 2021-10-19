# Pyroms

Welcome to Pyroms!

Pyroms is a collection of tools to process input and output files
from the Regional Ocean Modeling System, [ROMS](https://www.myroms.org/). It was originally
started by Rob Hetland as a googlecode project, then he morphed it
into octant, also at googlecode. Frederic Castruccio then created a
fork and renamed it back to pyroms.

Pyroms is now hosted on GitHub.com in the [ESMG/pyroms](https://github.com/ESMG/pyroms) repository. This version is on the [python3](https://github.com/ESMG/pyroms/tree/python3) branch. It requires Python 3.4 or later.

## Installation

Pyroms is still a bit rough around the edges, particularly with regard to installation. Recent development has been done in Python environments managed by [Conda](https://docs.conda.io/en/latest/). However Pyroms itself cannot yet be installed with Conda.

If you are starting from scratch, we recommend that you install
[Anaconda](https://www.anaconda.com/) or
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) and create a Python 3 environment (as of December 2020, version 3.8 is your best bet) for Pyroms and your other scientific software. You should also consider making conda-forge your default channel. See the [conda-forge tips and tricks page](https://conda-forge.org/docs/user/tipsandtricks.html).

NB Joe McGovern - make conda-forge your default channel on conda

If you don't want to use Conda, that's fine, but you will have to do more of the work yourself.

## Prerequisites

The following are required and are all available from [Conda-Forge](https://conda-forge.org/).

   * Python >= 3.4 (Python 3.8 currently recommended for new environments)
   * [numpy](https://numpy.org/)
   * [scipy](https://www.scipy.org/)
   * [matplotlib](https://matplotlib.org/)
   * [basemap](https://matplotlib.org/basemap/)
   * [netcdf4](https://unidata.github.io/netcdf4-python/netCDF4/index.html)
   * [cftime](https://unidata.github.io/cftime/)
   * [lpsolve55]
     * conda install -c conda-forge lpsolve55
   * [pip](https://pypi.org/project/pip/)

If you experience an issue with the following dedent error (https://github.com/matplotlib/basemap/issues/494):
   * from matplotlib.cbook import dedent ==> "ImportError: cannot import name 'dedent'" in matplotlib
   * Open .../site-packages/mpl_toolkits/basemap/__init__.py at relevant breakpoint
   * replace from ... import statement with
     * from inspect import cleandoc as dedent
     * THIS WILL NEED TO BE REDONE AT EACH IDE LAUNCH DUE TO LACK OF OWNERSHIP CREDENTIALS FOR MATPLOTLIB/BASEMAP
   * Alternatively, execute: pip install -U matplotlib==3.2
   * If you get KeyError: 'PROJ_LIB' as part of basemap running/installation, see the following links
     * set os.environ['PROJ_LIB'] = '/home/pete/anaconda3/envs/pyroms_MI/share/proj'
       * or os.environ['PROJ_LIB'] = anaconda location + '/pyroms_MI/share/proj'
     * https://stackoverflow.com/questions/56313027/keyerror-proj-lib-while-installing-basemap-using-conda-install
     * https://stackoverflow.com/questions/59376657/proj-lib-error-in-spyder-with-importing-basemap
     * https://github.com/conda-forge/basemap-feedstock/issues/30#issuecomment-423512069

If you encounter a memory issue with meshgrid, try the following
replace
lon, lat = np.meshgrid(long, lat)
with
lon, lat = np.meshgrid(long, lat, sparse=True)

If you get OSError: no file with expected extension:
(related to pyroms.grid.Gridgen), in hgrid.py (compiler will probably open in debud mode)
replace
        self._libgridgen = np.ctypeslib.load_library('libgridgen.', pyroms.__path__[0]) 
with
        self._libgridgen = np.ctypeslib.load_library('libgridgen.so', pyroms.__path__[0]) 

uninstall and reinstall pyroms, pyroms_toolbox and bathy_smoother

If this doesn't work, try the following solution from Mike Hadfield, NZ
https://github.com/ESMG/pyroms/issues/26
1) conda install -c conda-forge gridgen
2) then add the following 
   a) np.ctypeslib.load_library('libgridgen','/home/pete/anaconda3/envs/pyroms_MI/lib')


The following is optional: Pyroms can be built and run without it but some of the functionality will be missing.

   * scrip, a Python implementation of [SCRIP](https://github.com/SCRIP-Project/SCRIP),
     the Spherical Coordinate Remapping and Interpolation Package. This is used by the pyroms
     module. The Python scrip code (a rather old version) is
     [bundled in pyroms](https://github.com/ESMG/pyroms/tree/python3/pyroms/external/scrip)
     and can be built and installed separately as described below. In future we plan to
     move from the bundled scrip code to a stand-alone package like
     [ESMF/ESMPy](https://www.earthsystemcog.org/projects/esmpy/) or
     [PySCRIP](https://github.com/dchandan/PySCRIP).

The following is optional and provides high-resolution coastlines for basemap:

   * [basemap-data-hires](https://anaconda.org/conda-forge/basemap-data-hires/)

## Install from source

To clone a copy of the source and install the pyroms packages, you can use the following commands
```
# Cd to a convenient directory
$ git clone https://github.com/ESMG/pyroms.git
$ pip install -e pyroms/pyroms
$ pip install -e pyroms/pyroms_toolbox
$ pip install -e pyroms/bathy_smoother
```

This installs three PIP packages with the names pyroms, pyroms\_toolbox and bathy\_smoother,
each with an [eponymous](https://en.wiktionary.org/wiki/eponymous) module.

An [editable-mode](https://pip.pypa.io/en/stable/reference/pip_install/#editable-installs) installation is recommended becauses it means changes you make to your copy of the source code will take effect when you import the modules. If you don't want this you can omit the "-e" option

The "pip install" command runs "python setup.py install" (or "python setup.py develop" with the "-e" switch) in each of the subdirectories listed. The "pip install" form is recommended because it allow easy removal (see below)

The above should work on most Linuces and on OSX with the system gcc and gfortran compilers.
They have also been verified to work in a Conda environment on Windows,
provided you install the
[m2w64-gcc](https://anaconda.org/msys2/m2w64-gcc) and [m2w64-gfortran](https://anaconda.org/msys2/m2w64-gcc-fortran) compilers.

## Install scrip

If you install as above and try to import the three Pyroms modules without having installed
scrip you will get a warning like this:

```
$ python
Python 3.8.5 | packaged by conda-forge | (default, Aug 29 2020, 01:22:49)
[GCC 7.5.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import pyroms
WARNING:root: scrip could not be imported. Remapping functions will not be available
>>> import pyroms_toolbox
>>> import bathy_smoother
```

The scrip module is not available via Conda or any other package repository and we are looking at alternatives. In the meantime, scrip can be built and installed from source as follows

```
# Start in the directory into which you cloned pyroms and cd to the SCRIP
# source directory
$ cd pyroms/pyroms/external/scrip/source/

# Print the location of the active Conda environment (which is called "python38"
# in this case). The active environment location is used to find the netCDF and
# other libraries.
$ conda info | grep "active env location"
    active env location : /home/hadfield/miniconda3/envs/python38

# make sure that all requirements in the scrip makefile are installed
Consider this link, if issues with libhdf detection....
https://stackoverflow.com/questions/50593848/error-building-caffe-cannot-find-lhdf5-hl

also 

sudo apt-get install libhdf5-serial-dev netcdf-bin libnetcdf-dev

and this needed to be added to scrip makefile 

HDF_LIBDIR    ?= $(PREFIX)/lib/x86_64-linux-gnu/hdf5/serial

and change
LIB = -L$(SRCDIR) -L$(NETCDF_LIBDIR) $(NETCDF_LIBS)
to
LIB = -L$(SRCDIR) -L$(NETCDF_LIBDIR) -L$(HDF_LIBDIR) $(NETCDF_LIBS)

Needed to replace ?= with :=

Tried 
sudo apt-get install libtool

Needed to create a new symlink in /usr/lib using the following command
sudo ln -s /usr/lib/libdf.so.0.0.0 /usr/lib/libdf.so

# Run make to build the scrip Python extension and install it into the Conda
# environment. The makefile calculates a variable called SCRIP_EXT_DIR, into
# which it installs the scrip Python extension. If pyroms has been installed
# in editable (development) mode, set the DEVELOP variable to a non-empty value.
$ export PREFIX=/home/hadfield/miniconda3/envs/python38
$ make DEVELOP=1 PREFIX=$PREFIX install
$ mv -vf scrip*.so ../../../pyroms
‘scrip.cpython-38-x86_64-linux-gnu.so’ -> ‘../../../pyroms/scrip.cpython-38-x86_64-linux-gnu.so’
```

## Removal

To remove the three Pyroms packages you can use the "pip uninstall" command, referring to the packages by their package names

```
# Run from any directory in the same environment as you installed
# and use the package name
$ pip uninstall pyroms
$ pip uninstall pyroms_toolbox
$ pip uninstall bathy_smoother
```

If you have built and installed the scrip extension from the makefile as above, you can also uninstall it with the makefile. The PREFIX does not need to be set in this case.

```
# Start in the directory into which you cloned pyroms and cd to the SCRIP
# source directory
$ cd pyroms/pyroms/external/scrip/source/

# Remove with make.
$ make DEVELOP=1 uninstall
```

## Running Pyroms

We have a gridid.txt file that's pointed to by the PYROMS\_GRIDID\_FILE
environment variable. If you are operating on files containing
sufficient grid information already, you won't need to use this.
An example is provided in the examples directory.


## Doxygen

Running "doxygen .doxygen" in any of pyroms, pyroms\_toolbox or
bathy\_smoother will generate doxygen files. Edit the .doxygen files to
specify html vs. some other output format.


# Notes JVMCGOVERN
# To get lpsolve working, Try:
# 1: (didn't work)
#  conda install -c snorfalorpagus lpsolve
# 2:
# https://stackoverflow.com/questions/14906603/how-to-install-lpsolve-module-for-python-on-linux-ubuntu-10-04
# https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_Python_source.tar.gz/download
# https://sourceforge.net/projects/lpsolve/files/lpsolve/5.5.2.11/lp_solve_5.5.2.11_dev_win64.zip/download
# 3:
# with intel mkl 2019 installed
# conda install -c conda-forge lpsolve55

# NB use the following: to ensure .pyd modules are recognised by pycharm
# https://intellij-support.jetbrains.com/hc/en-us/community/posts/360008700780-Pycharm-does-not-recognize-pyd-file-as-a-module
# https://intellij-support.jetbrains.com/hc/en-us/community/posts/360009758419-failed-to-load-pyd-package
# 'C:\\Program Files\\JetBrains\\PyCharm Community Edition 2020.2.3\\plugins\\python-ce\\helpers\\pydev'
# 'C:\\Program Files\\JetBrains\\PyCharm Community Edition 2020.2.3\\plugins\\python-ce\\helpers\\third_party\\thriftpy'
# 'C:\\Program Files\\JetBrains\\PyCharm Community Edition 2020.2.3\\plugins\\python-ce\\helpers\\pydev'
# 'C:\\Users\\jmcgovern\\Anaconda3\\envs\\CELTICSEA\\python36.zip'
# 'C:\\Users\\jmcgovern\\Anaconda3\\envs\\CELTICSEA\\DLLs'
# 'C:\\Users\\jmcgovern\\Anaconda3\\envs\\CELTICSEA\\lib'
# 'C:\\Users\\jmcgovern\\Anaconda3\\envs\\CELTICSEA'
# 'C:\\Users\\jmcgovern\\Anaconda3\\envs\\CELTICSEA\\lib\\site-packages'
# 'c:\\users\\jmcgovern\\pycharmprojects\\pyroms_mi\\pyroms'
# 'c:\\users\\jmcgovern\\pycharmprojects\\pyroms_mi\\pyroms_toolbox'
# 'c:\\users\\jmcgovern\\pycharmprojects\\pyroms_mi\\bathy_smoother'
# 'C:\\Users\\jmcgovern\\PycharmProjects\\pyroms_MI'
# 'C:/Users/jmcgovern/PycharmProjects/pyroms_MI'
# import sys
#
# sys.path.append('C:\\Users\\jmcgovern\\PycharmProjects\\pyroms_MI\\pyroms')
# sys.path.append('C:\\Users\\jmcgovern\\PycharmProjects\\pyroms_MI\\pyroms_toolbox')
# sys.path.append('C:\\Users\\jmcgovern\\PycharmProjects\\pyroms_MI\\bathy_smoother')

# when handling issues with gdal installation, use the following command in the terminal window
# conda install -c anaconda gdal