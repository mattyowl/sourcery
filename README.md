# Sourcery

A web application for interacting with astronomical source catalogues.

Sourcery displays an object catalog as a table, with clickable links that open a page for each source. 
Each source information page contains an image plot (with support for multiple images from different 
surveys, and overplotting catalogs and contours), controls for editing user-defined database fields 
(e.g., comments on each object), a table of objects found in NED in the same region, and a full table
of the object properties.

The catalog can be queried to select subsets of objects, and the output can be exported as a .fits table, 
.csv, or DS9 .reg (region) file.

To get an idea of how Sourcery works and see some screenshots, take a look at the 
[example](examples/E-D56/README.md).


## Software needed

Sourcery itself is written in Python (2.7.x - currently) and uses MongoDB (tested on 2.6.10) for data 
storage. In addition, it requires the following Python modules to be installed (current versions used
by the author are given in brackets, earlier and later versions may also work):

* astLib (0.10.1)
* astropy (2.0.4)
* CherryPy (11.0.0)
* Cython (0.23.4)
* IPython (5.5.0; used for debugging only)
* matplotlib (2.1.0)
* NumPy (1.14.1)
* PIL (1.1.7; or Pillow 3.1.2)
* pymongo (3.2)
* requests (2.18.4)
* SciPy (0.17.0)
* pyvips (2.0.4)
* PyYAML (3.12)
* passlib (1.7.1)
* concurrent.futures

## Installation

As root:
    
```
sudo python setup.py install
```

Or, in your home directory:
    
```
python setup.py install --prefix=$HOME/local
```

Then add `$HOME/local/bin` to $PATH, and e.g., `$HOME/local/lib/python2.7/site-packages` to $PYTHONPATH.

```
export PATH=$HOME/local/bin:$PATH    
export PYTHONPATH=$HOME/local/lib/python2.7/site-packages:$PYTHONPATH
```


## Running Sourcery

See the [example](examples/E-D56/README.md) for a tutorial, some screenshots of what to expect, and an 
explanation of the configuration file.

See [here](APACHE_DEPLOYMENT.md) for a step-by-step guide to deploying Sourcery on the Apache webserver
running on Ubuntu 16.04.


## Comments, bug reports, help, suggestions etc..

Please contact Matt Hilton <matt.hilton@mykolab.com>.
