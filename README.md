# sourcery

A web application for interacting with astronomical source catalogues.

Displays the catalog as a table, with clickable links that open a page for each object. Each object page
contains an image plot (including support for multiple images from different surveys, and overplotting 
catalogs and contours), controls for editing user-defined database fields (e.g., comments on each object), 
a table of objects found in NED in the same region, and a full table of the object properties.

The catalog can be queried to select subsets of objects, and the output can be exported as a .fits table, 
.csv, or DS9 .reg (region) file.

## Software needed

Sourcery itself is written in python (2.7.x) and uses MongoDB (tested on 2.6.10) for data storage. In 
addition, it requires the following python modules to be installed (current versions used by the author are
given in brackets, earlier and later versions may also work):

* astLib (0.9.2)
* astropy (2.0.2)
* CherryPy (11.0.0)
* Cython (0.23.4)
* IPython (2.4.1; used for debugging only)
* matplotlib (1.5.1)
* NumPy (1.11)
* PIL (1.1.7; or Pillow 3.1.2)
* pyfits (3.4; will eventually be replace by astropy.io)
* pymongo (3.2)
* requests (2.18.4)
* SciPy (0.17.0)

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

## Running sourcery

Instructions for running an example to be added here...

## Deployment using Apache (on Ubuntu 16.04)

Note that this is quite clunky at present and may be tidied up...

1. Install all needed python packages (see above), MongoDB, and Apache 2.

2. Install mod_wsgi:

   ```
   sudo apt-get install libapache2-mod-wsgi
   ```

3. Install sourcery and build the MongoDB and cache

   ```
   sourcery_build_db webserver-sourcery.config
   sourcery_build_cache webserver-sourcery.config 
   ```

   An example .config and tutorial will be added later...

4. Copy the driver script `cgi-bin/webdb.wsgi` to a path where it is accessible by your webserver, and edit
it to point to your sourcery .config file.

5. Added a .conf file to your Apache 2 set-up that will enable cgi scripts to be executed from directory
where you installed the webdb.wsgi file. 

   For example:

   ```
   sudo vim /etc/apache2/conf-available/cgi-enabled.conf
   ```

   and paste into this file:

   ```
   <Directory "/home/yourusername/public_html/sourcery/cgi-bin">
   Options ExecCGI
   AddHandler wsgi-script .wsgi
   </Directory>
   ```

   and then:

   ```
   sudo a2enconf cgi-enabled
   sudo service apache2 reload
   ```

6. In your global site config (000-default-le-ssl.conf if using Let's Encrypt, otherwise 000-default.conf), 
above the `</VirtualHost>` line, add:

   ```
   WSGIDaemonProcess sourcery
   WSGIProcessGroup sourcery
   WSGIApplicationGroup %{GLOBAL}
   WSGIScriptAlias /sourcery /home/yourusername/public_html/sourcery/cgi-bin/webdb.wsgi
   ```

   The last line here isn't essential, but will enable your database to be accessed as, e.g.,
   http://www.example.com/sourcery (substituting your host name as appropriate).

7. Optionally (to stop astropy complaining about missing config files, although this doesn't seem to 
essential), edit /etc/apache2/envvars and add:

   ```
   ## For astropy
   export XDG_CONFIG_HOME=/var/www/astropyconfig
   export XDG_CACHE_HOME=/var/www/astropycache
   ```

   You will need to make these directories if they didn't exist, and change their ownership to www-data:www-data.

## Comments, bug reports, help, suggestions etc..

Please contact Matt Hilton <matt.hilton@mykolab.com>.
