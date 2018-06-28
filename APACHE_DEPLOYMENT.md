# Deployment using Apache (on Ubuntu 16.04)

Note that this is quite clunky at present and may be tidied up...

1. Install all needed python packages ([see README.md](README.md)), MongoDB, and Apache 2.

2. Install mod_wsgi:

   ```
   % sudo apt-get install libapache2-mod-wsgi
   ```

3. Install sourcery and build the MongoDB and cache

   ```
   % sourcery_build_db webserver-sourcery.yml
   % sourcery_build_cache webserver-sourcery.yml 
   ```

   See [examples/E-D56/README.md](examples/E-D56/README.md) for an example .yml file and a tutorial.

4. Copy the driver script `cgi-bin/example-sourcery.wsgi` to a path where it is accessible by your 
webserver (e.g., you could make a directory `/usr/local/www/wsgi-scripts/` to hold all of your .wsgi 
files), rename it to something appropriate, and edit it to point to your sourcery .yml file.

5. Add a .conf file to your Apache 2 set-up that will enable cgi scripts to be executed from directory
where you installed the .wsgi file. 

   For example:

   ```
   sudo vim /etc/apache2/conf-available/cgi-enabled.conf
   ```

   and paste into this file:

   ```
    <Directory /usr/local/www/wsgi-scripts>
    Require all granted
    Options ExecCGI
    AddHandler wsgi-script .wsgi
    <IfModule mod_rewrite.c>
            #Disable rewriting
            RewriteEngine Off
    </IfModule>
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
    # Essential for cherrypy apps to work
    WSGIApplicationGroup %{GLOBAL}

    # Mounting cherrypy apps
    WSGIScriptAlias /example-sourcery /usr/local/www/wsgi-scripts/example-sourcery.wsgi 
   ```

   The last line here isn't strictly essential, but will enable your database to be accessed as, e.g.,
   https://www.example.com/sourcery (substituting your host name as appropriate).

7. Optionally (to stop astropy complaining about missing config files, although this doesn't seem to 
essential), edit /etc/apache2/envvars and add:

   ```
   ## For astropy
   export XDG_CONFIG_HOME=/var/www/astropyconfig
   export XDG_CACHE_HOME=/var/www/astropycache
   ```

   You will need to make these directories if they didn't exist, and change their ownership to www-data:www-data.
