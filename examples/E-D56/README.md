# Running sourcery: ACTPol E-D56 cluster catalog example

Here is an example of how to use sourcery to serve a web database
containing the [two-season ACTPol E-D56 cluster catalog](http://adsabs.harvard.edu/abs/2017arXiv170905600H).
This tutorial covers running the database locally on your machine,
using cherrypy's built-in webserver - see ADD_LINK for a description
of how to deploy sourcery on Apache.

## The configuration .yml file

All of the options for sourcery are controlled from a YAML file - 
this directory contains an example called `E-D56Clusters.yml`. This
file has comments that describe what each of the options does. We
will go through switching these on in turn, but for now, we will 
start with simply setting up and serving a database that anyone can
access.

## Creating the database

Sourcery works by creating a MongoDB collection (which we will call
a "database") from a number of input .fits table files. The primary
.fits table is specified using the `catalogFileName` parameter in 
the .yml file. This may then be cross-matched against other 
.fits catalogs and a MongoDB collection of editable fields 
(see later). To create the database, run:

```
% sourcery_build_db E-D56Clusters.yml
```

Building the database itself (in this case) takes ~2 seconds, but 
this script will then also proceed to build an image cache. As 
specified in the .yml file, this will download images from SDSS,
unWISE, PS1, and fetch information from NED also. It may take a while
to run (~45 minutes, depending on your network speed). If there are 
any network problems while this is running (e.g., an error message 
like 'Connection reset by peer', then you can resume building the 
cache from where it left off by running:

```
% sourcery_build_cache E-D56Clusters.yml
```

Once this is complete, you can test the database by running:

```
% sourcery_test E-D56Clusters.yml 
> Point your browser at http://localhost:8080/sourcery
```

You should see a table page like Fig. 1 below when you navigate to
<http://localhost:8080/sourcery> using your web browser:

![alt text](figs/table.png "Fig. 1: The index page for the E-D56 example.")

If you click on one of the links, you should find yourself presented
with a page that contains all of the data on the object in the catalog,
including images we specified, as in Fig. 2 below:

![alt text](figs/sourcepage.png "Fig. 2: A source information page for the E-D56 example.")

Note that the images pulled from the SDSS DR13 webserver are pretty dark
by default, and you will probably need to crank up the brightness using the 
slider, as shown.


