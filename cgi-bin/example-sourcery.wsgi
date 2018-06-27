"""

Webserver script that runs sourcery.
Change the paths below to match your installation.

"""

import os
import sys
import numpy
import urllib
import matplotlib
import matplotlib.cbook
matplotlib.use("Agg")
import pylab
import cherrypy
from sourcery.sourceBrowser import SourceBrowser
import atexit

cherrypy.config.update({'environment': 'embedded', 'request.show_tracebacks': True,
                        'tools.sessions.on': True, 'tools.sessions.persistent': True,
                        'tools.sessions.storage_type': 'file',
                        'tools.sessions.storage_path': '/full/path/where/to/store/sessions',
                        'tools.sessions.name': 'sourcerydb_id',
                        'tools.sessions.clean_thread': True,
			'tools.auth.on': True})

if cherrypy.__version__.startswith('3.0') and cherrypy.engine.state == 0:
    cherrypy.engine.start(blocking=False)
    atexit.register(cherrypy.engine.stop)

conf=None
application=cherrypy.Application(SourceBrowser("/full/path/to/your/sourcery.yml", preprocess = False), script_name='/sourcery', config=conf)
