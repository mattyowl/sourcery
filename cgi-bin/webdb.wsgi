"""

Webserver script that runs sourcery.

"""

import os
import sys
import numpy
import urllib
#os.environ['MPLCONFIGDIR'] = "/home/xcs/matplotlib"
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
                        'tools.sessions.storage_path': '/home/xcs/sessions',
                        'tools.sessions.name': 'xcsdb_id',
                        'tools.sessions.clean_thread': True})

if cherrypy.__version__.startswith('3.0') and cherrypy.engine.state == 0:
    cherrypy.engine.start(blocking=False)
    atexit.register(cherrypy.engine.stop)

conf=None
application=cherrypy.Application(SourceBrowser("/home/xcs/sourcery.config", preprocess = False), script_name=None, config=conf)
