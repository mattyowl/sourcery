"""

Custom error pages etc. - need to be set in config - see e.g. sourcery_test

"""

import os
import cherrypy
import sourcery
from pkg_resources import resource_filename

templatesDir=resource_filename('sourcery','templates')
staticDir=resource_filename('sourcery','templates')

def error_page(status, message, traceback, version):
    # We clear all query stuff if something goes wrong - user can then just refresh page
    if not cherrypy.session.loaded: cherrypy.session.load()
    cherrypy.session['viewTopRow']=0
    cherrypy.session['queryRADeg']="0:360"
    cherrypy.session['queryDecDeg']="-90:90"
    cherrypy.session['querySearchBoxArcmin']=""
    cherrypy.session['queryOtherConstraints']=""
    cherrypy.session.save()
    with open(templatesDir+os.path.sep+"error.html", "r") as inFile:
        lines=inFile.readlines()
        html=""
        for line in lines:
            # line=line.replace("$STATICDIR", self.staticDir)
            line=line.replace("$VERSION", sourcery.__version__)
            line=line.replace("$SCRIPT_NAME", cherrypy.request.script_name)
            html=html+line
        html=html % {'message': message, 'traceback': traceback, 'status': status}
    return html
