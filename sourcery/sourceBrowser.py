"""

    Copyright 2014-2024 Matt Hilton (matt.hilton@mykolab.com)
    
    This file is part of Sourcery.

    Sourcery is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    sourcery is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sourcery.  If not, see <http://www.gnu.org/licenses/>.

"""

import os
import matplotlib
matplotlib.use('Agg')
import astropy
import astropy.table as atpy
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import operator
import urllib3
import urllib
try:
    from urllib import quote_plus  # Python 2.X
except ImportError:
    from urllib.parse import quote_plus  # Python 3+
import glob
import tarfile
from astLib import *
import astropy.io.fits as pyfits
import numpy as np
import pylab as plt
import matplotlib.patches as patches
from scipy import ndimage
import sourcery
from sourcery import catalogTools
#from sourcery import specFeatures
#import ConfigParser
import yaml
import requests
import sys
import time
import datetime
import string
import re
import base64
from PIL import Image
import io
Image.MAX_IMAGE_PIXELS=100000001 
import copy
import json
from io import BytesIO
from tempfile import TemporaryDirectory
import tempfile
import pymongo
from bson.son import SON
import pyximport; pyximport.install()
#import sourceryCython
import cherrypy
import pickle
#import pyvips
import IPython
from sourcery import sourceryAuth
from sourcery import tileDir
from passlib.hash import pbkdf2_sha256
import logging
from jinja2 import Environment, FileSystemLoader

# Logging
logger=logging.getLogger('sourcery')

_TEMPLATES_DIR = os.path.join(os.path.dirname(__file__), 'templates')

#-------------------------------------------------------------------------------------------------------------
def makeDegreesDistanceMap(degreesMap, wcs, RADeg, decDeg, maxDistDegrees):
    """Fills (in place) the 2d array degreesMap with distance in degrees from the given position,
    out to some user-specified maximum distance.

    Args:
        degreesMap (:obj:`np.ndarray`): Map (2d array) that will be filled with angular distance
            from the given coordinates. Probably you should feed in an array set to some extreme
            initial value (e.g., 1e6 everywhere) to make it easy to filter for pixels near the
            object coords afterwards.
        wcs (:obj:`astWCS.WCS`): WCS corresponding to degreesMap.
        RADeg (float): RA in decimal degrees of position of interest (e.g., object location).
        decDeg (float): Declination in decimal degrees of position of interest (e.g., object
            location).
        maxDistDegrees: The maximum radius out to which distance will be calculated.

    Returns:
        A map (2d array) of distance in degrees from the given position,
        (min x, max x) pixel coords corresponding to maxDistDegrees box,
        (min y, max y) pixel coords corresponding to maxDistDegrees box

    Note:
        This routine measures the pixel scale local to the given position, then assumes that it
        does not change. So, this routine may only be accurate close to the given position,
        depending upon the WCS projection used.

    """

    x0, y0=wcs.wcs2pix(RADeg, decDeg)
    ra0, dec0=RADeg, decDeg
    ra1, dec1=wcs.pix2wcs(x0+1, y0+1)
    xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)

    xDistPix=int(round((maxDistDegrees)/xPixScale))
    yDistPix=int(round((maxDistDegrees)/yPixScale))

    Y=degreesMap.shape[0]
    X=degreesMap.shape[1]

    minX=int(round(x0))-xDistPix
    maxX=int(round(x0))+xDistPix
    minY=int(round(y0))-yDistPix
    maxY=int(round(y0))+yDistPix
    if minX < 0:
        minX=0
    if maxX > X:
        maxX=X
    if minY < 0:
        minY=0
    if maxY > Y:
        maxY=Y

    xDeg=(np.arange(degreesMap.shape[1])-x0)*xPixScale
    yDeg=(np.arange(degreesMap.shape[0])-y0)*yPixScale
    for i in range(minY, maxY):
        degreesMap[i][minX:maxX]=np.sqrt(yDeg[i]**2+xDeg[minX:maxX]**2)

    return degreesMap, [minX, maxX], [minY, maxY]

#-------------------------------------------------------------------------------------------------------------
class SourceBrowser(object):
    
    def __init__(self, configFileName, preprocess = False, buildDatabase = False):
              
        # Parse config file
        self.parseConfig(configFileName)
        
        # Access control (optional)
        if 'userListFile' in self.configDict.keys():
            with open(self.configDict['userListFile'], "r") as stream:
                self.usersDict=yaml.safe_load(stream)
            for userName in self.usersDict.keys():
                if 'role' not in self.usersDict[userName].keys():
                    raise Exception("'role' must be defined in user list file %s" % (self.configDict['userListFile']))
                if self.usersDict[userName]['role'] not in ['editor', 'viewer']:
                    raise Exception("unknown user role - check user list file %s" % (self.configDict['userListFile']))
        else:
            self.usersDict={}
        
        # Displayed when failed login
        if 'contactInfo' in self.configDict.keys():
            self.contactInfo=self.configDict['contactInfo']
        else:
            self.contactInfo=""

        # Logo image encoded as data URL for embedding in all page headers
        self.logo_data_url = None
        logo_path = self.configDict.get('logoFileName')
        if logo_path and os.path.exists(logo_path):
            ext = os.path.splitext(logo_path)[1].lower()
            mime_map = {'.png': 'image/png', '.jpg': 'image/jpeg', '.jpeg': 'image/jpeg',
                        '.gif': 'image/gif', '.svg': 'image/svg+xml', '.webp': 'image/webp'}
            mime = mime_map.get(ext, 'image/png')
            with open(logo_path, 'rb') as f:
                self.logo_data_url = 'data:%s;base64,%s' % (mime, base64.b64encode(f.read()).decode('ascii'))
            
        # Image choices
        if 'defaultImageType' not in self.configDict.keys():
            self.configDict['defaultImageType']='best'
        if 'imagePrefs' not in self.configDict.keys():
            self.configDict['imagePrefs']=['DECaLS', 'DES', 'SDSS', 'unWISE']
        
        # Add news into self.configDict, if there is any...
        if 'newsFileName' in self.configDict.keys():
            self.addNews()
        
        if 'skyviewPath' in self.configDict.keys():
            self.skyviewPath=self.configDict['skyviewPath']
        else:
            # Default - if we have run the sourcery_fetch_skyview script
            # Disabled until find a better way of doing this for e.g., apache on webserver
            self.skyviewPath=None#os.environ['HOME']+os.path.sep+".sourcery"+os.path.sep+"skyview.jar"

        # For DES usage - load credentials, if they are there...
        #if 'DESServicesConfigPath' in self.configDict.keys():
            #configParser=ConfigParser.RawConfigParser()   
            #DESConfigPath=self.configDict['DESServicesConfigPath']#os.environ['HOME']+os.path.sep+".desservices.ini"
            #configParser.read(DESConfigPath)
            #self.DESUser=configParser.get('db-dessci', 'user')
            #self.DESPasswd=configParser.get('db-dessci', 'passwd')
        #else:
            #self.DESUser=None
            #self.DESPasswd=None
        #self.DESTokenID=None    # Used for interacting with DESCuts server
        
        # Below will be enabled if we have exactly one image in an imageDir
        self.mapPageEnabled=False

        # Set up storage dirs
        self.cacheDir=self.configDict['cacheDir']
        self.skyCacheDir=self.configDict['skyviewCacheDir']
        if os.path.exists(self.cacheDir) == False:
            os.makedirs(self.cacheDir)
        self.nedDir=self.cacheDir+os.path.sep+"NED"
        if os.path.exists(self.nedDir) == False:
            os.makedirs(self.nedDir)
        self.sdssRedshiftsDir=self.cacheDir+os.path.sep+"SDSSRedshifts"
        if os.path.exists(self.sdssRedshiftsDir) == False:
            os.makedirs(self.sdssRedshiftsDir)
        
        # More storage dirs... (we can't make these on-the-fly when running threaded)
        sdssCacheDir=self.cacheDir+os.path.sep+"SDSS"
        os.makedirs(sdssCacheDir, exist_ok = True)
        decalsCacheDir=self.cacheDir+os.path.sep+"DECaLS"
        os.makedirs(decalsCacheDir, exist_ok = True)
        ps1CacheDir=self.cacheDir+os.path.sep+"PS1"
        if os.path.exists(ps1CacheDir) == False:
            os.makedirs(ps1CacheDir)
        ps1CacheDir=self.cacheDir+os.path.sep+"PS1IR"
        if os.path.exists(ps1CacheDir) == False:
            os.makedirs(ps1CacheDir)
        wiseCacheDir=self.cacheDir+os.path.sep+"unWISE"
        os.makedirs(wiseCacheDir, exist_ok = True)

        # tileDirs set-up - KiDS, IAC-S82 etc..
        self.tileDirs={}
        if 'tileDirs' in self.configDict.keys():
            for tileDirDict in self.configDict['tileDirs']:
                if tileDirDict['label'] not in self.tileDirs.keys():
                    self.tileDirs[tileDirDict['label']]=tileDir.TileDir(tileDirDict['label'], tileDirDict['path'], 
                                                                        self.cacheDir, sizePix = tileDirDict['sizePix'])            
        
        # Big redshifts table - let's try and keep it in memory (may be a challenge on the webserver)
        # Really we should just put into a database table...
        if self.configDict['specRedshiftsTable'] is not None:
            self.specRedshiftsTab=atpy.Table().read(self.configDict['specRedshiftsTable'])  
        else:
            self.specRedshiftsTab=None
        
        # So we can display a status message on the index page in other processes if the database or cache is being rebuilt
        self.dbLockFileName=self.cacheDir+os.path.sep+"db.lock"
        self.cacheLockFileName=self.cacheDir+os.path.sep+"cache.lock"

        # Add descriptions of field (displayed on help page only)
        self.descriptionsDict=self.parseColumnDescriptionsFile()

        # MongoDB set up
        self.dbName=self.configDict['MongoDBName']
        if 'TagsDBName' not in self.configDict.keys():
            self.tagsDBName=self.dbName
        else:
            self.tagsDBName=self.configDict['TagsDBName']
        self.client=pymongo.MongoClient('localhost', 27017)
        self.mongoSess=self.client.start_session()
        self.db=self.client[self.dbName]
        self.sourceCollection=self.db['sourceCollection']
        self.sourceCollection.create_index([('loc', pymongo.GEOSPHERE)])
        self.fieldTypesCollection=self.db['fieldTypes']
        self.tagsDB=self.client[self.tagsDBName]
        self.tagsCollection=self.tagsDB['tagsCollection']
        self.tagsCollection.create_index([('loc', pymongo.GEOSPHERE)])
        if buildDatabase == True:
            self.buildDatabase()

        # Column to display info
        # Table pages
        self.tableDisplayColumns=[{'name': "name",
                                   'label': "Name",
                                   'fmt': "%s"},
                                  {'name': "RADeg",
                                   'label': "RA (deg)",
                                   'fmt': "%.6f"},
                                  {'name': "decDeg",
                                   'label': "Dec. (deg)",
                                   'fmt': "%.6f"}]
        for colDict in self.configDict['tableDisplayColumns']:
            for key in colDict.keys():
                dispDict={'name': key, 'label': key, 'fmt': colDict[key]}
                self.tableDisplayColumns.append(dispDict)

        # Support for tagging, classification etc. of candidates
        if 'fields' in self.configDict.keys():
            for fieldDict in self.configDict['fields']:
                dispDict={'name': fieldDict['name'], 'label': fieldDict['name']}
                if fieldDict['type'] == 'number':
                    dispDict['fmt']='%.3f'
                elif fieldDict['type'] == 'text':
                    dispDict['fmt']='%s'
                else:
                    raise Exception("only valid field types are 'number' and 'text'")
                if 'tableAlign' in fieldDict:
                    dispDict['tableAlign']=fieldDict['tableAlign']
                if 'displaySize' in fieldDict:
                    dispDict['displaySize']=fieldDict['displaySize']
                if 'showOnIndexPage' in fieldDict and fieldDict['showOnIndexPage'] == True:
                    self.tableDisplayColumns.append(dispDict)
                
        if 'classifications' in self.configDict.keys():
            dispDict={'name': "classification", 'label': "classification", 'fmt': "%s"}
            self.tableDisplayColumns.append(dispDict)
        
        # Now tracking when changes are made
        if 'fields' in self.configDict.keys():
            dispDict={'name': "lastUpdated", 'label': "lastUpdated", 'fmt': "%s"}
            dispDict={'name': "user", 'label': "user", 'fmt': "%s"}
            self.tableDisplayColumns.append(dispDict)
                                       
        # Full list of image directories - for adding image_ tags in buildCacheForObject
        # NOTE: handling of surveys (e.g., SDSS) is clunfky and getting unwieldy...
        # NOTE: made clunkier to allow different size images (imageDirs only for now)
        self.imDirLabelsList=[]
        self.imDirMaxSizeArcminList=[]
        if self.configDict['addDECaLSImage'] == True:
            self.imDirLabelsList.append("DECaLS")
            self.imDirMaxSizeArcminList.append(self.configDict['plotSizeArcmin'])
        if self.configDict['addSDSSImage'] == True:
            self.imDirLabelsList.append("SDSS")
            self.imDirMaxSizeArcminList.append(self.configDict['plotSizeArcmin'])
        if self.configDict['addPS1Image'] == True:
            self.imDirLabelsList.append("PS1")
            self.imDirMaxSizeArcminList.append(self.configDict['plotSizeArcmin'])
        if self.configDict['addPS1IRImage'] == True:
            self.imDirLabelsList.append("PS1IR")
            self.imDirMaxSizeArcminList.append(self.configDict['plotSizeArcmin'])
        if self.configDict['addUnWISEImage'] == True:
            self.imDirLabelsList.append("unWISE")
            self.imDirMaxSizeArcminList.append(self.configDict['plotSizeArcmin'])
        if self.configDict['addHSCImage'] == True:
            self.imDirLabelsList.append("HSC")
            self.imDirMaxSizeArcminList.append(self.configDict['plotSizeArcmin'])
        if 'tileDirs' in self.configDict.keys():
            for tileDirDict in self.configDict['tileDirs']:
                self.imDirLabelsList.append(tileDirDict['label'])
                self.imDirMaxSizeArcminList.append(self.configDict['plotSizeArcmin'])
        if 'imageDirs' in self.configDict.keys():
            for imDirDict in self.configDict['imageDirs']:
                self.imDirLabelsList.append(imDirDict['label'])
                if 'maxSizeArcmin' in imDirDict.keys():
                    self.imDirMaxSizeArcminList.append(imDirDict['maxSizeArcmin'])
                else:
                    self.imDirMaxSizeArcminList.append(self.configDict['plotSizeArcmin'])
                    
        # Pre-processing - building image cache etc.
        self.http=urllib3.PoolManager()
        if preprocess == True:
            self.preprocess()

        # This sets size of table view - view is controlled with session variables
        self.tableViewRows=40

        self.jinja_env = Environment(
            loader=FileSystemLoader(_TEMPLATES_DIR),
            autoescape=False
        )


    def matchTags(self, obj):
        """Find match in MongoDB to obj row from tab. If we don't find one, return a dictionary with blank
        values where fields would be. We now allow "overloading" - i.e., user can specify as editable fields
        columns which already exist in the database. In that case, we set as default values (if key not 
        found) the values in the input catalog.
        
        """
        
        if obj['RADeg'] > 180:
            lon=360.0-obj['RADeg']
        else:
            lon=obj['RADeg']
        matches=self.tagsCollection.find({'loc': SON({'$nearSphere': [lon, obj['decDeg']], '$maxDistance': np.radians(self.configDict['MongoDBCrossMatchRadiusArcmin']/60.0)})}).limit(1)
        matches=list(matches)

        if matches == []:
            newPost={'loc': {'type': 'Point', 'coordinates': [lon, obj['decDeg']]}}
            self.tagsCollection.insert_one(newPost)
            mongoDict={}
        else:
            mongoDict=matches[0]
        
        # Check we don't have a blank entry in terms of fields we expect
        if 'classifications' in self.configDict.keys() and 'classification' not in mongoDict.keys():
            mongoDict['classification']=""
        if 'fields' in self.configDict.keys():
            for fieldDict in self.configDict['fields']:
                if fieldDict['name'] not in mongoDict.keys():
                    # Allow overloading: first check if in obj keys
                    if fieldDict['name'] in obj.keys():
                        defValue=obj[fieldDict['name']]
                    elif fieldDict['type'] == 'number':
                        defValue=0.0
                    elif fieldDict['type'] == 'text':
                        defValue=""
                    mongoDict[fieldDict['name']]=defValue
        
        # Strip out _id and loc, because whatever is calling this routine won't want them
        if '_id' in mongoDict.keys():
            del mongoDict['_id']
        if 'loc' in mongoDict.keys():
            del mongoDict['loc']
    
        return mongoDict
    
    
    def addSourceryIDs(self, tab):
        """Adds a sourceryID column to the given astropy table object.
        
        Returns astropy table object
        
        """
        
        sourceryIDs=[]
        if 'sourceList' in tab.keys():
            # Takes < 1 sec for 36,000 sources
            for row in tab:
                # sourceryIDs.append(row['sourceList']+"_"+row['name'].replace(" ", "_"))
                sourceryIDs.append(row['sourceList']+"_%.6f_%.6f" % (row['RADeg'], row['decDeg']))
        else:
            for row in tab:
                sourceryIDs.append(row['name'].replace(" ", "_"))
        tab.add_column(atpy.Column(sourceryIDs, "sourceryID"))
        
        return tab
    
    
    def addFootprintColumns(self, tab):
        """Adds footprint_ columns to the table, for finding overlap with e.g. DES Y3.
        
        """
        
        for footprintDict in self.configDict['footprints']:
            colLabel='footprint_%s' % (footprintDict['label'])
            print("... adding %s ..." % (colLabel))
            tab.add_column(atpy.Column(np.zeros(len(tab), dtype = bool), colLabel))
            for maskPath in footprintDict['maskList']:
                with pyfits.open(maskPath) as img:
                    for extName in range(0, 2):
                        mask=img[extName].data
                        if mask is not None:
                            break
                    wcs=astWCS.WCS(img[extName].header, mode = 'pyfits')
                xy=np.array(wcs.wcs2pix(tab['RADeg'].data, tab['decDeg'].data))
                xs=np.array(xy[:, 0], dtype = int)
                ys=np.array(xy[:, 1], dtype = int)
                for i in range(len(tab)):
                    x=xs[i]
                    y=ys[i]
                    if x >= 0 and x < mask.shape[1] and y >= 0 and y < mask.shape[0]:
                        if mask[y, x] > 0:
                            tab[colLabel][i]=1
            del mask, wcs
        
        return tab
        

    def buildDatabase(self):
        """Import .fits table into MongoDB database as sourceCollection. Delete any pre-existing catalog
        there. Do all the cross matching at this stage also.
        
        We also cross match the tagsCollection onto sourceCollection, to save doing a full cross match again
        later. When we need to update, we'll update both tagsCollection and sourceCollection.
        
        We also store a list of fields and types in a collection, so that we can (a) use this for the 
        constraints help page; (b) keep fields in a sensible order (assuming input catalogs are in sensible 
        order)
        
        """
        
        print(">>> Building database ...")
        t0=time.time()
        with open(self.dbLockFileName, "w") as dbLockFile:
            pass

        # Delete any pre-existing entries
        self.db.drop_collection('sourceCollection')
        self.db.drop_collection('fieldTypes')
        
        # Table set up
        tab=atpy.Table().read(self.configDict['catalogFileName'])
        origLen=len(tab)
        tab.sort(["RADeg", "decDeg"])
        
        # Optionally zap all but first 10 rows (for testing)
        if 'quickTest' in self.configDict.keys() and self.configDict['quickTest'] == True:
            tab=tab[:10]
            print("... WARNING: quickTest mode enabled ...")
            
        # In case we don't have a 'name' column, we relabel a given column
        if 'nameColumn' in self.configDict.keys() and self.configDict['nameColumn'] != "":
            if 'name' not in tab.keys():
                tab.rename_column(self.configDict['nameColumn'], 'name')
        
        # NOTE: sourceList is now a special column: if present, we use that to make a hidden sourceryID column
        # We need this to ensure that on displaySourcePage, we show the right properties table for the selected object
        # However, we don't want to put this info in the tags table... as that needs to be based on positional matching
        tab=self.addSourceryIDs(tab)

        # Add footprints, if any
        if 'footprints' in self.configDict.keys():
            tab=self.addFootprintColumns(tab)
        
        # NOTE: another special column - this is for tracking whether the cache files (images, redshifts) have been
        # fetched or not, for a given object. We set this to 0 each time we rebuild the database, and set to 1 each
        # in preprocess after we process each object. This allows the cache building process to re-start from where
        # it left off without checking every single object again
        tab.add_column(atpy.Column(np.zeros(len(tab), dtype = bool), "cacheBuilt"))
        
        # Cross matching based on sourceryID
        # We do this first before position matching, because the join operation jumbles row order
        if 'crossMatchCatalogs' in self.configDict.keys():
            for xMatchDict in self.configDict['crossMatchCatalogs']:
                f=xMatchDict['fileName']
                label=xMatchDict['label']
                radiusArcmin=xMatchDict['crossMatchRadiusArcmin']
                xTab=atpy.Table().read(f)
                if 'nameCol' in xMatchDict.keys():
                    xTab.rename_column(xMatchDict['nameCol'], 'name')
                RAKeys=['ra', 'RA', 'Ra']
                decKeys=['dec', 'DEC', 'Dec']
                for r in RAKeys:
                    if r in xTab.keys():
                        xTab.rename_column(r, 'RADeg')
                for d in decKeys:
                    if d in xTab.keys():
                        xTab.rename_column(d, 'decDeg')
                # Cross match based on sourceryID
                if 'sourceList' in xTab.keys():
                    print("... matching %s based on sourceryID ..." % (label))
                    excludeKeys=['name', 'RADeg', 'decDeg', 'sourceList']
                    for key in xTab.keys():
                        if key not in excludeKeys:
                            xTab.rename_column(key, '%s_%s' % (label, key))
                    xTab=self.addSourceryIDs(xTab)
                    xTab.remove_columns(excludeKeys)
                    tab=atpy.join(tab, xTab, keys = ['sourceryID'], join_type = 'left')
                    for key in xTab.keys():
                        try:
                            tab[key][tab[key].mask]=-99
                        except:
                            continue
                    tab=atpy.Table(tab, masked = False)
                    tab.sort(["RADeg", "decDeg"]) # For some reason join jumbles row order
                    if len(tab) != origLen:
                        raise Exception("Length of table has changed - this probably means multiple objects with the same name but with the same sourceList. Fix table %s." % (f))
        
        # Cross matching based on positions
        if 'crossMatchCatalogs' in self.configDict.keys():
            tab.add_column(atpy.Column(np.arange(len(tab)), 'matchIndices'))
            origLen=len(tab)
            cat1=SkyCoord(ra = tab['RADeg'], dec = tab['decDeg'], unit = 'deg')
            for xMatchDict in self.configDict['crossMatchCatalogs']:
                f=xMatchDict['fileName']
                label=xMatchDict['label']
                radiusArcmin=xMatchDict['crossMatchRadiusArcmin']
                xTab=atpy.Table().read(f)
                # Repetition...
                if 'nameCol' in xMatchDict.keys():
                    xTab.rename_column(xMatchDict['nameCol'], 'name')
                RAKeys=['ra', 'RA', 'Ra']
                decKeys=['dec', 'DEC', 'Dec']
                for r in RAKeys:
                    if r in xTab.keys():
                        xTab.rename_column(r, 'RADeg')
                for d in decKeys:
                    if d in xTab.keys():
                        xTab.rename_column(d, 'decDeg')
                if 'sourceList' not in xTab.keys():
                    print("... matching %s based on position ..." % (label))
                    xMatchRadiusDeg=radiusArcmin/60.
                    cat2=SkyCoord(ra = xTab['RADeg'].data, dec = xTab['decDeg'].data, unit = 'deg')
                    xIndices, rDeg, sep3d = match_coordinates_sky(cat1, cat2, nthneighbor = 1)
                    mask=np.less(rDeg.value, xMatchRadiusDeg)
                    for key in xTab.keys():
                        xTab.rename_column(key, '%s_%s' % (label, key))
                    tab['matchIndices'][:]=-1
                    tab['matchIndices']=xIndices
                    # Could not get join to work
                    for key in xTab.keys():
                        if key not in tab.keys():
                            if xTab[key].dtype.kind == 'S':
                                tab.add_column(atpy.Column(np.array([""]*len(tab), dtype = xTab[key].dtype), key))
                            else:
                                dtype=xTab[key].dtype
                                if dtype == np.uint8:
                                    dtype=int
                                tab.add_column(atpy.Column(np.ones(len(tab), dtype = dtype)*-99, key))
                            tab[key][mask]=xTab[key][xIndices[mask]]
                    tab.add_column(atpy.Column(np.zeros(len(tab), dtype = int), '%s_match' % (label)))
                    tab.add_column(atpy.Column(np.ones(len(tab), dtype = float)*-99, '%s_distArcmin' % (label)))
                    tab['%s_match' % (label)][mask]=1
                    tab['%s_distArcmin' % (label)][mask]=rDeg.value[mask]*60.0
                del xTab
            tab.remove_column("matchIndices")
        assert(len(tab) == origLen)
        
        # Cache the result of the cross matches: we need this for speed later on when downloading catalogs
        # Otherwise, for large catalogs, we're hitting memory issues
        cachedTabFileName=self.cacheDir+os.path.sep+"%s_xMatchedTable.fits" % (self.configDict['catalogDownloadFileName'])
        #cachedTabFileName=self.cacheDir+os.path.sep+"%s_xMatchedTable.csv" % (self.configDict['catalogDownloadFileName'])
        if len(tab.columns) > 999:
            raise Exception("FITS format is limited to a maximum 1000 columns - so prune your cross-match tables.")
        tab.write(cachedTabFileName, overwrite = True)
        print("... written %s ..." % (cachedTabFileName))
        
        # Import each object into MongoDB - now doing this in bulk (slightly quicker)
        idCount=0
        fieldTypesList=[]   # Used for making sensible column order later
        fieldTypesDict={}   # Used for tracking types for help page
        postsList=[]
        for row in tab:
            #t0=time.time()
            # Need an id number for table display
            idCount=idCount+1
            newPost={'index': idCount}
            newPost['name']=row['name']
            newPost['RADeg']=row['RADeg']
            newPost['decDeg']=row['decDeg']
            # MongoDB coords for spherical geometry
            if row['RADeg'] > 180:
                lon=360.0-row['RADeg']
            else:
                lon=row['RADeg']
            newPost['loc']={'type': 'Point', 'coordinates': [lon, row['decDeg']]}
            
            # Properties in the table
            for key in tab.keys():
                # Just to make sure MongoDB happy with data types
                # e.g., redmapper .fits table doesn't play nicely by default
                if tab.columns[key].dtype.name.find("int") != -1:
                    # Changed for numpy 2+
                    val=row[key]
                    if type(val) == np.ma.core.MaskedConstant:
                        newPost[key]=-99
                    else:
                        newPost[key]=int(row[key])
                    if key not in fieldTypesList:
                        fieldTypesList.append(key)
                        fieldTypesDict[key]="number"
                elif tab.columns[key].dtype.name.find("str") != -1 or tab.columns[key].dtype.name.find("bytes") != -1:
                    newPost[key]=str(row[key])
                    if key not in fieldTypesList:
                        fieldTypesList.append(key)
                        fieldTypesDict[key]="text"
                elif tab.columns[key].dtype.name.find("float") != -1:
                    if key not in fieldTypesList:
                        fieldTypesList.append(key)
                        fieldTypesDict[key]="number"
                    newPost[key]=float(row[key])
                elif tab.columns[key].dtype.name.find("bool") != -1:
                    if key not in fieldTypesList:
                        fieldTypesList.append(key)
                        fieldTypesDict[key]="number"
                    newPost[key]=int(row[key]) # Was bool, but that causes some issues with queries that we'd need to fix properly
                else:
                    raise Exception("Unknown data type in column '%s'" % (key))
                        
            # NED cross match
            if 'addNEDMatches' in self.configDict.keys() and self.configDict['addNEDMatches'] == True:
                self.findNEDMatch(newPost, NEDObjTypes = self.configDict['NEDObjTypes'])
                stringKeys=['NED_name']
                numberKeys=['NED_z', 'NED_distArcmin', 'NED_RADeg', 'NED_decDeg']
                typesList=['text', 'number']
                for t, l in zip(typesList, [stringKeys, numberKeys]):
                    for key in l:
                        if key not in fieldTypesList:
                            fieldTypesList.append(key)
                            fieldTypesDict[key]=t
            
            # Match with tagsCollection
            tagsDict=self.matchTags(newPost)
            # if row['name'] == 'ACT-CL J1015.7+1813':
            #     print("zap")
            #     import IPython
            #     IPython.embed()
            #     sys.exit()
            for key in tagsDict:
                newPost[key]=tagsDict[key]
            if self.configDict['insertMode'] == 'single':
                print("... adding %s to database (%d/%d) ..." % (row['name'], idCount, len(tab)))
                self.sourceCollection.insert_one(newPost)
            else:
                postsList.append(newPost)
        
        # Insert all posts at once
        if self.configDict['insertMode'] == 'many':
            print("... inserting all posts into database ...")
            self.sourceCollection.insert_many(postsList)

        # Make collection of field types
        index=0
        for key in fieldTypesList:
            fieldDict={}
            fieldDict['name']=key
            fieldDict['type']=fieldTypesDict[key]
            fieldDict['index']=index
            if key in self.descriptionsDict.keys():
                fieldDict['description']=self.descriptionsDict[key]
            else:
                fieldDict['description']="-"
            self.fieldTypesCollection.insert_one(fieldDict)
            index=index+1

        t1=time.time()
        if os.path.exists(self.dbLockFileName) ==True:
            os.remove(self.dbLockFileName)

        print("... building database complete: took %.1f sec ..." % (t1-t0))
            

    def parseColumnDescriptionsFile(self):
        """Reads a text file containing column descriptions. The format of the file is:
        
        columnName: description
        
        Any white space between : and description is stripped.
        
        Returns a dictionary with keys {'columnName': 'description'}
        
        """
        
        descriptionsDict={}
        if 'descriptionsFileName' in self.configDict.keys():
            inFile=open(self.configDict['descriptionsFileName'], "r")
            lines=inFile.readlines()
            inFile.close()
            for line in lines:
                if line[0] != "#":
                    splitIndex=line.find(":")
                    key=line[:splitIndex]
                    desc=line[splitIndex+1:].lstrip().rstrip()
                    descriptionsDict[key]=desc

        return descriptionsDict
    
                            
    def parseConfig(self, configFileName):
        """Parse config file, unpacking parameters into the SourceBrowser object.
        
        NOTE: config file format is now yaml
        """
        with open(configFileName, "r") as stream:
            self.configDict=yaml.safe_load(stream)

        # Add root path where necessary in place
        if 'sourceryPath' in self.configDict.keys() and self.configDict['sourceryPath'] != "":
            rootDir=self.configDict['sourceryPath'].rstrip(os.path.sep)
            keysToFix=["userListFile", "cacheDir", "skyviewCacheDir", "newsFileName", "logoFileName", "crossMatchCatalogs", "tileDirs"]
            for k in keysToFix:
                if k in self.configDict.keys():
                    if type(self.configDict[k]) == list:
                        for i in range(len(self.configDict[k])):
                            for pathKey in ['fileName', 'path']:
                                if pathKey in self.configDict[k][i].keys():
                                    self.configDict[k][i][pathKey]=rootDir+os.path.sep+self.configDict[k][i][pathKey]                                
                    else:
                        self.configDict[k]=rootDir+os.path.sep+self.configDict[k]
        
        if 'defaultViewSizeArcmin' not in self.configDict.keys():
            self.configDict['defaultViewSizeArcmin']=self.configDict['plotSizeArcmin']
            
        if "NEDObjTypes" not in self.configDict.keys():
            self.configDict['NEDObjTypes']=['GClstr']

        if 'insertMode' not in self.configDict.keys():
            self.configDict['insertMode']='single'
        
        # We now support keeping a huge .fits table of spec-zs
        # We start with SDSS but this could (should?) be made generic
        if 'specRedshiftsTable' not in self.configDict.keys():
            self.configDict['specRedshiftsTable']=None
        #else:
        #    if 'sourceryPath' in self.configDict.keys() and self.configDict['sourceryPath'] != "":
        #        rootDir=self.configDict['sourceryPath'].rstrip(os.path.sep)
        #        self.configDict['specRedshiftsTable']=rootDir+os.path.sep+self.configDict['specRedshiftsTable']


    def addNews(self):
        """Parse news file, if there is one, filling up self.configDict['newsItems'].
        
        """
        
        if os.path.exists(self.configDict['newsFileName']) == True:
            inFile=open(self.configDict['newsFileName'], "r")
            lines=inFile.readlines()
            inFile.close()
        else:
            lines=[]
            
        newsItems=[]
        for line in lines:
            if line[0] != "#" and len(line) > 3:
                newsItems.append(line.replace("\n", ""))
        self.configDict['newsItems']=newsItems
        
        
    def fetchNEDInfo(self, name, RADeg, decDeg, retryFails = False):
        """Fetches NED info for given obj (which must have name, RADeg, decDeg keys) - just stores on disk 
        in cacheDir - we'll retrieve it later as needed.
        
        """
        halfMatchBoxLengthDeg=5.0/60.0
        RAMin=RADeg-halfMatchBoxLengthDeg
        RAMax=RADeg+halfMatchBoxLengthDeg
        decMin=decDeg-halfMatchBoxLengthDeg
        decMax=decDeg+halfMatchBoxLengthDeg
        outFileName=self.nedDir+os.path.sep+name.replace(" ", "_")+".txt"        
        if os.path.exists(outFileName) == False:
            # Open the file straightaway so we don't clash if running threaded
            with open(outFileName, 'wb') as f:
                print("... fetching NED info for %s ..." % (name))
                urlString="http://ned.ipac.caltech.edu/cgi-bin/objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon=%.6fd&lat=%.6fd&radius=%.2f&dot_include=ANY&in_objtypes1=GGroups&in_objtypes1=GClusters&in_objtypes1=QSO&in_objtypes2=Radio&in_objtypes2=SmmS&in_objtypes2=Infrared&in_objtypes2=Xray&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=ascii_tab&zv_breaker=30000.0&list_limit=5&img_stamp=YES" % (RADeg, decDeg, halfMatchBoxLengthDeg*60.0)
                resp=self.http.request('GET', urlString)
                f.write(resp.data)
                f.close()


    def findNEDMatch(self, obj, NEDObjTypes = ["GClstr"]):
        """Checks if there is a NED match for obj. Uses matching radius specified in config file by
        crossMatchRadiusArcmin. Only returns match if it matches the first entry in NEDObjTypes list.
        
        """
                    
        nedFileName=self.nedDir+os.path.sep+obj['name'].replace(" ", "_")+".txt"
        nedObjs=catalogTools.parseNEDResult(nedFileName, onlyObjTypes = NEDObjTypes)
            
        # Flag matches against clusters - choose nearest one
        rMin=10000
        crossMatchRadiusDeg=self.configDict['NEDCrossMatchRadiusArcmin']/60.0
        clusterMatch={}
        if len(nedObjs['RAs']) > 0:
            for i in range(len(nedObjs['RAs'])):
                ned=nedObjs
                if ned['sourceTypes'][i] == NEDObjTypes[0]:
                    r=astCoords.calcAngSepDeg(ned['RAs'][i], ned['decs'][i], obj['RADeg'], obj['decDeg'])
                    if r < rMin and r < crossMatchRadiusDeg:
                        keepName=False
                        if 'name' in clusterMatch:
                            if "ABELL" in clusterMatch['name']:
                                keepName=True
                        if keepName == False:
                            rMin=r
                            clusterMatch['name']=ned['names'][i]
                            if ned['redshifts'][i] != 'N/A':
                                clusterMatch['z']=float(ned['redshifts'][i])
                            else:
                                clusterMatch['z']=None
                            clusterMatch['rArcmin']=rMin*60.0
                            clusterMatch['NED_RADeg']=float(ned['RAs'][i])
                            clusterMatch['NED_decDeg']=float(ned['decs'][i])
        if clusterMatch != {}:
            obj['NED_name']=clusterMatch['name']  
            obj['NED_z']=clusterMatch['z']
            obj['NED_distArcmin']=clusterMatch['rArcmin']
            obj['NED_RADeg']=clusterMatch['NED_RADeg']
            obj['NED_decDeg']=clusterMatch['NED_decDeg']
        else:
            obj['NED_name']=None 
            obj['NED_z']=np.nan
            obj['NED_distArcmin']=np.nan
            obj['NED_RADeg']=np.nan
            obj['NED_decDeg']=np.nan        
        

    def fetchPS1Image(self, name, RADeg, decDeg, bands = "gri", refetch = False):
        """Fetches Pan-STARRS gri or izy .jpg using the cutout webservice.
        
        """
       
        # 1920 pixels is 8 arcmin on PS1 scale
        PS1PlotSizePix=int(round(self.configDict['plotSizeArcmin']*240))
            
        if bands == 'gri':
            cacheDirLabel="PS1"
            urlString="http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=%.6f+%.6f&filter=color&filter=g&filter=r&filter=i&filetypes=stack&auxiliary=data&size=%d&output_size=1024&verbose=0&autoscale=99.500000&catlist=" % (RADeg, decDeg, PS1PlotSizePix)
        elif bands == 'izy':
            cacheDirLabel="PS1IR"
            urlString="http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=%.6f+%.6f&filter=color&filter=i&filter=z&filter=y&filetypes=stack&auxiliary=data&size=%d&output_size=1024&verbose=0&autoscale=99.500000&catlist=" % (RADeg, decDeg, PS1PlotSizePix)
        else:
            raise Exception("PS1 bands should be 'gri' or 'izy'")

        ps1CacheDir=self.cacheDir+os.path.sep+cacheDirLabel
        
        if decDeg < -30:
            print("... outside PS1 area - skipping ...")
            return None
        
        subDir=str(RADeg).split(".")[0]
        os.makedirs(ps1CacheDir+os.path.sep+subDir, exist_ok = True)
        outFileName=ps1CacheDir+os.path.sep+subDir+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".jpg"
        tmpFile, tmpFileName=tempfile.mkstemp()
        
        if os.path.exists(outFileName) == False or refetch == True:
        
            if os.path.exists(tmpFileName) == True:
                os.remove(tmpFileName)

            resp=self.http.request('GET', urlString)
            with open(tmpFileName, 'wb') as f:
                f.write(resp.data)
                f.close()

            inFile=open(tmpFileName, 'r')
            lines=inFile.readlines()
            inFile.close()
            foundLine=False
            for line in lines:
                if line.find("fitscut.cgi") != -1 and line.find("green") != -1:
                    foundLine=True
                    break
            if foundLine == True:
                urlString='http://'+line.split('src="//')[-1].split('"')[0]
                resp=self.http.request('GET', urlString)
                with open(outFileName, 'wb') as f:
                    f.write(resp.data)
                    f.close()
            else:
                print("... WARNING: couldn't retrieve PS1 image ...")
        
        if os.path.exists(tmpFileName) == True:
            os.close(tmpFile)
            os.remove(tmpFileName)
                

    def fetchSDSSImage(self, name, RADeg, decDeg, refetch = False):
        """Fetches the SDSS .jpg for the given image size using the casjobs webservice.
        
        makeSDSSPlots loads these jpegs in, and use matplotlib to make them into plots with
        coord axes etc.
        
        The way we're handling image directories is a bit clunky at the moment...
        
        """
    
        sdssCacheDir=self.cacheDir+os.path.sep+"SDSS"

        subDir=str(RADeg).split(".")[0]
        os.makedirs(sdssCacheDir+os.path.sep+subDir, exist_ok = True)

        outFileName=sdssCacheDir+os.path.sep+subDir+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".jpg"
        SDSSWidth=int(round((1200.0/8.0)*self.configDict['plotSizeArcmin']))
        SDSSScale=(self.configDict['plotSizeArcmin']*60.0)/SDSSWidth # 0.396127
        if os.path.exists(outFileName) == False or refetch == True:
            #urlString="http://skyservice.pha.jhu.edu/DR10/ImgCutout/getjpeg.aspx?ra="+str(RADeg)+"&dec="+str(decDeg)
            urlString="http://skyserver.sdss.org/dr13/SkyServerWS/ImgCutout/getjpeg?TaskName=Skyserver.Chart.Image&ra="+str(RADeg)+"&dec="+str(decDeg)
            urlString=urlString+"&scale="+str(SDSSScale)+"&width="+str(int(SDSSWidth))+"&height="+str(int(SDSSWidth))
            resp=self.http.request('GET', urlString)
            with open(outFileName, 'wb') as f:
                f.write(resp.data)
                f.close()
            # Old
            #try:
                #urllib.urlretrieve(urlString, filename = outFileName)
            #except:
                #print "... WARNING: couldn't get SDSS image ..."
                #print urlString
                #outFileName=None


    def fetchLegacySurveyImage(self, name, RADeg, decDeg, sizePix = 800, refetch = False, layer = 'ls-dr10-early-grz',\
                               bands = None, cacheSubDir = None):
        """Fetches .jpg cut-out from legacysurvey.org sky viewer. Based on the code in sourcery.

        Valid layers include ls-dr9, ls-dr10-early-grz etc.

        """

        if cacheSubDir is None:
            decalsCacheDir=self.cacheDir+os.path.sep+"legacy_%s" % (layer)
        else:
            decalsCacheDir=self.cacheDir+os.path.sep+cacheSubDir

        os.makedirs(decalsCacheDir, exist_ok = True)

        subDir=str(RADeg).split(".")[0]
        os.makedirs(decalsCacheDir+os.path.sep+subDir, exist_ok = True)

        outFileName=decalsCacheDir+os.path.sep+subDir+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".jpg"

        decalsWidth=sizePix
        decalsPixScale=(self.configDict['plotSizeArcmin']*60.0)/float(decalsWidth)
        if os.path.exists(outFileName) == False or refetch == True:
            #http://legacysurvey.org/viewer/jpeg-cutout?ra=52.102810&dec=-21.670020&size=2048&layer=des-dr1&pixscale=0.3809&bands=grz
            urlString="http://legacysurvey.org/viewer/jpeg-cutout?ra=%.6f&dec=%.6f&size=%d&layer=%s&pixscale=%.4f" % (RADeg, decDeg, decalsWidth, layer, decalsPixScale)
            if bands is not None:
                urlString=urlString+"&bands=%s" % (bands)
            print(urlString)
            resp=self.http.request('GET', urlString)
            if resp.data.find(b'Server Error') != -1:
                return None
            with open(outFileName, 'wb') as f:
                f.write(resp.data)
                f.close()


    def fetchDECaLSImage(self, name, RADeg, decDeg, refetch = False):
        self.fetchLegacySurveyImage(name, RADeg, decDeg, refetch = refetch, layer = 'ls-dr10-grz',
                                    cacheSubDir = "DECaLS")


    def fetchUnWISEImage(self, name, RADeg, decDeg, refetch = False):
        self.fetchLegacySurveyImage(name, RADeg, decDeg, refetch = refetch, layer = 'unwise-neo6',
                                    cacheSubDir = "unWISE")


    def fetchHSCImage(self, name, RADeg, decDeg, refetch = False):
        self.fetchLegacySurveyImage(name, RADeg, decDeg, refetch = refetch, layer = 'hsc2',
                                    cacheSubDir = "HSC")
                
    @cherrypy.expose
    def makePlotFromJPEG(self, name, RADeg, decDeg, surveyLabel, plotNEDObjects = "false", plotSpecObjects = "false",\
                         plotSourcePos = "false", plotXMatch = "false", plotContours = "false", showAxes = "false",\
                         clipSizeArcmin = None, gamma = 1.0, redshift = "none", plotRedshift = "false"):
        """Makes plot of .jpg image with coordinate axes and NED, SDSS objects overlaid.
        
        To test this:
        
        http://localhost:8080/makeSDSSPlot?name=XMMXCS%20J001737.5-005234.2&RADeg=4.406325&decDeg=-0.876192
        
        To test zoom:
        
        http://localhost:8080/sourcery/makePlotFromJPEG?name=XMMXCS%20J001737.5-005234.2&RADeg=4.406325&decDeg=-0.876192&surveyLabel=SDSS&clipSizeArcmin=3.0
        
        """
        
        # Just in case they are passed as strings (e.g., if direct from the url)
        RADeg=float(RADeg)
        decDeg=float(decDeg)
        
        # This is only used for scaling the size of plotted points
        if clipSizeArcmin == None:
            sizeDeg=self.configDict['plotSizeArcmin']/60.0
        else:
            sizeDeg=float(clipSizeArcmin)/60.
        
        # Load data
        # inJPGPath=self.cacheDir+os.path.sep+surveyLabel+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".jpg"
        subDir=str(RADeg).split(".")[0]
        inJPGPath=self.cacheDir+os.path.sep+surveyLabel+os.path.sep+subDir+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".jpg"
        logger.info("Image path: %s" % (inJPGPath))
        if os.path.exists(inJPGPath) == False:
            # Live fetching of not yet cached images from web services, where we can
            if surveyLabel == 'unWISE':
                self.fetchUnWISEImage(name, RADeg, decDeg, refetch = False)
            elif surveyLabel == 'DECaLS':
                self.fetchDECaLSImage(name, RADeg, decDeg, refetch = False)
            elif surveyLabel == 'HSC':
                self.fetchHSCImage(name, RADeg, decDeg, refetch = False)
            elif surveyLabel == 'SDSS':
                self.fetchSDSSImage(name, RADeg, decDeg, refetch = False)
            elif surveyLabel == 'PS1':
                self.fetchPS1Image(name, RADeg, decDeg, bands = "gri", refetch = False)
            elif surveyLabel == 'PS1IR':
                self.fetchPS1Image(name, RADeg, decDeg, bands = "izy", refetch = False)
            if os.path.exists(inJPGPath) == False:
                inJPGPath=sourcery.__path__[0]+os.path.sep+"data"+os.path.sep+"noData.jpg"
        
        # Gets set to None if we got it from legacysurvey.org and we didn't write to disk
        if inJPGPath is not None:
            im=Image.open(inJPGPath)
            data=np.array(im)
            data=np.power(data, 1.0/float(gamma))
            try:
                data=np.flipud(data)
                #data=np.fliplr(data)
            except:
                #"... something odd about image (1d?) - aborting ..."
                return None
                    
        R=data[:, :, 0]
        G=data[:, :, 1]
        B=data[:, :, 2]
                
        # HACK: for ACT maps, with huge pixels, we can get offsets in .jpg relative to original
        # So, if we have a .fits image, load that and use to set centre coords
        #fitsFileName=inJPGPath.replace(".jpg", ".fits")
        #if os.path.exists(fitsFileName) == True:
            #hackWCS=astWCS.WCS(fitsFileName)
            #CRVAL1, CRVAL2=hackWCS.getCentreWCSCoords()
        #else:
            #CRVAL1, CRVAL2=RADeg, decDeg
        # Make a WCS
        CRVAL1, CRVAL2=RADeg, decDeg
        for imageType, maxSizeArcmin in zip(self.imDirLabelsList, self.imDirMaxSizeArcminList):
            if imageType == surveyLabel:
                sizeArcmin=maxSizeArcmin
                break
        xSizeDeg, ySizeDeg=sizeArcmin/60.0, sizeArcmin/60.0
        xSizePix=float(R.shape[1])
        ySizePix=float(R.shape[0])
        xRefPix=xSizePix/2.0
        yRefPix=ySizePix/2.0
        xOutPixScale=xSizeDeg/xSizePix
        yOutPixScale=ySizeDeg/ySizePix
        newHead=pyfits.Header()
        newHead['NAXIS']=2
        newHead['NAXIS1']=xSizePix
        newHead['NAXIS2']=ySizePix
        newHead['CTYPE1']='RA---TAN'
        newHead['CTYPE2']='DEC--TAN'
        newHead['CRVAL1']=CRVAL1
        newHead['CRVAL2']=CRVAL2
        newHead['CRPIX1']=xRefPix+1
        newHead['CRPIX2']=yRefPix+1
        newHead['CDELT1']=-xOutPixScale
        newHead['CDELT2']=xOutPixScale    # Makes more sense to use same pix scale
        newHead['CUNIT1']='DEG'
        newHead['CUNIT2']='DEG'
        wcs=astWCS.WCS(newHead, mode='pyfits')

        #cutLevels=[[R.min(), R.max()], [G.min(), G.max()], [B.min(), B.max()]]
        cutLevels=[[0, 255], [0, 255], [0, 255]]
        
        # Optional zoom
        if clipSizeArcmin != None:
            clipSizeArcmin=float(clipSizeArcmin)
            RClip=astImages.clipImageSectionWCS(R, wcs, RADeg, decDeg, clipSizeArcmin/60.0)
            GClip=astImages.clipImageSectionWCS(G, wcs, RADeg, decDeg, clipSizeArcmin/60.0)
            BClip=astImages.clipImageSectionWCS(B, wcs, RADeg, decDeg, clipSizeArcmin/60.0)
            R=RClip['data']
            G=GClip['data']
            B=BClip['data']
            wcs=RClip['wcs']
        #astImages.saveFITS("test.fits", R, wcs)
        
        # Make plot
        if showAxes == "true":
            axes=[0.1,0.085,0.9,0.85]
            axesLabels="sexagesimal"
            figSize=self.configDict['figSize']
        else:
            axes=[0, 0, 1, 1]
            axesLabels="sexagesimal"    # Avoid dealing with axis flips
            figSize=(max(self.configDict['figSize']), max(self.configDict['figSize']))
        fig=plt.figure(figsize = figSize)
        
        p=astPlots.ImagePlot([R, G, B], wcs, cutLevels = cutLevels, title = name.replace("_", " "), axes = axes, 
                            axesLabels = axesLabels)
        
        if showAxes != "true":
            scaleBarSizeArcmin=1.0
            p.addScaleBar('NW', scaleBarSizeArcmin*60.0, color='yellow', fontSize=20, width=2.0, label = "1'")
            plt.figtext(0.025, 0.95, name.replace("_", " "), ha = 'left', size = 24, color = 'yellow')
            #if plotTitle != None:
            #plt.figtext(0.965, 0.88, plotTitle, ha = 'right', size = 24)
        
        if plotSourcePos == "true":
            p.addPlotObjects([RADeg], [decDeg], 'clusterPos', symbol='cross', size=sizeDeg/20.0*3600.0, color='white')
                
        if plotNEDObjects == "true":
            # We should already have the files for this from doing addNEDInfo earlier
            nedFileName=self.nedDir+os.path.sep+name.replace(" ", "_")+".txt"
            nedObjs=catalogTools.parseNEDResult(nedFileName, onlyObjTypes = self.configDict['NEDObjTypes'])
            if len(nedObjs['RAs']) > 0:
                p.addPlotObjects(nedObjs['RAs'], nedObjs['decs'], 'nedObjects', objLabels = nedObjs['labels'],
                                    size = sizeDeg/40.0*3600.0, color = "#7cfc00")
        
        if plotRedshift == "true" and redshift != "none":
            plt.figtext(0.025, 0.03, "z = %.2f" % (float(redshift)), ha = 'left', size = 24, color = 'white')
    
        if plotSpecObjects == "true":
            specRedshifts=catalogTools.fetchSpecRedshifts(name, RADeg, decDeg, 
                                                          redshiftsTable = self.specRedshiftsTab)
            if specRedshifts is not None:
                specRAs=[]
                specDecs=[]
                specLabels=[]
                specCount=0
                for specObj in specRedshifts:
                    specCount=specCount+1
                    specRAs.append(specObj['RADeg'])
                    specDecs.append(specObj['decDeg'])
                    specLabels.append(str(specCount))
                if len(specRAs) > 0:
                    p.addPlotObjects(specRAs, specDecs, 'specObjects', objLabels = specLabels,
                                    size = sizeDeg/40.0*3600.0, symbol = 'box', color = "red")
                              
        if plotXMatch == "true":
            obj=self.sourceCollection.find_one({'name': name})
            xMatchRAs=[]
            xMatchDecs=[]
            xMatchLabels=[]
            if 'crossMatchCatalogs' in self.configDict.keys():
                for xMatchDict in self.configDict['crossMatchCatalogs']:
                    if "plotXMatch" in xMatchDict.keys() and xMatchDict['plotXMatch'] == False:
                        continue
                    label=xMatchDict['label']
                    RAKey='%s_RADeg' % (label)
                    decKey='%s_decDeg' % (label)
                    if RAKey in obj.keys() and decKey in obj.keys():
                        # We only want to show coords that are different, and not exactly cross-matched
                        # (e.g., cross matched zCluster results would have exact same RADeg, decDeg - useless to show)
                        if obj[RAKey] != obj['RADeg'] and obj[decKey] != obj['decDeg']:
                            xMatchRAs.append(obj[RAKey])
                            xMatchDecs.append(obj[decKey])
                            xMatchLabels.append(label)            
            # Other coords in obj dictionary that weren't yet picked up (useful for e.g. editable BCG coords)
            # Editable fields
            fieldsList=[]
            for fieldDict in self.configDict['fields']:
                fieldsList.append(str(fieldDict['name']))
            for key in fieldsList:
                if key.split("_")[-1] == 'RADeg':
                    if key.split("_")[0] not in xMatchLabels and key != "RADeg":
                        xMatchRAs.append(obj[key])
                        xMatchDecs.append(obj[key.replace("RADeg", "decDeg")])
                        xMatchLabels.append(str(key.split("_")[0])) 

            if len(xMatchRAs) > 0:
                p.addPlotObjects(xMatchRAs, xMatchDecs, 'xMatchObjects', objLabels = xMatchLabels,
                                 size = sizeDeg/40.0*3600.0, symbol = "diamond", color = 'cyan')
        
        if plotContours == "true":
            if 'contourImage' in self.configDict.keys() and self.configDict['contourImage'] != None:
                clipFileName=self.cacheDir+os.path.sep+self.configDict['contourImage']+os.path.sep+subDir+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".fits"
                if os.path.exists(clipFileName):
                    contourImg=pyfits.open(clipFileName)
                    contourWCS=astWCS.WCS(contourImg[0].header, mode = 'pyfits')
                    if self.configDict['contour1Sigma'] == "measureFromImage":
                        contourData=contourImg[0].data
                        # Choose level from clipped stdev
                        sigmaCut=2.0
                        mean=0
                        sigma=1e6
                        for i in range(20):
                            #nonZeroMask=np.not_equal(contourData, 0)
                            mask=np.less(abs(contourData-mean), sigmaCut*sigma)
                            #mask=np.logical_and(nonZeroMask, mask)
                            mean=np.mean(contourData[mask])
                            sigma=np.std(contourData[mask])
                    else:
                        sigma=self.configDict['contour1Sigma']
                    contourSigmaLevels=np.array(self.configDict['contourSigmaLevels'])
                    contourLevels=contourSigmaLevels*sigma
                    #contourLevels=[self.configDict['contour1Sigma'], 2*self.configDict['contour1Sigma'],
                                    #4*self.configDict['contour1Sigma'], 8*self.configDict['contour1Sigma'],
                                    #16*self.configDict['contour1Sigma']]
                    #contourLevels=np.linspace(self.configDict['contour1Sigma'], 
                                                #20*self.configDict['contour1Sigma'], 20)
                    p.addContourOverlay(contourImg[0].data, contourWCS, 'contour', levels = contourLevels, 
                                        width = self.configDict['contourWidth'],     
                                        color = self.configDict['contourColour'], 
                                        smooth = self.configDict['contourSmoothingArcsec'],
                                        highAccuracy = False)
                else:
                    plt.figtext(0.05, 0.05, "Adding contours failed - missing file: %s" % (clipFileName), color = 'red', backgroundcolor = 'black')
        
        # NOTE: Currently, this is not displaying images under python3 - something encoding related?
        #cherrypy.response.headers['Content-Type']="image/jpg"
        cherrypy.response.headers['Content-Type']="data:image/jpg;base64"
        buf=BytesIO()
        plt.savefig(buf, dpi = 96, format = 'jpg')
        plt.close()
        buf.seek(0)

        return base64.b64encode(buf.getvalue())
    
    
    @cherrypy.expose
    def makeSpectrumPlot(self, name, RADeg, decDeg):
        """Returns plot of spectrum that matches obj position (assuming it's within a few arcmin). 
        Marks on positions of spectral lines if redshift field is set appropriately.
        
        Relies on many parameters being set in the .config file.
        
        Returns None if no suitable matching spectrum found
        
        """

        # Just in case they are passed as strings (e.g., if direct from the url)
        RADeg=float(RADeg)
        decDeg=float(decDeg)
        name=self.URLToSourceName(name)
        
        obj=self.sourceCollection.find_one({'name': name})
        mongoDict=self.matchTags(obj)
        
        # Spin through to find matching file
        # If we have multiple files with the same RA, dec, then we'll take the last one we see
        fitsFiles=glob.glob(self.configDict['specDir']+os.path.sep+"*.fits")
        rDegMin=1e6
        matchFileName=None
        for f in fitsFiles:
            img=pyfits.open(f)
            headRA=img[0].header[self.configDict['specHeaderRAKey']]
            headDec=img[0].header[self.configDict['specHeaderDecKey']]
            if type(headRA) == str:
                headRA=astCoords.hms2decimal(headRA, ":")
            if type(headDec) == str:
                headDec=astCoords.dms2decimal(headDec, ":")
            rDeg=astCoords.calcAngSepDeg(RADeg, decDeg, headRA, headDec)
            if rDeg < rDegMin:
                rDegMin=rDeg
                matchFileName=f

        # Make plot
        if rDegMin > self.configDict['plotSizeArcmin']/60.:
            plt.figure(figsize=(10, 1))
            plt.figtext(0.5, 0.5, "No spectrum found", ha = 'center', va = 'center', size = 20)
        else:
            
            img=pyfits.open(matchFileName)
            wavelength=img[self.configDict['specExtName']].data[self.configDict['specLambdaKey']]
            flux=img[self.configDict['specExtName']].data[self.configDict['specFluxKey']]
            sky=img[self.configDict['specExtName']].data[self.configDict['specSkyKey']]
            flux=ndimage.uniform_filter1d(flux, int(self.configDict['specSmoothPix']))
            flux=flux/flux.max()
            
            plt.figure(figsize=(10,6))
            plt.title(matchFileName)
            plt.plot(wavelength, flux, 'k-')

            plt.xlabel("Wavelength (Angstroms)")
            plt.ylabel("Relative Flux")

            # Plots the spectral features in turn
            z=mongoDict[self.configDict['specRedshiftField']]
            yRange=np.linspace(0, flux.max()*1.2)
            for f in specFeatures.lineList:
                featureLabel=f[0]
                featureLambda=f[1]*(1+z)
                if featureLambda > wavelength.min() and featureLambda < wavelength.max():
                    # Greek letters? eta will cause a problem here!
                    featureLabel=featureLabel.replace("alpha", "$\\alpha$")
                    featureLabel=featureLabel.replace("beta", "$\\beta$")
                    featureLabel=featureLabel.replace("gamma", "$\gamma$")
                    featureLabel=featureLabel.replace("delta", "$\delta$")
                    featureLabel=featureLabel.replace("epsilon", "$\\epsilon$")
                    featureLabel=featureLabel.replace("zeta", "$\zeta$")
                    featureLabel=featureLabel.replace("theta", "$\\theta$")
                    plt.text(featureLambda, flux.max()*1.1, featureLabel, 
                            ha='center', va='top', size=12, rotation='vertical', color = 'red')
                    plt.plot([featureLambda]*len(yRange), yRange, 'k--')                
            
            # Sky, inc. main telluric absorption features
            plt.plot(wavelength, sky/sky.max()*0.3, 'b-', label='Sky')
            c=patches.Rectangle((6860, 0), (6930-6860), 1.2, fill=True, edgecolor=(0.8, 0.8, 0.8), 
                            facecolor=(0.8, 0.8, 0.8), linewidth=1)
            plt.gca().add_patch(c)
            c=patches.Rectangle((7590, 0), (7710-7590), 1.2, fill=True, edgecolor=(0.8, 0.8, 0.8), 
                            facecolor=(0.8, 0.8, 0.8), linewidth=1)
            plt.gca().add_patch(c)

            #pylab.legend(loc="upper right")
            #plt.savefig("test.jpg")

            plt.xlim(wavelength.min(), wavelength.max())
            plt.ylim(0, flux.max()*1.2)
            
        cherrypy.response.headers['Content-Type']="image/jpg"
        buf=StringIO.StringIO()
        plt.savefig(buf, dpi = 96, format = 'jpg')
        plt.close()
        buf.seek(0)
        
        return cherrypy.lib.file_generator(buf)
    
    
    def fetchSkyviewJPEG(self, name, RADeg, decDeg, RGBSurveysString, surveyLabel, lowSigmaCut = 2.0, highSigmaCut = 2.0,
                         refetch = False, cleanUp = True):
        """Fetch .fits images using skyview for the given survey - e.g., for 2MASS RGB = KHJ:
        
        RGBSurveysString = "2mass-k,2mass-h,2mass-j"
        
        Then makes an RGB .jpg.
        
        We will use surveyLabel to set where images are saved and as labels/links in the source browser.
        
        """
        
        # We need to fix up the skyview path in here later...
        if os.path.exists(self.cacheDir+os.path.sep+surveyLabel) == False:
            os.makedirs(self.cacheDir+os.path.sep+surveyLabel)
        
        rootFileName=self.cacheDir+os.path.sep+surveyLabel+os.path.sep+"%s.fits" % (name.replace(" ", "_"))
        outFileName=rootFileName.replace(".fits", ".jpg")
        
        if refetch == True or os.path.exists(outFileName) == False:
            #print "java -jar %s position=%.6f,%.6f output=%s size=%.3f pixels=%d survey=%s cache=%s/" % (self.skyviewPath, RADeg, decDeg, rootFileName, self.configDict['plotSizeArcmin']/60.0, self.configDict['plotSizePix'], RGBSurveysString, self.skyCacheDir)
            os.system("java -jar %s position=%.6f,%.6f output=%s size=%.3f pixels=%d survey=%s cache=%s/" % (self.skyviewPath, RADeg, decDeg, rootFileName, self.configDict['plotSizeArcmin']/60.0, self.configDict['skyviewPlotSizePix'], RGBSurveysString, self.skyCacheDir))
            
            # If skyview doesn't return something (e.g., FIRST in southern hemisphere) then copy the noData image
            if len(RGBSurveysString.split(",")) == 3:
                testFileName=rootFileName.replace(".fits", "_1.fits")
            else:
                testFileName=rootFileName
            if os.path.exists(testFileName) == False:
                noDataPath=sourcery.__path__[0]+os.path.sep+"data"+os.path.sep+"noData.jpg"
                os.system("cp %s %s" % (noDataPath, outFileName))
                return None
            
            # This handles either 3 surveys (r,g,b) or 1 survey (becomes greyscale)
            if len(RGBSurveysString.split(",")) == 3:
                RImg=pyfits.open(rootFileName.replace(".fits", "_1.fits"))
                GImg=pyfits.open(rootFileName.replace(".fits", "_2.fits"))
                BImg=pyfits.open(rootFileName.replace(".fits", "_3.fits"))
            else:
                RImg=pyfits.open(rootFileName)
                GImg=RImg
                BImg=RImg
                
            wcs=astWCS.WCS(RImg[0].header, mode = 'pyfits', zapKeywords = ['HISTORY', 'COMMENT'])
            
            R=RImg[0].data
            G=GImg[0].data
            B=BImg[0].data
            imData=np.array([RImg[0].data.transpose(), GImg[0].data.transpose(), BImg[0].data.transpose()])
            imData=imData.transpose()
            for i in range(imData.shape[2]):
                channel=imData[:, :, i]
                
                std=np.std(channel)
                med=np.median(channel)
                minLevel=med-lowSigmaCut*std
                
                # This is better
                freq, binEdges=np.histogram(channel, bins = int((channel.shape[0]*channel.shape[1])/100.0))
                binCentres=binEdges[:-1]+(binEdges[1]-binEdges[0])/2.0
                minLevel=binCentres[freq.tolist().index(freq.max())]
                
                lowMask=np.less(channel, minLevel)
                channel=channel-(minLevel)
                channel[lowMask]=0
                maxLevel=med+highSigmaCut*std
                if maxLevel > channel.max():
                    maxLevel=channel.max()
                highMask=np.greater(channel, maxLevel)
                channel=channel/maxLevel+0.001
                channel[highMask]=1.0
                channel=np.log10(channel)
                channel=channel-channel.min()
                channel=channel/channel.max()
                imData[:, :, i]=channel
            
            dpi=96.0
            plt.figure(figsize=(self.configDict['skyviewPlotSizePix']/dpi, self.configDict['skyviewPlotSizePix']/dpi), dpi = dpi)
            plt.axes([0, 0, 1, 1])
            plt.imshow(imData, interpolation="bilinear", origin='lower')
            plt.savefig(outFileName, dpi = dpi)
            plt.close()
            
            # Clean up .fits images
            if cleanUp == True:
                for i in range(1, 4):
                    toRemove=rootFileName.replace(".fits", "_%d.fits" % (i))
                    if os.path.exists(toRemove) == True:
                        os.remove(toRemove)
                    
        
    @cherrypy.expose
    def updateQueryParams(self, queryRADeg, queryDecDeg, querySearchBoxArcmin, queryOtherConstraints, 
                          queryApply = None, queryReset = None):
        """Updates query params in session, and then calls index again (which runs the query).
        
        """

        if not cherrypy.session.loaded: cherrypy.session.load()
        
        if queryReset:
            cherrypy.session['queryRADeg']="0:360"
            cherrypy.session['queryDecDeg']="-90:90"
            cherrypy.session['querySearchBoxArcmin']=""
            cherrypy.session['viewTopRow']=0
            cherrypy.session['queryOtherConstraints']=""
        
        if queryApply:
            cherrypy.session['queryRADeg']=queryRADeg
            cherrypy.session['queryDecDeg']=queryDecDeg
            cherrypy.session['querySearchBoxArcmin']=querySearchBoxArcmin
            cherrypy.session['viewTopRow']=0
            cherrypy.session['queryOtherConstraints']=queryOtherConstraints
        
        raise cherrypy.HTTPRedirect(cherrypy.request.script_name)


    def onLogin(self, username):
        """Called on successful login.
        
        Checks the permissions of the user and sets session variables accordingly
        
        """
        sourceryAuth.setEditPermissions(username, self.usersDict)
        
    
    def onLogout(self, username):
        """Called on logout"""
    
    
    def getLoginForm(self, username, msg="", from_page=cherrypy.request.script_name):
        return self.jinja_env.get_template('login.html').render(
            script_name=cherrypy.request.script_name,
            title=self.configDict.get('indexTitle', 'Sourcery Database'),
            from_page=from_page,
            msg=msg,
            username=username,
            login_message=self.configDict.get('loginMessage', ''),
            logo_data_url=self.logo_data_url,
        )
    
    
    @cherrypy.expose
    def login(self, username=None, password=None, from_page=cherrypy.request.script_name):
        if self.usersDict == {}:
            username="public"
            cherrypy.session[sourceryAuth.SESSION_KEY] = cherrypy.request.login = username
            self.onLogin(username) 
        else:
            if username is None or password is None:
                return self.getLoginForm("", from_page=from_page)
        
        error_msg = sourceryAuth.checkCredentials(username, password, self.usersDict, contactStr = self.contactInfo)
        if error_msg:
            return self.getLoginForm(username, error_msg, from_page)
        else:
            cherrypy.session[sourceryAuth.SESSION_KEY] = cherrypy.request.login = username
            self.onLogin(username)
        raise cherrypy.HTTPRedirect(from_page or cherrypy.request.script_name)
    
    
    @cherrypy.expose
    def logout(self, from_page=cherrypy.request.script_name):
        sess = cherrypy.session
        username = sess.get(sourceryAuth.SESSION_KEY, None)
        sess[sourceryAuth.SESSION_KEY] = None
        if username:
            cherrypy.request.login = None
            self.onLogout(username)
        raise cherrypy.HTTPRedirect(from_page or cherrypy.request.script_name)
    

    @cherrypy.expose
    @sourceryAuth.require()
    def index(self):
        """Shows the table page.
        
        """
        
        if os.path.exists(self.dbLockFileName) == True:
            return("""Apologies for the inconvenience - the database is being rebuilt. Please check again in about 10 minutes (probably much sooner).""")
        
        # Session variables: where in the table are we looking, query constraints
        if not cherrypy.session.loaded: cherrypy.session.load()
        if 'viewTopRow' not in cherrypy.session:
            cherrypy.session['viewTopRow']=0
        if 'queryRADeg' not in cherrypy.session:
            cherrypy.session['queryRADeg']="0:360"
        if 'queryDecDeg' not in cherrypy.session:
            cherrypy.session['queryDecDeg']="-90:90"
        if 'querySearchBoxArcmin' not in cherrypy.session:
            cherrypy.session['querySearchBoxArcmin']=""
        if 'queryOtherConstraints' not in cherrypy.session:
            cherrypy.session['queryOtherConstraints']=""
        queryRADeg=cherrypy.session.get('queryRADeg')
        queryDecDeg=cherrypy.session.get('queryDecDeg')
        querySearchBoxArcmin=cherrypy.session.get('querySearchBoxArcmin')
        queryOtherConstraints=cherrypy.session.get('queryOtherConstraints')
        
        queryPosts, numPosts=self.runQuery(queryRADeg, queryDecDeg, querySearchBoxArcmin, queryOtherConstraints)
        if 'numPosts' not in cherrypy.session:
            cherrypy.session['numPosts']=numPosts
        viewPosts=queryPosts[cherrypy.session['viewTopRow']:cherrypy.session['viewTopRow']+self.tableViewRows]

        # Censor queries on columns hidden from this user
        user=cherrypy.session['_sourcery_username']
        for userNameKey in self.usersDict.keys():
            userDict=self.usersDict[userNameKey]
            if userNameKey == user and 'hiddenXMatchTables' in userDict.keys():
                for prefix in userDict['hiddenXMatchTables']:
                    if queryOtherConstraints.find(prefix) != -1:
                        numPosts=0
                        viewPosts=[]

        # Quick links
        quick_links=[]
        for linkDict in self.configDict.get('quickLinks', []):
            url="updateQueryParams?queryRADeg=0%3A360&queryDecDeg=-90%3A90&querySearchBoxArcmin=&queryOtherConstraints="
            url=url+linkDict['constraints'].replace("+", "%2B")+"&queryApply=Apply"
            quick_links.append({'label': linkDict['label'], 'url': url})

        # Shareable query link
        share_url="updateQueryParams?queryRADeg=%s&queryDecDeg=%s&querySearchBoxArcmin=%s&queryOtherConstraints=%s&queryApply=Apply" % (
            self.sourceNameToURL(str(queryRADeg)), self.sourceNameToURL(str(queryDecDeg)),
            self.sourceNameToURL(str(querySearchBoxArcmin)), self.sourceNameToURL(queryOtherConstraints))
        share_query_url=cherrypy.request.base+cherrypy.request.script_name+"/"+share_url

        # Table columns — defaults plus any queried-on columns
        columnsShownList=[c['name'] for c in self.tableDisplayColumns]
        displayColumns=list(self.tableDisplayColumns)
        operators=["<", ">", "=", "!"]
        for logOp in [' and ', ' or ']:
            for c in queryOtherConstraints.split(logOp):
                for o in operators:
                    colName=c.split(o)[0].lstrip().rstrip()
                    if numPosts > 0 and colName in viewPosts[0].keys() and colName not in columnsShownList:
                        fieldTypeDict=self.fieldTypesCollection.find_one({'name': colName})
                        dispDict={'name': colName, 'label': colName}
                        if fieldTypeDict is None:
                            dispDict['fmt']='%s'
                        elif fieldTypeDict['type'] == 'number':
                            dispDict['fmt']='%.3f'
                        elif fieldTypeDict['type'] == 'text':
                            dispDict['fmt']='%s'
                        else:
                            raise Exception("unknown type for field '%s'" % (colName))
                        displayColumns.append(dispDict)
                        columnsShownList.append(colName)

        # Build table rows as lists of cell dicts
        table_rows=[]
        for obj in viewPosts:
            cells=[]
            for colDict in displayColumns:
                key=colDict['name']
                alignStr=colDict.get('tableAlign', 'center')
                if 'displaySize' in colDict:
                    widthStr='width: %dem;' % colDict['displaySize']
                    useDiv=True
                else:
                    useDiv=False
                    if colDict['fmt'] == '%s':
                        widthStr='width: %dem;' % len(obj[key]) if key in obj and type(obj[key]) == str else ''
                    else:
                        widthStr='width: %dem;' % len(colDict['fmt'] % 1)

                value=obj.get(key, 0.0 if colDict['fmt'] != '%s' else '')
                if key == 'name':
                    linkURL="displaySourcePage?sourceryID=%s&clipSizeArcmin=%.2f" % (
                        self.sourceNameToURL(obj['sourceryID']), self.configDict['defaultViewSizeArcmin'])
                    if 'defaultImageType' in self.configDict:
                        linkURL+="&imageType=%s" % self.configDict['defaultImageType']
                    value_html='<a href="%s" target=new>%s</a>' % (linkURL, obj['name'])
                elif key == 'NED_name' and obj.get('NED_name') != 'None':
                    nedName=obj[key]
                    nedLinkURL="http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES" % nedName.replace("+", "%2B").replace(" ", "+")
                    value_html='<a href=%s>%s</a>' % (nedLinkURL, nedName)
                elif key == 'NED_z' and obj.get('NED_z') == 'nan':
                    value_html='-'
                else:
                    try:
                        value_html=colDict['fmt'] % value
                    except Exception:
                        raise Exception("IndexError: check .config file tableDisplayColumns are actually in the .fits table, or for mixed '' \"\" inside []")
                cells.append({'value_html': value_html, 'align': alignStr, 'width_style': widthStr, 'use_div': useDiv})
            table_rows.append(cells)

        # Download links
        downloadLinkStr="downloadCatalog?queryRADeg=%s&queryDecDeg=%s&querySearchBoxArcmin=%s&queryOtherConstraints=%s&" % (
            queryRADeg, queryDecDeg, querySearchBoxArcmin, queryOtherConstraints)
        minimalDownloadLinkStr=downloadLinkStr+"minimalColumnSet=true&"
        download_link=quote_plus(downloadLinkStr, safe='&?=')
        minimal_download_link=quote_plus(minimalDownloadLinkStr, safe='&?=') if 'catalogDownloadMinimalColumns' in self.configDict else None

        # Meta data
        classificationDesc=self.configDict.get('classificationDescription', '')
        catalog_comments=self.configDict.get('catalogComments', '')
        if classificationDesc:
            catalog_comments=(catalog_comments+" <p>"+classificationDesc+"</p>") if catalog_comments else classificationDesc
        news_items=self.configDict.get('newsItems', [])[-5:] if 'newsItems' in self.configDict else []
        latest_news_str=("    –    Latest news: %s" % self.configDict['newsItems'][-1].split(":")[0]) if news_items else ''
        cache_rebuild_str="    –    [REBUILDING IMAGE CACHE]" if os.path.exists(self.cacheLockFileName) else ''
        hiddenConstraintsMessage=''
        for userNameKey in self.usersDict.keys():
            if userNameKey == user and 'hiddenConstraintsMessage' in self.usersDict[userNameKey]:
                hiddenConstraintsMessage="; %s" % self.usersDict[userNameKey]['hiddenConstraintsMessage']

        return self.jinja_env.get_template('index.html').render(
            script_name=cherrypy.request.script_name,
            title=self.configDict.get('indexTitle', 'Sourcery Database'),
            query_radeg=queryRADeg,
            query_decdeg=queryDecDeg,
            query_searchboxarcmin=querySearchBoxArcmin,
            query_otherconstraints=queryOtherConstraints,
            share_query_url=share_query_url,
            quick_links=quick_links,
            object_type_string=self.configDict['objectTypeString'],
            num_posts=numPosts,
            total_sources=self.sourceCollection.estimated_document_count(),
            hidden_constraints_message=hiddenConstraintsMessage,
            latest_news_str=latest_news_str,
            cache_rebuild_str=cache_rebuild_str,
            catalog_comments=catalog_comments,
            news_items=news_items,
            color_coding_rows=[],
            display_columns=displayColumns,
            table_rows=table_rows,
            download_link=download_link,
            minimal_download_link=minimal_download_link,
            catalog_download_name=self.configDict['catalogDownloadFileName'],
            hosted_by=self.configDict.get('hostedBy', ''),
            logo_data_url=self.logo_data_url,
        )


    @cherrypy.expose
    @sourceryAuth.require()
    def downloadCatalog(self, queryRADeg = "0:360", queryDecDeg = "-90:90", querySearchBoxArcmin = "",
                        queryOtherConstraints = "", fileFormat = "cat", minimalColumnSet = "false"):
        """Provide user with the current table view as a downloadable catalog.
        
        """

        # Fetch the cached table and update that with any changed classifications info
        cachedTabFileName=self.cacheDir+os.path.sep+"%s_xMatchedTable.fits" % (self.configDict['catalogDownloadFileName'])
        xTab=atpy.Table().read(cachedTabFileName)

        # Need this and change below for overloading with e.g. editable BCG coords
        editableFieldsList=[]
        for f in self.configDict['fields']:
            editableFieldsList.append(f['name'])
        for key in xTab.keys():
            if key in editableFieldsList:
                xTab.remove_column(key)
        
        # If minimal set of columns requested, delete all that aren't included
        minimalCols=xTab.keys()
        if minimalColumnSet == 'true':
            if 'catalogDownloadMinimalColumns' in self.configDict.keys():
                minimalCols=['sourceryID']
                for keyStr in self.configDict['catalogDownloadMinimalColumns']:
                    if keyStr[-1] != '_':
                        minimalCols.append(keyStr)
                    else:
                        for col in xTab.keys():
                            if col[:len(keyStr)] == keyStr:
                                minimalCols.append(col)
                for col in xTab.keys():
                    if col not in minimalCols:
                        xTab.remove_column(col)

        # NOTE: there may be fun unicode-related stuff here: e.g., u'BCG_RADeg' versus 'BCG_RADeg'
        keysList, typeNamesList, descriptionsList=self.getFieldNamesAndTypes(excludeKeys = [])
        keysToAdd=['sourceryID', 'RADeg', 'decDeg', 'classification']
        typeNamesToAdd=['text', 'number', 'number', 'text']
        for k, t in zip(keysList, typeNamesList):
            if k not in xTab.keys() and k not in keysToAdd and k in minimalCols:
                keysToAdd.append(k)
                typeNamesToAdd.append(t)
            if k not in keysToAdd and k in editableFieldsList:
                keysToAdd.append(k)
                typeNamesToAdd.append(t)

        # This part takes ~1 min for a catalog with 43k sources on laptop so is the only bit that needs speeding up
        # The query part is 0.05 sec - it's the big for loop that is slow [i.e., converting from posts -> astropy table]
        count=0
        keepIndices=[]
        posts, tabLength=self.runQuery(queryRADeg, queryDecDeg, querySearchBoxArcmin, queryOtherConstraints)
        tab=atpy.Table()
        tab.table_name=self.configDict['catalogDownloadFileName']
        for key, typeName in zip(keysToAdd, typeNamesToAdd):
            if typeName == 'number':
                tab.add_column(atpy.Column(np.zeros(tabLength, dtype = np.float64), str(key)))
            else:
                tab.add_column(atpy.Column(np.zeros(tabLength, dtype = 'S1000'), str(key)))
        for post in posts:
            # Trim xTab to match posts based on sourceryID
            if post['sourceryID'] in xTab['sourceryID']:
                keepIndices.append(np.where(xTab['sourceryID'] == post['sourceryID'])[0][0])
            for key in keysToAdd:
                if key in post.keys():           # NOTE: this handles image_ tags, which are 1 if present, and absent otherwise
                    tab[key][count]=post[key]
            count=count+1
        xTab=xTab[keepIndices]
        tab.rename_column('RADeg', 'tag_RADeg')
        tab.rename_column('decDeg', 'tag_decDeg')

        tab=atpy.join(xTab, tab, keys = ['sourceryID'])
        zapCols=['tag_RADeg', 'tag_decDeg', 'sourceryID', 'cacheBuilt']
        for z in zapCols:
            if z in tab.keys():
                tab.remove_column(z)

        # Zap any columns that this user shouldn't see
        user=cherrypy.session['_sourcery_username']
        for userNameKey in self.usersDict.keys():
            userDict=self.usersDict[userNameKey]
            if userNameKey == user and 'hiddenXMatchTables' in userDict.keys():
                for prefix in userDict['hiddenXMatchTables']:
                    colsToDelete=[]
                    for key in tab.keys():
                        if key.find(prefix) != -1 and key[:len(prefix)] == prefix:
                            colsToDelete.append(key)
                    tab.remove_columns(colsToDelete)

        # Zap any columns that contain only sentinels
        for key in tab.keys():
            try:
                if np.all(tab[key]) == -99:
                    tab.remove_column(key)
            except:
                pass

        tmpFile, tmpFileName=tempfile.mkstemp()
        if fileFormat == 'cat':
            tab.write(tmpFileName+".cat", format = 'ascii')
        elif fileFormat == 'fits':
            if len(tab.keys()) > 999:
                raise Exception("FITS table format does not support 1000+ columns")
            tab.write(tmpFileName+".fits", format = 'fits')
        elif fileFormat == 'reg':
            catalogTools.tab2DS9(tab, tmpFileName+".reg")
        
        cherrypy.response.headers['Content-Disposition']='attachment; filename="%s.%s"' % (self.configDict['catalogDownloadFileName'], fileFormat)
        
        # This may not be the nicest thing to do... but seems to work
        f=open(tmpFileName+"."+fileFormat, 'rb')
        
        return f
    
    
    @cherrypy.expose
    @sourceryAuth.require()
    def downloadThumbnailFITS(self, sourceryID):
        """Returns .fits image for download by the user.
        
        """
        obj=self.sourceCollection.find_one({'sourceryID': sourceryID})
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']
        subDir=str(RADeg).split(".")[0]
        imgPath=self.cacheDir+os.path.sep+self.configDict['downloadableFITS']+os.path.sep+subDir+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".fits"
        #imgPath=self.cacheDir+os.path.sep+self.configDict['downloadableFITS']+os.path.sep+catalogTools.makeRADecString(obj['RADeg'], obj['decDeg'])+".fits"
        cherrypy.response.headers['Content-Disposition']='attachment; filename="%s.%s"' % (obj['name'].replace(" ", "_"), 'fits')
        f=open(imgPath, 'rb')
        
        return f
        
        
    def sourceNameToURL(self, name):
        """Replaces + and spaces in source names so that they will be valid URLs.
        
        """
        return name.replace("+", "%2B").replace(" ", "%20").replace(":", "%3A")

        
    def URLToSourceName(self, url):
        """Replaces %20 and %2b in URLs with spaces and + signs.
        
        """
        return url.replace("%2b", "+").replace("%20", " ").replace("%3A", ":")


    def runQuery(self, queryRADeg, queryDecDeg, querySearchBoxArcmin, queryOtherConstraints, collection = 'source'):           
        """Runs a query, returns the posts found.
        
        If collection = 'source', runs on self.sourceCollection (default, whole catalog).
        
        If collection = 'tags', runs on self.tagsCollection
        
        """
        
        # Hidden constraints: for access control, e.g., show only DES users regions inside DES footprint
        user=cherrypy.session['_sourcery_username']
        hiddenConstraints=""
        for userNameKey in self.usersDict.keys():
            userDict=self.usersDict[userNameKey]
            if userNameKey == user and 'hiddenConstraints' in userDict.keys() and userDict['hiddenConstraints'] != None:
                hiddenConstraints=userDict['hiddenConstraints']
                if hiddenConstraints[0] == '"' or hiddenConstraints[0] == "'":
                    hiddenConstraints=hiddenConstraints[1:]
                if hiddenConstraints[-1] == '"' or hiddenConstraints[-1] == "'":
                    hiddenConstraints=hiddenConstraints[:-1]
        if hiddenConstraints != "":
            queryOtherConstraints=queryOtherConstraints+" and "+hiddenConstraints
                    
        # Build query document piece by piece...
        queryDict={}

        # Position
        if ":" not in queryRADeg and ":" not in queryDecDeg:
            queryRADeg=float(queryRADeg)
            queryDecDeg=float(queryDecDeg)
            querySearchBoxArcmin=float(querySearchBoxArcmin)
            RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(queryRADeg, queryDecDeg, querySearchBoxArcmin/60.0)
        else:
            RAMin, RAMax=queryRADeg.split(":")
            decMin, decMax=queryDecDeg.split(":")
            RAMin=float(RAMin)
            RAMax=float(RAMax)
            decMin=float(decMin)
            decMax=float(decMax)
        queryDict['decDeg']={'$lte': decMax, '$gte': decMin}
        if RAMin >= 0:
            queryDict['RADeg']={'$lte': RAMax, '$gte': RAMin}
        else:
            queryDict['$or']=[{'RADeg': {'$gte': 0, '$lte': RAMax}}, {'RADeg': {'$gte': 360+RAMin, '$lte': 360}}]

        # Other constraints
        constraintsDict=self.extractConstraintsDict(queryOtherConstraints)
        for key in constraintsDict:
            queryDict[key]=constraintsDict[key]

        # Execute query
        # NOTE: converting to list here is very slow
        if collection == 'source':
            self.sourceCollection.create_index([("RADeg", pymongo.ASCENDING)])
            #queryPosts=list(self.sourceCollection.find(queryDict).sort('decDeg').sort('RADeg'))
            numPosts=self.sourceCollection.count_documents(queryDict)
            queryPosts=self.sourceCollection.find(queryDict).sort('decDeg').sort('RADeg')  
        elif collection == 'tags':
            self.tagsCollection.create_index([("RADeg", pymongo.ASCENDING)])
            #queryPosts=list(self.sourceCollection.find(queryDict).sort('decDeg').sort('RADeg')) 
            numPosts=self.tagsCollection.count_documents(queryDict)
            queryPosts=self.tagsCollection.find(queryDict).sort('decDeg').sort('RADeg')
        else:
            raise Exception("collection should be 'source' or 'tags' only")
                        
        # If we wanted to store all this in its own collection
        #self.makeSessionCollection(queryPosts)
                
        return queryPosts, numPosts
        

    def makeSessionCollection(self, queryPosts):
        """Inserts all the posts from a query into a mongodb collection associated with this session. This
        is fairly slow for large databases (e.g., 10000+ objects). 
        
        Should not needed anymore, but left here just in case...
        
        """
        if not cherrypy.session.loaded: cherrypy.session.load()
            
        # Store results of query in another collection (empty it first if documents are in it)
        self.db.collection[cherrypy.session.id].remove({})
        
        # Add a date so we can expire the data after a couple of hours
        for q in queryPosts:
            q['lastModifiedDate']=datetime.datetime.utcnow()
            self.db.collection[cherrypy.session.id].insert_one(q)

        # This makes the session data self destruct after some time
        self.db.collection[cherrypy.session.id].create_index([('lastModifiedDate', 1)], expireAfterSeconds = 7200)

        
    def extractConstraintsDict(self, constraints):
        """Returns a dictionary of constraints extracted from string constraints, parsing all the operators
        >, <, =, ! etc.
        
        """
        
        # Order matters... equals should be last
        transDict={'<':  '$lt', 
                   '>':  '$gt',
                   '<=': '$lte',
                   '>=': '$gte',
                   '!=': '$ne',
                   '=': ''}

        # Newlines cause problems
        constraints=constraints.replace("\n", " ")
        
        # 'and' has precedence - which in practice means split on 'or' first, and then or all those together
        # Still need to handle () though, which means some recursion?
        orConstraints=constraints.split(' or ')
        allConstraintsDict={'$or': []}
        for orc in orConstraints:
            andConstraints=orc.split(" and ")
            constraintsDict={}
            for c in andConstraints:
                for op in transDict.keys():
                    bits=c.split(op)
                    if len(bits) == 2:
                        # Better way of checking for valid constraints than what we did before
                        key=bits[0].lstrip().rstrip()
                        value=bits[1].lstrip().rstrip()
                        before=c[c.find(op)-1]
                        after=c[c.find(op)+len(op)]
                        if before not in ['<', '>', '=', '!'] and after not in ['<', '>', '=', '!']:
                            validConstraint=True
                        else:
                            validConstraint=False
                        if validConstraint == True:    
                            # Strip " or ' from strings (saves confusion by user)
                            value=value.replace("'", "")
                            value=value.replace('"', '')
                            if key not in constraintsDict.keys():
                                constraintsDict[key]={}
                            # Queries won't work if we use strings instead of numbers when needed...
                            if op not in ['=', '!=']:
                                try:
                                    constraintsDict[key][transDict[op]]=float(value)
                                except:
                                    constraintsDict[key][transDict[op]]=value
                            else:
                                if op == '=':
                                    opStr="$in"
                                elif op == '!=':
                                    opStr="$nin"
                                if opStr not in constraintsDict[key].keys():
                                    constraintsDict[key][opStr]=[]
                                try:
                                    constraintsDict[key][opStr].append(float(value))
                                except:
                                    if '*' in value:
                                        regex='(?i)'    # make case insensitive
                                        if value[0] != '*':
                                            regexStr=regex+"^"+value
                                        else:
                                            regexStr=regex+value
                                        regexStr=regexStr.replace("*", ".*")
                                        constraintsDict[key][opStr].append(re.compile(regexStr))
                                    else:
                                        constraintsDict[key][opStr].append(value)
            allConstraintsDict['$or'].append(constraintsDict)
            #IPython.embed()
            #sys.exit()
        #print constraintsDict
        
        return allConstraintsDict
    
    
    @cherrypy.expose
    def changeTablePage(self, nextButton = None, prevButton = None):
        """Changes the viewed table page.
        
        """
        if not cherrypy.session.loaded: cherrypy.session.load()
        # Sometimes after a while we lose our session: need to check this out...
        if 'viewTopRow' not in cherrypy.session:
            cherrypy.session['viewTopRow']=0
        if 'queryRADeg' not in cherrypy.session:
            cherrypy.session['queryRADeg']="0:360"
        if 'queryDecDeg' not in cherrypy.session:
            cherrypy.session['queryDecDeg']="-90:90"
        if 'querySearchBoxArcmin' not in cherrypy.session:
            cherrypy.session['querySearchBoxArcmin']=""
        if 'queryOtherConstraints' not in cherrypy.session:
            cherrypy.session['queryOtherConstraints']=""
            
        if nextButton:
            viewTopRow=cherrypy.session.get('viewTopRow')
            viewTopRow=viewTopRow+self.tableViewRows
            if viewTopRow >= cherrypy.session.get('numPosts'):
                viewTopRow=cherrypy.session.get('numPosts')-self.tableViewRows
            cherrypy.session['viewTopRow']=viewTopRow
        if prevButton:
            viewTopRow=cherrypy.session.get('viewTopRow')
            viewTopRow=viewTopRow-self.tableViewRows
            if viewTopRow < 0:
                viewTopRow=0
            cherrypy.session['viewTopRow']=viewTopRow
            
        raise cherrypy.HTTPRedirect(cherrypy.request.script_name)
    
    
    @cherrypy.expose
    def displayConstraintsHelp(self):
        """Displays a page which lists all columns in the source list in a table and their types, together
        with a little blurb on how to define constraints.
        
        """
                
        excludeKeys=['RADeg', 'decDeg', 'sourceryID', 'cacheBuilt']
        keysList, typeNamesList, descriptionsList=self.getFieldNamesAndTypes(excludeKeys=excludeKeys)
        return self.jinja_env.get_template('constraints_help.html').render(
            script_name=cherrypy.request.script_name,
            hosted_by=self.configDict.get('hostedBy', ''),
            fields=list(zip(keysList, typeNamesList, descriptionsList)),
            logo_data_url=self.logo_data_url,
        )


    def getFieldNamesAndTypes(self, excludeKeys = []):
        """Fetches lists of field names, types and descriptions, for when displaying constraints help and 
        saving tables.
        
        """
        keysList=[]
        typeNamesList=[]
        descList=[]
        for post in self.fieldTypesCollection.find().sort('index'):
            if post['name'] not in excludeKeys:
                keysList.append(post['name'])
                typeNamesList.append(post['type'])
                descList.append(post['description'])
        # Add MongoDB object properties        
        if 'classifications' in self.configDict.keys():
            keysList.append('classification')
            typeNamesList.append("text")
            descList.append(self.configDict['classificationDescription'])
        if 'fields' in self.configDict.keys():
            for fieldDict in self.configDict['fields']:
                keysList.append(fieldDict['name'])
                typeNamesList.append(fieldDict['type'])
                descList.append(fieldDict['description'])

        # Override descriptions from the ones in the file
        for k in self.descriptionsDict.keys():
            if k in keysList:
                index=keysList.index(k)
                descList[index]=self.descriptionsDict[k]

        return keysList, typeNamesList, descList
                
    
    @cherrypy.expose
    def updateMongoDB(self, sourceryID, returnURL, **kwargs):
        """Update info on source in MongoDB.
        
        """
        
        if not cherrypy.session.loaded: cherrypy.session.load()
        
        obj=self.sourceCollection.find_one({'sourceryID': sourceryID})
        name=obj['name']

        # Bizarrely, legacy coordinates are given as degrees (lon, lat) but max distance has to be in radians...
        # Also, need lon between -180, +180
        if obj['RADeg'] > 180:
            lon=360.0-obj['RADeg']
        else:
            lon=obj['RADeg']
        matches=self.tagsCollection.find({'loc': SON({'$nearSphere': [lon, obj['decDeg']], '$maxDistance': np.radians(self.configDict['MongoDBCrossMatchRadiusArcmin']/60.0)})}).limit(1)
        mongoDict=matches.next()
        
        post={}
        for key in kwargs.keys():
            for fieldDict in self.configDict['fields']:
                if key == fieldDict['name']:
                    if fieldDict['type'] == 'number':
                        post[key]=float(kwargs[key])
                    else:
                        post[key]=kwargs[key]
            if key == 'classification':
                post[key]=kwargs[key]
        post['lastUpdated']=datetime.date.today().isoformat()
        post['user']=cherrypy.session['_sourcery_username']
        self.tagsCollection.update_one({'_id': mongoDict['_id']}, {'$set': post}, upsert = False)
                
        # Update source collection too - here we will do this for all sources that share the same name (we can have multiple source lists)
        # This could be done using cross matching based on coords instead, but this could cause confusion in the case of multiple sources
        # that are not the same object, located at separations < MongoDBCrossMatchRadiusArcmin
        # Example of this: erroneously deblended objects in the XCS list, where you only want to flag some objects as junk, and others not
        objs=self.sourceCollection.find({'name': name})
        for obj in objs:
            self.sourceCollection.update_one({'_id': obj['_id']}, {'$set': post}, upsert = False)
        
        # Would reset zoom level if changed
        if 'defaultImageType' in self.configDict.keys():
            imageType=self.configDict['defaultImageType']
        else:
            imageType="SDSS"
        
        raise cherrypy.HTTPRedirect(cherrypy.request.script_name+"/displaySourcePage?sourceryID=%s" % (self.sourceNameToURL(obj['sourceryID'])))
        #return self.displaySourcePage(obj['sourceryID'], clipSizeArcmin = self.configDict['plotSizeArcmin'],
                                      #imageType = imageType)


    def offlineUpdateTags(self, name, tagsToInsertDict):
        """Update tags for object matching name, offline version. Used by sourcery_fast_tag.
        
        """
               
        obj=self.sourceCollection.find_one({'name': name})
        
        # Bizarrely, legacy coordinates are given as degrees (lon, lat) but max distance has to be in radians...
        # Also, need lon between -180, +180
        if obj['RADeg'] > 180:
            lon=360.0-obj['RADeg']
        else:
            lon=obj['RADeg']
        matches=self.tagsCollection.find({'loc': SON({'$nearSphere': [lon, obj['decDeg']], '$maxDistance': np.radians(self.configDict['MongoDBCrossMatchRadiusArcmin']/60.0)})}).limit(1)
        mongoDict=matches.next()
        
        post={}
        for key in tagsToInsertDict.keys():
            for fieldDict in self.configDict['fields']:
                if key == fieldDict['name']:
                    if fieldDict['type'] == 'number':
                        post[key]=float(kwargs[key])
                    else:
                        post[key]=kwargs[key]
            if key == 'classification':
                post[key]=kwargs[key]
        post['lastUpdated']=datetime.date.today().isoformat()
        
        #print "How to avoid overwrites of things we shouldn't?"
        #IPython.embed()
        #sys.exit()
        
        self.tagsCollection.update_one({'_id': mongoDict['_id']}, {'$set': post}, upsert = False)
        
        # Update source collection too
        self.sourceCollection.update_one({'_id': obj['_id']}, {'$set': post}, upsert = False)
        
    
    @cherrypy.expose
    @sourceryAuth.require()
    def displaySourcePage(self, sourceryID, imageType = 'best', clipSizeArcmin = None, plotNEDObjects = "false",\
                          plotSpecObjects = "false", plotSourcePos = "false", plotXMatch = "false", plotContours = "false",\
                          showAxes = "false", gamma = 1.0, plotRedshift = "false", redshift = "none"):
        """Retrieve data on a source and display source page, showing given image plot.
        
        This should have form controls somewhere for editing the assigned redshift, redshift type, redshift 
        source and candidate status (e.g., confirmed cluster, junk etc.).
        
        imageType should match to an entry in self.imDirLabelsList
        
        """
        if not cherrypy.session.loaded: cherrypy.session.load()
        if 'viewTopRow' not in cherrypy.session:
            cherrypy.session['viewTopRow']=0
        if 'queryRADeg' not in cherrypy.session:
            cherrypy.session['queryRADeg']="0:360"
        if 'queryDecDeg' not in cherrypy.session:
            cherrypy.session['queryDecDeg']="-90:90"
        if 'querySearchBoxArcmin' not in cherrypy.session:
            cherrypy.session['querySearchBoxArcmin']=""
        if 'queryOtherConstraints' not in cherrypy.session:
            cherrypy.session['queryOtherConstraints']=""
        queryRADeg=cherrypy.session.get('queryRADeg')
        queryDecDeg=cherrypy.session.get('queryDecDeg')
        querySearchBoxArcmin=cherrypy.session.get('querySearchBoxArcmin')
        queryOtherConstraints=cherrypy.session.get('queryOtherConstraints')
        
        user=cherrypy.session['_sourcery_username']

        sourceryID=self.URLToSourceName(sourceryID)

        obj=self.sourceCollection.find_one({'sourceryID': sourceryID})
        mongoDict=self.matchTags(obj)
        name=obj['name']

        # For avoiding display of e.g. catalogs in which we don't have a cross match
        skipColumnPrefixList=[]
        for key in obj.keys():
            prefix=key.split("_")[0]
            matchKey="%s_match" % (prefix)
            if matchKey in obj.keys() and obj[matchKey] == 0 and prefix not in skipColumnPrefixList:
                skipColumnPrefixList.append(prefix)

        # For censoring some cross match tables from particular users
        for userNameKey in self.usersDict.keys():
            userDict=self.usersDict[userNameKey]
            if userNameKey == user and 'hiddenXMatchTables' in userDict.keys():
                for prefix in userDict['hiddenXMatchTables']:
                    matchKey="%s_match" % (prefix)
                    if matchKey in obj.keys() and obj[matchKey] == 1 and prefix not in skipColumnPrefixList:
                        skipColumnPrefixList.append(prefix)

        # Pick the best available image given the preference given in the config file
        if imageType == 'best':
            for key in self.configDict['imagePrefs']:
                if 'image_%s' % (key) in obj.keys() and obj['image_%s' % (key)] == 1:
                    imageType=key
                    break
            if imageType == 'best':
                imageType='DECaLS'

        # sourceryConfig JSON for JS
        sourcery_config_dict={
            'objectName': obj['name'],
            'RADeg': obj['RADeg'],
            'decDeg': obj['decDeg'],
            'redshift': obj.get('redshift', -99),
            'defaultClipSizeArcmin': self.configDict['defaultViewSizeArcmin'],
            'plotDisplayWidthPix': self.configDict['plotDisplayWidthPix'],
            'clickableRAField': self.configDict.get('clickableRAField', ''),
            'clickableDecField': self.configDict.get('clickableDecField', ''),
        }
        sourcery_config_json=json.dumps(sourcery_config_dict)

        # Image type radio buttons
        image_types=[{'label': label, 'max_size_arcmin': maxSize, 'selected': label == imageType}
                     for label, maxSize in zip(self.imDirLabelsList, self.imDirMaxSizeArcminList)]

        # Checkbox states
        checked_ned=plotNEDObjects == "true"
        checked_spec=plotSpecObjects == "true"
        checked_source_pos=plotSourcePos == "true"
        checked_xmatch=plotXMatch == "true"
        checked_contours=plotContours == "true"
        checked_show_axes=showAxes == "true"
        checked_redshift=plotRedshift == "true"

        has_contours='contourImage' in self.configDict and self.configDict['contourImage'] is not None
        contour_image=self.configDict.get('contourImage', '')

        has_download_fits='downloadableFITS' in self.configDict
        download_fits_label=self.configDict.get('downloadableFITS', '')

        imageMaxSizeArcmin=30.
        current_size_arcmin=clipSizeArcmin if clipSizeArcmin is not None else self.configDict['defaultViewSizeArcmin']
        plot_display_height_pix=self.configDict['plotDisplayWidthPix'] * 1.02

        # Editing controls
        has_edit_controls='fields' in self.configDict or 'classifications' in self.configDict
        edit_permission=cherrypy.session.get('editPermission', False)

        classifications=[]
        if 'classifications' in self.configDict:
            for i, c in enumerate(self.configDict['classifications']):
                classifications.append({
                    'value': c,
                    'checked': c == mongoDict.get('classification', ''),
                    'first': i == 0,
                })

        editable_fields=[]
        last_updated='-'
        editor_user='-'
        if 'fields' in self.configDict:
            for fieldDict in self.configDict['fields']:
                editable_fields.append({
                    'name': fieldDict['name'],
                    'value': str(mongoDict.get(fieldDict['name'], '')),
                    'display_size': fieldDict['displaySize'],
                })
            last_updated=mongoDict.get('lastUpdated', '-')
            editor_user=mongoDict.get('user', '-')

        # Optional spectrum
        if 'displaySpectra' in self.configDict and self.configDict['displaySpectra'] == True:
            spectrum_url="makeSpectrumPlot?name=%s&RADeg=%s&decDeg=%s" % (
                self.sourceNameToURL(name), obj['RADeg'], obj['decDeg'])
            spec_display_width_pix=self.configDict['specPlotDisplayWidthPix']
        else:
            spectrum_url=''
            spec_display_width_pix=0

        # NED matches
        self.fetchNEDInfo(obj['name'], obj['RADeg'], obj['decDeg'])
        self.findNEDMatch(obj, NEDObjTypes=self.configDict['NEDObjTypes'])
        nedFileName=self.nedDir+os.path.sep+obj['name'].replace(" ", "_")+".txt"
        nedObjs=catalogTools.parseNEDResult(nedFileName, onlyObjTypes=self.configDict['NEDObjTypes'])
        ned_rows=[]
        for i in range(len(nedObjs['RAs'])):
            nedLinkURL=("http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s"
                        "&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1"
                        "&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude"
                        "&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES"
                        % (nedObjs['names'][i].replace("+", "%2B").replace(" ", "+")))
            ned_rows.append({
                'label': nedObjs['labels'][i],
                'name': nedObjs['names'][i],
                'ned_url': nedLinkURL,
                'RA': nedObjs['RAs'][i],
                'dec': nedObjs['decs'][i],
                'source_type': nedObjs['sourceTypes'][i],
                'redshift': nedObjs['redshifts'][i],
            })

        # Spectroscopic redshifts
        spec_rows=[]
        if 'addSpecRedshifts' in self.configDict and self.configDict['addSpecRedshifts'] == True:
            specRedshifts=catalogTools.fetchSpecRedshifts(obj['name'], obj['RADeg'], obj['decDeg'],
                                                          redshiftsTable=self.specRedshiftsTab)
            for i, specObj in enumerate(specRedshifts, 1):
                spec_rows.append({
                    'id': i,
                    'RA': specObj['RADeg'],
                    'dec': specObj['decDeg'],
                    'redshift': specObj['z'],
                    'catalog': specObj['catalog'],
                })

        # Source properties
        properties=[]
        fieldTypes=self.fieldTypesCollection.find().sort('index')
        for f in fieldTypes:
            if f['name'] in obj.keys() and f['name'] not in ['sourceryID', 'cacheBuilt']:
                pkey=f['name']
                prefix=pkey.split("_")[0]
                if prefix in skipColumnPrefixList:
                    continue
                if obj[pkey] is None or obj[pkey] == "None":
                    value_html="-"
                elif pkey == "NED_name":
                    nedName=obj[pkey]
                    nedLinkURL=("http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s"
                                "&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1"
                                "&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude"
                                "&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES"
                                % (nedName.replace("+", "%2B").replace(" ", "+")))
                    value_html="<a href=%s>%s</a>" % (nedLinkURL, nedName)
                else:
                    value_html=str(obj[pkey])
                properties.append({'label': pkey, 'value_html': value_html})

        return self.jinja_env.get_template('source.html').render(
            script_name=cherrypy.request.script_name,
            source_name=name,
            hosted_by=self.configDict.get('hostedBy', ''),
            sourcery_config_json=sourcery_config_json,
            plot_display_height_pix=plot_display_height_pix,
            image_types=image_types,
            has_download_fits=has_download_fits,
            download_fits_label=download_fits_label,
            sourcery_id=sourceryID,
            object_name=obj['name'],
            return_url=cherrypy.url(),
            checked_show_axes=checked_show_axes,
            checked_ned=checked_ned,
            checked_spec=checked_spec,
            checked_source_pos=checked_source_pos,
            checked_xmatch=checked_xmatch,
            checked_contours=checked_contours,
            checked_redshift=checked_redshift,
            has_contours=has_contours,
            contour_image=contour_image,
            max_size_arcmin=imageMaxSizeArcmin,
            current_size_arcmin=current_size_arcmin,
            current_gamma=gamma,
            has_edit_controls=has_edit_controls,
            edit_permission=edit_permission,
            classifications=classifications,
            editable_fields=editable_fields,
            last_updated=last_updated,
            editor_user=editor_user,
            spectrum_url=spectrum_url,
            spec_display_width_pix=spec_display_width_pix,
            ned_rows=ned_rows,
            spec_rows=spec_rows,
            properties=properties,
            logo_data_url=self.logo_data_url,
        )
    
    
    def addImageDirTags(self):
        """Adds image_ fields to the database for image dirs where matching images are found in the cache.
        This takes about 1-2 min to run (serial) for databases with ~30,000 objects.
        
        """
        print(">>> Adding image_* tags ...")
        for label in self.imDirLabelsList:
            if self.fieldTypesCollection.find_one({'name': 'image_%s' % (label)}) == None:
                keysList, typeNamesList, descriptionsList=self.getFieldNamesAndTypes()
                fieldDict={}
                fieldDict['name']='image_%s' % (label)
                fieldDict['type']='number'
                fieldDict['description']='1 if object has image in the database; 0 otherwise'
                fieldDict['index']=len(keysList)+1
                self.fieldTypesCollection.insert_one(fieldDict)
        t0=time.time()
        # We need to do this to avoid hitting 32 Mb limit below when using large databases
        self.sourceCollection.create_index([("RADeg", pymongo.ASCENDING)])
        posts=self.sourceCollection.find({}).sort('decDeg').sort('RADeg')
        for obj in posts:
            name=obj['name']
            RADeg=obj['RADeg']
            decDeg=obj['decDeg']
            #print("... %s ..." % (name))
            fileNameLabel=catalogTools.makeRADecString(RADeg, decDeg)
            # Add imageDir tags
            # NOTE: we look in other survey dirs (e.g., SDSS) while we're here
            minSizeBytes=40000
            for label in self.imDirLabelsList:
                f=self.configDict['cacheDir']+os.path.sep+label+os.path.sep+fileNameLabel+".jpg"
                if os.path.exists(f) == True:
                    # image size check: don't include SDSS if image size is tiny as no data
                    skipImage=False
                    if os.stat(f).st_size < minSizeBytes and label == 'SDSS':
                        skipImage=True
                    if skipImage == False:
                        post={'image_%s' % (label): 1}
                        self.sourceCollection.update_one({'_id': obj['_id']}, {'$set': post}, upsert = False)
        t1=time.time()
        print("... took %.3f sec ..." % (t1-t0))
        
        
    def preprocess(self):
        """This re-runs pre-processing steps (e.g., NED matching, SDSS image fetching etc.).
        
        If the user specified their own imageDirs, then the .jpg images from these are constructed here
        
        Directories containing ready-made .jpgs can also be directly added into the cacheDir folder. So long
        as these have a corresponding entry in the .config file they will be picked up. We spin through
        those folders and also add 'image_<imageDirLabel>' tags in the MongoDB too. 
        
        """
        
        # So we can display a status message on the index page in other processes if the cache is being rebuilt
        outFile=open(self.cacheLockFileName, "wb")
        outFile.close()

        # Add image_* tags first - so that queries that need this still work while cache rebuilding
        self.addImageDirTags()
        
        # Make .jpg images from local, user-supplied .fits images
        if 'imageDirs' in self.configDict.keys():
            self.makeImageDirJPEGs()
            
        # tileDirs set-up - DES, KiDS, IAC-S82 etc..
        if 'tileDirs' in self.configDict.keys():
            print(">>> Setting up tileDir WCS info ...")
            for tileDirDict in self.configDict['tileDirs']:
                if tileDirDict['label'] not in self.tileDirs.keys():
                    self.tileDirs[tileDirDict['label']]=tileDir.TileDir(tileDirDict['label'], tileDirDict['path'], self.cacheDir,
                                                                        sizePix = tileDirDict['sizePix'])
                self.tileDirs[tileDirDict['label']].setUpWCSDict()
                                              
        # We need to do this to avoid hitting 32 Mb limit below when using large databases
        self.sourceCollection.create_index([("RADeg", pymongo.ASCENDING)])
        
        # Threaded
        # NOTE: Threads have occassionally given weird issues (e.g., mismatched WISE images)
        # Check that when thread write to disk they don't clash with each other
        if 'threadedCacheBuild' in self.configDict.keys() and self.configDict['threadedCacheBuild'] == True:
            cursor=self.sourceCollection.find({'cacheBuilt': 0}, no_cursor_timeout = True, session = self.mongoSess).sort('decDeg').sort('RADeg')
            sourceryIDs=[]
            for obj in cursor:
                sourceryIDs.append(obj['sourceryID'])
            cursor.close()
            import multiprocessing
            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor(max_workers = multiprocessing.cpu_count()) as executor:
                executor.map(self.buildCacheForObject, sourceryIDs)
        else:
            # Serial - still useful if debugging
            cursor=self.sourceCollection.find({'cacheBuilt': 0}, no_cursor_timeout = True, session = self.mongoSess).sort('decDeg').sort('RADeg')
            lastRefreshTime=time.time()
            for obj in cursor:
                # Keeping the session alive - refresh every 10 min or so
                if time.time()-lastRefreshTime > 600:
                    self.db.command({"refreshSessions" : [self.mongoSess.session_id]})
                    lastRefreshTime=time.time()
                self.buildCacheForObject(obj['sourceryID'], refetch = False)
            cursor.close()
                
        self.addImageDirTags()
            
        # This will stop index displaying "cache rebuilding" message
        if os.path.exists(self.cacheLockFileName) == True:
            os.remove(self.cacheLockFileName)

    
    @cherrypy.expose
    def buildCacheForObject(self, sourceryID, refetch = False, from_page = None):
        """Given a sourceryID (unique ID number), (re)fetch all the available imaging.
        As well as allowing threading, this also enables 'spot fixes' by clicking a button on the candidate 
        page (so if user spots image coords off, they can fix rather than manually deleting / re-running
        build cache). 
        
        NOTE: We may remove the Rebuild Cache button if the change to coord-based file names is accurate
        enough.
        
        """
        
        # If called via web...
        if refetch == "true":
            refetch=True
        
        obj=self.sourceCollection.find_one({'sourceryID': sourceryID})

        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']
        fileNameLabel=catalogTools.makeRADecString(RADeg, decDeg)
                    
        print(">>> Fetching data to cache for object %s" % (name))
        self.fetchNEDInfo(name, RADeg, decDeg)
        # Web services
        #if self.configDict['addSpecRedshifts'] == True:
            #catalogTools.fetchSDSSRedshifts(self.sdssRedshiftsDir, name, RADeg, decDeg,
                                            #redshiftsTable = self.specRedshiftsTab)
        if self.configDict['addSDSSImage'] == True:
            self.fetchSDSSImage(name, RADeg, decDeg, refetch = refetch)
        if self.configDict['addDECaLSImage'] == True:
            self.fetchDECaLSImage(name, RADeg, decDeg, refetch = refetch)
        if self.configDict['addHSCImage'] == True:
            self.fetchHSCImage(name, RADeg, decDeg, refetch = refetch)
        if self.configDict['addPS1Image'] == True:
            self.fetchPS1Image(name, RADeg, decDeg, refetch = refetch, bands = 'gri')
        if self.configDict['addPS1IRImage'] == True:
            self.fetchPS1Image(name, RADeg, decDeg, refetch = refetch, bands = 'izy')
        if self.configDict['addUnWISEImage'] == True:
            self.fetchUnWISEImage(name, RADeg, decDeg, refetch = refetch)
        # Tile dirs
        for key in self.tileDirs.keys():
            self.tileDirs[key].fetchImage(name, RADeg, decDeg, self.configDict['plotSizeArcmin'], refetch = refetch)
        # Skyview (not maintained or tested recently as very slow)
        if 'skyviewLabels' in self.configDict.keys():
            for surveyString, label in zip(self.configDict['skyviewSurveyStrings'], self.configDict['skyviewLabels']):
                self.fetchSkyviewJPEG(name, RADeg, decDeg, surveyString, label, refetch = refetch)
        
        # Add imageDir tags
        # NOTE: we look in other survey dirs (e.g., SDSS) while we're here
        minSizeBytes=40000
        for label in self.imDirLabelsList:
            f=self.configDict['cacheDir']+os.path.sep+label+os.path.sep+fileNameLabel+".jpg"
            if os.path.exists(f) == True:
                # image size check: don't include SDSS if image size is tiny as no data
                skipImage=False
                if os.stat(f).st_size < minSizeBytes and label == 'SDSS':
                    skipImage=True
                if skipImage == False:
                    post={'image_%s' % (label): 1}
                    self.sourceCollection.update_one({'_id': obj['_id']}, {'$set': post}, upsert = False)
        
        # Flag this as done
        self.sourceCollection.update_one({'_id': obj['_id']}, {'$set': {'cacheBuilt': 1}}, upsert = False)
        
        # If using for spot fixes in web interface, need to refresh page when done
        if from_page != None:
            raise cherrypy.HTTPRedirect(cherrypy.request.script_name+"/"+from_page)
    

    def makeImageDirJPEGs(self):
        """Actually makes .jpg images from .fits images in given directories. We figure out which image to use
        from spinning through the headers. 
        
        For cases where there may be more than one suitable image, you can use matchKey in the config file to
        specify a field in the sources database to match on (e.g., obsID in the case of XCS). Otherwise, the
        first image in which the object is found will be taken.
        
        If we were going to add a clickable map image (e.g., for ACT), this would be the place to put it
        in...
                        
        """
        print(">>> Making imageDir .jpgs ...")
        for imDirDict in self.configDict['imageDirs']:
            
            imageDir=imDirDict['path']
            label=imDirDict['label']
            colorMap=imDirDict['colorMap']
            sizePix=imDirDict['sizePix']
            if 'minMaxRadiusArcmin' in imDirDict.keys():
                minMaxRadiusArcmin=imDirDict['minMaxRadiusArcmin']
            else:
                minMaxRadiusArcmin=None
            scaling=imDirDict['scaling']
            matchKey=imDirDict['matchKey']
            if 'maxSizeArcmin' in imDirDict.keys():
                maxSizeArcmin=imDirDict['maxSizeArcmin']
            else:
                maxSizeArcmin=self.configDict['plotSizeArcmin']
            
            print("... %s ..." % (label))

            if 'skipMakingNewImages' in self.configDict.keys() and label in self.configDict['skipMakingNewImages']:
                print("... WARNING: skipMakingNewImages enabled for %s ..." % (label))
                continue
                    
            # NOTE: Need to worry at some point about labels with spaces...
            outDir=self.configDict['cacheDir']+os.path.sep+label
            if os.path.exists(outDir) == False:
                os.makedirs(outDir)
                            
            imgList=glob.glob(imageDir+os.path.sep+"*.fits")
            imgList=imgList+glob.glob(imageDir+os.path.sep+"*.fits.gz")

            # Pickled dictionary of WCS headers, for speed
            pickleFileName=imageDir+os.path.sep+"headerDict.pickled"
            if os.path.exists(pickleFileName) == True:
                pickleFile=open(pickleFileName, "rb")
                unpickler=pickle.Unpickler(pickleFile)
                headerDict=unpickler.load()
                pickleFile.close()
            else:
                headerDict={}
                # Takes ~160 sec to build the first time for ~10,000 images, ~2 sec for subsequent runs
                print("... building headerDict pickle for .fits images under %s/ ..." % (imageDir))
                t0=time.time()
                origLength=len(headerDict.keys())
                for imgFileName in imgList:
                    if imgFileName not in headerDict.keys():
                        with pyfits.open(imgFileName) as img:
                            for ext in img:
                                if ext.data is not None:
                                    break
                            wcs=astWCS.WCS(ext.header, mode = 'pyfits')
                        headerDict[imgFileName]=wcs.header.copy()
                t1=time.time()
                # Write pickled headerDict, in case it was updated
                if len(headerDict.keys()) > origLength:
                    print("... writing updated headerDict pickle to %s/ ..." % (imageDir))
                    pickleFile=open(pickleFileName, "wb")
                    pickler=pickle.Pickler(pickleFile)
                    pickler.dump(headerDict)
                    pickleFile.close()
            
            # Convert headerDict to wcsDict - takes ~30 sec for ~10,000 images
            t0=time.time()
            wcsDict={}
            for key in headerDict.keys():
                wcsDict[key]=astWCS.WCS(headerDict[key], mode = 'pyfits')
            t1=time.time()
                        
            # If only one image, set flag so that we will enable map web page
            # This needs to go somewhere else...
            if len(imgList) == 1:
                self.mapPageEnabled=True
            else:
                self.mapPageEnabled=False
            
            self.sourceCollection.create_index([("RADeg", pymongo.ASCENDING)])  # Avoid 32 Mb limit
            objList=self.sourceCollection.find(no_cursor_timeout = True, session = self.mongoSess).sort('decDeg').sort('RADeg')
            lastRefreshTime=time.time()

            for obj in objList:
                # Keeping the session alive - refresh every 10 min or so
                if time.time()-lastRefreshTime > 600:
                    self.db.command({"refreshSessions" : [self.mongoSess.session_id]})
                    lastRefreshTime=time.time()

                subDir=str(obj['RADeg']).split(".")[0]
                os.makedirs(outDir+os.path.sep+subDir, exist_ok = True)
                outFileName=outDir+os.path.sep+subDir+os.path.sep+catalogTools.makeRADecString(obj['RADeg'], obj['decDeg'])+".jpg"
                if os.path.exists(outFileName) == False:
                    print("... making image for %s ..." % (obj['name']))
                    for imgFileName in imgList:
                        wcs=wcsDict[imgFileName]
                        useThisImage=False
                        if matchKey != None:
                            if imgFileName.find(obj[matchKey]) != -1:
                                useThisImage=True
                        else:
                            pixCoords=wcs.wcs2pix(obj['RADeg'], obj['decDeg'])
                            if pixCoords[0] >= 0 and pixCoords[0] < wcs.header['NAXIS1'] and pixCoords[1] >= 0 and pixCoords[1] < wcs.header['NAXIS2']:
                                useThisImage=True
                        if useThisImage == False:
                            continue

                        # Clip image
                        with pyfits.open(imgFileName) as img:
                            for ext in img:
                                if ext.data is not None:
                                    break
                            data=ext.data
                        clip=astImages.clipImageSectionWCS(data, wcs, obj['RADeg'], obj['decDeg'],
                                                           maxSizeArcmin/60.0)
                        
                        # Optional check - did we pick an image where object is in zeroed border?
                        pixCoords=clip['wcs'].wcs2pix(obj['RADeg'], obj['decDeg'])
                        if 'checkAllZeros' in imDirDict.keys() and imDirDict['checkAllZeros'] == True:
                            if clip['data'][int(pixCoords[1]), int(pixCoords[0])] == 0:
                                continue
                            if clip['data'].sum() == 0:
                                continue
                        
                        # Save .fits if needed
                        if 'contourImage' in self.configDict.keys() and self.configDict['contourImage'] == label:
                            fitsOutFileName=outFileName.replace(".jpg", ".fits")
                            astImages.saveFITS(fitsOutFileName, clip['data'], clip['wcs'])

                        # Save .jpg                            
                        # Try to pick sensible cut levels
                        # Min-Max scaling
                        # Should probably stick with this, but also add log option for optical
                        if scaling == 'auto' and minMaxRadiusArcmin is not None:
                            rMap=np.zeros(clip['data'].shape, dtype = float)
                            rMap, blah1, blah2=makeDegreesDistanceMap(rMap, clip['wcs'], obj['RADeg'], obj['decDeg'], 100.0)
                            minMaxData=clip['data'][np.less(rMap, minMaxRadiusArcmin/60.0)]
                            cuts=[clip['data'].min(), clip['data'].max()]
                        elif scaling == 'log':
                            logClip=np.array(clip['data'])
                            logClip[logClip > 0]=np.log10(logClip[logClip > 0])
                            logCutMin=logClip[clip['data'] < imDirDict['cuts'][0]].max()
                            logCutMax=logClip[clip['data'] > imDirDict['cuts'][1]].min()
                            logClip=(logClip-logCutMin)/logCutMax
                            logClip[clip['data'] < imDirDict['cuts'][0]]=0.0
                            logClip[clip['data'] > imDirDict['cuts'][1]]=1.0
                            clip['data']=logClip
                            cuts=[0, 1]
                        else:
                            scaleMin, scaleMax=scaling.split(":")
                            scaleMin=float(scaleMin)
                            scaleMax=float(scaleMax)
                            cuts=[scaleMin, scaleMax]
                        
                        dpi=96.0
                        f=plt.figure(figsize=(sizePix/dpi, sizePix/dpi), dpi = dpi)
                        plt.axes([0, 0, 1, 1])
                        plt.imshow(clip['data'], interpolation = "none", origin = 'lower', 
                                    cmap = colorMap, norm = plt.Normalize(cuts[0], cuts[1]))
                        try:
                            plt.savefig(outFileName, dpi = dpi)
                        except:
                            raise Exception("if you see this, you probably need to update PIL/Pillow")
                        plt.close()

        
        
        
        
