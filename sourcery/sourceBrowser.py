"""

    Copyright 2014-2018 Matt Hilton (matt.hilton@mykolab.com)
    
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
from io import BytesIO
from tempfile import TemporaryDirectory
import tempfile
import pymongo
from bson.son import SON
import pyximport; pyximport.install()
import sourceryCython
import cherrypy
import pickle
#import pyvips
import IPython
from sourcery import sourceryAuth
from sourcery import tileDir
from passlib.hash import pbkdf2_sha256

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
            
        # Image choices
        if 'defaultImageType' not in self.configDict.keys():
            self.configDict['defaultImageType']='best'
        if 'imagePrefs' not in self.configDict.keys():
            self.configDict['imagePrefs']=['DES', 'SDSS', 'unWISE']
        
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
        if os.path.exists(sdssCacheDir) == False:
            os.makedirs(sdssCacheDir)
        decalsCacheDir=self.cacheDir+os.path.sep+"DECaLS"
        if os.path.exists(decalsCacheDir) == False:
            os.makedirs(decalsCacheDir)
        ps1CacheDir=self.cacheDir+os.path.sep+"PS1"
        if os.path.exists(ps1CacheDir) == False:
            os.makedirs(ps1CacheDir)
        ps1CacheDir=self.cacheDir+os.path.sep+"PS1IR"
        if os.path.exists(ps1CacheDir) == False:
            os.makedirs(ps1CacheDir)
        wiseCacheDir=self.cacheDir+os.path.sep+"unWISE"
        if os.path.exists(wiseCacheDir) == False:
            os.makedirs(wiseCacheDir)
            
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
        
        # MongoDB set up
        self.dbName=self.configDict['MongoDBName']
        self.client=pymongo.MongoClient('localhost', 27017)
        self.db=self.client[self.dbName]
        self.sourceCollection=self.db['sourceCollection']
        self.sourceCollection.ensure_index([('loc', pymongo.GEOSPHERE)])
        self.fieldTypesCollection=self.db['fieldTypes']
        self.tagsCollection=self.db['tagsCollection']
        self.tagsCollection.ensure_index([('loc', pymongo.GEOSPHERE)])
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
        # NOTE: handling of surveys (e.g., SDSS) is clunky and getting unwieldy...
        # NOTE: made clunkier to allow different size images (imageDirs only for now)
        self.imDirLabelsList=[]
        self.imDirMaxSizeArcminList=[]
        if self.configDict['addSDSSImage'] == True:
            self.imDirLabelsList.append("SDSS")
            self.imDirMaxSizeArcminList.append(self.configDict['plotSizeArcmin'])
        if self.configDict['addDECaLSImage'] == True:
            self.imDirLabelsList.append("DECaLS")
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
        if matches.count() == 0:
            newPost={'loc': {'type': 'Point', 'coordinates': [lon, obj['decDeg']]}}
            self.tagsCollection.insert(newPost)
            mongoDict={}
        else:
            mongoDict=matches.next()
        
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
                sourceryIDs.append(row['sourceList']+"_"+row['name'].replace(" ", "_"))           
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
            tab.add_column(atpy.Column(np.zeros(len(tab)), colLabel))
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
        tab.add_column(atpy.Column(np.zeros(len(tab)), "cacheBuilt"))
        
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
                        tab[key][tab[key].mask]=-99
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
                                tab.add_column(atpy.Column(np.ones(len(tab), dtype = xTab[key].dtype)*-99, key))
                            tab[key][mask]=xTab[key][xIndices[mask]]
                    tab.add_column(atpy.Column(np.zeros(len(tab), dtype = int), '%s_match' % (label)))
                    tab.add_column(atpy.Column(np.ones(len(tab), dtype = float)*-99, '%s_distArcmin' % (label)))
                    tab['%s_match' % (label)][mask]=1
                    tab['%s_distArcmin' % (label)][mask]=rDeg.value[mask]*60.0
            tab.remove_column("matchIndices")
        assert(len(tab) == origLen)
        
        # Cache the result of the cross matches: we need this for speed later on when downloading catalogs
        # Otherwise, for large catalogs, we're hitting memory issues
        cachedTabFileName=self.cacheDir+os.path.sep+"%s_xMatchedTable.fits" % (self.configDict['catalogDownloadFileName'])
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
            
            #print "... adding %s to database (%d/%d) ..." % (row['name'], idCount, len(tab))
            
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
                    newPost[key]=bool(row[key])
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
            for key in tagsDict:
                newPost[key]=tagsDict[key]
                
            postsList.append(newPost)
        
        # Insert all posts at once
        self.sourceCollection.insert_many(postsList)

        # Add descriptions of field (displayed on help page only)
        descriptionsDict=self.parseColumnDescriptionsFile()
        
        # Make collection of field types
        index=0
        for key in fieldTypesList:
            fieldDict={}
            fieldDict['name']=key
            fieldDict['type']=fieldTypesDict[key]
            fieldDict['index']=index
            if key in descriptionsDict.keys():
                fieldDict['description']=descriptionsDict[key]
            else:
                fieldDict['description']="-"
            self.fieldTypesCollection.insert(fieldDict)
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
            keysToFix=["userListFile", "cacheDir", "skyviewCacheDir", "newsFileName", "crossMatchCatalogs", "tileDirs"]
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
        
        outFileName=ps1CacheDir+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".jpg"
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
                          
        outFileName=sdssCacheDir+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".jpg"
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


    def fetchDECaLSImage(self, name, RADeg, decDeg, refetch = False):
        """Fetches DECaLS .jpg cut-out. Unfortunatly at the moment, these are limited to 512 pixels
        maximum at the moment (DR7).
        
        """
    
        decalsCacheDir=self.cacheDir+os.path.sep+"DECaLS"
                          
        outFileName=decalsCacheDir+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".jpg"
        decalsWidth=512 # Max set by DECaLS server
        decalsPixScale=(self.configDict['plotSizeArcmin']*60.0)/float(decalsWidth)
        if os.path.exists(outFileName) == False or refetch == True:
            urlString="http://legacysurvey.org/viewer/jpeg-cutout?ra=%.6f&dec=%.6f&size=%d&layer=decals-dr7&pixscale=%.4f&bands=grz" % (RADeg, decDeg, decalsWidth, decalsPixScale)
            resp=self.http.request('GET', urlString)
            with open(outFileName, 'wb') as f:
                f.write(resp.data)
                f.close()
                

    def fetchUnWISEImage(self, name, RADeg, decDeg, refetch = False):
        """Retrieves unWISE W1, W2 .fits images and makes a colour .jpg.
        
        """
        
        wiseCacheDir=self.cacheDir+os.path.sep+"unWISE"
                
        # 2.75" pixels in the unWISE images (max for query is 250 pixels though)
        sizePix=int(round(self.configDict['plotSizeArcmin']*60.0/2.75))  
        if sizePix > 250:
            raise Exception("WISE web service does not work for images > 250 WISE pixels - use tileDir set-up instead or use smaller plotSizeArcmin")
        
        outFileName=wiseCacheDir+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".jpg"
        if os.path.exists(outFileName) == False or refetch == True:
            with TemporaryDirectory() as tmpDirPath:
                targzPath=tmpDirPath+os.path.sep+"wise.tar.gz"            
                print("... fetching unWISE data for %s ..." % (name))
                resp=self.http.request('GET', "http://unwise.me/cutout_fits?version=neo1&ra=%.6f&dec=%.6f&size=%d&bands=12" % (RADeg, decDeg, sizePix))
                with open(targzPath, 'wb') as f:
                    f.write(resp.data)
                    f.close()
                with tarfile.open(targzPath) as t:
                    t.extractall(path = os.path.split(targzPath)[0])
                    t.close()
                wiseFiles=glob.glob(tmpDirPath+os.path.sep+"unwise-*-img-m.fits")
                w1FileName=None
                w2FileName=None
                for w in wiseFiles:
                    # Weirdly, the archives sometimes have multiple images, some of odd dimensions
                    if w.find("-img-m.fits") != -1:
                        img=pyfits.open(w)
                        if img[0].data.shape == (sizePix, sizePix):
                            if w.find("-w1-img") != -1:
                                w1FileName=w
                            if w.find("-w2-img") != -1:
                                w2FileName=w
                if w1FileName == None or w2FileName == None:
                    # In this case, we have to stitch all the images together
                    # Takes ~1.6 sec per image
                    for band in ['w1', 'w2']:
                        # Make a WCS
                        CRVAL1, CRVAL2=RADeg, decDeg
                        sizeArcmin=self.configDict['plotSizeArcmin']
                        xSizeDeg, ySizeDeg=sizeArcmin/60.0, sizeArcmin/60.0
                        xSizePix=sizePix
                        ySizePix=sizePix
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
                        newHead['CDELT1']=xOutPixScale
                        newHead['CDELT2']=xOutPixScale    # Makes more sense to use same pix scale
                        newHead['CUNIT1']='DEG'
                        newHead['CUNIT2']='DEG'
                        wcs=astWCS.WCS(newHead, mode='pyfits')
                        outData=np.zeros([sizePix, sizePix])
                        imgFileNames=glob.glob(("*-%s-img-m.fits" % (band)))
                        for fileName in imgFileNames:
                            img=pyfits.open(fileName)
                            imgWCS=astWCS.WCS(img[0].header, mode = 'pyfits')
                            for y in range(sizePix):
                                for x in range(sizePix):
                                    outRADeg, outDecDeg=wcs.pix2wcs(x, y)
                                    inX, inY=imgWCS.wcs2pix(outRADeg, outDecDeg)
                                    # Once, this returned infinity...
                                    try:
                                        inX=int(round(inX))
                                        inY=int(round(inY))
                                    except:
                                        continue
                                    if inX >= 0 and inX < img[0].data.shape[1]-1 and inY >= 0 and inY < img[0].data.shape[0]-1:
                                        outData[y, x]=img[0].data[inY, inX]
                        if band == 'w1':
                            bClip={'wcs': wcs, 'data': outData}
                        elif band == 'w2':
                            rClip={'wcs': wcs, 'data': outData}
                else:
                    wcs=astWCS.WCS(w1FileName)
                    with pyfits.open(w1FileName) as bImg:
                        bClip={'wcs': wcs, 'data': bImg[0].data}
                    with pyfits.open(w2FileName) as rImg:
                        rClip={'wcs': wcs, 'data': rImg[0].data}
            
                try:
                    gClip={'wcs': rClip['wcs'], 'data': (rClip['data']+bClip['data'])/2.0}
                except:
                    raise Exception("W1, W2 images not same dimensions")

            ## Clean up
            #fileList=os.listdir(tmpDirPath)
            #for f in fileList:
                #os.remove(tmpDirPath+os.path.sep+f)
            #os.removedirs(tmpDirPath)
            
            # Make colour .jpg
            # Nicer log scaling - twiddle with the min, max levels here and cuts below as you like
            dpi=96.0
            bData=bClip['data']
            gData=gClip['data']
            rData=rClip['data']
            rData[np.less(rData, 1e-5)]=1e-5
            rData[np.greater(rData, 1000)]=1000.0
            rData=np.log10(rData)
            gData[np.less(gData, 1e-5)]=1e-5
            gData[np.greater(gData, 1000)]=1000.0
            gData=np.log10(gData)
            bData[np.less(bData, 1e-5)]=1e-5
            bData[np.greater(bData, 1000)]=1000.0
            bData=np.log10(bData)

            cuts=[0, 3]

            sizePix=1024
            f=plt.figure(figsize=(sizePix/dpi, sizePix/dpi), dpi = dpi)
            axes=[0., 0., 1.0, 1.0]
            plot=astPlots.ImagePlot([rData, gData, bData], bClip['wcs'], axes = axes, 
                                cutLevels = [cuts, cuts, cuts], axesFontSize = 18.0)
            plt.savefig(outFileName, dpi = dpi)
            plt.close()

                
    @cherrypy.expose
    def makePlotFromJPEG(self, name, RADeg, decDeg, surveyLabel, plotNEDObjects = "false", plotSpecObjects = "false", 
                         plotSourcePos = "false", plotXMatch = "false", plotContours = "false", showAxes = "false", clipSizeArcmin = None, gamma = 1.0):
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
        inJPGPath=self.cacheDir+os.path.sep+surveyLabel+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".jpg"
        if os.path.exists(inJPGPath) == False:
            # Testing fall back option of live fetch from legacysurvey.org
            fetchWidth=512 # Max set by DECaLS server
            fetchPixScale=(self.configDict['plotSizeArcmin']*60.0)/float(fetchWidth)
            if surveyLabel == 'SDSS':
                layer='sdss'
            elif surveyLabel == 'DECaLS':
                layer='decals-dr7'
            elif surveyLabel == 'unWISE':
                layer='unwise-w1w2'
            elif surveyLabel == 'DES':
                layer='des-dr1'
            else:
                layer=None 
            if layer is not None:
                urlString="http://legacysurvey.org/viewer/jpeg-cutout?ra=%.6f&dec=%.6f&size=%d&layer=%s&pixscale=%.4f&bands=grz" % (RADeg, decDeg, fetchWidth, layer, fetchPixScale)
                resp=self.http.request('GET', urlString)
                #with open('test.jpg', 'wb') as f:
                    #f.write(resp.data)
                im=Image.open(io.BytesIO(resp.data))
                data=np.array(im)
                data=np.power(data, 1.0/float(gamma))
                if data.shape != (fetchWidth, fetchWidth, 3):
                    layer=None  # We get (256, 256) if not in footprint (I think)
            if layer is None:
                inJPGPath=sourcery.__path__[0]+os.path.sep+"data"+os.path.sep+"noData.jpg"
            else:
                inJPGPath=None # success in this case
        
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

        cutLevels=[[R.min(), R.max()], [G.min(), G.max()], [B.min(), B.max()]]
        
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
            p.addScaleBar('NW', scaleBarSizeArcmin*60.0, color='yellow', fontSize=16, width=2.0, label = "1'")
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
                clipFileName=self.cacheDir+os.path.sep+self.configDict['contourImage']+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".fits"
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
        html="""<html><body style="font-family: sans-serif; vertical align: top; justify: full;">
        <table cellpadding="4" cellspacing="0" border="0" style="text-align: left; width: 100%;">
            <tbody>
                <tr>
                    <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                        text-align: center; vertical-align: middle; font-size: 150%;">
                        <b>$TITLE</b>
                    </td>
                </tr>
            </tbody>
        </table>
        <br>
        <fieldset>
        <legend><b>Enter Login Information</b></legend>
        <p>
        $LOGIN_MESSAGE
            <form method="post" action="$SCRIPT_NAME/login">
            <input type="hidden" name="from_page" value="$FROM_PAGE" />
            <label for="username"><b>Username:</b></label>
            <input type="text" name="username" value="$USERNAME" /><br />
            <br>
            <label for="username"><b>Password:</b></label>            
            <input type="password" name="password" /><br />
            <br>
            <input type="submit" style="font-size: 1.05em;" value="Log in" />
        </p>
        <p>$MSG</p>
        </fieldset>
        </body></html>"""
        html=html.replace("$SCRIPT_NAME", cherrypy.request.script_name)
        html=html.replace("$FROM_PAGE", from_page)
        html=html.replace("$MSG", msg)
        html=html.replace("$USERNAME", username)
        if 'indexTitle' in self.configDict.keys():
            html=html.replace("$TITLE", self.configDict['indexTitle'])
        else:
            html=html.replace("$TITLE", "Sourcery Database")
        if 'loginMessage' in self.configDict.keys():
            html=html.replace("$LOGIN_MESSAGE", "<p>%s</p>" % (self.configDict['loginMessage']))
        else:
            html=html.replace("$LOGIN_MESSAGE", "")
            
        return html
    
    
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
        
        templatePage="""<html>
        <head>
            <meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
            <style>
            .f {
                    float: right;
            }
            </style>
            <style>
            .tablefont, .tablefont TD, .tablefont TH {
                font-family: monospace;
            }
            </style>
            <title>$TITLE</title>
        </head>
        
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js"></script>

        <script type="text/javascript">
        
        $(function(){
           // Set cursor to pointer and add click function
           $("legend").css("cursor","pointer").click(function(){
               var legend = $(this);
               var value = $(this).children("span").html();
               if(value == "hide")
                   value="expand";
               else
                   value="hide";
               $(this).siblings().slideToggle(0, function() { legend.children("span").html(value); } );
           });    
        });
        
        $(document).ready(function() {
            var legends = document.getElementsByTagName("legend");
            for (var i=0; i<legends.length; i++)
            {
                var spans = legends[i].getElementsByTagName("span");
                var value=spans[0].innerHTML;
                if (value == "expand") {
                    spans[0].innerHTML="hide";
                    legends[i].click();
                }
            }
        });
            
        </script>
        
        <style>
        #logout a:link, #logout a:visited {
            background-color: white;
            color: black;
            padding: 5px 5px;
            text-align: center;
            text-decoration: none;
            display: inline-block;
        }

        #logout a:hover, #logout a:active {
            background-color: white;
        }
        </style>
        
        <body style="font-family: sans-serif; vertical align: top; justify: full;" onload="hideShow()">
        <table cellpadding="4" cellspacing="0" border="0" style="text-align: left; width: 100%;">
            <tbody>
                <tr>
                    <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                        text-align: center; vertical-align: middle; font-size: 150%;">
                        <b>$TITLE</b>
                    </td>
                    <div id='logout' style="text-align: right; vertical-align: middle; font-size: 100%;  
                        position: absolute; right: 20px; top: 20px;">
                        <a href='$SCRIPT_NAME/logout' color>logout</a>
                    </div>
                </tr>
            </tbody>
        </table>
        
        $META_DATA
        $COLOR_CODING
        
        <br>
        <form method="get" action="updateQueryParams">
        <fieldset>
        $QUICK_QUERY_LINKS
        <legend><span style='border: black 1px solid; color: gray; padding: 2px'>hide</span><b>Constraints</b></legend>
        <p>Enter coordinate ranges (e.g., 120:220) or set the search box length. Use negative RA values to wrap around 0 degrees (e.g., -60:60).</p>
        <p>
        <label for="queryRADeg">RA (degrees)</label>
        <input type="text" value="$QUERY_RADEG" name="queryRADeg"/>
        <label for="queryDecDeg">Dec. (degrees)</label>
        <input type="text" value="$QUERY_DECDEG" name="queryDecDeg"/>
        <label for="querySearchBoxArcmin">Search box length (arcmin)</label>
        <input type="text" value="$QUERY_SEARCHBOXARCMIN" name="querySearchBoxArcmin"/>
        </p>
        </p>
        <label for="queryOtherConstraints">Other constraints <a href=$CONSTRAINTS_HELP_LINK target=new>(help)</a></label>
        <textarea style="width:100%" name="queryOtherConstraints">$QUERY_OTHERCONSTRAINTS</textarea>
        <input type="submit" class="f" style="font-size: 1.05em;" name="queryApply" value="Apply">
        <input type="submit" class="f" style="font-size: 1.05em;" name="queryReset" value="Reset"><br>
        <div style="display: table-cell; vertical-align: middle;"><i>$SHARE_QUERY_LINK</i></div>
        </p>
        </fieldset>
        </form>
            
        <form method="post" action="changeTablePage" style="border">
        <div id="buttons">
            <input type="submit" class="f" style="font-size: 1.05em;" name="nextButton" value=">">
            <input type="submit" class="f" style="font-size: 1.05em;" name="prevButton" value="<">
            <div style="clear:both"></div><!-- Need this to have the buttons actually inside div#buttons -->
        </div>
                
        <table frame=border cellspacing=0 cols=$TABLE_COLS rules=all border=2 width=100% align=center class=tablefont>
        <tbody>
            <tr style="background-color: rgb(0, 0, 0); font-family: monospace; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 125%;">
            $TABLE_COL_NAMES
            </tr>
            <font size="1">
            $TABLE_DATA
            </font>
        </tbody>
        </table>

        <div id="buttons">
            <input type="submit" class="f" style="font-size: 1.05em;" name="nextButton" value=">">
            <input type="submit" class="f" style="font-size: 1.05em;" name="prevButton" value="<">
            <div style="clear:both"></div><!-- Need this to have the buttons actually inside div#buttons -->
        </div>
        </form>
        
        $DOWNLOAD_LINKS

        </tbody>
        </table>
        <hr>
        <a href="https://github.com/mattyowl/sourcery"><i>Sourcery</i></a> - $HOSTED_STR
        <br>
        <br>
        </body>
        </html>
        """
                
        html=templatePage
        html=html.replace("$SCRIPT_NAME", cherrypy.request.script_name)

        if 'indexTitle' in self.configDict.keys():
            html=html.replace("$TITLE", self.configDict['indexTitle'])
        else:
            html=html.replace("$TITLE", "Sourcery Database")
                
        # First need to apply query parameters here
        queryPosts=self.runQuery(queryRADeg, queryDecDeg, querySearchBoxArcmin, queryOtherConstraints)        
        numPosts=queryPosts.count()
        if 'numPosts' not in cherrypy.session:
            cherrypy.session['numPosts']=numPosts
        
        # Then cut to number of rows to view as below
        viewPosts=queryPosts[cherrypy.session['viewTopRow']:cherrypy.session['viewTopRow']+self.tableViewRows]
        
        # Quick query link(s) - at top of 'constraints' box
        if 'quickLinks' in self.configDict.keys():            
            quickLinkStr="<p><i>Quick links:</i> "
            for linkDict in self.configDict['quickLinks']:
                url="updateQueryParams?queryRADeg=0%3A360&queryDecDeg=-90%3A90&querySearchBoxArcmin=&queryOtherConstraints="
                url=url+linkDict['constraints'].replace("+", "%2B") # quick and dirty fix for + in query constraints text
                url=url+"&queryApply=Apply"
                quickLinkStr=quickLinkStr+'<a href="%s">%s</a>' % (url, linkDict['label'])
                quickLinkStr=quickLinkStr+" - "
            quickLinkStr=quickLinkStr[:-3]+"</p>"
        else:
            quickLinkStr=""
        html=html.replace("$QUICK_QUERY_LINKS", quickLinkStr)
            
        # Fill in query params
        html=html.replace("$QUERY_RADEG", queryRADeg)
        html=html.replace("$QUERY_DECDEG", queryDecDeg)
        html=html.replace("$QUERY_SEARCHBOXARCMIN", querySearchBoxArcmin)
        html=html.replace("$QUERY_OTHERCONSTRAINTS", queryOtherConstraints)
        html=html.replace("$OBJECT_TYPE_STRING", self.configDict['objectTypeString'])
        html=html.replace("$NUMBER_SOURCES", str(numPosts))#str(len(queryPosts)))
        if 'hostedBy' in self.configDict.keys():
            html=html.replace("$HOSTED_STR", self.configDict['hostedBy'])
        else:
            html=html.replace("$HOSTER_STR", "")
        html=html.replace("$CONSTRAINTS_HELP_LINK", "displayConstraintsHelp?")
        
        # Shareable query link
        # NOTE: we should do this properly really...
        shareQueryURL=cherrypy.request.base+cherrypy.request.script_name+"/"
        url="updateQueryParams?queryRADeg=%s&queryDecDeg=%s&querySearchBoxArcmin=%s&queryOtherConstraints=" % (self.sourceNameToURL(str(queryRADeg)), 
         self.sourceNameToURL(str(queryDecDeg)), 
         self.sourceNameToURL(str(querySearchBoxArcmin)))
        url=url+self.sourceNameToURL(queryOtherConstraints)+"&queryApply=Apply"
        shareQueryURL=shareQueryURL+url
        html=html.replace("$SHARE_QUERY_LINK", "<a href='%s'>Shareable link for this query</a>" %  (shareQueryURL))
        
        # This would hide censored columns, but still permit queries based on them to be shown
        #columnsHiddenList=[]
        #for userNameKey in self.usersDict.keys():
            #userDict=self.usersDict[userNameKey]
            #if userNameKey == user and 'hiddenXMatchTables' in userDict.keys():
                #for prefix in userDict['hiddenXMatchTables']:
                    #for key in viewPosts[0].keys():
                        #if key.find(prefix) != -1 and key[:len(prefix)] == prefix:
                            #columnsHiddenList.append(key)
        # This is a crude way of breaking queries on censored columns by users who don't have permission to see them...
        user=cherrypy.session['_sourcery_username']
        for userNameKey in self.usersDict.keys():
            userDict=self.usersDict[userNameKey]
            if userNameKey == user and 'hiddenXMatchTables' in userDict.keys():
                for prefix in userDict['hiddenXMatchTables']:
                    if queryOtherConstraints.find(prefix) != -1:
                        numPosts=0
                        viewPosts=[]
      
        # Table columns - as well as defaults, add ones we query on
        columnsShownList=[]
        for colDict in self.tableDisplayColumns:
            if colDict['name'] not in columnsShownList:
                columnsShownList.append(colDict['name'])
        displayColumns=[]+self.tableDisplayColumns
        operators=["<", ">", "=", "!"]
        logicalOps=[' and ', ' or ']
        for logOp in logicalOps:
            constraints=queryOtherConstraints.split(logOp)
            for c in constraints:
                for o in operators:
                    colName=c.split(o)[0].lstrip().rstrip()
                    if numPosts > 0 and colName in viewPosts[0].keys() and colName not in columnsShownList:# and colName not in columnsHiddenList:
                        fieldTypeDict=self.fieldTypesCollection.find_one({'name': colName})
                        dispDict={'name': colName, 'label': colName}
                        if fieldTypeDict['type'] == 'number':
                            dispDict['fmt']='%.3f'
                        elif fieldTypeDict['type'] == 'text':
                            dispDict['fmt']='%s'
                        else:
                            raise Exception("unknown type for field '%s'" % (colName))
                        displayColumns.append(dispDict)
                        
        columnHeadings=""
        for colDict in displayColumns:
            columnHeadings=columnHeadings+"\n           <td><b>%s</b></td>" % (colDict['label'])
        html=html.replace("$TABLE_COL_NAMES", columnHeadings)
        html=html.replace("$TABLE_COLS", str(len(displayColumns)))
        
        # Meta data
        if 'classificationDescription' in self.configDict.keys():
            classificationDesc=self.configDict['classificationDescription']
        else:
            classificationDesc=""
        if 'catalogComments' not in self.configDict.keys():
            commentsString=classificationDesc
        else:
            commentsString=self.configDict['catalogComments']+" <p>"+classificationDesc+"</p>"
        
        if 'newsItems' in self.configDict.keys() and len(self.configDict['newsItems']) > 0:
            latestNewsStr="    &#8211;    Latest news: %s" % (self.configDict['newsItems'][-1].split(":")[0])
        else:
            latestNewsStr=""
        
        # We now display a message if cache rebuild is in progress
        if os.path.exists(self.cacheLockFileName) == True:
            cacheRebuildStr="    &#8211;    [REBUILDING IMAGE CACHE]"
        else:
            cacheRebuildStr=""
        
        # For users with hiddenConstraints set
        hiddenConstraintsMessage=""
        user=cherrypy.session['_sourcery_username']
        for userNameKey in self.usersDict.keys():
            userDict=self.usersDict[userNameKey]
            if userNameKey == user and 'hiddenConstraintsMessage' in userDict.keys():
                hiddenConstraintsMessage="; %s" % (userDict['hiddenConstraintsMessage'])
        
        metaData="""<br><fieldset>
        <legend><span style='border: black 1px solid; color: gray; padding: 2px'>hide</span><b>Source List Information</b></legend>
        Total number of %s: %d (original source list: %d%s) %s %s
        <p>%s</p>
        $NEWS
        </fieldset>""" % (self.configDict['objectTypeString'], numPosts, self.sourceCollection.count(), hiddenConstraintsMessage, latestNewsStr,
                          cacheRebuildStr, commentsString)
        if 'newsItems' in self.configDict.keys():
            newsStr="<p>News:<ul>\n"
            for item in self.configDict['newsItems'][-5:]:
                newsStr=newsStr+"<li>%s</li>\n" % (item)
            newsStr=newsStr+"</ul></p>\n"
            metaData=metaData.replace("$NEWS", newsStr)
        else:
            metaData=metaData.replace("$NEWS", "")
            
        html=html.replace("$META_DATA", metaData)        
        
        # Catalog download links
        #http://localhost:8080/downloadCatalog?queryRADeg=0%3A360&queryDecDeg=-90%3A90&querySearchBoxArcmin=&queryOtherConstraints=softCts+%3E+300
        shortCatalogName=self.configDict['catalogDownloadFileName']+".cat"
        shortFITSName=shortCatalogName.replace(".cat", ".fits")
        shortRegName=shortCatalogName.replace(".cat", ".reg")
        downloadLinkStr="downloadCatalog?queryRADeg=%s&queryDecDeg=%s&querySearchBoxArcmin=%s&queryOtherConstraints=%s&" % (queryRADeg, queryDecDeg, querySearchBoxArcmin, queryOtherConstraints)
        downloadLinkStr=quote_plus(downloadLinkStr, safe='&?=')
        downloadLinks="""<fieldset>
        <legend><span style='border: black 1px solid; color: gray; padding: 2px'>hide</span><b>Download Catalog</b></legend>
        <ul>
        <li><a href=%sfileFormat=cat>%s</a>   (plain text)</li>
        <li><a href=%sfileFormat=fits>%s</a>   (FITS table format)</li>
        <li><a href=%sfileFormat=reg>%s</a>   (DS9 region file)</li></ul>
        <p>Note that current constraints are applied to downloaded catalogs.</p>
        </fieldset><br>
        """ % (downloadLinkStr, shortCatalogName, downloadLinkStr, shortFITSName, downloadLinkStr, shortRegName)
        html=html.replace("$DOWNLOAD_LINKS", downloadLinks)
        
        tableData=""
        usedBckColors=[]
        usedBckKeys=[]

        for obj in viewPosts:
            
            bckColor="white"
                
            # Row for each object in table
            rowString="<tr>\n"
            for colDict in displayColumns:
                # Python 2/3 thing
                try:
                    htmlKey="$"+string.upper(colDict['name'])+"_KEY"
                except:
                    htmlKey="$"+colDict['name'].upper()+"_KEY"                    
                if 'tableAlign' in colDict.keys():
                    alignStr=colDict['tableAlign']
                else:
                    alignStr="center"
                if 'displaySize' in colDict.keys():
                    widthStr='width: %dem;' % (colDict['displaySize'])
                    useDiv=True
                else:
                    useDiv=False
                    if colDict['fmt'] == '%s':
                        if colDict['name'] in obj.keys():
                            widthStr='width: %dem;' % (len(obj[colDict['name']]))
                        else:
                            widthStr=""
                    else:
                        widthStr='width: %dem;' % (len(colDict['fmt'] % (1)))
                rowString=rowString+"   <td style='background-color: "+bckColor+"; "+widthStr+"' align="+alignStr+">"
                if useDiv == True:
                    rowString=rowString+'<div style="text-overflow: ellipsis; white-space: nowrap; overflow: hidden; '+widthStr+'">'
                rowString=rowString+htmlKey
                if useDiv == True:
                    rowString=rowString+"</div>\n"
                rowString=rowString+"</td>\n"
            rowString=rowString+"</tr>\n"
            
            # Insert values - note name is special
            for colDict in displayColumns:
                key=colDict['name']
                if key in obj.keys():
                    try:
                        value=obj[key]
                    except:
                        raise Exception("missing key %s" % (key))
                else:
                    # No entry in MongoDB tags yet
                    if colDict['fmt'] != "%s":
                        value=0.0
                    else:
                        value=""
                # Python 2/3 thing
                try:
                    htmlKey="$"+string.upper(key)+"_KEY"
                except:
                    htmlKey="$"+key.upper()+"_KEY"
                if key == "name":
                    #linksDir="dummy"
                    linkURL="displaySourcePage?sourceryID=%s&clipSizeArcmin=%.2f" % (self.sourceNameToURL(obj['sourceryID']), self.configDict['defaultViewSizeArcmin'])
                    if 'defaultImageType' in self.configDict.keys():
                        linkURL=linkURL+"&imageType=%s" % (self.configDict['defaultImageType'])
                    nameLink="<a href=\"%s\" target=new>%s</a>" % \
                        (linkURL, obj['name'])
                    rowString=rowString.replace(htmlKey, "%s" % (nameLink))
                elif key == "NED_name" and obj['NED_name'] != "None":
                    nedName=obj[key]
                    nedLinkURL="http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES" % (nedName.replace("+", "%2B").replace(" ", "+"))
                    rowString=rowString.replace(htmlKey, "<a href=%s>%s</a>" % (nedLinkURL, nedName))
                elif key == "NED_z" and obj['NED_z'] == "nan":
                    rowString=rowString.replace(htmlKey, "-")
                else:
                    try:
                        rowString=rowString.replace(htmlKey, colDict['fmt'] % (value))
                    except:
                        IPython.embed()
                        sys.exit()
                        raise Exception("""IndexError: check .config file tableDisplayColumns are actually in the .fits table, or for mixed '' "" inside [] """ )
                           
            tableData=tableData+rowString
            
        html=html.replace("$TABLE_DATA", tableData)
        
        # Colour coding table key
        if len(usedBckColors) > 0:
            colorCoding="""<table frame=border cellspacing=0 cols=2 rules=all border=2 width=60% align=center>
                        <tbody>
                            <tr>
                                <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                                    text-align: center; vertical-align: middle; font-size: 110%;" colspan=$NUM_COLORS>Color coding</td>
                            </tr>
                            <tr style="background-color: rgb(0, 0, 0); font-family: sans-serif; 
                                    text-align: center; vertical-align: middle; font-size: 110%;">
                            $COLOR_KEY
                            </tr>
                        </tbody>
                        </table>
                        <br><br>
                        """
            colorCoding=colorCoding.replace("$NUM_COLORS", str(len(usedBckColors)))
            keyString=""
            for bckColor, bckKey in zip(usedBckColors, usedBckKeys):
                keyString=keyString+'<td style="background-color:'+bckColor+';" width='+str(100.0/len(usedBckColors))+'%>'+bckKey+'</td>\n'
            colorCoding=colorCoding.replace("$COLOR_KEY", keyString)
            html=html.replace("$COLOR_CODING", colorCoding)
        else:
            html=html.replace("$COLOR_CODING", "")
        
        return html    


    @cherrypy.expose
    @sourceryAuth.require()
    def downloadCatalog(self, queryRADeg = "0:360", queryDecDeg = "-90:90", querySearchBoxArcmin = "",
                        queryOtherConstraints = "", fileFormat = "cat"):
        """Provide user with the current table view as a downloadable catalog.
        
        """
                
        # Fetch the cached table and update that with any changed classifications info
        cachedTabFileName=self.cacheDir+os.path.sep+"%s_xMatchedTable.fits" % (self.configDict['catalogDownloadFileName'])
        xTab=atpy.Table().read(cachedTabFileName)
                
        posts=list(self.runQuery(queryRADeg, queryDecDeg, querySearchBoxArcmin, queryOtherConstraints))
        tabLength=len(posts)#posts.count()
        
        # Trim xTab to match posts based on sourceryID
        keepIndices=[]
        for p in posts:
            if p['sourceryID'] in xTab['sourceryID']:
                keepIndices.append(np.where(xTab['sourceryID'] == p['sourceryID'])[0][0])
        xTab=xTab[keepIndices]

        # Need this and change below for overloading with e.g. editable BCG coords
        editableFieldsList=[]
        for f in self.configDict['fields']:
            editableFieldsList.append(f['name'])
        for key in xTab.keys():
            if key in editableFieldsList:
                xTab.remove_column(key)
        
        # NOTE: there may be fun unicode-related stuff here: e.g., u'BCG_RADeg' versus 'BCG_RADeg'
        keysList, typeNamesList, descriptionsList=self.getFieldNamesAndTypes(excludeKeys = [])
        keysToAdd=['sourceryID', 'RADeg', 'decDeg']
        typeNamesToAdd=['text', 'number', 'number']
        for k, t in zip(keysList, typeNamesList):
            if k not in xTab.keys() or k in editableFieldsList:
                if k not in keysToAdd:
                    keysToAdd.append(k)
                    typeNamesToAdd.append(t)
                
        tab=atpy.Table()
        tab.table_name=self.configDict['catalogDownloadFileName']
        for key, typeName in zip(keysToAdd, typeNamesToAdd):
            if typeName == 'number':
                tab.add_column(atpy.Column(np.zeros(tabLength, dtype = float), str(key)))
            else:
                tab.add_column(atpy.Column(np.zeros(tabLength, dtype = 'S1000'), str(key)))
        
        count=0
        for post in posts:
            for key in keysToAdd:
                if key in post.keys():           # NOTE: this handles image_ tags, which are 1 if present, and absent otherwise
                    tab[key][count]=post[key]
            count=count+1
        tab.rename_column('RADeg', 'tag_RADeg')
        tab.rename_column('decDeg', 'tag_decDeg')
        
        tab=atpy.join(xTab, tab, keys = ['sourceryID'])
        tab.remove_columns(['tag_RADeg', 'tag_decDeg', 'sourceryID', 'cacheBuilt'])
        
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
        
        tmpFile, tmpFileName=tempfile.mkstemp()
        if fileFormat == 'cat':
            tab.write(tmpFileName+".cat", format = 'ascii')
        elif fileFormat == 'fits':
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
        imgPath=self.cacheDir+os.path.sep+self.configDict['downloadableFITS']+os.path.sep+catalogTools.makeRADecString(obj['RADeg'], obj['decDeg'])+".fits"
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
            self.sourceCollection.ensure_index([("RADeg", pymongo.ASCENDING)])
            #queryPosts=list(self.sourceCollection.find(queryDict).sort('decDeg').sort('RADeg'))        
            queryPosts=self.sourceCollection.find(queryDict).sort('decDeg').sort('RADeg')  
        elif collection == 'tags':
            self.tagsCollection.ensure_index([("RADeg", pymongo.ASCENDING)])
            #queryPosts=list(self.sourceCollection.find(queryDict).sort('decDeg').sort('RADeg')) 
            queryPosts=self.sourceCollection.find(queryDict).sort('decDeg').sort('RADeg')
        else:
            raise Exception("collection should be 'source' or 'tags' only")
                        
        # If we wanted to store all this in its own collection
        #self.makeSessionCollection(queryPosts)
                
        return queryPosts
        

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
            self.db.collection[cherrypy.session.id].insert(q)

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
                
        templatePage="""<html>
        <head>
            <meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
            <title>Constraints Help</title>
        </head>
        <body style="font-family: sans-serif; vertical align: top; justify: full;">
        <table cellpadding="4" cellspacing="0" border="0" style="text-align: left; width: 100%;">
            <tbody>
                <tr>
                    <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                        text-align: center; vertical-align: middle; font-size: 125%;">
                        Constraints Help
                    </td>
                </tr>
            </tbody>
        </table>
        
        <br>
        <p>
        Constraints can be placed on the source list columns listed below. Comparison operators which are understood are <, >, >=, <=, =, !=. Logical operators which are understood are 'and', 'or'.</p>
        <p>Each constraint should be
        separated by 'and' or 'or', e.g.,</p>
        <p><tt>redshift >= 0 and redshift < 0.4</tt></p>
        <p><tt>RM_match = 1 or PSZ2_match = 1</tt></p>
        <p><b>Note that when querying for multiple values in the same column using '=' or '!=', 'and' acts like a delimiter, rather than in a strictly logical sense</b>. For example, to fetch all objects with classification of 'cluster' and 'not cluster', one can write</p>
        <tt>classification = 'cluster' and classification = 'not cluster'</tt> 
        <p>This will leave out all table rows which have classification set to some other value (e.g., 'probable cluster'
        or 'possible cluster'). This query is also equivalent to</p>
        <tt>classification = 'cluster' or classification = 'not cluster'</tt>
        <p>The wildcard '*' is supported in text searches, e.g.,</p> 
        <tt>classification = '* cluster'</tt> 
        <p>will return all objects flagged as 'probable cluster', 'possible cluster', or 'not cluster', but not objects with classification = 'cluster'.</p>
        <tt>notes = '*high-z*'</tt> 
        <p>will return all objects where the string 'high-z' appears in the notes field somewhere. <b>Note that wildcard text searches are case insensitive</b>.
        </p>
        <br>
        <table frame=border cellspacing=0 cols=3 rules=all border=2 width=80% align=center>
        <tbody>
            <tr style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 110%;">
            <td><b>Column</td>
            <td><b>Type</b></td>
            <td><b>Description</b></td>
            </tr>
            $TABLE_DATA
        </tbody>
        </table>        

        <hr>
        <i>Sourcery</i> - $HOSTED_STR
        <br>
        <br>
        </body>
        </html>
        """
        html=templatePage

        html=html.replace("$HOSTED_STR", self.configDict['hostedBy'])        
        
        # Fill in table of data types
        bckColor="white"
        tableData=""
        excludeKeys=['RADeg', 'decDeg', 'sourceryID', 'cacheBuilt'] # because we handle differently
        keysList, typeNamesList, descriptionsList=self.getFieldNamesAndTypes(excludeKeys = excludeKeys)
        for key, typeName, description in zip(keysList, typeNamesList, descriptionsList):
            # Row for each column in table
            rowString="<tr>\n"                
            rowString=rowString+"   <td style='background-color: "+bckColor+";' align=left width=10%><b>"+key+"</b></td>\n"
            rowString=rowString+"   <td style='background-color: "+bckColor+";' align=left width=10%>"+typeName+"</td>\n"
            rowString=rowString+"   <td style='background-color: "+bckColor+";' align=left width=80%>"+description+"</td>\n"
            rowString=rowString+"</tr>\n"                           
            tableData=tableData+rowString
        html=html.replace("$TABLE_DATA", tableData)
        
        return html   


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
        
        return keysList, typeNamesList, descList
                
    
    @cherrypy.expose
    def updateMongoDB(self, name, returnURL, **kwargs):
        """Update info on source in MongoDB.
        
        """
        
        if not cherrypy.session.loaded: cherrypy.session.load()
        
        # To start with, match ONLY on object name... this assumes that objects with the same name will have sufficiently
        # similar coordinates (we handle updating for objects in multiple source lists below)
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
        self.tagsCollection.update({'_id': mongoDict['_id']}, {'$set': post}, upsert = False)
        
        # Update source collection too - here we will do this for all sources that share the same name (we can have multiple source lists)
        # This could be done using cross matching based on coords instead, but this could cause confusion in the case of multiple sources
        # that are not the same object, located at separations < MongoDBCrossMatchRadiusArcmin
        # Example of this: erroneously deblended objects in the XCS list, where you only want to flag some objects as junk, and others not
        objs=self.sourceCollection.find({'name': name})
        for obj in objs:
            self.sourceCollection.update({'_id': obj['_id']}, {'$set': post}, upsert = False)
        
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
        
        self.tagsCollection.update({'_id': mongoDict['_id']}, {'$set': post}, upsert = False)
        
        # Update source collection too
        self.sourceCollection.update({'_id': obj['_id']}, {'$set': post}, upsert = False)
        
    
    @cherrypy.expose
    @sourceryAuth.require()
    def displaySourcePage(self, sourceryID, imageType = 'best', clipSizeArcmin = None, plotNEDObjects = "false", plotSpecObjects = "false", plotSourcePos = "false", plotXMatch = "false", plotContours = "false", showAxes = "false", gamma = 1.0):
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
                
        templatePage="""<html>
        <head>
            <meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
            <title>$SOURCE_NAME</title>
        </head>
        <body style="font-family: sans-serif; vertical align: top; justify: full;">
        <table cellpadding="4" cellspacing="0" border="0" style="text-align: left; width: 100%;">
            <tbody>
                <tr>
                    <!-- $PREV_LINK_CODE -->
                    <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                        text-align: center; vertical-align: middle; font-size: 150%;">
                        <b>$SOURCE_NAME</b>
                    </td>
                    <!-- $NEXT_LINK_CODE -->
                </tr>
            </tbody>
        </table>
        
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/2.1.1/jquery.min.js"></script>

        <br>

        <div style="display: table; width: 85%; margin:0 auto;">
            <div style="display: table-row; height: auto; margin:0 auto;">
                <div style="display: table-cell;">
                    <div style="display: block;">
                        $PLOT_CONTROLS
                    </div>
                    <div style="display: block;">
                        $TAG_CONTROLS
                    </div>
                </div>
                <div style="display: table-cell;">
                    <fieldset style="height: 100%;">
                    <legend><b>Source Image</b></legend>
                    <div id="imagePlot"></div>
                    </fieldset>
                </div>
            </div>
        </div>

        <table frame=border cellspacing=0 cols=1 rules=None border=0 width=100%>
        <tbody>        

        $SPECTRUM_PLOT

        <tr>
            <td align=center>$NED_MATCHES_TABLE</td>
        </tr>

        <tr>
            <td align=center>$SPEC_MATCHES_TABLE</td>
        </tr>
        
        <tr>
            <td align=center>$PROPERTIES_TABLE</td>
        </tr>
        
        </tbody>
        </table>
        <br>
        <hr>
        <i>Sourcery</i> - $HOSTED_STR
        <br>
        <br>
        </body>
        </html>
        """
        templatePage=templatePage.replace("$PLOT_DISPLAY_WIDTH_PIX", str(self.configDict['plotDisplayWidthPix']))

        # Taken this out from above the caption line
        #<tr><td align=center><b>$SIZE_ARC_MIN' x $SIZE_ARC_MIN'</b></td></tr>

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
            # Fall back option
            if imageType == 'best':
                imageType='SDSS'
        
        # Controls for image zoom, plotting NED, SDSS, etc.       
        plotFormCode="""
        <script>
        function printValue(sliderID, textbox) {
            var x = document.getElementById(textbox);
            var y = document.getElementById(sliderID);
            x.value = y.value;
        }
        function updateSlider(value) {
            slider = document.getElementById('sizeSlider');
            slider.value = value;
            printValue('sizeSlider','sizeSliderValue');
        }
        function parseImageTypeValue(imageType, index) {
            var res = imageType.split(";")[index];
            return res;
        }
        window.onload = function() { printValue('sizeSlider', 'sizeSliderValue'); printValue('gammaSlider', 'gammaSliderValue');}
        </script>

        <script type="text/javascript">
        
            $(document).ready(function() {
                    $.post('makePlotFromJPEG', 
                           {name: '$OBJECT_NAME',
                            RADeg: $OBJECT_RADEG,
                            decDeg: $OBJECT_DECDEG,
                            surveyLabel: parseImageTypeValue($('input:radio[name=imageType]:checked').val(), 0),
                            plotNEDObjects: $('input:checkbox[name=plotNEDObjects]').prop('checked'),
                            plotSpecObjects: $('input:checkbox[name=plotSpecObjects]').prop('checked'),
                            plotSourcePos: $('input:checkbox[name=plotSourcePos]').prop('checked'),
                            plotXMatch: $('input:checkbox[name=plotXMatch]').prop('checked'),
                            plotContours: $('input:checkbox[name=plotContours]').prop('checked'),
                            showAxes: $('input:checkbox[name=showAxes]').prop('checked'),
                            clipSizeArcmin: $DEFAULT_CLIP_SIZE_ARCMIN,
                            gamma: $("#gamma").val()}, 
                            function(data) {
                                // directly insert the image
                                $("#imagePlot").html('<img src="data:image/jpg;base64,' + data + '" align="middle" border=0 width="$PLOT_DISPLAY_WIDTH_PIX"/>') ;
                           });
                    return false;
            });
            
            $(document).ready(function() {
                $(':checkbox').change(function(){
                    $( "#imageForm" ).submit();
                });
                $('input:radio[name=imageType]').change(function(){
                    $( "#imageForm" ).submit();
                });
            });
            
            $CLICKABLE_RADEC_SCRIPT
            
            $(function() {
                // When the form is submitted...
                $("#imageForm").submit(function() {   
                    var maxSizeArcmin = parseImageTypeValue($('input:radio[name=imageType]:checked').val(), 1);
                    if (maxSizeArcmin > $DEFAULT_CLIP_SIZE_ARCMIN) {
                        updateSlider(parseImageTypeValue($('input:radio[name=imageType]:checked').val(), 1));
                    }
                    if (parseFloat($("#sizeSliderValue").val()) > maxSizeArcmin) {
                        updateSlider(maxSizeArcmin);
                    }
                    $.post('makePlotFromJPEG', 
                           {name: '$OBJECT_NAME',
                            RADeg: $OBJECT_RADEG,
                            decDeg: $OBJECT_DECDEG,
                            surveyLabel: parseImageTypeValue($('input:radio[name=imageType]:checked').val(), 0),
                            plotNEDObjects: $('input:checkbox[name=plotNEDObjects]').prop('checked'),
                            plotSpecObjects: $('input:checkbox[name=plotSpecObjects]').prop('checked'),
                            plotSourcePos: $('input:checkbox[name=plotSourcePos]').prop('checked'),
                            plotXMatch: $('input:checkbox[name=plotXMatch]').prop('checked'),
                            plotContours: $('input:checkbox[name=plotContours]').prop('checked'),
                            showAxes: $('input:checkbox[name=showAxes]').prop('checked'),
                            clipSizeArcmin: $("#sizeSliderValue").val(),
                            gamma: $("#gammaSliderValue").val()}, 
                            function(data) {
                                // directly insert the image
                                $("#imagePlot").html('<img src="data:image/jpg;base64,' + data + '" align="middle" border=0 width="$PLOT_DISPLAY_WIDTH_PIX"/>') ;
                                //alert($('input:radio[name=imageType]:checked').val());
                           });
                    return false;
                });                
            });

        </script>

        <fieldset style="height: 100%;">
        <legend><b>Image Controls</b></legend>
        
        <form id="buildCache" method="get" action="buildCacheForObject"></form>   
        $THUMB_FORM_DECLARED
        <div style="display: inline-block;">
        <input form="buildCache" type="hidden" value="$SOURCERY_ID" name="sourceryID"/>
        <input form="buildCache" type="hidden" value="true" name="refetch"/>
        <input form="buildCache" type="hidden" value="displaySourcePage?sourceryID=$SOURCERY_URL_ID" name="from_page"/>
        <input form="buildCache" type="submit" style="display: inline-block;" value="Update Cache [NB: Slow]">
        $THUMB_FORM_CONTROLS
        </div>
        
        <form action="#" id="imageForm" method="post">        
        <input name="name" value="$OBJECT_NAME" type="hidden">
        <p><b>Survey:</b></p> 
        <p>
        $IMAGE_TYPES
        </p>      
        <p><b>Plot:</b></p>
        <p>
        <span style="margin-left: 1.2em; display: inline-block">
        <input type="checkbox" name="showAxes" value=1 $CHECKED_SHOWAXES>
        <label for="showAxes">Show coordinate axes</label>
        </span>
        <span style="margin-left: 1.2em; display: inline-block">
        $CONTOUR_CODE
        </span>
        <span style="margin-left: 1.2em; display: inline-block">
        <input type="checkbox" name="plotSourcePos" value=1 $CHECKED_SOURCEPOS>
        <label for="plotSourcePos">Source position</label>
        </span>
        <span style="margin-left: 1.2em; display: inline-block">
        <input type="checkbox" name="plotNEDObjects" value=1 $CHECKED_NED>
        <label for="plotNEDObjects">NED objects</label>
        </span>
        <span style="margin-left: 1.2em; display: inline-block">
        <input type="checkbox" name="plotSpecObjects" value=1 $CHECKED_SDSS>
        <label for="plotSpecObjects">Spec redshifts</label>
        </span>
        <span style="margin-left: 1.2em; display: inline-block">
        <input type="checkbox" name="plotXMatch" value=1 $CHECKED_XMATCH>
        <label for="plotXMatch">Cross match objects</label>
        </span>
        </p>

        <p align="right">
        <span style="margin-left: 1.2em; display: inline-block">
        <label for="clipSizeArcmin">Image Size (arcmin)</label>
        <input id="sizeSlider" name="clipSizeArcmin" type="range" min="1.0" max="$MAX_SIZE_ARCMIN" step="0.5" value=$CURRENT_SIZE_ARCMIN onchange="printValue('sizeSlider','sizeSliderValue')">
        <input id="sizeSliderValue" type="text" size="2"/>
        </span>
        <span style="margin-left: 1.2em; display: inline-block">
        <label for="gamma">Brightness (&gamma;)</label>
        <input id="gammaSlider" name="gamma" type="range" min="0.2" max="3.0" step="0.2" value=$CURRENT_GAMMA onchange="printValue('gammaSlider','gammaSliderValue')">
        <input id="gammaSliderValue" type="text" size="2"/>
        </span>
        </p>
        <p align="right">
        <input type="submit" style="font-size: 1.05em;" value="Apply">
        </p>
        </form>
                
        </fieldset>
     
        """ 
        
        # For cache build and downloadable thumbnail buttons
        plotFormCode=plotFormCode.replace("$SOURCERY_ID", obj['sourceryID'])
        plotFormCode=plotFormCode.replace("$DEFAULT_CLIP_SIZE_ARCMIN", str(self.configDict['defaultViewSizeArcmin']))
        plotFormCode=plotFormCode.replace("$SOURCERY_URL_ID", self.sourceNameToURL(obj['sourceryID']))        
        if 'downloadableFITS' in self.configDict.keys():
            plotFormCode=plotFormCode.replace('$THUMB_FORM_DECLARED', '<form id="downloadThumbnail" method="get" action="downloadThumbnailFITS">')
            thumbForm="""<div style="display: inline-block;">
            <input form="downloadThumbnail" type="hidden" value="$SOURCERY_ID" name="sourceryID"/>
            <input form="downloadThumbnail" type="submit" style="display: inline-block;" value="Download $IMGDIRLABEL FITS">
            </form>"""
            thumbForm=thumbForm.replace("$SOURCERY_ID", obj['sourceryID'])
            thumbForm=thumbForm.replace("$IMGDIRLABEL", self.configDict['downloadableFITS'])
            plotFormCode=plotFormCode.replace("$THUMB_FORM_CONTROLS", thumbForm)
        else:
            plotFormCode=plotFormCode.replace("$THUMB_FORM_DECLARED", "")
            plotFormCode=plotFormCode.replace("$THUMB_FORM_CONTROLS", "")
        
        # Clickable coordinate script
        clickCoordScript="""$(document).ready(function() {
                $("#imagePlot").on("click", function(event) {
                
                    var x = event.pageX - this.offsetLeft;
                    var y = event.pageY - this.offsetTop;
                    var width = $("#imagePlot").width();
                    var height = $("#imagePlot").height();
                    var CDELT1 = ($("#sizeSliderValue").val()/60.0) / width;
                    var CDELT2 = ($("#sizeSliderValue").val()/60.0) / height;
                    var CRPIX1 = (width/2.0)+1;
                    var CRPIX2 = (height/2.0)+1;
                    var CRVAL1 = $OBJECT_RADEG;
                    var CRVAL2 = $OBJECT_DECDEG; 
                    
                    // debugging
                    //alert("CDELT1=" + CDELT1 + " CDELT2=" + CDELT2 + " CRPIX1=" + CRPIX1 + " CRPIX2=" + CRPIX2 + " CRVAL1=" + CRVAL1 + " CRVAL2=" + CRVAL2);
                    
                    // Adjust for different origin to FITS
                    x=width-x;
                    y=height-y;
                    
                    // Intermediate coords
                    var deg2Rad = Math.PI / 180;
                    var xint = deg2Rad*(CDELT1*(x-CRPIX1));
                    var yint = deg2Rad*(CDELT2*(y-CRPIX2));

                    var Rtheta = Math.sqrt(Math.pow(xint, 2) + Math.pow(yint, 2));
                    var phi = Math.atan2(xint, -yint);
                    var theta = Math.atan(1./Rtheta);

                    var alpha_p = deg2Rad*(CRVAL1);
                    var delta_p = deg2Rad*(CRVAL2);
                    var phi_p = deg2Rad*(180.);

                    var ra = alpha_p + Math.atan2(-Math.cos(theta)*Math.sin(phi-phi_p), Math.sin(theta)*Math.cos(delta_p) - Math.cos(theta)*Math.sin(delta_p)*Math.cos(phi-phi_p));
                    var dec = Math.asin(Math.sin(theta)*Math.sin(delta_p) + Math.cos(theta)*Math.cos(delta_p)*Math.cos(phi-phi_p));
                    
                    ra=(180/Math.PI)*ra;
                    dec=(180/Math.PI)*dec;
                    
                    // debugging
                    //alert("X Coordinate / RA: " + x + " / " + ra + " Y Coordinate / dec: " + y + " / " + dec);
                    
                    if ($('input:checkbox[name=showAxes]').prop('checked') == true)  {
                        alert("Coords not recorded when showing coordinate axes - deselect 'Show coordinate axes'");
                    } else {
                        $('[name=$CLICKABLE_RA]').val(ra);
                        $('[name=$CLICKABLE_DEC]').val(dec);
                    };
                    
                });
            });
        """
        if 'clickableRAField' in self.configDict.keys() and 'clickableDecField' in self.configDict.keys():
            clickCoordScript=clickCoordScript.replace("$CLICKABLE_RA", self.configDict['clickableRAField'])
            clickCoordScript=clickCoordScript.replace("$CLICKABLE_DEC", self.configDict['clickableDecField'])
        else:
            clickCoordScript=""
        plotFormCode=plotFormCode.replace("$CLICKABLE_RADEC_SCRIPT", clickCoordScript)
            
        # Taken out: onChange="this.form.submit();" from all checkboxes ^^^
        plotFormCode=plotFormCode.replace("$PLOT_DISPLAY_WIDTH_PIX", str(self.configDict['plotDisplayWidthPix']))
        plotFormCode=plotFormCode.replace("$OBJECT_NAME", obj['name'])
        plotFormCode=plotFormCode.replace("$OBJECT_RADEG", str(obj['RADeg']))
        plotFormCode=plotFormCode.replace("$OBJECT_DECDEG", str(obj['decDeg']))
        plotFormCode=plotFormCode.replace("$OBJECT_SURVEY", imageType) 
        if 'contourImage' in self.configDict.keys() and self.configDict['contourImage'] != None:
            contourCode='<input type="checkbox" name="plotContours" value=1 $CHECKED_CONTOURS>\n'
            contourCode=contourCode+'<label for="plotContours">Contours ($CONTOUR_IMAGE)</label>'
            contourCode=contourCode.replace("$CONTOUR_IMAGE", self.configDict['contourImage'])
            plotFormCode=plotFormCode.replace("$CONTOUR_CODE", contourCode)
        else:
            plotFormCode=plotFormCode.replace("$CONTOUR_CODE", "")
        
        # NOTE: a bit of hack here: we're incorporated max size arcmin into imageType value
        # We use javascript ^^^ above to split that around ; delimiter
        imageTypesCode=""            
        for label, maxSizeArcmin in zip(self.imDirLabelsList, self.imDirMaxSizeArcminList):
            if label == imageType:
                imageTypesCode=imageTypesCode+'<span style="margin-left: 1.2em; display: inline-block"><input type="radio" name="imageType" value="%s;%.1f" checked>%s</span>\n' % (label, maxSizeArcmin, label)
            else:
                imageTypesCode=imageTypesCode+'<span style="margin-left: 1.2em; display: inline-block"><input type="radio" name="imageType" value="%s;%.1f">%s</span>\n' % (label, maxSizeArcmin, label)
        plotFormCode=plotFormCode.replace("$IMAGE_TYPES", imageTypesCode)
        
        if plotNEDObjects == "true":
            plotFormCode=plotFormCode.replace("$CHECKED_NED", " checked")
        else:
            plotFormCode=plotFormCode.replace("$CHECKED_NED", "")
        if plotSpecObjects == "true":
            plotFormCode=plotFormCode.replace("$CHECKED_SPEC", " checked")
        else:
            plotFormCode=plotFormCode.replace("$CHECKED_SPEC", "")
        if plotSourcePos == "true":
            plotFormCode=plotFormCode.replace("$CHECKED_SOURCEPOS", " checked")
        else:
            plotFormCode=plotFormCode.replace("$CHECKED_SOURCEPOS", "")
        if plotXMatch == "true":
            plotFormCode=plotFormCode.replace("$CHECKED_XMATCH", " checked")
        else:
            plotFormCode=plotFormCode.replace("$CHECKED_XMATCH", "")
        if plotContours == "true":
            plotFormCode=plotFormCode.replace("$CHECKED_CONTOURS", " checked")
        else:
            plotFormCode=plotFormCode.replace("$CHECKED_CONTOURS", "")
        if showAxes == "true":
            plotFormCode=plotFormCode.replace("$CHECKED_SHOWAXES", " checked")
        else:
            plotFormCode=plotFormCode.replace("$CHECKED_SHOWAXES", "")
        
        # This block here probably is hardly useful
        imageMaxSizeArcmin=30.
        #for imType, maxSizeArcmin in zip(self.imDirLabelsList, self.imDirMaxSizeArcminList):
            #if imType == imageType:
                #imageMaxSizeArcmin=mdaxSizeArcmin
                #break
        
        plotFormCode=plotFormCode.replace("$MAX_SIZE_ARCMIN", str(imageMaxSizeArcmin))        
        if clipSizeArcmin == None:
            plotFormCode=plotFormCode.replace("$CURRENT_SIZE_ARCMIN", str(self.configDict['defaultViewSizeArcmin']))
        else:
            plotFormCode=plotFormCode.replace("$CURRENT_SIZE_ARCMIN", str(clipSizeArcmin))
        
        plotFormCode=plotFormCode.replace("$CURRENT_GAMMA", str(gamma))
                
        # Tagging controls (including editable properties of catalog, e.g., for assigning classification or redshifts)
        tagFormCode="""
        <form method="post" action="updateMongoDB">    
        <input name="name" value="$OBJECT_NAME" type="hidden">
        <input name="returnURL" value=$RETURN_URL" type="hidden">
        <fieldset style="height: 100%;">
        <legend><b>Editing Controls</b></legend>
        $CLASSIFICATION_CONTROLS
        $FIELD_CONTROLS
        <p align="right">
        <input type="submit" class="f" style="font-size: 1.05em;" value="Update" $DISABLED_STR>
        </p>
        </fieldset>
        </form>
        """
        if 'fields' in self.configDict.keys():
            if cherrypy.session['editPermission'] == False:
                readOnlyStr="readonly"
                tagFormCode=tagFormCode.replace("$DISABLED_STR", "disabled")
            else:
                readOnlyStr=""
                tagFormCode=tagFormCode.replace("$DISABLED_STR", "")
            tagFormCode=tagFormCode.replace("$PLOT_DISPLAY_WIDTH_PIX", str(self.configDict['plotDisplayWidthPix']))
            tagFormCode=tagFormCode.replace("$OBJECT_NAME", name)
            tagFormCode=tagFormCode.replace("$RETURN_URL", cherrypy.url())
            if 'fields' in self.configDict.keys():
                #fieldsCode="<p><b>Fields:</b>"
                fieldsCode='<p align="left">'
                for fieldDict in self.configDict['fields']:
                    fieldsCode=fieldsCode+'<span style="display: inline-block; margin-bottom: 6pt; margin-right: 8pt">'
                    fieldsCode=fieldsCode+'<label for="%s"><b>%s: </b></label>\n' % (fieldDict['name'], fieldDict['name'])
                    fieldsCode=fieldsCode+'<input type="text" value="%s" name="%s" size=%d %s/>\n' % (str(mongoDict[fieldDict['name']]), 
                                                                                                   fieldDict['name'], 
                                                                                                   fieldDict['displaySize'],
                                                                                                   readOnlyStr)
                    fieldsCode=fieldsCode+"</span>"
                fieldsCode=fieldsCode+"</p>"
                if 'lastUpdated' in mongoDict.keys():
                    lastUpdated=mongoDict['lastUpdated']
                else:
                    lastUpdated='-'
                if 'user' in mongoDict.keys():
                    userName=mongoDict['user']
                else:
                    userName='-'
                fieldsCode=fieldsCode+'<p><label for = "lastUpdated"><b>Last Updated:</b></label>\n'
                fieldsCode=fieldsCode+'<input type="text" value="%s" name="lastUpdated" size=10 readonly/>\n' % (lastUpdated)
                fieldsCode=fieldsCode+'<label for = "user"><b>By User:</b></label>\n'
                fieldsCode=fieldsCode+'<input type="text" value="%s" name="userName" size=10 readonly/></p>\n' % (userName)
                fieldsCode=fieldsCode+"</p>"
            tagFormCode=tagFormCode.replace('$FIELD_CONTROLS', fieldsCode)
        else:
            tagFormCode=tagFormCode.replace('$FIELD_CONTROLS', "")
        
        if 'classifications' in self.configDict.keys():
            classificationsCode="<p><b>Classification:</b></p>\n<p>"
            if cherrypy.session['editPermission'] == False:
                readOnlyStr="disabled"
            else:
                readOnlyStr=""
            for c in self.configDict['classifications']:
                if c == mongoDict['classification']:
                    classificationsCode=classificationsCode+'<span style="margin-left: 1.2em; display: inline-block;"><input type="radio" name="classification" value="%s" checked %s>%s</span>\n' % (c, readOnlyStr, c)
                else:
                    classificationsCode=classificationsCode+'<span style="margin-left: 1.2em; display: inline-block;"><input type="radio" onChange="this.form.submit();" name="classification" value="%s" %s>%s</span>\n' % (c, readOnlyStr, c)
            classificationsCode=classificationsCode+"</p>"
            tagFormCode=tagFormCode.replace("$CLASSIFICATION_CONTROLS", classificationsCode)
        else:
            tagFormCode=tagFormCode.replace("$CLASSIFICATION_CONTROLS", "")
        
        # Optional spectrum plot
        if 'displaySpectra' in self.configDict.keys() and self.configDict['displaySpectra'] == True:
            spectrumCode="""<br><table frame=border cellspacing=0 cols=1 rules=all border=2 width=80% align=center>
            <tbody>
            <tr>
                <th style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 110%;" colspan=6>
                    <b>Spectrum</b>
                </th>
            </tr>
            <tr>
            <td align=center><img src="makeSpectrumPlot?name=$NAME&RADeg=$RA&decDeg=$DEC" align="middle" border=2 width=$SPEC_DISPLAY_WIDTH_PIX/></td>
            </tr>
            </tbody>
            </table>
            """
            spectrumCode=spectrumCode.replace("$NAME", self.sourceNameToURL(name))
            spectrumCode=spectrumCode.replace("$RA", str(obj['RADeg']))
            spectrumCode=spectrumCode.replace("$DEC", str(obj['decDeg']))
            spectrumCode=spectrumCode.replace("$SPEC_DISPLAY_WIDTH_PIX", str(self.configDict['specPlotDisplayWidthPix']))
        else:
            spectrumCode=""
            
        # Put it all together...
        html=templatePage
        html=html.replace("$SPECTRUM_PLOT", spectrumCode)
        html=html.replace("$PLOT_CONTROLS", plotFormCode)
        if 'fields' in self.configDict.keys() or 'classifications' in self.configDict.keys():
            html=html.replace("$TAG_CONTROLS", tagFormCode)
        else:
            html=html.replace("$TAG_CONTROLS", "")
        html=html.replace("$SOURCE_NAME", name)
        html=html.replace("$SIZE_ARC_MIN", "%.1f" % (self.configDict['plotSizeArcmin']))
        html=html.replace("$HOSTED_STR", self.configDict['hostedBy'])

        #for label, caption in zip(self.imDirLabelsList, self.imageCaptions):
            #if label == imageType:
                #html=html.replace("$CAPTION", "%s" % (caption))
                            
        # NED matches table
        self.fetchNEDInfo(obj['name'], obj['RADeg'], obj['decDeg'])
        self.findNEDMatch(obj, NEDObjTypes = self.configDict['NEDObjTypes'])
        nedFileName=self.nedDir+os.path.sep+obj['name'].replace(" ", "_")+".txt"
        nedObjs=catalogTools.parseNEDResult(nedFileName, onlyObjTypes = self.configDict['NEDObjTypes'])
        if len(nedObjs['RAs']) > 0:
            nedTable="""<br><table frame=border cellspacing=0 cols=6 rules=all border=2 width=80% align=center>
            <tbody>
            <tr>
                <th style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 110%;" colspan=6>
                    <b>NED Matches</b>
                </th>
            </tr>
            <tr>
                <td><b>ID</b></td>
                <td><b>Name</b></td>
                <td><b>RA</b></td>
                <td><b>Dec.</b></td>
                <td><b>Object Type</b></td>
                <td><b>Redshift</b></td>
            </tr>
            """                
            for i in range(len(nedObjs['RAs'])):
                rowString="""<tr>
                    <td align=center width=10%>$ID</td>
                    <td align=center width=10%>$NAME</td>
                    <td align=center width=10%>$RA</td>
                    <td align=center width=10%>$DEC</td>
                    <td align=center width=10%>$TYPE</td>
                    <td align=center width=10%>$REDSHIFT</td>
                </tr>
                """
                nedLinkURL="http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES" % (nedObjs['names'][i].replace("+", "%2B").replace(" ", "+"))
                rowString=rowString.replace("$ID", "%s" % (nedObjs['labels'][i]))
                rowString=rowString.replace("$NAME", "<a href =%s>%s</a>" % (nedLinkURL, nedObjs['names'][i]))
                rowString=rowString.replace("$RA", "%.5f" % (nedObjs['RAs'][i]))
                rowString=rowString.replace("$DEC", "%.5f" % (nedObjs['decs'][i]))
                rowString=rowString.replace("$TYPE", "%s" % (nedObjs['sourceTypes'][i]))
                rowString=rowString.replace("$REDSHIFT", "%s" % (nedObjs['redshifts'][i]))
                nedTable=nedTable+rowString
            nedTable=nedTable+"</tbody></table>"
        else:
            nedTable=""
        html=html.replace("$NED_MATCHES_TABLE", nedTable)
        
        # SDSS matches table
        if 'addSpecRedshifts' in self.configDict.keys() and self.configDict['addSpecRedshifts'] == True:
            specRedshifts=catalogTools.fetchSpecRedshifts(obj['name'], obj['RADeg'], obj['decDeg'], 
                                                          redshiftsTable = self.specRedshiftsTab)
            specTable="""<br><table frame=border cellspacing=0 cols=7 rules=all border=2 width=80% align=center>
            <tbody>
            <tr>
                <th style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 110%;" colspan=6>
                    <b>Spectroscopic Redshifts</b>
                </th>
            </tr>
            <tr>
                <td><b>ID</b></td>
                <td><b>RA</b></td>
                <td><b>Dec.</b></td>
                <td><b>z</b></td>
                <td><b>zWarning</b></td>
                <td><b>catalog</b></td> 
            </tr>
            """              
            specCount=0
            for specObj in specRedshifts:
                specCount=specCount+1
                rowString="""<tr>
                    <td align=center width=10%>$ID</td>
                    <td align=center width=10%>$RA</td>
                    <td align=center width=10%>$DEC</td>
                    <td align=center width=10%>$REDSHIFT</td>
                    <td align=center width=10%>$Z_WARNING</td>
                    <td align=center width=10%>$Z_CATALOG</td>
                </tr>
                """
                rowString=rowString.replace("$ID", "%d" % (specCount))
                rowString=rowString.replace("$RA", "%.5f" % (specObj['RADeg']))
                rowString=rowString.replace("$DEC", "%.5f" % (specObj['decDeg']))
                rowString=rowString.replace("$REDSHIFT", "%.3f" % (specObj['z']))
                rowString=rowString.replace("$Z_WARNING", "%s" % (specObj['zWarning']))
                rowString=rowString.replace("$Z_CATALOG", "%s" % (specObj['catalog']))
                specTable=specTable+rowString
            specTable=specTable+"</tbody></table>"
        else:
            specTable=""
        html=html.replace("$SPEC_MATCHES_TABLE", specTable)

        # Source properties table
        propTable="""<br><table frame=border cellspacing=0 cols=2 rules=all border=2 width=80% align=center>
        <tbody>
        <tr>
            <th style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                text-align: center; vertical-align: middle; font-size: 110%;" colspan=2>
                <b>Source Properties</b>
            </th>
        </tr>
        """
        fieldTypes=self.fieldTypesCollection.find().sort('index')
        for f in fieldTypes:
            if f['name'] in obj.keys() and f['name'] not in ['sourceryID', 'cacheBuilt']:
                pkey=f['name']
                rowString="""<tr><td align=left width=50%><b>$KEY_LABEL</b></td>
                <td align=center width=50%>$KEY_VALUE</td></tr>
                """
                rowString=rowString.replace("$KEY_LABEL", pkey)
                if obj[pkey] == "None" or obj[pkey] == None:
                    rowString=rowString.replace("$KEY_VALUE", "-")
                else:
                    if pkey == "NED_name":
                        nedName=obj[pkey]
                        nedLinkURL="http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES" % (nedName.replace("+", "%2B").replace(" ", "+"))
                        rowString=rowString.replace("$KEY_VALUE", "<a href=%s>%s</a>" % (nedLinkURL, nedName))
                    else:
                        rowString=rowString.replace("$KEY_VALUE", str(obj[pkey]))
                # Skip over cross-matched tables with no matches
                prefix=pkey.split("_")[0]
                if prefix not in skipColumnPrefixList:
                    propTable=propTable+rowString
        propTable=propTable+"</td></tr></tbody></table>"
        html=html.replace("$PROPERTIES_TABLE", propTable)
                
        return html   
    
    
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
                self.fieldTypesCollection.insert(fieldDict)
        t0=time.time()
        # We need to do this to avoid hitting 32 Mb limit below when using large databases
        self.sourceCollection.ensure_index([("RADeg", pymongo.ASCENDING)])
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
                        self.sourceCollection.update({'_id': obj['_id']}, {'$set': post}, upsert = False)      
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
        self.sourceCollection.ensure_index([("RADeg", pymongo.ASCENDING)])
        
        # Threaded
        # NOTE: Threads have occassionally given weird issues (e.g., mismatched WISE images)
        # Check that when thread write to disk they don't clash with each other
        if 'threadedCacheBuild' in self.configDict.keys() and self.configDict['threadedCacheBuild'] == True:
            cursor=self.sourceCollection.find({'cacheBuilt': 0}, no_cursor_timeout = True).sort('decDeg').sort('RADeg')
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
            cursor=self.sourceCollection.find({'cacheBuilt': 0}, no_cursor_timeout = True).sort('decDeg').sort('RADeg')
            for obj in cursor:
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
        if self.configDict['addSpecRedshifts'] == True:
            catalogTools.fetchSDSSRedshifts(self.sdssRedshiftsDir, name, RADeg, decDeg,
                                            redshiftsTable = self.specRedshiftsTab)
        if self.configDict['addSDSSImage'] == True:
            self.fetchSDSSImage(name, RADeg, decDeg, refetch = refetch)
        if self.configDict['addDECaLSImage'] == True:
            self.fetchDECaLSImage(name, RADeg, decDeg, refetch = refetch)
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
                    self.sourceCollection.update({'_id': obj['_id']}, {'$set': post}, upsert = False)
        
        # Flag this as done
        self.sourceCollection.update({'_id': obj['_id']}, {'$set': {'cacheBuilt': 1}}, upsert = False)
        
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
                        wcs=astWCS.WCS(imgFileName)
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
            
            self.sourceCollection.ensure_index([("RADeg", pymongo.ASCENDING)])  # Avoid 32 Mb limit
            objList=self.sourceCollection.find(no_cursor_timeout = True).sort('decDeg').sort('RADeg')

            for obj in objList:
                outFileName=outDir+os.path.sep+catalogTools.makeRADecString(obj['RADeg'], obj['decDeg'])+".jpg"
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
                            data=img[0].data                     
                        clip=astImages.clipImageSectionWCS(data, wcs, obj['RADeg'], obj['decDeg'],
                                                           maxSizeArcmin/60.0)
                        
                        # Sanity check - did we pick an image where object is in zeroed border?
                        pixCoords=clip['wcs'].wcs2pix(obj['RADeg'], obj['decDeg'])
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
                            clip['data']=catalogTools.byteSwapArr(clip['data'])
                            # Avoid cython type troubles
                            if clip['data'].dtype != np.float32:
                                clip['data']=np.array(clip['data'], dtype = np.float32)
                            rMap=sourceryCython.makeDegreesDistanceMap(clip['data'], clip['wcs'], obj['RADeg'], obj['decDeg'], 100.0)
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

        
        
        
        
