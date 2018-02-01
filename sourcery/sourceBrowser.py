"""

    Copyright 2014 Matt Hilton (matt.hilton@mykolab.com)
    
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
import astropy.table as atpy
import operator
import urllib
import urllib2
import glob
from astLib import *
import pyfits
import numpy as np
import pylab as plt
import matplotlib.patches as patches
from scipy import ndimage
import sourcery
from sourcery import catalogTools
from sourcery import specFeatures
import ConfigParser
import requests
import sys
import time
import datetime
import string
import re
import base64
from PIL import Image
import copy
import StringIO
import tempfile
import pymongo
from bson.son import SON
import pyximport; pyximport.install()
import sourceryCython
import cherrypy
import pickle
import pyvips
import IPython
 
#-------------------------------------------------------------------------------------------------------------
class SourceBrowser(object):
    
    def __init__(self, configFileName, preprocess = False, buildDatabase = False):
        
        # Parse config file
        self.parseConfig(configFileName)
        
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
        if 'DESServicesConfigPath' in self.configDict.keys():
            configParser=ConfigParser.RawConfigParser()   
            DESConfigPath=self.configDict['DESServicesConfigPath']#os.environ['HOME']+os.path.sep+".desservices.ini"
            configParser.read(DESConfigPath)
            self.DESUser=configParser.get('db-dessci', 'user')
            self.DESPasswd=configParser.get('db-dessci', 'passwd')
        else:
            self.DESUser=None
            self.DESPasswd=None
        self.DESTokenID=None    # Used for interacting with DESCuts server
        
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
        self.tableDisplayColumns=['name', 'RADeg', 'decDeg']+self.configDict['tableDisplayColumns']
        self.tableDisplayColumnLabels=['Name', 'RA (degrees)', 'Dec. (degrees)']+self.configDict['tableDisplayColumns']
        self.tableDisplayColumnFormats=['%s', '%.6f', '%.6f']+self.configDict['tableDisplayColumnFormats']
        # Source pages
        self.sourceDisplayColumns=[]

        # Support for tagging, classification etc. of candidates
        # We have three things here:
        #   tags: these work as flag columns - eventually, users can add any tag, we can then search for objects matching that tag
        #   classification: what kind of object is in the catalog? List of types defined in config file
        #   fields: these are editable, used to, e.g., assign a redshift and source
        # We create empty entries in the MongoDB and corresponding database columns if we cannot find an existing entry
        # We will only add/populate columns of atpy table when writing output
        #if 'tags' in self.configDict.keys():
            #for t in self.configDict['tags']:
                #self.tab.add_column(t, np.zeros(len(self.tab), dtype = bool))
            #self.sourceDisplayColumns=self.sourceDisplayColumns+self.configDict['tags']
        if 'fields' in self.configDict.keys():
            formatsList=[]
            for f, t in zip(self.configDict['fields'], self.configDict['fieldTypes']):
                if t == 'number':
                    formatsList.append('%.3f')
                elif t == 'text':
                    formatsList.append('%s')
            self.sourceDisplayColumns=self.sourceDisplayColumns+self.configDict['fields']
            self.tableDisplayColumns=self.tableDisplayColumns+self.configDict['fields']
            self.tableDisplayColumnLabels=self.tableDisplayColumnLabels+self.configDict['fields']
            self.tableDisplayColumnFormats=self.tableDisplayColumnFormats+formatsList
        if 'classifications' in self.configDict.keys():
            self.sourceDisplayColumns=self.sourceDisplayColumns+["classification"]
            self.tableDisplayColumns=self.tableDisplayColumns+["classification"]
            self.tableDisplayColumnLabels=self.tableDisplayColumnLabels+["classification"]
            self.tableDisplayColumnFormats=self.tableDisplayColumnFormats+["%s"]
        
        # Now tracking when changes are made
        self.sourceDisplayColumns.append('lastUpdated')
        self.tableDisplayColumns.append('lastUpdated')
        self.tableDisplayColumnLabels.append('lastUpdated')
        self.tableDisplayColumnFormats.append('%s')
                               
        # We will generate images dynamically... here we set up info like labels and captions
        self.imageLabels=[]      # labels at top of each source page that allow us to select image to view
        self.imageCaptions=[]    # caption that goes under image shown on the source pages
        if "imageDirsLabels" in self.configDict.keys():
            for label in self.configDict['imageDirsLabels']:
                self.imageLabels.append(label)
                self.imageCaptions.append("I don't think we're using this any more")
        # SDSS colour .jpgs
        if "addSDSSImage" in self.configDict.keys() and self.configDict['addSDSSImage'] == True:
            label="SDSS"
            self.imageLabels.append(label)
            self.imageCaptions.append("%.1f' x %.1f' false color (g,r,i) SDSS DR13 image. The source position is marked with the white cross.<br>Objects marked with green circles are in NED; objects marked with red squares have SDSS DR14 spectroscopic redshifts." % (self.configDict['plotSizeArcmin'], self.configDict['plotSizeArcmin']))
        # DES colour .jpgs
        if "addDESImage" in self.configDict.keys() and self.configDict['addDESImage'] == True:
            label="DES"
            self.imageLabels.append(label)
            self.imageCaptions.append("%.1f' x %.1f' false color (g,r,i) DES DR1 image. The source position is marked with the white cross.<br>Objects marked with green circles are in NED; objects marked with red squares have SDSS DR14 spectroscopic redshifts." % (self.configDict['plotSizeArcmin'], self.configDict['plotSizeArcmin']))
        # PS1 colour .jpgs
        if "addPS1Image" in self.configDict.keys() and self.configDict['addPS1Image'] == True:
            label="PS1"
            self.imageLabels.append(label)
            self.imageCaptions.append("%.1f' x %.1f' false color (g,r,i) Pan-STARSS PS1 3pi image. The source position is marked with the white cross.<br>Objects marked with green circles are in NED; objects marked with red squares have SDSS DR14 spectroscopic redshifts." % (self.configDict['plotSizeArcmin'], self.configDict['plotSizeArcmin']))
        # PS1IR colour .jpgs
        if "addPS1IRImage" in self.configDict.keys() and self.configDict['addPS1IRImage'] == True:
            label="PS1IR"
            self.imageLabels.append(label)
            self.imageCaptions.append("%.1f' x %.1f' false color (i,z,y) Pan-STARSS PS1 3pi image. The source position is marked with the white cross.<br>Objects marked with green circles are in NED; objects marked with red squares have SDSS DR14 spectroscopic redshifts." % (self.configDict['plotSizeArcmin'], self.configDict['plotSizeArcmin']))
        # CFHTLS colour .jpgs
        if "addCFHTLSImage" in self.configDict.keys() and self.configDict['addCFHTLSImage'] == True:
            label="CFHTLS"
            self.imageLabels.append(label)
            self.imageCaptions.append("%.1f' x %.1f' false color (g,r,i) CFHT Legacy Survey image. The source position is marked with the white cross.<br>Objects marked with green circles are in NED; objects marked with red squares have SDSS DR14 spectroscopic redshifts." % (self.configDict['plotSizeArcmin'], self.configDict['plotSizeArcmin']))
        # unWISE colour .jpgs
        if "addUnWISEImage" in self.configDict.keys() and self.configDict['addUnWISEImage'] == True:
            label="unWISE"
            self.imageLabels.append(label)
            self.imageCaptions.append("%.1f' x %.1f' false color (W1, W2) unWISE image. The source position is marked with the white cross.<br>Objects marked with green circles are in NED; objects marked with red squares have SDSS DR14 spectroscopic redshifts." % (self.configDict['plotSizeArcmin'], self.configDict['plotSizeArcmin']))
        # Skyview images
        if 'skyviewLabels' in self.configDict.keys():
            for label in self.configDict['skyviewLabels']:
                self.imageLabels.append(label)
                self.imageCaptions.append("%.1f' x %.1f' false color %s image. The source position is marked with the white cross.<br>Objects marked with green circles are in NED; objects marked with red squares have SDSS DR14 spectroscopic redshifts." % (self.configDict['plotSizeArcmin'], self.configDict['plotSizeArcmin'], label))
       
        # Pre-processing
        # NOTE: this includes generating .jpgs from user-specified, probably proprietary, image dirs
        # We can run this on the webserver by going to localhost:8080/preprocess
        # We might want to prevent that and force it to run manually only...
        if preprocess == True:
            self.preprocess()

        # This sets size of table view - view is controlled with session variables
        self.tableViewRows=40


    def matchTags(self, obj):
        """Find match in MongoDB to obj row from tab. If we don't find one, return a dictionary with blank
        values where fields would be.
        
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
            for f, t in zip(self.configDict['fields'], self.configDict['fieldTypes']):
                if f not in mongoDict.keys():
                    if t == 'number':
                        mongoDict[f]=0.0
                    elif t == 'text':
                        mongoDict[f]=""
        
        # Strip out _id and loc, because whatever is calling this routine won't want them
        if '_id' in mongoDict.keys():
            del mongoDict['_id']
        if 'loc' in mongoDict.keys():
            del mongoDict['loc']
    
        return mongoDict
    

    def buildDatabase(self):
        """Import .fits table into MongoDB database as sourceCollection. Delete any pre-existing catalog
        there. Do all the cross matching at this stage also.
        
        We also cross match the tagsCollection onto sourceCollection, to save doing a full cross match again
        later. When we need to update, we'll update both tagsCollection and sourceCollection.
        
        We also store a list of fields and types in a collection, so that we can (a) use this for the 
        constraints help page; (b) keep fields in a sensible order (assuming input catalogs are in sensible 
        order)
        
        """
        
        print ">>> Building database ..."
        t0=time.time()
        
        # Delete any pre-existing entries
        self.db.drop_collection('sourceCollection')
        self.db.drop_collection('fieldTypes')
        
        # Table set up
        tab=atpy.Table().read(self.configDict['catalogFileName'])
        tab.sort(["RADeg", "decDeg"])
        
        # In case we don't have a 'name' column, we relabel a given column
        if 'nameColumn' in self.configDict.keys() and self.configDict['nameColumn'] != "":
            if 'name' not in tab.keys():
                tab.rename_column(self.configDict['nameColumn'], 'name')

        # Only load tables once for cross matching
        xTabsDict={}
        xMatchRadiusDeg={}
        if 'crossMatchCatalogFileNames' in self.configDict.keys():
            for f, label, radiusArcmin in zip(self.configDict['crossMatchCatalogFileNames'], 
                                              self.configDict['crossMatchCatalogLabels'], 
                                              self.configDict['crossMatchRadiusArcmin']):
                xTabsDict[label]=atpy.Table().read(f)
                xMatchRadiusDeg[label]=radiusArcmin/60.
        #if 'crossMatchRadiusArcmin' in self.configDict.keys():
            #crossMatchRadiusDeg=self.configDict['crossMatchRadiusArcmin']/60.0
                
        # Import each object into MongoDB, cross matching as we go...
        idCount=0
        fieldTypesList=[]   # Used for making sensible column order later
        fieldTypesDict={}   # Used for tracking types for help page
        for row in tab:
            # Need an id number for table display
            idCount=idCount+1
            newPost={'index': idCount}
            newPost['name']=row['name']
            newPost['RADeg']=row['RADeg']
            newPost['decDeg']=row['decDeg']
            
            print "... adding %s to database (%d/%d) ..." % (row['name'], idCount, len(tab))
            
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
                elif tab.columns[key].dtype.name.find("string") != -1:
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
                    raise Exception, "Unknown data type in column '%s' of table cross match table '%s'" % (key, label)            
                        
            # NED cross match
            if 'addNEDMatches' in self.configDict.keys() and self.configDict['addNEDMatches'] == True:
                self.findNEDMatch(newPost)
                stringKeys=['NED_name']
                numberKeys=['NED_z', 'NED_distArcmin', 'NED_RADeg', 'NED_decDeg']
                typesList=['text', 'number']
                for t, l in zip(typesList, [stringKeys, numberKeys]):
                    for key in l:
                        if key not in fieldTypesList:
                            fieldTypesList.append(key)
                            fieldTypesDict[key]=t
            
            # All other cross matches
            # We won't add name, RADeg, decDeg, if name matches our source catalog name
            # This is an easy/lazy way to add extra columns like SZ properties or photo-zs
            if 'crossMatchCatalogFileNames' in self.configDict.keys():
                for label in self.configDict['crossMatchCatalogLabels']:
                    xTab=xTabsDict[label]
                    crossMatchRadiusDeg=xMatchRadiusDeg[label]
                    r=astCoords.calcAngSepDeg(row['RADeg'], row['decDeg'], xTab['RADeg'], xTab['decDeg'])
                    if r.min() < crossMatchRadiusDeg:
                        xMatch=xTab[np.where(r == r.min())][0]
                        xKeysList=list(xTab.keys())
                        # Undocumentated feature - for the BCG position table, which only has positions
                        if xTab.keys() == ['name', 'RADeg', 'decDeg']:  # for atpy, this was a tuple - for astropy, a list?
                            zapPosKeys=False
                        else:
                            zapPosKeys=True
                        if 'name' in xKeysList and xMatch['name'] == newPost['name']:
                            del xKeysList[xKeysList.index('name')]
                            if zapPosKeys == True:
                                del xKeysList[xKeysList.index('RADeg')]
                                del xKeysList[xKeysList.index('decDeg')]
                        for key in xKeysList:
                            newKey='%s_%s' % (label, key)
                            # Just to make sure MongoDB happy with data types
                            # e.g., redmapper .fits table doesn't play nicely by default
                            if xTab.columns[key].dtype.name.find("int") != -1:
                                newPost[newKey]=int(xMatch[key])
                                if newKey not in fieldTypesList:
                                    fieldTypesList.append(newKey)
                                    fieldTypesDict[newKey]="number"
                            elif xTab.columns[key].dtype.name.find("string") != -1:
                                newPost[newKey]=str(xMatch[key])
                                if newKey not in fieldTypesList:
                                    fieldTypesList.append(newKey)
                                    fieldTypesDict[newKey]="text"
                            elif xTab.columns[key].dtype.name.find("float") != -1:
                                newPost[newKey]=float(xMatch[key])
                                if newKey not in fieldTypesList:
                                    fieldTypesList.append(newKey)
                                    fieldTypesDict[newKey]="number"
                            elif xTab.columns[key].dtype.name.find("bool") != -1:
                                newPost[newKey]=bool(xMatch[key])
                                if newKey not in fieldTypesList:
                                    fieldTypesList.append(newKey)
                                    fieldTypesDict[newKey]="number"
                            else:
                                raise Exception, "Unknown data type in column '%s' of table cross match table '%s'" % (key, label)
                        skipMatchKeys=False
                        newPost['%s_match' % (label)]=1
                        numberKeys=['match']
                        if 'name' in xTab.keys() and xMatch['name'] == newPost['name']:
                            skipMatchKeys=True    
                        if skipMatchKeys == False:
                            newPost['%s_RADeg' % (label)]=float(xMatch['RADeg'])
                            newPost['%s_decDeg' % (label)]=float(xMatch['decDeg'])
                            newPost['%s_distArcmin' % (label)]=r.min()*60.0
                            numberKeys=numberKeys+['RADeg', 'decDeg', 'distArcmin']
                        for key in numberKeys:
                            if '%s_%s' % (label, key) not in fieldTypesList:
                                fieldTypesList.append('%s_%s' % (label, key))
                                fieldTypesDict['%s_%s' % (label, key)]="number"                        
                    else:
                        newPost['%s_match' % (label)]=0

            # Match with tagsCollection
            tagsDict=self.matchTags(newPost)
            for key in tagsDict:
                newPost[key]=tagsDict[key]
            
            self.sourceCollection.insert(newPost)

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
        print "... building database complete: took %.1f sec ..." % (t1-t0)
            

    def parseColumnDescriptionsFile(self):
        """Reads a text file containing column descriptions. The format of the file is:
        
        columnName: description
        
        Any white space between : and description is stripped.
        
        Returns a dictionary with keys {'columnName': 'description'}
        
        """
        
        descriptionsDict={}
        if 'descriptionsFileName' in self.configDict.keys():
            inFile=file(self.configDict['descriptionsFileName'], "r")
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
        
        """
        self.configDict={}
        inFile=file(configFileName, "r")
        lines=inFile.readlines()
        inFile.close()
        for line in lines:
            if line[0] != "#" and len(line) > 3 and line.find("=") != -1:
                # We do things this way to allow = in values
                equalIndex=line.find('=')
                value=line[equalIndex+1:]
                key=line[:equalIndex]
                value=value.rstrip().lstrip()
                key=key.rstrip()
                if value == 'True':
                    value=True
                elif value == 'False':
                    value=False
                elif value[0] == '[':
                    # This is now complicated by having to handle , inside "" or ''
                    lst=[]
                    items=value.replace("[", "").replace("]", "")#.split(",")
                    if len(items) == 0:
                        continue
                    if items[0] == '"':
                        delim='"'
                    elif items[0] == "'":
                        delim="'"
                    else:
                        delim=""
                    if delim != "":
                        delimIndices=[]
                        for i in range(len(items)):
                            if items[i] == delim:
                                delimIndices.append(i)
                        extractedItems=[]
                        for i in range(len(delimIndices)-1):
                            extractedItems.append(items[delimIndices[i]:delimIndices[i+1]+1])                   
                        validItems=[]
                        for b in extractedItems:
                            if b[0] == delim and b[-1] == delim:
                                candidateItem=b.replace(delim, "").lstrip().rstrip()
                                if candidateItem != ',':
                                    validItems.append(candidateItem)
                        value=validItems
                    else:
                        # In this case, a list of numbers or True, False
                        value=[]
                        extractedItems=items.split(",")
                        for b in extractedItems:
                            if b.lstrip().rstrip() == 'True':
                                value.append(True)
                            elif b.lstrip().rstrip() == 'False':
                                value.append(False)
                            elif b.lstrip().rstrip() == 'None':
                                value.append(None)
                            else:
                                value.append(float(b))
                elif value[0] == '(':
                    items=value.replace("(", "").replace(")", "").split(",")
                    lst=[]
                    for i in items:
                        lst.append(float(i))
                    value=tuple(lst)
                elif value[0] == "'" or value[0] == '"':
                    value=str(value.replace("'", "").replace('"', ""))
                self.configDict[key]=value
                
        # Not sure of a better way to force things which look like numbers to be floats, not strings
        for key in self.configDict.keys():
            if type(self.configDict[key]) == str:
                try:
                    self.configDict[key]=float(self.configDict[key])
                except:
                    continue
        
        # Add root path where necessary in place
        if 'sourceryPath' in self.configDict.keys() and self.configDict['sourceryPath'] != "":
            rootDir=self.configDict['sourceryPath'].rstrip(os.path.sep)
            keysToFix=["cacheDir", "skyviewCacheDir", "newsFileName", "crossMatchCatalogFileNames"]
            for k in keysToFix:
                if type(self.configDict[k]) == list:
                    for i in range(len(self.configDict[k])):
                        self.configDict[k][i]=rootDir+os.path.sep+self.configDict[k][i]
                else:
                    self.configDict[k]=rootDir+os.path.sep+self.configDict[k]
        

    def addNews(self):
        """Parse news file, if there is one, filling up self.configDict['newsItems'].
        
        """
        
        if os.path.exists(self.configDict['newsFileName']) == True:
            inFile=file(self.configDict['newsFileName'], "r")
            lines=inFile.readlines()
            inFile.close()
        else:
            lines=[]
            
        newsItems=[]
        for line in lines:
            if line[0] != "#" and len(line) > 3:
                newsItems.append(line.replace("\n", ""))
        self.configDict['newsItems']=newsItems
        
        
    def fetchNEDInfo(self, obj, retryFails = False):
        """Fetches NED info for given obj (which must have name, RADeg, decDeg keys) - just stores on disk 
        in cacheDir - we'll retrieve it later as needed.
        
        """
        halfMatchBoxLengthDeg=5.0/60.0
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']
        RAMin=RADeg-halfMatchBoxLengthDeg
        RAMax=RADeg+halfMatchBoxLengthDeg
        decMin=decDeg-halfMatchBoxLengthDeg
        decMax=decDeg+halfMatchBoxLengthDeg
        outFileName=self.nedDir+os.path.sep+name.replace(" ", "_")+".txt"        
        if os.path.exists(outFileName) == False:
            print "... fetching NED info for %s ..." % (name)
            try:                
                urllib.urlretrieve("http://ned.ipac.caltech.edu/cgi-bin/objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon=%.6fd&lat=%.6fd&radius=%.2f&dot_include=ANY&in_objtypes1=GGroups&in_objtypes1=GClusters&in_objtypes1=QSO&in_objtypes2=Radio&in_objtypes2=SmmS&in_objtypes2=Infrared&in_objtypes2=Xray&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=ascii_tab&zv_breaker=30000.0&list_limit=5&img_stamp=YES" % (RADeg, decDeg, halfMatchBoxLengthDeg*60.0), filename = outFileName)
            except:
                # This will block if our server can't see NED - so, we'll write something in the file
                # We can test for this above and re-do the query if needed
                print "WARNING: couldn't get NED info"
                outFileName=None


    def findNEDMatch(self, obj, NEDObjType = "GClstr"):
        """Checks if there is a NED match for obj. Uses matching radius specified in config file by
        crossMatchRadiusArcmin.
        
        """
        
        if "NEDObjType" in self.configDict.keys():
            NEDObjType=self.configDict['NEDObjType']
            
        nedFileName=self.nedDir+os.path.sep+obj['name'].replace(" ", "_")+".txt"
        nedObjs=catalogTools.parseNEDResult(nedFileName)
            
        # Flag matches against clusters - choose nearest one
        rMin=10000
        crossMatchRadiusDeg=self.configDict['NEDCrossMatchRadiusArcmin']/60.0
        clusterMatch={}
        if len(nedObjs['RAs']) > 0:
            for i in range(len(nedObjs['RAs'])):
                ned=nedObjs
                if ned['sourceTypes'][i] == NEDObjType:
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
            
            
    #def addCrossMatchTabs(self):
        #"""Cross matches external catalog crossMatchTab to self.tab, adding matches in place.
        #If there is a column called 'redshift', we include that
        
        #"""
        #if 'crossMatchCatalogFileNames' in self.configDict.keys():
            #for f, label in zip(self.configDict['crossMatchCatalogFileNames'], self.configDict['crossMatchCatalogLabels']):
                #xTab=atpy.Table().read(f)
                #self.tab.add_column('%s_name' % (label), ['__________________________']*len(self.tab))
                #self.tab['%s_name' % (label)]=None
                #self.tab.add_column('%s_z' % (label), [np.nan]*len(self.tab))
                #self.tab.add_column('%s_distArcmin' % (label), [np.nan]*len(self.tab))
                #self.tab.add_column('%s_RADeg' % (label), [np.nan]*len(self.tab))
                #self.tab.add_column('%s_decDeg' % (label), [np.nan]*len(self.tab))
                #self.tab.add_column('%s_match' % (label), np.zeros(len(self.tab)))
                #self.tableDisplayColumns=self.tableDisplayColumns+["%s_name" % (label)]
                #self.tableDisplayColumnLabels=self.tableDisplayColumnLabels+["%s" % (label)]
                #self.tableDisplayColumnFormats=self.tableDisplayColumnFormats+["%s"]
                #self.sourceDisplayColumns=self.sourceDisplayColumns+["%s_name" % (label), "%s_z" % (label), "%s_RADeg" % (label), "%s_decDeg" % (label), "%s_distArcmin" % (label)]

                ## Flag matches against clusters - choose nearest one
                #zKeys=['z', 'redshift', 'Z', 'REDSHIFT']
                #nameKeys=['name', 'id', 'NAME', 'ID']
                #crossMatchRadiusDeg=self.configDict['crossMatchRadiusArcmin']/60.0
                #for row in self.tab:
                    #r=astCoords.calcAngSepDeg(row['RADeg'], row['decDeg'], xTab['RADeg'], xTab['decDeg'])
                    #if r.min() < crossMatchRadiusDeg:
                        #xMatch=xTab[np.where(r == r.min())][0]
                        #for zKey in zKeys:
                            #if zKey in xTab.keys():
                                #row['%s_z' % (label)]=float(xMatch[zKey])
                        #for nameKey in nameKeys:
                            #if nameKey in xTab.keys():
                                #row['%s_name' % (label)]=str(xMatch[nameKey])
                        #row['%s_RADeg' % (label)]=float(xMatch['RADeg'])
                        #row['%s_decDeg' % (label)]=float(xMatch['decDeg'])
                        #row['%s_distArcmin' % (label)]=r.min()*60.0
                        #row['%s_match' % (label)]=1
            

    def fetchPS1Image(self, obj, refetch = False):
        """Fetches Pan-STARRS gri .jpg using the cutout webservice.
        
        """
       
        ps1CacheDir=self.cacheDir+os.path.sep+"PS1"
        if os.path.exists(ps1CacheDir) == False:
            os.makedirs(ps1CacheDir)
        
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']  
        if decDeg < -30:
            print "... outside PS1 area - skipping ..."
            return None
        
        outFileName=ps1CacheDir+os.path.sep+name.replace(" ", "_")+".jpg"
        
        if os.path.exists(outFileName) == False or refetch == True:
        
            if os.path.exists('ps1tmp.html') == True:
                os.remove('ps1tmp.html')
                
            # 1920 pixels is 8 arcmin on PS1 scale
            PS1PlotSizePix=int(round(self.configDict['plotSizeArcmin']*240))
            
            urlString="http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=%.6f+%.6f&filter=color&filter=g&filter=r&filter=i&filetypes=stack&auxiliary=data&size=%d&output_size=1024&verbose=0&autoscale=99.500000&catlist=" % (RADeg, decDeg, PS1PlotSizePix)
            urllib.urlretrieve(urlString, filename = 'ps1tmp.html')
            
            inFile=file('ps1tmp.html', 'r')
            lines=inFile.readlines()
            inFile.close()
            for line in lines:
                if line.find("fitscut.cgi") != -1 and line.find("green") != -1:
                    break
            urlString='http://'+line.split('src="//')[-1].split('"')[0]
            try:
                urllib.urlretrieve(urlString, filename = outFileName)
            except:
                print "... WARNING: couldn't get PS1 image ..."
                print urlString
                outFileName=None


    def fetchPS1IRImage(self, obj, refetch = False):
        """Fetches Pan-STARRS izy .jpg using the cutout webservice.
        
        """
       
        ps1CacheDir=self.cacheDir+os.path.sep+"PS1IR"
        if os.path.exists(ps1CacheDir) == False:
            os.makedirs(ps1CacheDir)
        
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']  
        if decDeg < -30:
            print "... outside PS1 area - skipping ..."
            return None
        
        outFileName=ps1CacheDir+os.path.sep+name.replace(" ", "_")+".jpg"
        
        if os.path.exists(outFileName) == False or refetch == True:
        
            if os.path.exists('ps1tmp.html') == True:
                os.remove('ps1tmp.html')
                
            # 1920 pixels is 8 arcmin on PS1 scale
            PS1PlotSizePix=int(round(self.configDict['plotSizeArcmin']*240))
            
            urlString="http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=%.6f+%.6f&filter=color&filter=i&filter=z&filter=y&filetypes=stack&auxiliary=data&size=%d&output_size=1024&verbose=0&autoscale=99.500000&catlist=" % (RADeg, decDeg, PS1PlotSizePix)
            urllib.urlretrieve(urlString, filename = 'ps1tmp.html')
            
            inFile=file('ps1tmp.html', 'r')
            lines=inFile.readlines()
            inFile.close()
            for line in lines:
                if line.find("fitscut.cgi") != -1 and line.find("green") != -1:
                    break
            urlString='http://'+line.split('src="//')[-1].split('"')[0]
            try:
                urllib.urlretrieve(urlString, filename = outFileName)
            except:
                print "... WARNING: couldn't get PS1 image ..."
                print urlString
                outFileName=None
                

    def fetchSDSSImage(self, obj, refetch = False):
        """Fetches the SDSS .jpg for the given image size using the casjobs webservice.
        
        makeSDSSPlots loads these jpegs in, and use matplotlib to make them into plots with
        coord axes etc.
        
        The way we're handling image directories is a bit clunky at the moment...
        
        """
    
        sdssCacheDir=self.cacheDir+os.path.sep+"SDSS"
        if os.path.exists(sdssCacheDir) == False:
            os.makedirs(sdssCacheDir)
            
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']                
        outFileName=sdssCacheDir+os.path.sep+name.replace(" ", "_")+".jpg"
        SDSSWidth=1200.0
        SDSSScale=(self.configDict['plotSizeArcmin']*60.0)/SDSSWidth # 0.396127
        if os.path.exists(outFileName) == False or refetch == True:
            #urlString="http://skyservice.pha.jhu.edu/DR10/ImgCutout/getjpeg.aspx?ra="+str(RADeg)+"&dec="+str(decDeg)
            urlString="http://skyserver.sdss.org/dr13/SkyServerWS/ImgCutout/getjpeg?TaskName=Skyserver.Chart.Image&ra="+str(RADeg)+"&dec="+str(decDeg)
            urlString=urlString+"&scale="+str(SDSSScale)+"&width="+str(int(SDSSWidth))+"&height="+str(int(SDSSWidth))
            try:
                urllib.urlretrieve(urlString, filename = outFileName)
            except:
                print "... WARNING: couldn't get SDSS image ..."
                print urlString
                outFileName=None
    
    
    def fetchDESImage(self, obj, refetch = False, numRetries = 5):
        """Make DES co-add .jpg using the publicly available DR1 .tiff files. 
        
        We avoid the descuts server, as currently that serves images which don't correspond to what we ask for.
        
        This will be monstrously slow because each tile image is ~160 Mb. But the code that takes the .tiff
        and spits out a clipped .jpg. is fast enough.
        
        NOTE: The below is no longer true... but left here in case we need to write it again...
        
        This will only work if you have DES data access rights - 'DESServicesConfigPath' in the
        .config file should specify the path to the file containing these 
        (e.g., $HOME/.desservices.ini if you are using easyaccess).
        
        """
        
        desCacheDir=self.cacheDir+os.path.sep+"DES"
        if os.path.exists(desCacheDir) == False:
            os.makedirs(desCacheDir)
        
        if obj['decDeg'] > 5:
            print "... outside DES dec range - skipping ..."
            return None
        
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']                
        outFileName=desCacheDir+os.path.sep+name.replace(" ", "_")+".jpg"
        
        # Procedure: spin through tile WCSs, find which tiles we need, download if necessary, paste pixels into low-res image (unWISE style)
        if os.path.exists(outFileName) == False or refetch == True:
                       
            # Blank WCS
            CRVAL1, CRVAL2=obj['RADeg'], obj['decDeg']
            sizePix=1024
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
            outWCS=astWCS.WCS(newHead, mode='pyfits')
            outData=np.zeros([sizePix, sizePix, 3], dtype = np.uint8)
            RAMin, RAMax, decMin, decMax=outWCS.getImageMinMaxWCSCoords()
            if RAMax-RAMin > 1.0:   # simple check for 0h crossing... assuming no-one wants images > a degree across
                RAMax=-(360-RAMax)
                temp=RAMin
                RAMin=RAMax
                RAMax=temp
            checkCoordsList=[[RAMin, decMin], [RAMin, decMax], [RAMax, decMin], [RAMax, decMax]]
            
            # Spin though all DES tile WCSs and identify which tiles contain our image (use all four corners; takes 0.4 sec)
            matchTilesList=[]
            for tileName in self.DESWCSDict.keys():
                wcs=self.DESWCSDict[tileName]
                for c in checkCoordsList:
                    pixCoords=wcs.wcs2pix(c[0], c[1])
                    if pixCoords[0] >= 0 and pixCoords[0] < wcs.header['NAXIS1'] and pixCoords[1] >= 0 and pixCoords[1] < wcs.header['NAXIS2']: 
                        if tileName not in matchTilesList:
                            matchTilesList.append(tileName)

            if matchTilesList == []:
                print "... object not in any DES tiles ..."
                return None
            else:

                # We work with the .tiff files... downloading .fits images could be done similarly
                for tileName in matchTilesList:
                    matchTab=self.DESTileTab[np.where(self.DESTileTab['TILENAME'] == tileName)][0]
                    tiffFileName=desCacheDir+os.path.sep+tileName+".tiff"
                    tileJPGFileName=tiffFileName.replace(".tiff", ".jpg")
                    if os.path.exists(tileJPGFileName) == False:
                        if os.path.exists(tiffFileName) == False:
                            print "... downloading .tiff image for tileName = %s ..." % (tileName)
                            try:
                                urllib.urlretrieve(str(matchTab['TIFF_COLOR_IMAGE']), tiffFileName)
                            except:
                                os.remove(tiffFileName)
                                raise Exception, "downloading DES .tiff image failed"
                        # NOTE: we use pyvips, because images are too big for PIL
                        # We save disk space by caching a lower quality version of the entire tile
                        print "... converting .tiff for tileName = %s to .jpg ..." % (tileName)
                        im=pyvips.Image.new_from_file(tiffFileName, access = 'sequential')
                        im.write_to_file(tileJPGFileName+'[Q=80]')
                        os.remove(tiffFileName)
                    else:
                        im=pyvips.Image.new_from_file(tileJPGFileName, access = 'sequential')
                    
                    # Splat pixels from the .tiff into our small image WCS, from which we'll make the .jpg
                    d=np.ndarray(buffer = im.write_to_memory(), dtype = np.uint8, shape = [im.height, im.width, im.bands])
                    d=np.flipud(d)
                    inWCS=self.DESWCSDict[tileName]
                    for y in range(sizePix):
                        for x in range(sizePix):
                            outRADeg, outDecDeg=outWCS.pix2wcs(x, y)
                            inX, inY=inWCS.wcs2pix(outRADeg, outDecDeg)
                            inX=int(round(inX))
                            inY=int(round(inY))
                            if inX >= 0 and inX < d.shape[1]-1 and inY >= 0 and inY < d.shape[0]-1:
                                outData[y, x]=d[inY, inX]
                
                # Flips needed to get N at top, E at left
                outData=np.flipud(outData)
                outData=np.fliplr(outData)
                
                # We could do this with vips... but lazy...
                outIm=Image.fromarray(outData)
                outIm.save(outFileName)
                print "... made cut-out .jpg ..."
                
    
    def fetchCFHTLSImage(self, obj, refetch = False):
        """Retrieves colour .jpg from CFHT legacy survey. Returns True if successful, False if not
        
        NOTE: Broken since CADC removed this service (I can't find it any more)
        """

        cfhtCacheDir=self.cacheDir+os.path.sep+"CFHTLS"
        if os.path.exists(cfhtCacheDir) == False:
            os.makedirs(cfhtCacheDir)
        
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']
                
        outFileName=cfhtCacheDir+os.path.sep+name.replace(" ", "_")+".jpg"
        
        #http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/access/cut.html
        
        url="http://www1.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/community/CFHTLS-SG/cgi/cutcfhtls.pl?ra=$RA&dec=$DEC&size=$SIZE_ARCMIN&units=arcminutes&wide=true&deep=true&preview=colour"
        url=url.replace("$RA", str(RADeg))
        url=url.replace("$DEC", str(decDeg))
        url=url.replace("$SIZE_ARCMIN", str(self.configDict['plotSizeArcmin']))
        print url
        foundCutoutInfo=False
        if os.path.exists(outFileName) == False or refetch == True:
            response=urllib2.urlopen(url)
            lines=response.read()
            lines=lines.split("\n")
            for line in lines:
                if line.find("cutout preview") != -1:
                    foundCutoutInfo=True
                    break
            # This happens if outside of footprint
            if line == '<img src="/community/CFHTLS-SG/cgi/CFHTLScolcut.pl?field=&amp;section=" alt="cutout preview">':
                foundCutoutInfo=False
            if foundCutoutInfo == True:
                imageURL="http://www1.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/"+line[line.find("src=")+4:].split('"')[1]
                urllib.urlretrieve(imageURL, outFileName)
                #except:
                    #noDataPath=sourcery.__path__[0]+os.path.sep+"data"+os.path.sep+"noData.jpg"
                    #os.system("cp %s %s" % (noDataPath, outFileName))
        
        return foundCutoutInfo


    def fetchUnWISEImage(self, obj, refetch = False):
        """Retrieves unWISE W1, W2 .fits images and makes a colour .jpg.
        
        """
        
        wiseCacheDir=self.cacheDir+os.path.sep+"unWISE"
        if os.path.exists(wiseCacheDir) == False:
            os.makedirs(wiseCacheDir)
        
        name=obj['name']
        RADeg=obj['RADeg']
        decDeg=obj['decDeg']
        
        # 2.75" pixels in the unWISE images (max for query is 250 pixels though)
        sizePix=int(round(self.configDict['plotSizeArcmin']*60.0/2.75))
        
        outFileName=wiseCacheDir+os.path.sep+name.replace(" ", "_")+".jpg"
        targzPath=wiseCacheDir+os.path.sep+"wise.tar.gz"
        if os.path.exists(outFileName) == False or refetch == True:
            print "... fetching unWISE data for %s ..." % (name) 
            
            urllib.urlretrieve("http://unwise.me/cutout_fits?version=neo1&ra=%.6f&dec=%.6f&size=%d&bands=12" % (RADeg, decDeg, sizePix), targzPath)
            
            os.system("tar -zxvf %s" % (targzPath))
            wiseFiles=glob.glob("unwise-*-img-m.fits")
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
                                if inX >= 0 and inX < img[0].data.shape[1]-1 and inY >= 0 and inY < img[0].data.shape[0]-1:
                                    outData[y, x]=img[0].data[inY, inX]
                    if band == 'w1':
                        bClip={'wcs': wcs, 'data': outData}
                    elif band == 'w2':
                        rClip={'wcs': wcs, 'data': outData}
            else:
                wcs=astWCS.WCS(w1FileName)
                bImg=pyfits.open(w1FileName)
                rImg=pyfits.open(w2FileName)
                bClip={'wcs': wcs, 'data': bImg[0].data}
                rClip={'wcs': wcs, 'data': rImg[0].data}
        
            try:
                gClip={'wcs': rClip['wcs'], 'data': (rClip['data']+bClip['data'])/2.0}
            except:
                raise Exception, "W1, W2 images not same dimensions"

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
            plt.axes([0, 0, 1, 1])
            axes=[0., 0., 1.0, 1.0]
            plot=astPlots.ImagePlot([rData, gData, bData], bClip['wcs'], axes = axes, 
                                cutLevels = [cuts, cuts, cuts], axesFontSize = 18.0)
            plt.savefig(outFileName, dpi = dpi)
            plt.close()
            
            # Clean up
            os.remove(targzPath)
            os.system("rm unwise-*.fits*")

                
    @cherrypy.expose
    def makePlotFromJPEG(self, name, RADeg, decDeg, surveyLabel, plotNEDObjects = "false", plotSDSSObjects = "false", plotSourcePos = "false", plotXMatch = "false", plotContours = "false", noAxes = "false", clipSizeArcmin = None, gamma = 1.0):
        """Makes plot of .jpg image with coordinate axes and NED, SDSS objects overlaid.
        
        To test this:
        
        http://localhost:8080/makeSDSSPlot?name=XMMXCS%20J001737.5-005234.2&RADeg=4.406325&decDeg=-0.876192
        
        To test zoom:
        
        http://localhost:8080/makePlotFromJPEG?name=XMMXCS%20J001737.5-005234.2&RADeg=4.406325&decDeg=-0.876192&surveyLabel=SDSS&clipSizeArcmin=3.0
        
        """
        
        # Just in case they are passed as strings (e.g., if direct from the url)
        RADeg=float(RADeg)
        decDeg=float(decDeg)
        
        sizeDeg=self.configDict['plotSizeArcmin']/60.0
        
        # Load data
        inJPGPath=self.cacheDir+os.path.sep+surveyLabel+os.path.sep+name.replace(" ", "_")+".jpg"
        if os.path.exists(inJPGPath) == False:
            inJPGPath=sourcery.__path__[0]+os.path.sep+"data"+os.path.sep+"noData.jpg"
        
        im=Image.open(inJPGPath)
        data=np.array(im)
        data=np.power(data, 1.0/float(gamma))
        try:
            data=np.flipud(data)
            data=np.fliplr(data)
        except:
            #"... something odd about image (1d?) - aborting ..."
            return None
        
        R=data[:, :, 0]
        G=data[:, :, 1]
        B=data[:, :, 2]
        
        # HACK: for ACT maps, with huge pixels, we can get offsets in .jpg relative to original
        # So, if we have a .fits image, load that and use to set centre coords
        fitsFileName=inJPGPath.replace(".jpg", ".fits")
        if os.path.exists(fitsFileName) == True:
            hackWCS=astWCS.WCS(fitsFileName)
            CRVAL1, CRVAL2=hackWCS.getCentreWCSCoords()
        else:
            CRVAL1, CRVAL2=RADeg, decDeg
        # Make a WCS
        sizeArcmin=self.configDict['plotSizeArcmin']
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
        newHead['CDELT1']=xOutPixScale
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
        if noAxes == "false":
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
        
        if noAxes != "false":
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
            nedObjs=catalogTools.parseNEDResult(nedFileName)
            if len(nedObjs['RAs']) > 0:
                p.addPlotObjects(nedObjs['RAs'], nedObjs['decs'], 'nedObjects', objLabels = nedObjs['labels'],
                                    size = sizeDeg/40.0*3600.0, color = "#7cfc00")
    
        if plotSDSSObjects == "true":
            SDSSRedshifts=catalogTools.fetchSDSSRedshifts(self.sdssRedshiftsDir, name, RADeg, decDeg)
            if SDSSRedshifts != None:
                sdssRAs=[]
                sdssDecs=[]
                sdssLabels=[]
                sdssCount=0
                for sdssObj in SDSSRedshifts:
                    sdssCount=sdssCount+1
                    sdssRAs.append(sdssObj['RADeg'])
                    sdssDecs.append(sdssObj['decDeg'])
                    sdssLabels.append(str(sdssCount))
                if len(sdssRAs) > 0:
                    p.addPlotObjects(sdssRAs, sdssDecs, 'sdssObjects', objLabels = sdssLabels,
                                    size = sizeDeg/40.0*3600.0, symbol = 'box', color = "red")
                              
        if plotXMatch == "true":
            obj=self.sourceCollection.find_one({'name': name})
            xMatchRAs=[]
            xMatchDecs=[]
            xMatchLabels=[]
            for label in self.configDict['crossMatchCatalogLabels']:
                RAKey='%s_RADeg' % (label)
                decKey='%s_decDeg' % (label)
                if RAKey in obj.keys() and decKey in obj.keys():
                    xMatchRAs.append(obj[RAKey])
                    xMatchDecs.append(obj[decKey])
                    xMatchLabels.append(label)
            if len(xMatchRAs) > 0:
                p.addPlotObjects(xMatchRAs, xMatchDecs, 'xMatchObjects', objLabels = xMatchLabels,
                                 size = sizeDeg/40.0*3600.0, symbol = "diamond", color = 'cyan')
        
        if plotContours == "true":
            if 'contourImage' in self.configDict.keys() and self.configDict['contourImage'] != None:
                contourImg=pyfits.open(self.cacheDir+os.path.sep+self.configDict['contourImage']+os.path.sep+name.replace(" ", "_")+".fits")
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
        
        cherrypy.response.headers['Content-Type']="image/jpg"
        buf=StringIO.StringIO()
        plt.savefig(buf, dpi = 96, format = 'jpg')
        plt.close()
        
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
        
        return self.index()
    
        
    @cherrypy.expose
    def index(self):
        """Shows the table page.
        
        """
        
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
                font-family:sans-serif;
                font-size:small;
            }
            </style>
            <title>$TITLE</title>
        </head>
        
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.2/jquery.min.js"></script>

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
        
        <body style="font-family: sans-serif; vertical align: top; justify: full;" onload="hideShow()">
        <table cellpadding="4" cellspacing="0" border="0" style="text-align: left; width: 100%;">
            <tbody>
                <tr>
                    <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                        text-align: center; vertical-align: middle; font-size: 125%;">
                        $TITLE
                    </td>
                </tr>
            </tbody>
        </table>
        
        $META_DATA
        $COLOR_CODING
        
        <br>
        <form method="get" action="updateQueryParams">
        <fieldset>
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
        <input type="submit" class="f" name="queryApply" value="Apply">
        <input type="submit" class="f" name="queryReset" value="Reset"><br>
        </p>
        </fieldset>
        </form>
            
        <form method="post" action="changeTablePage">
        <div id="buttons">
            <input type="submit" class="f" name="nextButton" value=">">
            <input type="submit" class="f" name="prevButton" value="<">
            <div style="clear:both"></div><!-- Need this to have the buttons actually inside div#buttons -->
        </div>
        </form>
                
        <table frame=border cellspacing=0 cols=$TABLE_COLS rules=all border=2 width=100% align=center class=tablefont>
        <tbody>
            <tr style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 110%;">
            $TABLE_COL_NAMES
            </tr>
            <font size="1">
            $TABLE_DATA
            </font>
        </tbody>
        </table>

        <form method="post" action="changeTablePage">
        <div id="buttons">
            <input type="submit" class="f" name="nextButton" value=">">
            <input type="submit" class="f" name="prevButton" value="<">
            <div style="clear:both"></div><!-- Need this to have the buttons actually inside div#buttons -->
        </div>
        
        $DOWNLOAD_LINKS

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
        
        if 'indexTitle' in self.configDict.keys():
            html=html.replace("$TITLE", self.configDict['indexTitle'])
        else:
            html=html.replace("$TITLE", "Sourcery Database")
        
        # First need to apply query parameters here
        queryPosts=self.runQuery(queryRADeg, queryDecDeg, querySearchBoxArcmin, queryOtherConstraints)        

        # Then cut to number of rows to view as below
        viewPosts=queryPosts[cherrypy.session['viewTopRow']:cherrypy.session['viewTopRow']+self.tableViewRows]
        
        # Fill in query params
        html=html.replace("$QUERY_RADEG", queryRADeg)
        html=html.replace("$QUERY_DECDEG", queryDecDeg)
        html=html.replace("$QUERY_SEARCHBOXARCMIN", querySearchBoxArcmin)
        html=html.replace("$QUERY_OTHERCONSTRAINTS", queryOtherConstraints)
        html=html.replace("$OBJECT_TYPE_STRING", self.configDict['objectTypeString'])
        html=html.replace("$NUMBER_SOURCES", str(len(queryPosts)))
        html=html.replace("$HOSTED_STR", self.configDict['hostedBy'])
        html=html.replace("$CONSTRAINTS_HELP_LINK", "displayConstraintsHelp?")
        
        # Table columns - as well as defaults, add ones we query on
        displayColumns=[]+self.tableDisplayColumns
        displayColumnLabels=[]+self.tableDisplayColumnLabels
        displayColumnFormats=[]+self.tableDisplayColumnFormats
        operators=["<", ">", "=", "!"]
        logicalOps=[' and ', ' or ']
        for logOp in logicalOps:
            constraints=queryOtherConstraints.split(logOp)
            for c in constraints:
                for o in operators:
                    colName=c.split(o)[0].lstrip().rstrip()
                    if len(viewPosts) > 0 and colName in viewPosts[0].keys() and colName not in displayColumns:
                        displayColumns.append(colName)
                        displayColumnLabels.append(colName)
                        fieldTypeDict=self.fieldTypesCollection.find_one({'name': colName})
                        if fieldTypeDict['type'] == 'number':
                            displayColumnFormats.append('%.3f')
                        elif fieldTypeDict['type'] == 'text':
                            displayColumnFormats.append('%s')
                        else:
                            raise Exception, "unknown type for field '%s'" % (colName)
        
        columnHeadings=""
        for key in displayColumnLabels:
            columnHeadings=columnHeadings+"\n           <td><b>%s</b></td>" % (key)
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

        metaData="""<br><fieldset>
        <legend><span style='border: black 1px solid; color: gray; padding: 2px'>expand</span><b>Source List Information</b></legend>
        Total number of %s: %d (original source list: %d) %s
        <p>Original source list = %s</p>
        <p>%s</p>
        $NEWS
        </fieldset>""" % (self.configDict['objectTypeString'], len(queryPosts), self.sourceCollection.count(), latestNewsStr, 
                          os.path.split(self.configDict['catalogFileName'])[-1], commentsString)
        if 'newsItems' in self.configDict.keys():
            newsStr="<p>News:<ul>\n"
            for item in self.configDict['newsItems']:
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
        downloadLinkStr=urllib.quote_plus(downloadLinkStr, safe='&?=')
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
            
            # On the fly MongoDB matching
            #tagsDict=self.matchTags(obj)
                
            # Highlighting of rows - obviously, order matters here!
            bckColor="white"
            #if 'observed 2009B' in obj.keys() and obj['observed 2009B'] == True:
                #bckColor="darkgray"
                #bckKey='observed 2009B'
            #if 'SPT cluster' in obj.keys() and obj['SPT cluster'] == True:
                #bckColor="gold"
                #bckKey='SPT cluster'
            #if 'ACT 2008 cluster' in obj.keys() and obj['ACT 2008 cluster'] == True:
                #bckColor="deeppink"
                #bckKey='ACT 2008 cluster'
            #if bckColor not in usedBckColors and bckColor != "white":
                #usedBckColors.append(bckColor)
                #usedBckKeys.append(bckKey)
                
            # Row for each object in table
            rowString="<tr>\n"
            for key in displayColumns:
                htmlKey="$"+string.upper(key)+"_KEY"
                rowString=rowString+"   <td style='background-color: "+bckColor+";' align=center width=10%>"+htmlKey+"</td>\n"
            rowString=rowString+"</tr>\n"
            
            # Insert values - note name is special
            for key, fmt in zip(displayColumns, displayColumnFormats):
                if key in obj.keys():
                    try:
                        value=obj[key]
                    except:
                        raise Exception, "missing key %s" % (key)
                else:
                    # No entry in MongoDB tags yet
                    if fmt != "%s":
                        value=0.0
                    else:
                        value=""
                htmlKey="$"+string.upper(key)+"_KEY"
                if key == "name":
                    #linksDir="dummy"
                    linkURL="displaySourcePage?name=%s&clipSizeArcmin=%.2f" % (self.sourceNameToURL(obj['name']), self.configDict['plotSizeArcmin'])
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
                        rowString=rowString.replace(htmlKey, fmt % (value))
                    except:
                        raise Exception, """IndexError: check .config file tableDisplayColumns are actually in the .fits table, or for mixed '' "" inside [] """ 
                           
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
    def downloadCatalog(self, queryRADeg = "0:360", queryDecDeg = "-90:90", querySearchBoxArcmin = "",
                        queryOtherConstraints = "", fileFormat = "cat"):
        """Provide user with the current table view as a downloadable catalog.
        
        """
        
        posts=self.runQuery(queryRADeg, queryDecDeg, querySearchBoxArcmin, queryOtherConstraints)        

        #if not cherrypy.session.loaded: cherrypy.session.load()
    
        keysList, typeNamesList, descriptionsList=self.getFieldNamesAndTypes(excludeKeys = [])
        
        #posts=list(self.db.collection[cherrypy.session.id].find().sort('decDeg').sort('RADeg'))    # Current view, including classification info
        tabLength=len(posts)
        
        tab=atpy.Table()
        for key, typeName in zip(keysList, typeNamesList):
            if typeName == 'number':
                tab.add_column(atpy.Column(np.zeros(tabLength, dtype = float), str(key)))
            else:
                tab.add_column(atpy.Column(np.zeros(tabLength, dtype = 'S1000'), str(key)))
        
        for i in range(tabLength):
            row=tab[i]
            post=posts[i]
            for key in keysList:
                if key in post.keys():
                    row[key]=post[key]

        tab.table_name=self.configDict['catalogDownloadFileName']
        
        tmpFileName=tempfile.mktemp()
        if fileFormat == 'cat':
            tab.write(tmpFileName+".cat", format = 'ascii')
        elif fileFormat == 'fits':
            tab.write(tmpFileName+".fits", format = 'fits')
        elif fileFormat == 'reg':
            catalogTools.tab2DS9(tab, tmpFileName+".reg")
        
        cherrypy.response.headers['Content-Disposition']='attachment; filename="%s.%s"' % (self.configDict['catalogDownloadFileName'], fileFormat)
        
        # This may not be the nicest thing to do... but seems to work
        f=file(tmpFileName+"."+fileFormat, 'rb')
        
        return f
    
        
    def sourceNameToURL(self, name):
        """Replaces + and spaces in source names so that they will be valid urls.
        
        """
        return name.replace("+", "%2b").replace(" ", "%20")

        
    def URLToSourceName(self, url):
        """Replaces %20 and %2b in URLs with spaces and + signs.
        
        """
        return url.replace("%2b", "+").replace("%20", " ")


    def runQuery(self, queryRADeg, queryDecDeg, querySearchBoxArcmin, queryOtherConstraints):           
        """Runs a query, returns the posts found.
        
        """

        #if not cherrypy.session.loaded: cherrypy.session.load()
            
        ## Store results of query in another collection (empty it first if documents are in it)
        #self.db.collection[cherrypy.session.id].remove({})

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
        self.sourceCollection.ensure_index([("RADeg", pymongo.ASCENDING)])
        queryPosts=list(self.sourceCollection.find(queryDict).sort('decDeg').sort('RADeg'))        
        
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
            cherrypy.session['viewTopRow']=viewTopRow
        if prevButton:
            viewTopRow=cherrypy.session.get('viewTopRow')
            viewTopRow=viewTopRow-self.tableViewRows
            if viewTopRow < 0:
                viewTopRow=0
            cherrypy.session['viewTopRow']=viewTopRow
            
        return self.index()        
    
    
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
        excludeKeys=['RADeg', 'decDeg'] # because we handle differently
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
            for f, t, d in zip(self.configDict['fields'], self.configDict['fieldTypes'], 
                               self.configDict['fieldDescriptions']):
                keysList.append(f)
                typeNamesList.append(t)
                descList.append(d)
        
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
            if key in self.configDict['fields']:
                if self.configDict['fieldTypes'][self.configDict['fields'].index(key)] == 'number':
                    post[key]=float(kwargs[key])
                else:
                    post[key]=kwargs[key]
            else:
                post[key]=kwargs[key]
        post['lastUpdated']=datetime.date.today().isoformat()
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
                        
        return self.displaySourcePage(name, clipSizeArcmin = self.configDict['plotSizeArcmin'],
                                      imageType = imageType)


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
            if key in self.configDict['fields']:
                if self.configDict['fieldTypes'][self.configDict['fields'].index(key)] == 'number':
                    post[key]=float(tagsToInsertDict[key])
                else:
                    post[key]=tagsToInsertDict[key]
            else:
                post[key]=tagsToInsertDict[key]
        
        #print "How to avoid overwrites of things we shouldn't?"
        #IPython.embed()
        #sys.exit()
        
        self.tagsCollection.update({'_id': mongoDict['_id']}, {'$set': post}, upsert = False)
        
        # Update source collection too
        self.sourceCollection.update({'_id': obj['_id']}, {'$set': post}, upsert = False)
        
    
    @cherrypy.expose
    def displaySourcePage(self, name, imageType = 'SDSS', clipSizeArcmin = None, plotNEDObjects = "false", plotSDSSObjects = "false", plotSourcePos = "false", plotXMatch = "false", plotContours = "false", noAxes = "false", gamma = 1.0):
        """Retrieve data on a source and display source page, showing given image plot.
        
        This should have form controls somewhere for editing the assigned redshift, redshift type, redshift 
        source and candidate status (e.g., confirmed cluster, junk etc.).
        
        imageType should match to an entry in self.imageLabels
        
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

        name=self.URLToSourceName(name)
        
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
                        text-align: center; vertical-align: middle; font-size: 125%;">
                        $SOURCE_NAME
                    </td>
                    <!-- $NEXT_LINK_CODE -->
                </tr>
            </tbody>
        </table>
        
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.2/jquery.min.js"></script>

        <table frame=border cellspacing=0 cols=1 rules=None border=0 width=100%>
        <tbody>
        
        <tr>
            <td align=center id="imagePlot">
            </td>
        </tr>

        <tr><td align=center>$PLOT_CONTROLS</td></tr>

        <tr><td align=center>$TAG_CONTROLS</td></tr>

        $SPECTRUM_PLOT
        
        <tr>
            <td align=center>$NED_MATCHES_TABLE</td>
        </tr>

        <tr>
            <td align=center>$SDSS_MATCHES_TABLE</td>
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

        obj=self.sourceCollection.find_one({'name': name})
        mongoDict=self.matchTags(obj)
        
        # Controls for image zoom, plotting NED, SDSS, etc.       
        plotFormCode="""
        <script>
        function printValue(sliderID, textbox) {
            var x = document.getElementById(textbox);
            var y = document.getElementById(sliderID);
            x.value = y.value;
        }
        window.onload = function() { printValue('sizeSlider', 'sizeSliderValue'); printValue('gammaSlider', 'gammaSliderValue');}
        </script>

        <script type="text/javascript">
        
            $(document).ready(function() {
                    $.post('makePlotFromJPEG', 
                           {name: '$OBJECT_NAME',
                            RADeg: $OBJECT_RADEG,
                            decDeg: $OBJECT_DECDEG,
                            surveyLabel: $('input:radio[name=imageType]:checked').val(),
                            plotNEDObjects: $('input:checkbox[name=plotNEDObjects]').prop('checked'),
                            plotSDSSObjects: $('input:checkbox[name=plotSDSSObjects]').prop('checked'),
                            plotSourcePos: $('input:checkbox[name=plotSourcePos]').prop('checked'),
                            plotXMatch: $('input:checkbox[name=plotXMatch]').prop('checked'),
                            plotContours: $('input:checkbox[name=plotContours]').prop('checked'),
                            noAxes: $('input:checkbox[name=noAxes]').prop('checked'),
                            clipSizeArcmin: $("#clipSizeArcmin").val(),
                            gamma: $("#gamma").val()}, 
                            function(data) {
                                // directly insert the image
                                $("#imagePlot").html('<img src="data:image/jpg;base64,' + data + '" align="middle" border=2 width="$PLOT_DISPLAY_WIDTH_PIX"/>') ;
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
            
            $(function() {
                // When the form is submitted...
                $("#imageForm").submit(function() {   
                    $.post('makePlotFromJPEG', 
                           {name: '$OBJECT_NAME',
                            RADeg: $OBJECT_RADEG,
                            decDeg: $OBJECT_DECDEG,
                            surveyLabel: $('input:radio[name=imageType]:checked').val(),
                            plotNEDObjects: $('input:checkbox[name=plotNEDObjects]').prop('checked'),
                            plotSDSSObjects: $('input:checkbox[name=plotSDSSObjects]').prop('checked'),
                            plotSourcePos: $('input:checkbox[name=plotSourcePos]').prop('checked'),
                            plotXMatch: $('input:checkbox[name=plotXMatch]').prop('checked'),
                            plotContours: $('input:checkbox[name=plotContours]').prop('checked'),
                            noAxes: $('input:checkbox[name=noAxes]').prop('checked'),
                            clipSizeArcmin: $("#sizeSliderValue").val(),
                            gamma: $("#gammaSliderValue").val()}, 
                            function(data) {
                                // directly insert the image
                                $("#imagePlot").html('<img src="data:image/jpg;base64,' + data + '" align="middle" border=2 width="$PLOT_DISPLAY_WIDTH_PIX"/>') ;
                                //alert($('input:radio[name=imageType]:checked').val());
                           });
                    return false;
                });
            });

        </script>

        <form action="#" id="imageForm" method="post">        
        <fieldset style="width:80%">
        <legend><b>Image Controls</b></legend>
        <input name="name" value="$OBJECT_NAME" type="hidden">
        <p><b>Survey:</b> $IMAGE_TYPES</p>      
        <p><b>Plot:</b>
        <input type="checkbox" name="noAxes" value=1 $CHECKED_NOAXES>Remove coordinate axes
        <input type="checkbox" name="plotContours" value=1 $CHECKED_CONTOURS>Contours ($CONTOUR_IMAGE)
        <input type="checkbox" name="plotSourcePos" value=1 $CHECKED_SOURCEPOS>Source position
        <input type="checkbox" name="plotNEDObjects" value=1 $CHECKED_NED>NED objects
        <input type="checkbox" name="plotSDSSObjects" value=1 $CHECKED_SDSS>SDSS DR14 objects
        <input type="checkbox" name="plotXMatch" value=1 $CHECKED_XMATCH>Cross match objects
        </p>
        
        <label for="clipSizeArcmin">Image Size (arcmin)</label>
        <input id="sizeSlider" name="clipSizeArcmin" type="range" min="1.0" max="$MAX_SIZE_ARCMIN" step="0.5" value=$CURRENT_SIZE_ARCMIN onchange="printValue('sizeSlider','sizeSliderValue')">
        <input id="sizeSliderValue" type="text" size="2"/>
        
        <label for="gamma">Gamma</label>
        <input id="gammaSlider" name="gamma" type="range" min="0.2" max="3.0" step="0.2" value=$CURRENT_GAMMA onchange="printValue('gammaSlider','gammaSliderValue')">
        <input id="gammaSliderValue" type="text" size="2"/>
        
        <input type="submit" value="Apply">
        </fieldset>
        </form>
     
        """ 
        # Taken out: onChange="this.form.submit();" from all checkboxes ^^^
        plotFormCode=plotFormCode.replace("$PLOT_DISPLAY_WIDTH_PIX", str(self.configDict['plotDisplayWidthPix']))
        plotFormCode=plotFormCode.replace("$OBJECT_NAME", name)
        plotFormCode=plotFormCode.replace("$OBJECT_RADEG", str(obj['RADeg']))
        plotFormCode=plotFormCode.replace("$OBJECT_DECDEG", str(obj['decDeg']))
        plotFormCode=plotFormCode.replace("$OBJECT_SURVEY", imageType) 
        if 'contourImage' in self.configDict.keys() and self.configDict['contourImage'] != None:
            plotFormCode=plotFormCode.replace("$CONTOUR_IMAGE", self.configDict['contourImage'])
        else:
            plotFormCode=plotFormCode.replace("$CONTOUR_IMAGE", "None")
        
        imageTypesCode=""            
        for label in self.imageLabels:
            if label == imageType:
                imageTypesCode=imageTypesCode+'<input type="radio" name="imageType" value="%s" checked>%s\n' % (label, label)
            else:
                imageTypesCode=imageTypesCode+'<input type="radio" name="imageType" value="%s">%s\n' % (label, label)
        plotFormCode=plotFormCode.replace("$IMAGE_TYPES", imageTypesCode)
        
        if plotNEDObjects == "true":
            plotFormCode=plotFormCode.replace("$CHECKED_NED", " checked")
        else:
            plotFormCode=plotFormCode.replace("$CHECKED_NED", "")
        if plotSDSSObjects == "true":
            plotFormCode=plotFormCode.replace("$CHECKED_SDSS", " checked")
        else:
            plotFormCode=plotFormCode.replace("$CHECKED_SDSS", "")
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
        if noAxes == "true":
            plotFormCode=plotFormCode.replace("$CHECKED_NOAXES", " checked")
        else:
            plotFormCode=plotFormCode.replace("$CHECKED_NOAXES", "")
            
        plotFormCode=plotFormCode.replace("$MAX_SIZE_ARCMIN", str(self.configDict['plotSizeArcmin']))        
        if clipSizeArcmin == None:
            plotFormCode=plotFormCode.replace("$CURRENT_SIZE_ARCMIN", str(self.configDict['plotSizeArcmin']))
        else:
            plotFormCode=plotFormCode.replace("$CURRENT_SIZE_ARCMIN", str(clipSizeArcmin))
        
        plotFormCode=plotFormCode.replace("$CURRENT_GAMMA", str(gamma))
                
        # Tagging controls (including editable properties of catalog, e.g., for assigning classification or redshifts)
        if 'enableEditing' in self.configDict.keys() and self.configDict['enableEditing'] == True:
            tagFormCode="""
            <form method="post" action="updateMongoDB">    
            <input name="name" value="$OBJECT_NAME" type="hidden">
            <input name="returnURL" value=$RETURN_URL" type="hidden">
            <fieldset style="width:80%">
            <legend><b>Editing Controls</b></legend>
            <p><b>Classification:</b>
            $CLASSIFICATION_CONTROLS
            </p>
            <p><b>Fields:</b>
            $FIELD_CONTROLS
            </p>
            <input type="submit" class="f" value="Update">
            </fieldset>
            </form>
            """
            tagFormCode=tagFormCode.replace("$PLOT_DISPLAY_WIDTH_PIX", str(self.configDict['plotDisplayWidthPix']))
            tagFormCode=tagFormCode.replace("$OBJECT_NAME", name)
            tagFormCode=tagFormCode.replace("$RETURN_URL", cherrypy.url())
            if 'fields' in self.configDict.keys():
                fieldsCode=""
                for f, s in zip(self.configDict['fields'], self.configDict['fieldDisplaySizes']):
                    fieldsCode=fieldsCode+'<label for="%s">%s</label>\n' % (f, f)
                    fieldsCode=fieldsCode+'<input type="text" value="%s" name="%s" size=%d/>\n' % (str(mongoDict[f]), f, s)
                if 'lastUpdated' in mongoDict.keys():
                    lastUpdated=mongoDict['lastUpdated']
                else:
                    lastUpdated='-'
                fieldsCode=fieldsCode+'<p><label for = "lastUpdated"><b>Last Updated:</b></label>\n'
                fieldsCode=fieldsCode+'<input type="text" value="%s" name="lastUpdated" size=10 readonly/></p>\n' % (lastUpdated)
            tagFormCode=tagFormCode.replace('$FIELD_CONTROLS', fieldsCode)
            classificationsCode=""
            for c in self.configDict['classifications']:
                if c == mongoDict['classification']:
                    classificationsCode=classificationsCode+'<input type="radio" name="classification" value="%s" checked>%s\n' % (c, c)
                else:
                    classificationsCode=classificationsCode+'<input type="radio" onChange="this.form.submit();" name="classification" value="%s">%s\n' % (c, c)
            tagFormCode=tagFormCode.replace("$CLASSIFICATION_CONTROLS", classificationsCode)
        
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
        if 'enableEditing' in self.configDict.keys() and self.configDict['enableEditing'] == True:
            html=html.replace("$TAG_CONTROLS", tagFormCode)
        else:
            html=html.replace("$TAG_CONTROLS", "")
        html=html.replace("$SOURCE_NAME", name)
        html=html.replace("$SIZE_ARC_MIN", "%.1f" % (self.configDict['plotSizeArcmin']))
        html=html.replace("$HOSTED_STR", self.configDict['hostedBy'])

        #for label, caption in zip(self.imageLabels, self.imageCaptions):
            #if label == imageType:
                #html=html.replace("$CAPTION", "%s" % (caption))
                            
        # NED matches table
        self.fetchNEDInfo(obj)
        self.findNEDMatch(obj)
        nedFileName=self.nedDir+os.path.sep+obj['name'].replace(" ", "_")+".txt"
        nedObjs=catalogTools.parseNEDResult(nedFileName)
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
        if 'addSDSSRedshifts' in self.configDict.keys() and self.configDict['addSDSSRedshifts'] == True:
            SDSSRedshifts=catalogTools.fetchSDSSRedshifts(self.sdssRedshiftsDir, obj['name'], obj['RADeg'],
                                                          obj['decDeg'])
            sdssTable="""<br><table frame=border cellspacing=0 cols=7 rules=all border=2 width=80% align=center>
            <tbody>
            <tr>
                <th style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 110%;" colspan=5>
                    <b>SDSS Redshifts</b>
                </th>
            </tr>
            <tr>
                <td><b>ID</b></td>
                <td><b>RA</b></td>
                <td><b>Dec.</b></td>
                <td><b>z</b></td>
                <td><b>zWarning</b></td>                    
            </tr>
            """              
            sdssCount=0
            for sdssObj in SDSSRedshifts:
                sdssCount=sdssCount+1
                rowString="""<tr>
                    <td align=center width=10%>$ID</td>
                    <td align=center width=10%>$RA</td>
                    <td align=center width=10%>$DEC</td>
                    <td align=center width=10%>$REDSHIFT</td>
                    <td align=center width=10%>$Z_WARNING</td>
                </tr>
                """
                rowString=rowString.replace("$ID", "%d" % (sdssCount))
                rowString=rowString.replace("$RA", "%.5f" % (sdssObj['RADeg']))
                rowString=rowString.replace("$DEC", "%.5f" % (sdssObj['decDeg']))
                rowString=rowString.replace("$REDSHIFT", "%.3f" % (sdssObj['z']))
                rowString=rowString.replace("$Z_WARNING", "%s" % (sdssObj['zWarning']))
                sdssTable=sdssTable+rowString
            sdssTable=sdssTable+"</tbody></table>"
        else:
            sdssTable=""
        html=html.replace("$SDSS_MATCHES_TABLE", sdssTable)

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
            if f['name'] in obj.keys():
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
                propTable=propTable+rowString
        propTable=propTable+"</td></tr></tbody></table>"
        html=html.replace("$PROPERTIES_TABLE", propTable)
                
        return html   
    
    
    @cherrypy.expose
    def preprocess(self):
        """This re-runs pre-processing steps (e.g., NED matching, SDSS image fetching etc.).
        
        If the use specified their own imageDirs, then the .jpg images from these are constructed here
        
        Directories containing ready-made .jpgs can also be directly added into the cacheDir folder. So long
        as these have a corresponding entry in the .config file they will be picked up. We spin through
        those folders and also add 'image_<imageDirLabel>' tags in the MongoDB too. 
        
        """
        
        # For DES public DR1 images access (we have to stitch together tiles if necessary anyway, as DESCuts has problems with objects near edge)
        # This is the result of SELECT * FROM DR1_TILE_INFO and contains WCS info for all the tiles
        if 'addDESImage' in self.configDict.keys() and self.configDict['addDESImage'] == True:
            self.DESTileTab=atpy.Table().read(sourcery.__path__[0]+os.path.sep+"data"+os.path.sep+"DES_DR1_TILE_INFO.csv")
            # Building the WCS dict here takes ~21 sec
            self.DESWCSDict={}
            for row in self.DESTileTab:
                newHead=pyfits.Header()
                newHead['NAXIS']=2
                newHead['NAXIS1']=row['NAXIS1']
                newHead['NAXIS2']=row['NAXIS2']
                newHead['CTYPE1']=row['CTYPE1']
                newHead['CTYPE2']=row['CTYPE2']
                newHead['CRVAL1']=row['CRVAL1']
                newHead['CRVAL2']=row['CRVAL2']
                newHead['CRPIX1']=row['CRPIX1']
                newHead['CRPIX2']=row['CRPIX2']
                newHead['CD1_1']=row['CD1_1']
                newHead['CD1_2']=row['CD1_2']    
                newHead['CD2_1']=row['CD2_1']    
                newHead['CD2_2']=row['CD2_2']    
                newHead['CUNIT1']='DEG'
                newHead['CUNIT2']='DEG'
                self.DESWCSDict[row['TILENAME']]=astWCS.WCS(newHead.copy(), mode = 'pyfits')  
        
        # In the CFHT dir, we keep a file that lists objects that don't have data
        # Saves us pinging CFHT servers again if we rerun
        cfhtCacheDir=self.cacheDir+os.path.sep+"CFHTLS"
        if os.path.exists(cfhtCacheDir) == False:
            os.makedirs(cfhtCacheDir)
        failsFileName=cfhtCacheDir+os.path.sep+"objectsNotFound.txt"
        CFHTFailsList=[]
        if os.path.exists(failsFileName) == True:
            inFile=file(failsFileName, "r")
            lines=inFile.readlines()
            inFile.close()
            for line in lines:
                CFHTFailsList.append(line.replace("\n", ""))
        
        # We need to do this to avoid hitting 32 Mb limit below when using large databases
        self.sourceCollection.ensure_index([("RADeg", pymongo.ASCENDING)])

        cursor=self.sourceCollection.find(no_cursor_timeout = True).sort('decDeg').sort('RADeg')
        for obj in cursor:

            print ">>> Fetching data to cache for object %s" % (obj['name'])            
            self.fetchNEDInfo(obj)
            if self.configDict['addSDSSRedshifts'] == True:
                catalogTools.fetchSDSSRedshifts(self.sdssRedshiftsDir, obj['name'], obj['RADeg'], obj['decDeg'])
            if self.configDict['addSDSSImage'] == True:
                self.fetchSDSSImage(obj)
            if self.configDict['addDESImage'] == True:
                self.fetchDESImage(obj)
            if self.configDict['addPS1Image'] == True:
                self.fetchPS1Image(obj)
            if self.configDict['addPS1IRImage'] == True:
                self.fetchPS1IRImage(obj)
            if self.configDict['addCFHTLSImage'] == True:
                if obj['name'] not in CFHTFailsList:
                    CFHTResult=self.fetchCFHTLSImage(obj)
                    if CFHTResult == False:
                        CFHTFailsList.append(obj['name'])
            if self.configDict['addUnWISEImage'] == True:
                try:
                    self.fetchUnWISEImage(obj)
                except:
                    print("... problem with UnWISE image for %s - skipping ..." % (obj['name']))
            if 'skyviewLabels' in self.configDict.keys():
                for surveyString, label in zip(self.configDict['skyviewSurveyStrings'], self.configDict['skyviewLabels']):
                    self.fetchSkyviewJPEG(obj['name'], obj['RADeg'], obj['decDeg'], surveyString, label)         
        cursor.close()
        
        # Update CFHT fails list
        outFile=file(failsFileName, "w")
        for objName in CFHTFailsList:
            outFile.write(objName+"\n")
        outFile.close()
        
        # Make .jpg images from user-supplied .fits images
        if 'imageDirs' in self.configDict.keys():
            self.makeImageDirJPEGs()
            
        # Now spin through cache imageDirs and add 'image_<imageDirLabel>' tags
        print ">>> Adding image_<imageDirLabel> tags to MongoDB ..."
        minSizeBytes=40000
        imageDirs=glob.glob(self.cacheDir+os.path.sep+"*")
        for imageDir in imageDirs:
            label=os.path.split(imageDir)[-1]
            print "... %s ..." % (label)
            fileNames=glob.glob(imageDir+os.path.sep+"*.jpg")
            for f in fileNames:
                
                # image size check: don't include SDSS if image size is tiny as no data
                skipImage=False
                if os.stat(f).st_size < minSizeBytes and label == 'SDSS':
                    skipImage=True
                
                objName=os.path.split(f)[-1].replace(".jpg", "")
                namesToTry=[objName, objName.replace("_", " ")]
                obj=None
                for n in namesToTry:
                    obj=self.sourceCollection.find_one({'name': n})
                    if obj != None:
                        break
                
                if obj != None and skipImage == False:
                    post={'image_%s' % (label): 1}
                    self.sourceCollection.update({'_id': obj['_id']}, {'$set': post}, upsert = False)
                    if self.fieldTypesCollection.find_one({'name': 'image_%s' % (label)}) == None:
                        keysList, typeNamesList, descriptionsList=self.getFieldNamesAndTypes()
                        fieldDict={}
                        fieldDict['name']='image_%s' % (label)
                        fieldDict['type']='number'
                        fieldDict['description']='1 if object has image in the database; 0 otherwise'
                        fieldDict['index']=len(keysList)+1
                        self.fieldTypesCollection.insert(fieldDict)


    def makeImageDirJPEGs(self):
        """Actual makes .jpg images from .fits images in given directories. We figure out which image to use
        from spinning through the headers. 
        
        For XCS, there may be more than one image... will need to think how to handle that. For now we will 
        take the first match.
        
        If there is only one image in a directory (like an ACT map say), and all objects fall in it, we will
        flag to make a clickable map page.
        
        For tracking e.g. follow-up (and using it), we add a key to object if it has an imageDir image,
        with name image_<imageDirLabel> (e.g., 'image_NICFPS-Ks'). So would be able to search on
        'image_NICFPS-Ks = 1'
                
        """
        print ">>> Making imageDir .jpgs ..."
        for imageDir, label, colourMap, sizePix, minMaxRadiusArcmin, scaling in zip(
                 self.configDict['imageDirs'], 
                 self.configDict['imageDirsLabels'],
                 self.configDict['imageDirsColourMaps'],
                 self.configDict['imageDirsSizesPix'],
                 self.configDict['imageDirsMinMaxRadiusArcmin'],
                 self.configDict['imageDirsScaling']):
            
            print "... %s ..." % (label)

            if 'skipMakingNewImages' in self.configDict.keys() and label in self.configDict['skipMakingNewImages']:
                print "... WARNING: skipMakingNewImages enabled for %s ..." % (label)
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
                pickleFile=file(pickleFileName, "rb")
                unpickler=pickle.Unpickler(pickleFile)
                headerDict=unpickler.load()
                pickleFile.close()
            else:
                headerDict={}
            # Takes ~160 sec to build the first time for ~10,000 images, ~2 sec for subsequent runs
            print "... building headerDict pickle for .fits images under %s/ ..." % (imageDir)
            t0=time.time()
            origLength=len(headerDict.keys())
            for imgFileName in imgList:
                if imgFileName not in headerDict.keys():
                    wcs=astWCS.WCS(imgFileName)
                    headerDict[imgFileName]=wcs.header.copy()
            t1=time.time()
            # Write pickled headerDict, in case it was updated
            if len(headerDict.keys()) > origLength:
                print "... writing updated headerDict pickle to %s/ ..." % (imageDir)
                pickleFile=file(pickleFileName, "wb")
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
            objList=list(self.sourceCollection.find().sort('decDeg').sort('RADeg'))

            for obj in objList:
                outFileName=outDir+os.path.sep+obj['name'].replace(" ", "_")+".jpg"
                
                if os.path.exists(outFileName) == True:
                    print "... image for %s exists..." % (obj['name'])
                else:
                    print "... making image for %s ..." % (obj['name'])
                                        
                    for imgFileName in imgList:
                        
                        wcs=wcsDict[imgFileName]
                                                
                        data=None
                        # coordsAreInImage sometimes gives spurious results, not clear why...
                        # Replacement with below works - need to check and fix in astWCS
                        pixCoords=wcs.wcs2pix(obj['RADeg'], obj['decDeg'])
                        if pixCoords[0] >= 0 and pixCoords[0] < wcs.header['NAXIS1'] and pixCoords[1] >= 0 and pixCoords[1] < wcs.header['NAXIS2']:                           
                                
                            # Add to mongodb - we now do this in preprocess
                            #post={'image_%s' % (label): 1}
                            #self.sourceCollection.update({'_id': obj['_id']}, {'$set': post}, upsert = False)
                            #if self.fieldTypesCollection.find_one({'name': 'image_%s' % (label)}) == None:
                                #keysList, typeNamesList=self.getFieldNamesAndTypes()
                                #fieldDict={}
                                #fieldDict['name']='image_%s' % (label)
                                #fieldDict['type']='number'
                                #fieldDict['index']=len(keysList)+1
                                #self.fieldTypesCollection.insert(fieldDict)

                            if 'contourImage' in self.configDict.keys() and self.configDict['contourImage'] == label:
                                fitsOutFileName=outFileName.replace(".jpg", ".fits")
                            else:
                                fitsOutFileName=None
                            if fitsOutFileName != None and os.path.exists(fitsOutFileName) == False:
                                if data == None:
                                    img=pyfits.open(imgFileName)
                                    data=img[0].data
                                clip=astImages.clipImageSectionWCS(data, wcs, obj['RADeg'], obj['decDeg'],
                                                                self.configDict['plotSizeArcmin']/60.0)
                                astImages.saveFITS(fitsOutFileName, clip['data'], clip['wcs'])

                            if os.path.exists(outFileName) == False:
                                
                                if np.any(data) == None:
                                    img=pyfits.open(imgFileName)
                                    data=img[0].data
                                clip=astImages.clipImageSectionWCS(data, wcs, obj['RADeg'], obj['decDeg'],
                                                                self.configDict['plotSizeArcmin']/60.0)

                                # Try to pick sensible cut levels
                                # Min-Max scaling
                                # Should probably stick with this, but also add log option for optical
                                if scaling == 'auto' and minMaxRadiusArcmin != None:
                                    clip['data']=catalogTools.byteSwapArr(clip['data'])
                                    # Avoid cython type troubles
                                    if clip['data'].dtype != np.float32:
                                        clip['data']=np.array(clip['data'], dtype = np.float32)
                                    rMap=sourceryCython.makeDegreesDistanceMap(clip['data'], clip['wcs'], obj['RADeg'], obj['decDeg'], 100.0)
                                    minMaxData=clip['data'][np.less(rMap, minMaxRadiusArcmin/60.0)]
                                    cuts=[clip['data'].min(), clip['data'].max()]
                                else:
                                    scaleMin, scaleMax=scaling.split(":")
                                    scaleMin=float(scaleMin)
                                    scaleMax=float(scaleMax)
                                    cuts=[scaleMin, scaleMax]
                                
                                # This should guard against picking up edges of images, if source position is not actually visible
                                # (e.g., XMM images)
                                if cuts[0] == 0 and cuts[1] == 0:
                                    continue
                                
                                dpi=96.0
                                f=plt.figure(figsize=(sizePix/dpi, sizePix/dpi), dpi = dpi)
                                plt.axes([0, 0, 1, 1])
                                #p=astPlots.ImagePlot(clip['data'], clip['wcs'], cutLevels = [cuts[0], cuts[1]], axesLabels = None, axes = [0., 0., 1.0, 1.0], interpolation = "none")
                                plt.imshow(clip['data'], interpolation = "none", origin = 'lower', 
                                            cmap = colourMap, norm = plt.Normalize(cuts[0], cuts[1]))
                                try:
                                    plt.savefig(outFileName, dpi = dpi)
                                except:
                                    raise Exception, "if you see this, you probably need to update PIL/Pillow"
                                plt.close()
                                
                                #IPython.embed()
                                #sys.exit()

        
        
        
        
        
