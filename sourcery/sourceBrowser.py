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
import atpy
import operator
import urllib
import urllib2
from astLib import *
import pyfits
import numpy
import pylab
import catalogTools
import sourcery
import sys
import time
import string
from PIL import Image
import copy
import StringIO
import tempfile
import cherrypy
import IPython
 
#-------------------------------------------------------------------------------------------------------------
class SourceBrowser(object):
    
    def __init__(self, configFileName, preprocess = False):
        
        # Parse config file
        self.parseConfig(configFileName)
        
        if 'skyviewPath' in self.configDict.keys():
            self.skyviewPath=self.configDict['skyviewPath']
        else:
            # Default - if we have run the sourcery_fetch_skyview script
            self.skyviewPath=os.environ['HOME']+os.path.sep+".sourcery"+os.path.sep+"skyview.jar"

        # Parse catalog
        self.tab=atpy.Table(self.configDict['catalogFileName'])
        self.tab.sort(["RADeg", "decDeg"])
        
        # In case we don't have a 'name' column, we relabel a given column
        if 'nameColumn' in self.configDict.keys() and self.configDict['nameColumn'] != "":
            self.tab.rename_column(self.configDict['nameColumn'], 'name')
                    
        # Column to display info
        # Table pages
        self.tableDisplayColumns=['name', 'RADeg', 'decDeg']+self.configDict['tableDisplayColumns']
        self.tableDisplayColumnLabels=['Name', 'R.A. (degrees)', 'Dec. (degrees)']+self.configDict['tableDisplayColumns']
        self.tableDisplayColumnFormats=['%s', '%.6f', '%.6f']+self.configDict['tableDisplayColumnFormats']
        # Source pages
        self.sourceDisplayColumns=list(self.tab.keys())
        #self.sourceDisplayColumnLabels=['Name', 'R.A. (degrees)', 'Dec. (degrees)']+self.configDict['sourceDisplayColumnLabels']
        #self.sourceDisplayColumnFormats=['%s', '%.6f', '%.6f']+self.configDict['sourceDisplayColumnFormats']
        
        # Add NED match columns - we will fill on the fly...
        if 'addNEDMatches' in self.configDict.keys() and self.configDict['addNEDMatches'] == True:
            self.tab.add_column('NED_name', ['_______________']*len(self.tab))
            self.tab['NED_name']=None
            self.tab.add_column('NED_z', [numpy.nan]*len(self.tab))
            self.tab.add_column('NED_distArcmin', [numpy.nan]*len(self.tab))
            self.tab.add_column('NED_RADeg', [numpy.nan]*len(self.tab))
            self.tab.add_column('NED_decDeg', [numpy.nan]*len(self.tab))
            self.tableDisplayColumns=self.tableDisplayColumns+["NED_name", "NED_z"]
            self.tableDisplayColumnLabels=self.tableDisplayColumnLabels+["NED Name", "NED z"]
            self.tableDisplayColumnFormats=self.tableDisplayColumnFormats+["%s", "%.3f"]
            self.sourceDisplayColumns=self.sourceDisplayColumns+["NED_name", "NED_z", "NED_RADeg", "NED_decDeg", "NED_distArcmin"]
            #self.sourceDisplayColumnLabels=self.sourceDisplayColumnLabels+["NED Name", "NED z", "NED R.A. (degrees)", "NED Dec. (degrees)", "Distance to NED object (arcmin)"]
            #self.sourceDisplayColumnFormats=self.sourceDisplayColumnFormats+["%s", "%.3f", "%.6f", "%.6f", "%.1f"]

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
            
        # We will generate images dynamically... here we set up info like labels and captions
        self.imageLabels=[]      # labels at top of each source page that allow us to select image to view
        self.imageCaptions=[]    # caption that goes under image shown on the source pages
        # SDSS colour .jpgs
        if "addSDSSImage" in self.configDict.keys() and self.configDict['addSDSSImage'] == True:
            label="SDSS"
            self.imageLabels.append(label)
            self.imageCaptions.append("%.1f' x %.1f' false color (g,r,i) SDSS DR8 image. The source position is marked with the white cross.<br>Objects marked with green circles are in NED; objects marked with red squares have SDSS DR12 spectroscopic redshifts." % (self.configDict['plotSizeArcmin'], self.configDict['plotSizeArcmin']))
        # Skyview images
        for label in self.configDict['skyviewLabels']:
            self.imageLabels.append(label)
            self.imageCaptions.append("%.1f' x %.1f' false color %s image. The source position is marked with the white cross.<br>Objects marked with green circles are in NED; objects marked with red squares have SDSS DR12 spectroscopic redshifts." % (self.configDict['plotSizeArcmin'], self.configDict['plotSizeArcmin'], label))
       
        # Pre-processing? 
        # We can run this on the webserver by going to localhost:8080/preprocess
        # We might want to prevent that and force it to run manually only...
        if preprocess == True:
            self.preprocess()

        # This sets size of table view - view is controlled with session variables
        self.tableViewRows=40
        
    
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
                    if items[0] == '"':
                        delim='"'
                    elif items[0] == "'":
                        delim="'"
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
                elif value[0] == '(':
                    items=value.replace("(", "").replace(")", "").split(",")
                    lst=[]
                    for i in items:
                        lst.append(float(i))
                    value=tuple(lst)
                elif value[0] == "'" or value[0] == '"':
                    value=str(value.replace("'", "").replace('"', ""))
                self.configDict[key]=value
        # Ugly, ugly hack
        self.configDict['plotSizeArcmin']=float(self.configDict['plotSizeArcmin'])
        self.configDict['plotSizePix']=1200.0   # Same as SDSS 
        
    
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


    def findNEDMatch(self, obj):
        """Checks if there is a NED match for obj. Uses a 2.5' matching radius.
        
        """
        nedFileName=self.nedDir+os.path.sep+obj['name'].replace(" ", "_")+".txt"
        nedObjs=catalogTools.parseNEDResult(nedFileName)
            
        # Flag matches against clusters - choose nearest one
        rMin=10000
        crossMatchRadiusDeg=2.5/60.0
        clusterMatch={}
        if len(nedObjs['RAs']) > 0:
            for i in range(len(nedObjs['RAs'])):
                ned=nedObjs
                if ned['sourceTypes'][i] == 'GClstr':
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
            obj['NED_z']=numpy.nan
            obj['NED_distArcmin']=numpy.nan
            obj['NED_RADeg']=numpy.nan
            obj['NED_decDeg']=numpy.nan
            
        
    def fetchSDSSRedshifts(self, name, RADeg, decDeg):
        """Queries SDSS for redshifts. 
        
        """
        url='http://skyserver.sdss3.org/dr10/en/tools/search/x_sql.aspx'
        outFileName=self.sdssRedshiftsDir+os.path.sep+"%s.csv" % (name.replace(" ", "_"))
        if os.path.exists(outFileName) == False:
            sql="""SELECT
            p.objid,p.ra,p.dec,p.r,
            s.specobjid,s.z, 
            dbo.fSpecZWarningN(s.zWarning) as warning,
            s.plate, s.mjd, s.fiberid
            FROM PhotoObj AS p
            JOIN SpecObj AS s ON s.bestobjid = p.objid
            WHERE 
            p.ra < %.6f+0.1 and p.ra > %.6f-0.1
            AND p.dec < %.6f+0.1 and p.dec > %.6f-0.1
            """ % (RADeg, RADeg, decDeg, decDeg)
            # Filter SQL so that it'll work
            fsql = ''
            for line in sql.split('\n'):
                fsql += line.split('--')[0] + ' ' + os.linesep;
            params=urllib.urlencode({'cmd': fsql, 'format': "csv"})
            response=urllib2.urlopen(url+'?%s' % (params))
            lines=response.read()
            lines=lines.split("\n")
            outFile=file(outFileName, "w")
            for line in lines:
                outFile.write(line+"\n")
            outFile.close()        
        else:
            inFile=file(outFileName, "r")
            lines=inFile.readlines()
            inFile.close()
            
        # Parse .csv into catalog
        SDSSRedshifts=[]
        if lines[0] == "No objects have been found\n":
            SDSSRedshifts=[]
        elif len(lines) > 1 and lines[1] == '"ERROR: Maximum 60 queries allowed per minute. Rejected query: SELECT \n':
            os.remove(outFileName)
            raise Exception, "Exceeded 60 queries/min on SDSS server. Take a breather and rerun (previous queries cached)."
        else:
            SDSSRedshifts=[]
            for line in lines[2:]: # first line (DR7) always heading, first two lines (DR10) always heading
                if len(line) > 3:
                    zDict={}
                    bits=line.replace("\n", "").split(",")
                    zDict['objID']=bits[0]
                    try:
                        zDict['RADeg']=float(bits[1])
                        zDict['decDeg']=float(bits[2])
                    except:
                        if len(lines) > 1 and lines[1].find('"ERROR: Maximum 60 queries allowed per minute. Rejected query: SELECT') != -1:
                            raise Exception, "Exceeded 60 queries/min on SDSS server. Take a breather and rerun nemo (previous queries cached)."
                        else:
                            print "Hmm. Not able to parse SDSS redshifts"
                            IPython.embed()
                            sys.exit()
                    zDict['rMag']=float(bits[3])
                    zDict['specObjID']=bits[4]
                    zDict['z']=float(bits[5])
                    zDict['zWarning']=bits[6]
                    zDict['plate']=bits[7]
                    zDict['mjd']=bits[8]
                    zDict['fiberID']=bits[9]
                    # Throw out stars/junk
                    if zDict['z'] > 0.02:
                        SDSSRedshifts.append(zDict)
        
        return SDSSRedshifts
                 

    def fetchSDSSDR8Image(self, obj, refetch = False):
        """Fetches the SDSS DR8 .jpg for the given image size using the casjobs webservice.
        
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
            try:
                urlString="http://skyservice.pha.jhu.edu/DR8/ImgCutout/getjpeg.aspx?ra="+str(RADeg)+"&dec="+str(decDeg)
                urlString=urlString+"&scale="+str(SDSSScale)+"&width="+str(int(SDSSWidth))+"&height="+str(int(SDSSWidth))
                urllib.urlretrieve(urlString, filename = outFileName)
            except:
                print "... WARNING: couldn't get SDSS DR8 image ..."
                outFileName=None
        
    
    @cherrypy.expose
    def makePlotFromJPEG(self, name, RADeg, decDeg, surveyLabel, plotNEDObjects = "False", plotSDSSObjects = "False", plotSourcePos = "False", clipSizeArcmin = None):
        """Makes plot of SDSS image with coordinate axes and NED, SDSS objects overlaid.
        
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
            return None
        
        im=Image.open(inJPGPath)
        data=numpy.array(im)
        try:
            data=numpy.flipud(data)
            data=numpy.fliplr(data)
        except:
            #"... something odd about image (1d?) - aborting ..."
            return None
        
        R=data[:, :, 0]
        G=data[:, :, 1]
        B=data[:, :, 2]
        
        # Make a WCS
        sizeArcmin=self.configDict['plotSizeArcmin']
        xSizeDeg, ySizeDeg=sizeArcmin/60.0, sizeArcmin/60.0
        xSizePix=R.shape[1]
        ySizePix=R.shape[0]
        xRefPix=xSizePix/2.0
        yRefPix=ySizePix/2.0
        xOutPixScale=xSizeDeg/xSizePix
        yOutPixScale=ySizeDeg/ySizePix
        cardList=pyfits.CardList()
        cardList.append(pyfits.Card('NAXIS', 2))
        cardList.append(pyfits.Card('NAXIS1', xSizePix))
        cardList.append(pyfits.Card('NAXIS2', ySizePix))
        cardList.append(pyfits.Card('CTYPE1', 'RA---TAN'))
        cardList.append(pyfits.Card('CTYPE2', 'DEC--TAN'))
        cardList.append(pyfits.Card('CRVAL1', RADeg))
        cardList.append(pyfits.Card('CRVAL2', decDeg))
        cardList.append(pyfits.Card('CRPIX1', xRefPix+1))
        cardList.append(pyfits.Card('CRPIX2', yRefPix+1))
        cardList.append(pyfits.Card('CDELT1', xOutPixScale))
        cardList.append(pyfits.Card('CDELT2', xOutPixScale))    # Makes more sense to use same pix scale
        cardList.append(pyfits.Card('CUNIT1', 'DEG'))
        cardList.append(pyfits.Card('CUNIT2', 'DEG'))
        newHead=pyfits.Header(cards=cardList)
        wcs=astWCS.WCS(newHead, mode='pyfits')
        
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

        cutLevels=[[R.min(), R.max()], [G.min(), G.max()], [B.min(), B.max()]]
                                               
        # Make plot
        fig=pylab.figure(figsize = self.configDict['figSize'])
        axes=[0.1,0.085,0.9,0.85]
        axesLabels="sexagesimal"
        p=astPlots.ImagePlot([R, G, B], wcs, cutLevels = cutLevels, title = name, axes = axes, 
                            axesLabels = axesLabels)
        if plotSourcePos == "True":
            p.addPlotObjects([RADeg], [decDeg], 'clusterPos', symbol='cross', size=sizeDeg/20.0*3600.0, color='white')
                
        if plotNEDObjects == "True":
            # We should already have the files for this from doing addNEDInfo earlier
            nedFileName=self.nedDir+os.path.sep+name.replace(" ", "_")+".txt"
            nedObjs=catalogTools.parseNEDResult(nedFileName)
            if len(nedObjs['RAs']) > 0:
                p.addPlotObjects(nedObjs['RAs'], nedObjs['decs'], 'nedObjects', objLabels = nedObjs['labels'],
                                    size = sizeDeg/40.0*3600.0, color = "#7cfc00")
    
        if plotSDSSObjects == "True":
            SDSSRedshifts=self.fetchSDSSRedshifts(name, RADeg, decDeg)
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
                                    
        #if contourImg != None:
            #p.addContourOverlay(contourImg[0].data, contourWCS, 'actContour', levels = contourLevels, width = 2,     
                                    #color = 'yellow', highAccuracy = False)
        
        #outFileName=plotsDir+os.path.sep+name.replace(" ", "_")+".jpg"
        #pylab.savefig(outFileName)
        #pylab.close()
        cherrypy.response.headers['Content-Type']="image/jpg"
        buf=StringIO.StringIO()
        pylab.savefig(buf, dpi = 96, format = 'jpg')
        pylab.close()
        
        return buf.getvalue()
     
     
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
            os.system("java -jar %s position=%.6f,%.6f output=%s size=%.3f pixels=%d survey=%s cache=%s/" % (self.skyviewPath, RADeg, decDeg, rootFileName, self.configDict['plotSizeArcmin']/60.0, self.configDict['plotSizePix'], RGBSurveysString, self.skyCacheDir))
            
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
            imData=numpy.array([RImg[0].data.transpose(), GImg[0].data.transpose(), BImg[0].data.transpose()])
            imData=imData.transpose()
            for i in range(imData.shape[2]):
                channel=imData[:, :, i]
                
                std=numpy.std(channel)
                med=numpy.median(channel)
                minLevel=med-lowSigmaCut*std
                
                # This is better
                freq, binEdges=numpy.histogram(channel, bins = (channel.shape[0]*channel.shape[1])/100.0)
                binCentres=binEdges[:-1]+(binEdges[1]-binEdges[0])/2.0
                minLevel=binCentres[freq.tolist().index(freq.max())]
                
                lowMask=numpy.less(channel, minLevel)
                channel=channel-(minLevel)
                channel[lowMask]=0
                maxLevel=med+highSigmaCut*std
                if maxLevel > channel.max():
                    maxLevel=channel.max()
                highMask=numpy.greater(channel, maxLevel)
                channel=channel/maxLevel+0.001
                channel[highMask]=1.0
                channel=numpy.log10(channel)
                channel=channel-channel.min()
                channel=channel/channel.max()
                imData[:, :, i]=channel
            
            dpi=96.0
            pylab.figure(figsize=(self.configDict['plotSizePix']/dpi, self.configDict['plotSizePix']/dpi), dpi = dpi)
            pylab.axes([0, 0, 1, 1])
            pylab.imshow(imData, interpolation="bilinear", origin='lower')
            pylab.savefig(outFileName, dpi = dpi)
            pylab.close()
            
            # Clean up .fits images
            if cleanUp == True:
                for i in range(1, 4):
                    toRemove=rootFileName.replace(".fits", "_%d.fits" % (i))
                    if os.path.exists(toRemove) == True:
                        os.remove(toRemove)
                    
        
    @cherrypy.expose
    def index(self):
        """Shows the table page.
        
        """
        
        # Session variables: where in the table are we looking, query constraints
        if not cherrypy.session.loaded: cherrypy.session.load()
        if 'viewTopRow' not in cherrypy.session:
            cherrypy.session['viewTopRow']=0
        if 'viewTab' not in cherrypy.session:
            cherrypy.session['viewTab']=copy.deepcopy(self.tab)
        if 'queryRADeg' not in cherrypy.session:
            cherrypy.session['queryRADeg']="0:360"
        if 'queryDecDeg' not in cherrypy.session:
            cherrypy.session['queryDecDeg']="-90:90"
        if 'querySearchBoxArcmin' not in cherrypy.session:
            cherrypy.session['querySearchBoxArcmin']=""
        if 'queryOtherConstraints' not in cherrypy.session:
            cherrypy.session['queryOtherConstraints']=""
        viewTab=cherrypy.session.get('viewTab')
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
            <title>Source Browser</title>
        </head>
        <body style="font-family: sans-serif; vertical align: top; justify: full;">
        <table cellpadding="4" cellspacing="0" border="0" style="text-align: left; width: 100%;">
            <tbody>
                <tr>
                    <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                        text-align: center; vertical-align: middle; font-size: 125%;">
                        Source Browser
                    </td>
                </tr>
            </tbody>
        </table>
        
        $META_DATA
        $COLOR_CODING
        
        <br>
        <form method="get" action="runQuery">
        <fieldset>
        <legend><b>Constraints</b></legend>
        <p>Enter coordinate ranges (e.g., 120:220) or set the search box length. Use negative R.A. values to wrap around 0 degrees (e.g., -60:60).</p>
        <p>
        <label for="queryRADeg">R.A. (degrees)</label>
        <input type="text" value="$QUERY_RADEG" name="queryRADeg"/>
        <label for="queryDecDeg">Dec. (degrees)</label>
        <input type="text" value="$QUERY_DECDEG" name="queryDecDeg"/>
        <label for="querySearchBoxArcmin">Search box length (arcmin)</label>
        <input type="text" value="$QUERY_SEARCHBOXARCMIN" name="querySearchBoxArcmin"/>
        </p>
        </p>
        <label for="queryOtherConstraints"><a href=$CONSTRAINTS_HELP_LINK target=new>Other constraints</a></label>
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
                
        <table frame=border cellspacing=0 cols=6 rules=all border=2 width=100% align=center>
        <tbody>
            <tr style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 110%;">"""
        for key in self.tableDisplayColumnLabels:
            templatePage=templatePage+"\n           <td><b>%s</b></td>" % (key)
        templatePage=templatePage+"""
            </tr>
            $TABLE_DATA
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
        
        # Fill in query params
        html=html.replace("$QUERY_RADEG", queryRADeg)
        html=html.replace("$QUERY_DECDEG", queryDecDeg)
        html=html.replace("$QUERY_SEARCHBOXARCMIN", querySearchBoxArcmin)
        html=html.replace("$QUERY_OTHERCONSTRAINTS", queryOtherConstraints)
        html=html.replace("$OBJECT_TYPE_STRING", self.configDict['objectTypeString'])
        html=html.replace("$NUMBER_SOURCES", str(len(viewTab)))
        html=html.replace("$HOSTED_STR", self.configDict['hostedBy'])
        html=html.replace("$CONSTRAINTS_HELP_LINK", "displayConstraintsHelp?")
        
        # Meta data
        READMEComment="""Matches to other catalogs (e.g. NED) listed on this page are within 2.5' radius of the 
        candidate position."""
        if 'catalogComments' not in self.configDict.keys():
            commentsString=READMEComment
        else:
            commentsString=self.configDict['catalogComments']+" "+READMEComment
        metaData="""<br><fieldset>
        <legend><b>Source List Information</b></legend>
        <p>Original source list = %s</p>
        <p>Total number of %s = %d (original source list: %d)</p>
        <p>%s</p>
        </fieldset>""" % (os.path.split(self.configDict['catalogFileName'])[-1], self.configDict['objectTypeString'], 
                          len(viewTab), len(self.tab), commentsString)
        html=html.replace("$META_DATA", metaData)        
        
        # Catalog download links
        shortCatalogName=self.configDict['catalogDownloadFileName']+".cat"
        shortFITSName=shortCatalogName.replace(".cat", ".fits")
        shortRegName=shortCatalogName.replace(".cat", ".reg")
        html=html.replace("$DOWNLOAD_LINKS", """<fieldset>
        <legend><b>Download Catalog</b></legend>
        <ul>
        <li><a href=downloadCatalog?fileFormat=cat>%s</a>   (plain text)</li>
        <li><a href=downloadCatalog?fileFormat=fits>%s</a>   (FITS table format)</li>
        <li><a href=downloadCatalog?fileFormat=reg>%s</a>   (DS9 region file)</li></ul>
        <p>Note that current constraints are applied to downloaded catalogs.</p>
        </fieldset><br>
        """ % (shortCatalogName, shortFITSName, shortRegName))
                
        tableData=""
        usedBckColors=[]
        usedBckKeys=[]
        viewTopRow=cherrypy.session.get('viewTopRow')
        viewBottomRow=viewTopRow+self.tableViewRows
        if viewBottomRow > len(viewTab):
            viewBottomRow=len(viewTab)
        for obj in viewTab[viewTopRow:viewBottomRow]:#viewTopRow:viewTopRow+self.tableViewRows]:
            
            # On the fly NED matching (we can cache these beforehand by running preprocess()
            self.fetchNEDInfo(obj)
            self.findNEDMatch(obj)

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
                
            # Row for each cluster in table
            rowString="<tr>\n"
            for key in self.tableDisplayColumns:
                htmlKey="$"+string.upper(key)+"_KEY"
                rowString=rowString+"   <td style='background-color: "+bckColor+";' align=center width=10%>"+htmlKey+"</td>\n"
            rowString=rowString+"</tr>\n"
            
            # Insert values - note name is special
            # Need to work out how to handle the clusterMatch column with Toby's stuff
            for key, fmt in zip(self.tableDisplayColumns, self.tableDisplayColumnFormats):
                htmlKey="$"+string.upper(key)+"_KEY"
                if key == "name":
                    #linksDir="dummy"
                    nameLink="<a href=\"%s\" target=new>%s</a>" % \
                        ("displaySourcePage?name=%s&clipSizeArcmin=%.2f" % (self.sourceNameToURL(obj['name']), self.configDict['plotSizeArcmin']), obj['name'])
                    rowString=rowString.replace(htmlKey, "%s" % (nameLink))
                elif key == "NED_name" and obj['NED_name'] != "None":
                    nedName=obj[key]
                    nedLinkURL="http://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES" % (nedName.replace("+", "%2B").replace(" ", "+"))
                    rowString=rowString.replace(htmlKey, "<a href=%s>%s</a>" % (nedLinkURL, nedName))
                elif key == "NED_z" and numpy.isnan(obj['NED_z']) == True:
                    rowString=rowString.replace(htmlKey, "-")
                else:
                    rowString=rowString.replace(htmlKey, fmt % (obj[key]))
                           
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
    def downloadCatalog(self, fileFormat = "cat"):
        """Provide user with the current table view as a downloadable catalog.
        
        """
        if not cherrypy.session.loaded: cherrypy.session.load()
        if 'viewTab' not in cherrypy.session:
            cherrypy.session['viewTab']=copy.deepcopy(self.tab)
        viewTab=cherrypy.session.get('viewTab')
        
        tmpFileName=tempfile.mktemp()
        if fileFormat == 'cat':
            viewTab.write(tmpFileName+".cat", type = 'ascii')
        elif fileFormat == 'fits':
            viewTab.write(tmpFileName+".fits", type = 'fits')
        elif fileFormat == 'reg':
            catalogTools.tab2DS9(viewTab, tmpFileName+".reg")
        
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


    @cherrypy.expose
    def runQuery(self, queryRADeg, queryDecDeg, querySearchBoxArcmin, queryOtherConstraints, queryApply = None, queryReset = None):            
        if not cherrypy.session.loaded: cherrypy.session.load()
        if queryApply:
            cherrypy.session['viewTab']=copy.deepcopy(self.tab)
            cherrypy.session['queryRADeg']=queryRADeg
            cherrypy.session['queryDecDeg']=queryDecDeg
            cherrypy.session['querySearchBoxArcmin']=querySearchBoxArcmin
            cherrypy.session['queryOtherConstraints']=queryOtherConstraints
            viewTab=cherrypy.session.get('viewTab')
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
            if RAMin >= 0:
                viewTab=viewTab.where(viewTab['RADeg'] > RAMin)
                viewTab=viewTab.where(viewTab['RADeg'] < RAMax)
            else:
                mask1=numpy.less(viewTab['RADeg'],  RAMax)
                mask2=numpy.greater(viewTab['RADeg']-360.0, RAMin) 
                viewTab=viewTab.where(numpy.logical_or(mask1, mask2))
            viewTab=viewTab.where(viewTab['decDeg'] > decMin)
            viewTab=viewTab.where(viewTab['decDeg'] < decMax)
            # Other constraints
            if queryOtherConstraints != "":
                result=self.applyOtherConstraints(viewTab, queryOtherConstraints)
                if type(result) == str:
                    return result   # Error message
                else:
                    viewTab=result
            cherrypy.session['viewTopRow']=0
            cherrypy.session['viewTab']=viewTab
        
        if queryReset:
            cherrypy.session['queryRADeg']="0:360"
            cherrypy.session['queryDecDeg']="-90:90"
            cherrypy.session['querySearchBoxArcmin']=""
            cherrypy.session['viewTab']=self.tab
            cherrypy.session['viewTopRow']=0
            cherrypy.session['queryOtherConstraints']=""
            
        return self.index()


    def applyOtherConstraints(self, tab, queryStr):
        """Returns a table with queryStr constraints applied, or returns an error message.
        
        """

        # Firstly, do nothing if the queryStr is empty
        if len(queryStr) == 0:
            return tab
        
        # Order matters for matching operators below (for... break)
        operators=['>=', '<=', '<', '>', '=']
        
        # Maybe we will only support 'and'
        constraints=queryStr.split(" and ")
        newTab=copy.deepcopy(tab)
        for c in constraints:
            foundOp=False
            for op in operators:
                if c.find(op) != -1:
                    foundOp=True
                    break
            if foundOp == True:
                bits=c.split(op)
                key=bits[0].lstrip().rstrip()
                value=bits[1].lstrip().rstrip()
                # All greater/less only make sense with numbers anyway, fail gracelessly if not floats at present
                if op == '>':
                    newTab=newTab.where(numpy.greater(newTab[key], float(value)))
                elif op == '>=':
                    newTab=newTab.where(numpy.greater_equal(newTab[key], float(value)))
                elif op == '<':
                    newTab=newTab.where(numpy.less(newTab[key], float(value)))
                elif op == '<=':
                    newTab=newTab.where(numpy.less_equal(newTab[key], float(value)))
                elif op == '=':
                    if newTab[key].dtype.name.find("string") != -1:
                        # Strip out "' if put there by the user
                        value=str(value.replace("'", "").replace('"', ''))
                    else:
                        value=float(value)
                    newTab=newTab.where(newTab[key] == value)
            else:
                return "Error: unrecognised operator in constraint '%s'" % (c)
        
        if newTab != None:
            return newTab
        else:
            # Or we should return error?
            return tab
    
    
    @cherrypy.expose
    def changeTablePage(self, nextButton = None, prevButton = None):
        """Changes the viewed table page.
        
        """
        if not cherrypy.session.loaded: cherrypy.session.load()
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
        if not cherrypy.session.loaded: cherrypy.session.load()
        if 'viewTab' not in cherrypy.session:
            cherrypy.session['viewTab']=copy.deepcopy(self.tab)
        viewTab=cherrypy.session.get('viewTab')
                
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
        Constraints can be placed on the source list columns listed below. Each constraint should be
        separated by 'and', e.g., <i>"redshift >= 0 and redshift < 0.4"</i> ('or' is not supported yet; comparison
        operators which are understood are '<', '>', '>=', '<=' and '=').
        </p>
        <br>
        <table frame=border cellspacing=0 cols=2 rules=all border=2 width=80% align=center>
        <tbody>
            <tr style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                    text-align: center; vertical-align: middle; font-size: 110%;">
            <td><b>Column</td>
            <td><b>Type</b></td>
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
        for key in viewTab.keys():
            # Row for each column in table
            a=viewTab[key]
            rowString="<tr>\n"                
            rowString=rowString+"   <td style='background-color: "+bckColor+";' align=center width=10%><b>"+key+"</b></td>\n"
            rowString=rowString+"   <td style='background-color: "+bckColor+";' align=center width=10%>"+a.dtype.str+" ("+a.dtype.name+")"+"</td>\n"
            rowString=rowString+"</tr>\n"                           
            tableData=tableData+rowString
        html=html.replace("$TABLE_DATA", tableData)
        
        return html   

    
    @cherrypy.expose
    def displaySourcePage(self, name, imageType = 'SDSS', clipSizeArcmin = None, plotNEDObjects = "False", plotSDSSObjects = "False", plotSourcePos = "False"):
        """Retrieve data on a source and display source page, showing given image plot.
        
        This should have form controls somewhere for editing the assigned redshift, redshift type, redshift 
        source and candidate status (e.g., confirmed cluster, junk etc.).
        
        imageType should match to an entry in self.imageLabels
        
        """
        if not cherrypy.session.loaded: cherrypy.session.load()
        if 'viewTab' not in cherrypy.session:
            cherrypy.session['viewTab']=copy.deepcopy(self.tab)
        viewTab=cherrypy.session.get('viewTab')
        
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
                    $PREV_LINK_CODE
                    <td style="background-color: rgb(0, 0, 0); font-family: sans-serif; color: rgb(255, 255, 255); 
                        text-align: center; vertical-align: middle; font-size: 125%;">
                        $SOURCE_NAME
                    </td>
                    $NEXT_LINK_CODE
                </tr>
            </tbody>
        </table>
        
        <table frame=border cellspacing=0 cols=1 rules=None border=0 width=100%>
        <tbody>
        
        <tr>
            <td align=center><a href="$IMAGE_PATH"><img src="$IMAGE_PATH" align="middle" border=2 width="$PLOT_DISPLAY_WIDTH_PIX"></a>
            </td>
        </tr>
        <tr><td align=center>$PLOT_CONTROLS</td></tr>

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

        objTabIndex=numpy.where(viewTab['name'] == name)[0][0]
        obj=viewTab[objTabIndex]
        
        # Controls for image zoom, plotting NED, SDSS, etc.       
        plotFormCode="""
        <script>
        function printValue(sliderID, textbox) {
            var x = document.getElementById(textbox);
            var y = document.getElementById(sliderID);
            x.value = y.value;
        }
        window.onload = function() { printValue('sizeSlider', 'sizeSliderValue');}
        </script>

        <form method="get" action="displaySourcePage">        
        <fieldset style="width:$PLOT_DISPLAY_WIDTH_PIXpx">
        <legend><b>Image Controls</b></legend>
        <input name="name" value="$OBJECT_NAME" type="hidden">
        <p>Survey: $IMAGE_TYPES</p>      
        <p>Show:
        <input type="checkbox" onChange="this.form.submit();" name="plotSourcePos" value="True"$CHECKED_SOURCEPOS>Source position
        <input type="checkbox" onChange="this.form.submit();" name="plotNEDObjects" value="True"$CHECKED_NED>NED objects
        <input type="checkbox" onChange="this.form.submit();" name="plotSDSSObjects" value="True"$CHECKED_SDSS>SDSSDR10 objects
        </p>
        <label for="clipSizeArcmin">Image Size (arcmin)</label>
        <input id="sizeSlider" name="clipSizeArcmin" type="range" min="1.0" max="$MAX_SIZE_ARCMIN" step="0.5" value=$CURRENT_SIZE_ARCMIN onchange="printValue('sizeSlider','sizeSliderValue')">
        <input id="sizeSliderValue" type="text" size="2"/>
        <input type="submit" class="f" value="Apply">
        </fieldset>
        </form>
        """ 
        plotFormCode=plotFormCode.replace("$PLOT_DISPLAY_WIDTH_PIX", str(self.configDict['plotDisplayWidthPix']))
        plotFormCode=plotFormCode.replace("$OBJECT_NAME", name)
        
        imageTypesCode=""
        for label in self.imageLabels:
            if label == imageType:
                imageTypesCode=imageTypesCode+'<input type="radio" name="imageType" value="%s" checked>%s\n' % (label, label)
            else:
                imageTypesCode=imageTypesCode+'<input type="radio" onChange="this.form.submit();" name="imageType" value="%s">%s\n' % (label, label)
        plotFormCode=plotFormCode.replace("$IMAGE_TYPES", imageTypesCode)
        
        if plotNEDObjects == "True":
            plotFormCode=plotFormCode.replace("$CHECKED_NED", " checked")
        else:
            plotFormCode=plotFormCode.replace("$CHECKED_NED", "")
        if plotSDSSObjects == "True":
            plotFormCode=plotFormCode.replace("$CHECKED_SDSS", " checked")
        else:
            plotFormCode=plotFormCode.replace("$CHECKED_SDSS", "")
        if plotSourcePos == "True":
            plotFormCode=plotFormCode.replace("$CHECKED_SOURCEPOS", " checked")
        else:
            plotFormCode=plotFormCode.replace("$CHECKED_SOURCEPOS", "")
            
        plotFormCode=plotFormCode.replace("$MAX_SIZE_ARCMIN", str(self.configDict['plotSizeArcmin']))        
        if clipSizeArcmin == None:
            plotFormCode=plotFormCode.replace("$CURRENT_SIZE_ARCMIN", str(self.configDict['plotSizeArcmin']))
        else:
            plotFormCode=plotFormCode.replace("$CURRENT_SIZE_ARCMIN", str(clipSizeArcmin))
        
        # Directly serving .jpg image
        if imageType == 'SDSS':
            self.fetchSDSSDR8Image(obj)
        else:
            skyviewIndex=self.configDict['skyviewLabels'].index(imageType)
            self.fetchSkyviewJPEG(obj['name'], obj['RADeg'], obj['decDeg'], self.configDict['skyviewSurveyStrings'][skyviewIndex], imageType)
        if clipSizeArcmin == None:
            imagePath="makePlotFromJPEG?name=%s&RADeg=%.6f&decDeg=%.6f&surveyLabel=%s&plotNEDObjects=%s&plotSDSSObjects=%s&plotSourcePos=%s" % (self.sourceNameToURL(obj['name']), obj['RADeg'], obj['decDeg'], imageType, plotNEDObjects, plotSDSSObjects, plotSourcePos)
        else:
            imagePath="makePlotFromJPEG?name=%s&RADeg=%.6f&decDeg=%.6f&surveyLabel=%s&clipSizeArcmin=%.3f&plotNEDObjects=%s&plotSDSSObjects=%s&plotSourcePos=%s" % (self.sourceNameToURL(obj['name']), obj['RADeg'], obj['decDeg'], imageType, float(clipSizeArcmin), plotNEDObjects, plotSDSSObjects, plotSourcePos)
            
        html=templatePage
        html=html.replace("$PLOT_CONTROLS", plotFormCode)
        html=html.replace("$SOURCE_NAME", name)
        html=html.replace("$IMAGE_PATH", imagePath)        
        html=html.replace("$SIZE_ARC_MIN", "%.1f" % (self.configDict['plotSizeArcmin']))
        html=html.replace("$HOSTED_STR", self.configDict['hostedBy'])

        #for label, caption in zip(self.imageLabels, self.imageCaptions):
            #if label == imageType:
                #html=html.replace("$CAPTION", "%s" % (caption))
                        
        # Previous and next object page links
        if objTabIndex > 0:
            prevObjIndex=objTabIndex-1
        else:
            prevObjIndex=None
        if objTabIndex < len(viewTab)-1:
            nextObjIndex=objTabIndex+1
        else:
            nextObjIndex=None
        if prevObjIndex != None:
            prevObj=viewTab[prevObjIndex]            
            prevLinkCode="""<td style="background-color: rgb(0, 0, 0); font-family: sans-serif; 
            color: rgb(255, 255, 255); text-align: center; vertical-align: middle; font-size: 125%;">
            <a href=$PREV_LINK><b><<</b></a></td>
            """ 
            prevLinkCode=prevLinkCode.replace("$PREV_LINK", "displaySourcePage?name=%s&imageType=%sclipSizeArcmin=%.2f&plotNEDObjects=%s&plotSDSSObjects=%s&plotSourcePos=%s" % (self.sourceNameToURL(prevObj['name']), self.sourceNameToURL(imageType), float(clipSizeArcmin), plotNEDObjects, plotSDSSObjects, plotSourcePos))
        else:
            prevLinkCode=""
        if nextObjIndex != None:
            nextObj=viewTab[nextObjIndex]
            nextLinkCode="""<td style="background-color: rgb(0, 0, 0); font-family: sans-serif; 
            color: rgb(255, 255, 255); text-align: center; vertical-align: middle; font-size: 125%;">
            <a href=$NEXT_LINK><b>>></b></a></td>
            """
            nextLinkCode=nextLinkCode.replace("$NEXT_LINK", "displaySourcePage?name=%s&imageType=%s&clipSizeArcmin=%.2f&plotNEDObjects=%s&plotSDSSObjects=%s&plotSourcePos=%s" % (self.sourceNameToURL(nextObj['name']), self.sourceNameToURL(imageType), float(clipSizeArcmin), plotNEDObjects, plotSDSSObjects, plotSourcePos))
        else:
            nextLinkCode=""
        html=html.replace("$PREV_LINK_CODE", prevLinkCode)
        html=html.replace("$NEXT_LINK_CODE", nextLinkCode)
    
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
                <td><b>R.A.</b></td>
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
            SDSSRedshifts=self.fetchSDSSRedshifts(obj['name'], obj['RADeg'], obj['decDeg'])
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
                <td><b>R.A.</b></td>
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
        for pkey in self.sourceDisplayColumns:
            if pkey in self.tab.keys():
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
        """This re-runs pre-processing steps (e.g., NED matching, SDSS image fetching etc.)
        
        """
        for obj in self.tab:
            print ">>> Fetching data to cache for object %s" % (obj['name'])
            self.fetchNEDInfo(obj)
            self.fetchSDSSRedshifts(obj['name'], obj['RADeg'], obj['decDeg'])
            if self.configDict['addSDSSImage'] == True:
                self.fetchSDSSDR8Image(obj)
            for surveyString, label in zip(self.configDict['skyviewSurveyStrings'], self.configDict['skyviewLabels']):
                self.fetchSkyviewJPEG(obj['name'], obj['RADeg'], obj['decDeg'], surveyString, label) 
        
