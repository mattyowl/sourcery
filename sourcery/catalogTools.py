"""

    Copyright 2014 Matt Hilton (matt.hilton@mykolab.com)
    
    This file is part of Sourcery.

    Sourcery is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sourcery is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sourcery.  If not, see <http://www.gnu.org/licenses/>.

"""

from astLib import *
from scipy import ndimage
import numpy
import operator
import os
import urllib
import urllib2
import sys
import time
import atpy
import datetime
import IPython

#-------------------------------------------------------------------------------------------------------------
XMATCH_RADIUS_DEG=1.4/60.0  # catalog matching radius, for sim comparisons

#-------------------------------------------------------------------------------------------------------------
def clipSmoothedTanResampledImage(obj, mapData, mapWCS, sizeDeg, gaussSmoothArcSecRadius, 
                                  outFileName = None, sizePix = 200):
    """Clips a tan resampled, (optionally smoothed) section around an object in an image, writes it out
    to outFileName, and returns a dictionary containing the clipped map data and WCS. 
        
    """
    
    RADeg=obj['RADeg']
    decDeg=obj['decDeg']
    
    # This solves for the RA, dec coords we need to clip to get exactly the right dimensions and dodge
    # the cea projection shenanigans
    tolerance=1e-8  # in degrees on sky
    targetHalfSizeSkyDeg=(sizeDeg*1.1)/2.0  # slightly bigger, trim down afterwards
    funcCalls=["astCoords.calcAngSepDeg(RADeg, decDeg, guess, decDeg)",
               "astCoords.calcAngSepDeg(RADeg, decDeg, guess, decDeg)",
               "astCoords.calcAngSepDeg(RADeg, decDeg, RADeg, guess)",
               "astCoords.calcAngSepDeg(RADeg, decDeg, RADeg, guess)"]
    coords=[RADeg, RADeg, decDeg, decDeg]
    signs=[1.0, -1.0, 1.0, -1.0]
    results=[]
    for f, c, sign in zip(funcCalls, coords, signs):
        # Initial guess range
        maxGuess=sign*targetHalfSizeSkyDeg*2.0
        minGuess=sign*targetHalfSizeSkyDeg/10.0
        guessStep=(maxGuess-minGuess)/10.0
        guesses=numpy.arange(minGuess+c, maxGuess+c, guessStep)
        for i in range(60):
            minSizeDiff=1e6
            bestGuess=None
            for guess in guesses:
                sizeDiff=abs(eval(f)-targetHalfSizeSkyDeg)
                if sizeDiff < minSizeDiff:
                    minSizeDiff=sizeDiff
                    bestGuess=guess
            if minSizeDiff < tolerance:
                break
            else:
                guessRange=abs((maxGuess-minGuess))
                maxGuess=bestGuess+guessRange/4.0
                minGuess=bestGuess-guessRange/4.0
                guessStep=(maxGuess-minGuess)/10.0
                guesses=numpy.arange(minGuess, maxGuess, guessStep)
        results.append(bestGuess)
    RAMax=results[0]
    RAMin=results[1]
    decMax=results[2]
    decMin=results[3]

    # To tan proj, scale, smooth, save
    # We'll make these 1% bigger than 12 arcmin actually to see if we can get rid of annoying edge effect
    # in contour plots
    tanClip=astImages.clipUsingRADecCoords(mapData, mapWCS, RAMin, RAMax, decMin, decMax)
    try:
        tanClip=astImages.resampleToTanProjection(tanClip['data'], tanClip['wcs'], outputPixDimensions = [sizePix, sizePix])
    except:
        print "Hmm - tan reprojection failed"
        IPython.embed()
        sys.exit()
    #tanClip=astImages.clipImageSectionWCS(tanClip['data'], tanClip['wcs'], RADeg, decDeg, sizeDeg*1.01)
    dataClip=tanClip['data']
    scaleFactor=float(sizePix)/float(tanClip['data'].shape[1])
    tanClip=astImages.scaleImage(tanClip['data'], tanClip['wcs'], scaleFactor)
    if gaussSmoothArcSecRadius != None:
        radPix=(gaussSmoothArcSecRadius/3600.0)/tanClip['wcs'].getPixelSizeDeg()
        tanClip['data']=ndimage.gaussian_filter(tanClip['data'], radPix)                        
    
    if outFileName != None:
        astImages.saveFITS(outFileName, tanClip['data'], tanClip['wcs']) 
    
    return tanClip
        
#-------------------------------------------------------------------------------------------------------------
def tab2DS9(tab, outFileName, color = "cyan"):
    """Writes atpy Table into a ds9 region file.
    
    If color == 'key', then use 'color' key in object dictionary to set color.
    
    """

    outFile=file(outFileName, "w")
    timeStamp=datetime.datetime.today().date().isoformat()
    comment="# DS9 region file"
    outFile.write(comment)
    outFile.write('global dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    for obj in tab:
        outFile.write("fk5;point(%.6f,%.6f) # point=boxcircle color={%s} text={%s}\n" \
                    % (obj['RADeg'], obj['decDeg'], color, obj['name']))
    outFile.close()
             
#-------------------------------------------------------------------------------------------------------------
def parseNEDResult(inFileName, onlyObjTypes = None):
    """Parses NED tab-delimited text file query result, returns dictionary.
    
    onlyObjTypes can be a string indicating types of objects only to include e.g. GClstr
    
    """
    
    if inFileName != None and os.path.exists(inFileName):
        inFile=file(inFileName, "r")
        lines=inFile.readlines()
        inFile.close()
    else:
        # Fail safe in case we couldn't contact NED
        lines=[]

    dataStarted=False
    labels=[]
    names=[]
    RAs=[]
    decs=[]
    sourceTypes=[]
    redshifts=[]
    for line in lines:
        bits=line.split("\t")
        if bits[0] == "1":
            dataStarted=True
        if dataStarted == True:
            if onlyObjTypes == str(bits[4]) or onlyObjTypes == None:
                labels.append(bits[0])
                names.append(bits[1])
                RAs.append(float(bits[2]))
                decs.append(float(bits[3]))
                sourceTypes.append(str(bits[4]))
                if bits[6] == '':
                    redshifts.append('N/A')
                else:
                    redshifts.append(str(bits[6]))

    return {'labels': labels, 'names': names, 'RAs': RAs, 'decs': decs, 'sourceTypes': sourceTypes, 'redshifts': redshifts}

#-------------------------------------------------------------------------------------------------------------
def addSDSSRedshifts(catalog, cacheDir = "SDSSQueryResults"):
    """Queries SDSS for redshifts. 
    
    """
    
    print ">>> Adding spec zs from SDSS ..."
    #url = 'http://cas.sdss.org/astrodr7/en/tools/search/x_sql.asp'
    #url = 'http://skyserver.sdss3.org/dr10/en/tools/search/x_sql.aspx'
    url = 'http://skyserver.sdss.org/dr12/en/tools/search/x_sql.aspx'

    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir)
    
    count=0
    consecutiveQueryCount=0
    for obj in catalog:

        count=count+1
        print "... %s (%d/%d) ..." % (obj['name'], count, len(catalog))
        
        outFileName=cacheDir+os.path.sep+"%s.csv" % (obj['name'].replace(" ", "_"))
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
            """ % (obj['RADeg'], obj['RADeg'], obj['decDeg'], obj['decDeg'])

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
            
            consecutiveQueryCount=consecutiveQueryCount+1
            if consecutiveQueryCount > 50:
                print "... sleeping to give SDSS server a break ..."
                time.sleep(60)
                consecutiveQueryCount=0
        
        else:
            
            inFile=file(outFileName, "r")
            lines=inFile.readlines()
            inFile.close()
        
        # Parse .csv into catalog
        if lines[0] == "No objects have been found\n":
            obj['SDSSRedshifts']=None
        elif len(lines) > 1 and lines[1] == '"ERROR: Maximum 60 queries allowed per minute. Rejected query: SELECT \n':
            os.remove(outFileName)
            raise Exception, "Exceeded 60 queries/min on SDSS server. Take a breather and rerun nemo (previous queries cached)."
        else:
            obj['SDSSRedshifts']=[]
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
                            raise Exception, "Exceeded 60 queries/min on SDSS server. Take a breather and rerun (previous queries cached)."
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
                    obj['SDSSRedshifts'].append(zDict)

    