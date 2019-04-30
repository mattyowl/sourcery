"""

    Copyright 2014-2018 Matt Hilton (matt.hilton@mykolab.com)
    
    This file is part of Sourcery.

    Sourcery is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sourcery is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FORf A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sourcery.  If not, see <http://www.gnu.org/licenses/>.

"""

from astLib import *
from scipy import ndimage
import numpy
import operator
import os
import urllib3
import sys
import time
import astropy.table as atpy
import datetime
import IPython

#-------------------------------------------------------------------------------------------------------------
XMATCH_RADIUS_DEG=1.4/60.0  # catalog matching radius, for sim comparisons

#-------------------------------------------------------------------------------------------------------------
def makeRADecString(RADeg, decDeg):
    """Switched to using %.5f_%.5f as part of image file names.
    
    """
    return "%.5f_%.5f" % (RADeg, decDeg)
    
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

    # To Tan proj, scale, smooth, save
    # We'll make these 1% bigger than 12 arcmin actually to see if we can get rid of annoying edge effect
    # in contour plots
    # NOTE: This has problems if crossing RA = 0, get RAMin spanning 0.4 deg to 359 deg say
    # Doesn't help if use other clipping routines either!
    tanClip=astImages.clipUsingRADecCoords(mapData, mapWCS, RAMin, RAMax, decMin, decMax)
    
    try:
        tanClip=astImages.resampleToTanProjection(tanClip['data'], tanClip['wcs'], outputPixDimensions = [sizePix, sizePix])
    except:
        print("WARNING: image needs clipping over 0h RA? Fix later")
        return None
    
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

    outFile=open(outFileName, "w")
    timeStamp=datetime.datetime.today().date().isoformat()
    comment="# DS9 region file"
    outFile.write(comment)
    outFile.write('global dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    for obj in tab:
        outFile.write("fk5;point(%.6f,%.6f) # point=boxcircle color={%s} text={%s}\n" \
                    % (obj['RADeg'], obj['decDeg'], color, obj['name']))
    outFile.close()
            
#-------------------------------------------------------------------------------------------------------------
def parseNEDResult(inFileName, onlyObjTypes = []):
    """Parses NED tab-delimited text file query result, returns dictionary.
    
    onlyObjTypes can be a string indicating types of objects only to include e.g. GClstr
    
    """
    
    if inFileName != None and os.path.exists(inFileName):
        inFile=open(inFileName, "r")
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
            try:
                if str(bits[4]) in onlyObjTypes:
                    labels.append(bits[0])
                    names.append(bits[1])
                    RAs.append(float(bits[2]))
                    decs.append(float(bits[3]))
                    sourceTypes.append(str(bits[4]))
                    if bits[6] == '':
                        redshifts.append('N/A')
                    else:
                        redshifts.append(str(bits[6]))
            except:
                continue # Probably a mal-formed line (we have seen things like '\t0\n' in the middle of the data)
                
    return {'labels': labels, 'names': names, 'RAs': RAs, 'decs': decs, 'sourceTypes': sourceTypes, 'redshifts': redshifts}

#-------------------------------------------------------------------------------------------------------------
def byteSwapArr(arr):
    """FITS is big-endian, but cython likes native-endian arrays (little-endian for x86)... so, byteswap
    if we need.
    
    """

    if arr.dtype.byteorder == '>':
        arr=arr.byteswap().newbyteorder('=')

    return arr

#-------------------------------------------------------------------------------------------------------------
def fetchSDSSRedshifts(cacheDir, name, RADeg, decDeg):
    """Queries SDSS DR13 for redshifts, writing output into cacheDir.
    
    Returns a list of dictionaries containing the redshift catalog.
    
    """
    
    if os.path.exists(cacheDir) == False:
        os.makedirs(cacheDir)
        
    if decDeg > -20:
        #url='http://skyserver.sdss3.org/dr10/en/tools/search/x_sql.aspx'
        url='http://skyserver.sdss.org/dr14/en/tools/search/x_results.aspx'
        outFileName=cacheDir+os.path.sep+"%s.csv" % (name.replace(" ", "_"))
        if os.path.exists(outFileName) == False:
            # Open file straightaway in case running threaded
            with open(outFileName, 'wb') as f:
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
                http=urllib3.PoolManager()
                resp=http.request('GET', url, fields = {'searchtool': 'SQL', 
                                                        'TaskName': 'Skyserver.Search.SQL', 
                                                        'cmd': fsql, 
                                                        'format': "csv"})
                f.write(resp.data)
      
        with open(outFileName, "r") as inFile:
            lines=inFile.readlines()
        if len(lines) == 0:
            print("... file %s is empty - retrying ..." % (outFileName))
            os.remove(outFileName)
            fetchSDSSRedshifts(cacheDir, name, RADeg, decDeg)

    else:
        return []
        
    # Parse .csv into catalog
    SDSSRedshifts=[]
    if lines[0] == "No objects have been found\n":
        SDSSRedshifts=[]
    elif len(lines) > 1 and lines[1] == '"ERROR: Maximum 60 queries allowed per minute. Rejected query: SELECT \n':
        os.remove(outFileName)
        raise Exception("Exceeded 60 queries/min on SDSS server. Take a breather and rerun (previous queries cached).")
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
                        raise Exception("Exceeded 60 queries/min on SDSS server. Take a breather and rerun nemo (previous queries cached).")
                    else:
                        print("Probably asking for too many queries from SDSS... waiting then trying again.")
                        time.sleep(60)
                        os.remove(outFileName)
                        SDSSRedshifts=fetchSDSSRedshifts(cacheDir, name, RADeg, decDeg)
                        break
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
    
