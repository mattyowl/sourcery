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

import astropy.io.fits as pyfits
import astropy.table as atpy
from astLib import astWCS
import numpy as np
from scipy import interpolate
try:
    import pyvips
except:
    print("WARNNG: couldn't import pyvips - won't be able to handle large DES .tiffs")
from PIL import Image
Image.MAX_IMAGE_PIXELS=100000001 
import os
import sourcery
from sourcery import catalogTools
import urllib3
import time
import IPython

class TileDir:
    """The TileDir class handles directories that contain entire surveys (e.g., DES, KiDS, 
    S82) which have been broken up into tiles. We extract images for sources in catalogs
    from these tiles, taking into account overlaps etc. (i.e., if you ask for an 8' clip,
    you will get one).
    
    """
    
    def __init__(self, label, tileDir, sourceryCacheDir, WCSTabPath = None, sizePix = 1024):
        """Initialise a TileDir object. TileDirs handle directories that contain preview .jpg
        images of an entire survey, broken into tiles. We follow how this was done for DES 
        DR1, and adapt it to other surveys.
        
        label specifies how this survey will be referred to within Sourcery - e.g., 'DES',
        'KiDS', 'IAC-S82'. Output images for sources will be stored under 
        sourceryCacheDir/label/
        
        tileDir gives the location of the preview .jpg images of each tile, and should also
        contain a file called WCSTab.fits, which is a .fits table containing WCS header 
        keywords corresponding to each tile image. 
        
        WCSTab.fits should have the following columns:
        TILENAME NAXIS NAXIS1 NAXIS2 CTYPE1... etc.
        
        The TILENAME column should correspond to the .jpg file names in tileDir (minus the
        .jpg suffix). The other columns should be the WCS-related FITS header keywords from 
        the original images that were used to make the .jpg preview images.
        
        WCSTabPath can be given if this file is stored in a different location to the default
        (tileDir/WCSTab.fits).
        
        NOTE: we handle DES itself slightly differently to this, as we can also fetch DES tile
        images as needed over the internet. So label = 'DES' is special (see below...)
        
        """
        
        self.label=label
        self.tileDir=tileDir
        self.outputCacheDir=sourceryCacheDir+os.path.sep+label
        
        self.sizePix=sizePix
        
        # For DES, we can fetch .tif images of tiles and convert to .jpg as needed
        # If doing that, we need to make sure this directory exists...
        if os.path.exists(self.tileDir) == False:
            os.makedirs(self.tileDir)
        
        # For fetching DES .tif images if needed
        self.http=urllib3.PoolManager()

        if os.path.exists(self.outputCacheDir) == False:
            os.makedirs(self.outputCacheDir)
        
        if WCSTabPath == None:
            self.WCSTabPath=self.tileDir+os.path.sep+"WCSTab.fits"
        else:
            self.WCSTabPath=WCSTabPath
        
        # Special treatment for DES
        if self.label == 'DES':
            self.WCSTabPath=sourcery.__path__[0]+os.path.sep+"data"+os.path.sep+"DES_DR1_TILE_INFO.csv"
            
        self.WCSDict=None
        
    
    def setUpWCSDict(self):
        """Sets-up WCS info, needed for fetching images (can take ~30 sec or so, don't do this lightly).
        
        """
        
        # Add some extra columns to speed up searching
        self.tileTab=atpy.Table().read(self.WCSTabPath)        
        self.tileTab.add_column(atpy.Column(np.zeros(len(self.tileTab)), 'RAMin'))
        self.tileTab.add_column(atpy.Column(np.zeros(len(self.tileTab)), 'RAMax'))
        self.tileTab.add_column(atpy.Column(np.zeros(len(self.tileTab)), 'decMin'))
        self.tileTab.add_column(atpy.Column(np.zeros(len(self.tileTab)), 'decMax'))
        self.WCSDict={}
        keyWordsToGet=['NAXIS', 'NAXIS1', 'NAXIS2', 'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 
                        'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CDELT1', 'CDELT2', 'CUNIT1', 'CUNIT2']
        for row in self.tileTab:
            newHead=pyfits.Header()
            for key in keyWordsToGet:
                if key in self.tileTab.keys():
                    newHead[key]=row[key]
            # Defaults if missing (needed for e.g. DES)
            if 'NAXIS' not in newHead.keys():
                newHead['NAXIS']=2
            if 'CUNIT1' not in newHead.keys():
                newHead['CUNIT1']='DEG'
            if 'CUNIT2' not in newHead.keys():
                newHead['CUNIT2']='DEG'
            self.WCSDict[row['TILENAME']]=astWCS.WCS(newHead.copy(), mode = 'pyfits')  
            ra0, dec0=self.WCSDict[row['TILENAME']].pix2wcs(0, 0)
            ra1, dec1=self.WCSDict[row['TILENAME']].pix2wcs(row['NAXIS1'], row['NAXIS2'])
            if ra1 > ra0:
                ra1=-(360-ra1)
            row['RAMin']=min([ra0, ra1])
            row['RAMax']=max([ra0, ra1])
            row['decMin']=min([dec0, dec1])
            row['decMax']=max([dec0, dec1])
        
            
    def fetchImage(self, name, RADeg, decDeg, sizeArcmin, refetch = False):
        """Make .jpg image of a source of sizeArcmin, using preview .jpg tiles covering a whole survey. 
                
        """

        if self.WCSDict == None:
            self.setUpWCSDict()        
        
        # Inside footprint check
        raMask=np.logical_and(np.greater_equal(RADeg, self.tileTab['RAMin']), np.less(RADeg, self.tileTab['RAMax']))
        decMask=np.logical_and(np.greater_equal(decDeg, self.tileTab['decMin']), np.less(decDeg, self.tileTab['decMax']))
        tileMask=np.logical_and(raMask, decMask)
        if tileMask.sum() == 0:
            print("... object not in any %s tiles ..." % (self.label))
            return None
                       
        outFileName=self.outputCacheDir+os.path.sep+catalogTools.makeRADecString(RADeg, decDeg)+".jpg"
        
        # Procedure: spin through tile WCSs, find which tiles we need, paste pixels into low-res image
        if os.path.exists(outFileName) == False or refetch == True:
            
            # Blank WCS
            CRVAL1, CRVAL2=RADeg, decDeg
            sizePix=self.sizePix
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
            newHead['CDELT1']=-xOutPixScale
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
    
            matchTilesList=self.tileTab['TILENAME'][tileMask].tolist()                        
            if matchTilesList == []:
                print("... object not in any %s tiles ..." % (self.label))
                return None
            else:

                # We work with .jpg preview files that we made with STIFF
                for tileName in matchTilesList:

                    matchTab=self.tileTab[np.where(self.tileTab['TILENAME'] == tileName)][0]
                    
                    # Special treament for DES - can fetch .tiff previews for tiles over network
                    # (i.e., don't have to make them ourselves)
                    # We then convert them to .jpg
                    if self.label == 'DES':
                        tiffFileName=self.tileDir+os.path.sep+tileName+".tiff"
                        tileJPGFileName=tiffFileName.replace(".tiff", ".jpg")
                        if os.path.exists(tileJPGFileName) == False:
                            if os.path.exists(tiffFileName) == False:
                                print("... downloading .tiff image for tileName = %s ..." % (tileName))
                                resp=self.http.request('GET', str(matchTab['TIFF_COLOR_IMAGE']))
                                with open(tiffFileName, 'wb') as f:
                                    f.write(resp.data)
                                    f.close()
                                # Old
                                #try:
                                    #urllib.urlretrieve(str(matchTab['TIFF_COLOR_IMAGE']), tiffFileName)
                                #except:
                                    #raise Exception, "downloading DES .tiff image failed"
                            # NOTE: we use pyvips, because images are too big for PIL
                            # We save disk space by caching a lower quality version of the entire tile
                            print("... converting .tiff for tileName = %s to .jpg ..." % (tileName))
                            im=pyvips.Image.new_from_file(tiffFileName, access = 'sequential')
                            im.write_to_file(tileJPGFileName+'[Q=80]')
                            os.remove(tiffFileName)

                    # Everything...
                    tileJPGFileName=self.tileDir+os.path.sep+tileName+".jpg"
                    if os.path.exists(tileJPGFileName) == True:
                        #try:
                        #im=pyvips.Image.new_from_file(tileJPGFileName, access = 'sequential')
                        #except:
                        im=Image.open(tileJPGFileName)
                    else:
                        print("... tile %s missing from %s tiles .jpg preview directory (probably missing i or g-band coverage) ..." % (tileName, self.label))
                        continue
                    
                    # New - several orders of magnitude quicker
                    # Assumes images are aligned N vertically, E at left
                    #d=np.ndarray(buffer = im.write_to_memory(), dtype = np.uint8, shape = [im.height, im.width, im.bands])
                    d=np.array(im)
                    d=np.flipud(d)
                    inWCS=self.WCSDict[tileName]

                    # NOTE: Linear interpolation like this is v. quick but wrong for TAN at large dec.
                    # Since we're only using this for display purposes, the trick we use here should be okay
                    # (well, introduces a little rotation on DES images)
                    RAc, decc=outWCS.getCentreWCSCoords()
                    xc, yc=inWCS.wcs2pix(RAc, decc)
                    xc, yc=int(xc), int(yc)
                    xIn=np.arange(d.shape[1])
                    yIn=np.arange(d.shape[0])
                    inRACoords=np.array(inWCS.pix2wcs(xIn, [yc]*len(xIn)))
                    inDecCoords=np.array(inWCS.pix2wcs([xc]*len(yIn), yIn))
                    inRA=inRACoords[:, 0]
                    inDec=inDecCoords[:, 1]
                    RAToX=interpolate.interp1d(inRA, xIn, fill_value = 'extrapolate')
                    DecToY=interpolate.interp1d(inDec, yIn, fill_value = 'extrapolate')
                    outRACoords=np.array(outWCS.pix2wcs(np.arange(outData.shape[1]), [0]*outData.shape[1]))
                    outDecCoords=np.array(outWCS.pix2wcs([0]*np.arange(outData.shape[0]), np.arange(outData.shape[0])))
                    outRA=outRACoords[:, 0]
                    outDec=outDecCoords[:, 1]
                    xIn=np.array(RAToX(outRA), dtype = int)
                    yIn=np.array(DecToY(outDec), dtype = int)
                    xMask=np.logical_and(xIn >= 0, xIn < im.width)
                    yMask=np.logical_and(yIn >= 0, yIn < im.height)
                    xOut=np.arange(outData.shape[1])
                    yOut=np.arange(outData.shape[0])
                    for i in yOut[yMask]:
                        outData[i][xMask]=d[yIn[i], xIn[xMask]]

                # Flips needed to get N at top, E at left
                outData=np.flipud(outData)
                #outData=np.fliplr(outData)
                
                # We could do this with vips... but lazy...
                outIm=Image.fromarray(outData)
                outIm.save(outFileName)
                outIm.close()
                print("... made %s cut-out .jpg ..." % (self.label))
                
