# Cython routines for sourcery
#

#include "python.pxi"
#include "numpy.pxi"

from astLib import *
import numpy as np
cimport numpy as np
import math
import time
import sys

#-------------------------------------------------------------------------------------------------------------
def makeDegreesDistanceMap(np.ndarray[np.float32_t, ndim=2] data, wcs, RADeg, decDeg, maxDistDegrees):
    """Returns an array of distance in degrees from given position. maxDistDegrees sets the limit around
    the given RADeg, decDeg position to which the distance is calculated.
    
    """
    
    # Pixel distance grid            
    cdef float x0, y0, ra0, dec0, ra1, dec1, xPixScale, yPixScale
    cdef Py_ssize_t x, y, X, Y, minX, maxX, minY, maxY, xDistPix, yDistPix
    cdef np.ndarray[np.float32_t, ndim=2] degreesMap
    
    x0, y0=wcs.wcs2pix(RADeg, decDeg)
    ra0, dec0=RADeg, decDeg
    ra1, dec1=wcs.pix2wcs(x0+1, y0+1)
    xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
    
    xDistPix=int(round((maxDistDegrees)/xPixScale))
    yDistPix=int(round((maxDistDegrees)/yPixScale))

    Y=data.shape[0]
    X=data.shape[1]
    
    # Real space map of angular distance in degrees, but only consider values near x0, y0
    degreesMap=np.ones([Y, X], dtype=np.float32)*1e6
    minX=int(round(x0))-xDistPix
    maxX=int(round(x0))+xDistPix
    minY=int(round(y0))-yDistPix
    maxY=int(round(y0))+yDistPix
    if minX < 0:
        minX=0
    if maxX >= X:
        maxX=X-1
    if minY < 0:
        minY=0
    if maxY >= Y:
        maxY=Y-1
    for y in range(minY, maxY):
        for x in range(minX, maxX):
            yDeg=(y-y0)*yPixScale
            xDeg=(x-x0)*xPixScale
            degreesMap[y, x]=math.sqrt(xDeg*xDeg+yDeg*yDeg)
    
    return degreesMap

#-------------------------------------------------------------------------------------------------------------
def makeXYDegreesDistanceMaps(np.ndarray[np.float32_t, ndim=2] data, wcs, RADeg, decDeg, maxDistDegrees):
    """Returns an array of distance along x, y axes in degrees from given position. maxDistDegrees sets the 
    limit around the given RADeg, decDeg position to which the distance is calculated.
    
    """
    
    # Pixel distance grid            
    cdef float x0, y0, ra0, dec0, ra1, dec1, xPixScale, yPixScale
    cdef Py_ssize_t x, y, X, Y, minX, maxX, minY, maxY, xDistPix, yDistPix
    cdef np.ndarray[np.float32_t, ndim=2] xDegreesMap
    cdef np.ndarray[np.float32_t, ndim=2] yDegreesMap

    x0, y0=wcs.wcs2pix(RADeg, decDeg)
    ra0, dec0=RADeg, decDeg
    ra1, dec1=wcs.pix2wcs(x0+1, y0+1)
    xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
    
    xDistPix=int(round((maxDistDegrees)/xPixScale))
    yDistPix=int(round((maxDistDegrees)/yPixScale))

    Y=data.shape[0]
    X=data.shape[1]
    
    # Real space map of angular distance in degrees, but only consider values near x0, y0
    xDegreesMap=np.ones([Y, X], dtype=np.float32)*1e6
    yDegreesMap=np.ones([Y, X], dtype=np.float32)*1e6
    minX=int(round(x0))-xDistPix
    maxX=int(round(x0))+xDistPix
    minY=int(round(y0))-yDistPix
    maxY=int(round(y0))+yDistPix
    if minX < 0:
        minX=0
    if maxX >= X:
        maxX=X-1
    if minY < 0:
        minY=0
    if maxY >= Y:
        maxY=Y-1
    for y in range(minY, maxY):
        for x in range(minX, maxX):
            yDeg=(y-y0)*yPixScale
            xDeg=(x-x0)*xPixScale
            xDegreesMap[y, x]=xDeg
            yDegreesMap[y, x]=yDeg
    
    return [xDegreesMap, yDegreesMap]
