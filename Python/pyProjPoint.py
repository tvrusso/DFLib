"""
   This simple Python script is meant to demonstrate how to use the DFLib
   python bindings with Point and Report classes implemented in Python.

   You do not have to build the C++ library version of DFLib to use this, only
   the Python library:

   python setup.py build_ext --inplace [--swig <path to swig>]

   This assumes you actually have installed swig, and are set up to compile
   C++ codes.  This should create the DFLib python module.

   This script also assumes you have installed the "pyproj" module from
   http://jswhit.github.io/pyproj/
   or via your system's package management system.

"""

import DFLib
import copy
import math
from pyproj import Proj

"""
  The two classes below constitute a simplistic implementation of the minimum
  work you'd have to do to use DFLib in Python.

  The "pyProjPoint" class derives from the DFLib::Abstract::Point interface,
  and implements its interface methods:
    __init__ : the constructor for the class
    setXY:     Sets coordinates in "XY" space, the one that is actually used
               by the ReportCollection class for computation.
    getXY:     Returns an STL vector containing the XY position.
    Clone:     Makes a new pyProjPoint and populates it with the values in the
               current object, returning a pointer to it.  This method is ONLY
               used by the "Fix Cut Average" computation.
    setUserCoords:  Allows you to specify the point location in some coordinate
                    system other than the XY coordinates used for computation.
                    In this implementation, XY is Mercator Projection 
                    coordinates, and UserCoords are lat/lon.
    getUserCoords:  Returns the point location in the user coordinates.
"""

class pyProjPoint(DFLib.Point):
   def __init__(self):
      super(pyProjPoint,self).__init__()
      self.myXY=DFLib.vectord(2)
      self.myXY[0]=0
      self.myXY[1]=0
      self.myMercProj=Proj('+proj=merc +datum=WGS84 +lat_ts=0')   
      self.myUserCoords=DFLib.vectord(2)
   def setXY(self,aPosition):
      self.myXY[0]=aPosition[0]
      self.myXY[1]=aPosition[1]
   def getXY(self):
      return self.myXY
  # This is NOT really a "clone" operation, but it serves the purpose enough
  # for DFLib
   def Clone(self):
      foo=pyProjPoint().__disown__()
      foo.setXY(self.getXY())
      return foo
   def getUserCoords(self):
      self.myUserCoords[0],self.myUserCoords[1]=self.myMercProj(self.myXY[0],self.myXY[1],inverse=True)
      return self.myUserCoords
   def setUserCoords(self,uPosition):
      tempVect=DFLib.vectord(2)
      tempVect[0],tempVect[1]=self.myMercProj(uPosition[0],uPosition[1])
      self.setXY(tempVect)

"""
  The pyProjReport class implements the DFLib::Abstract::Report interface.

  setReceiverLocation:  Sets the receiver location in XY coordinates.
  getReceiverLocation:  Returns the receiver location in XY coordinates (in 
                        an STL vector)
  setBearing:           Sets the report's bearing (in the XY space,
                        which means it needs to have grid convergence and
                        magnetic declination taken into account!) This is in
                        degrees.
  setSigma:             Sets the standard deviation of the distribution of
                        random errors expected from this type of receiver.
  getReportBearingRadians:  Returns the bearing converted to radians.
                        (0<=bearing<=2pi)
  getBearing:           Returns the bearing in degrees
  getBearingStandardDeviationRadians:  Returns sigma in units of radians
                        (0<=sigma<=2pi)
  getSigma:             Returns sigma in degrees.
"""

class pyProjReport(DFLib.Report):
   def __init__(self,name,valid,Point,bearing,sigma):
      super(pyProjReport,self).__init__(name,valid)
      self.myPoint=pyProjPoint()
      self.myPoint.setXY(Point.getXY())
      self.setBearing(bearing)
      self.setSigma(sigma)
   def setReceiverLocation(self,theLocation):
      self.myPoint.setXY(theLocation)
   def getReceiverLocation(self):
      return self.myPoint.getXY()
   def setBearing(self,theBearing):
      self.bearing=theBearing*math.pi/180.0
      while (self.bearing < 0):
         self.bearing = self.bearing + 2*math.pi
      while (self.bearing >= 2*math.pi):
         self.bearing = self.bearing - 2*math.pi
   def setSigma(self,theSigma):
      self.sigma = theSigma*math.pi/180.0
   def getReportBearingRadians(self):
      return self.bearing
   def getBearing(self):
      return self.getReportBearingRadians()*180.0/math.pi
   def getBearingStandardDeviationRadians(self):
      return self.sigma
   def getSigma(self):
      return self.getBearingStandardDeviationRadians()*180.0/math.pi

#Create 3 reports.  These correspond exactly to the "receivers3" points
#in the main DFLib directory.
# The bearings are the product of a single run of testlsDFfix, randomizing
# around the true bearings to a transmitter at 106d35'W 34d55'N
# according to the sigmas associated with the transmitters
#

v=DFLib.vectord(2)
v[0],v[1]=-(106+38./60.+15.4/3600.) , (34+56./60.+48.7/3600)
r1Point=pyProjPoint()
r1Point.setUserCoords(v)
r1=pyProjReport("r1",True,r1Point,123.866,1.0)
r1p=r1.getReceiverLocation()

v[0],v[1]=-(106+18./60.+17.0/3600.) , (34+58./60.+12.5/3600)
r2Point=pyProjPoint()
r2Point.setUserCoords(v)
r2=pyProjReport("r2",True,r2Point,-104.265,3.0)
r2p=r2.getReceiverLocation()

v[0],v[1]=-(106+28./60.+53.8/3600.) , (35+10./60.+33.9/3600)
r3Point=pyProjPoint()
r3Point.setUserCoords(v)
r3=pyProjReport("r3",True,r3Point,-160.354,2.0)
r3p=r3.getReceiverLocation()

print 'r1: (x,y)=(%f,%f)  bearing=%f sigma=%f'%(r1p[0],r1p[1],r1.getBearing(),r1.getSigma())
print 'r1: (x,y)=(%f,%f)  bearing(rad)=%f sigma(rad)=%f'%(r1p[0],r1p[1],r1.getReportBearingRadians(),r1.getBearingStandardDeviationRadians())
print 'r2: (x,y)=(%f,%f)  bearing=%f sigma=%f'%(r2p[0],r2p[1],r2.getBearing(),r2.getSigma())
print 'r2: (x,y)=(%f,%f)  bearing(rad)=%f sigma(rad)=%f'%(r2p[0],r2p[1],r2.getReportBearingRadians(),r2.getBearingStandardDeviationRadians())

print 'r3: (x,y)=(%f,%f)  bearing=%f sigma=%f'%(r3p[0],r3p[1],r3.getBearing(),r3.getSigma())
print 'r3: (x,y)=(%f,%f)  bearing(rad)=%f sigma(rad)=%f'%(r3p[0],r3p[1],r3.getReportBearingRadians(),r3.getBearingStandardDeviationRadians())

rc=DFLib.ReportCollection()
rc.addReport(r1)
rc.addReport(r2)
rc.addReport(r3)

LSfix=pyProjPoint()
rc.computeLeastSquaresFix(LSfix)

LSv=LSfix.getXY()
LSu=LSfix.getUserCoords()
print 'ls fix: x=%f y=%f'%(LSv[0],LSv[1])
print 'ls fix(LL): x=%f y=%f'%(LSu[0],LSu[1])

FCAfix=pyProjPoint()
FCA_stddev=DFLib.vectord(2)
rc.computeFixCutAverage(FCAfix,FCA_stddev)
   
FCAv=FCAfix.getXY()
FCAu=FCAfix.getUserCoords()

print 'FCA fix: x=%f y=%f'%(FCAv[0],FCAv[1])
print 'FCA fix(LL): x=%f y=%f'%(FCAu[0],FCAu[1])
print 'FCA stddevs(LL) x=%f y=%f'%(FCA_stddev[0],FCA_stddev[1])

MLfix = LSfix.Clone()
rc.aggressiveComputeMLFix(MLfix)


MLv=MLfix.getXY()
MLu=MLfix.getUserCoords()
print 'ML fix: x=%f y=%f'%(MLv[0],MLv[1])
print 'ML fix(LL): x=%f y=%f'%(MLu[0],MLu[1])
