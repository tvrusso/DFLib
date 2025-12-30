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

#  The two classes below constitute a simplistic implementation of the minimum
#  work you'd have to do to use DFLib in Python.
#
#  The "pyProjPoint" class derives from the DFLib::Abstract::Point interface,
#  and implements its interface methods:
#    __init__ : the constructor for the class
#    setXY:     Sets coordinates in "XY" space, the one that is actually used
#               by the ReportCollection class for computation.
#    getXY:     Returns an STL vector containing the XY position.
#    Clone:     Makes a new pyProjPoint and populates it with the values in the
#               current object, returning a pointer to it.  This method is ONLY
#               used by the "Fix Cut Average" computation.
#    setUserCoords:  Allows you to specify the point location in some coordinate
#                    system other than the XY coordinates used for computation.
#                    In this implementation, XY is Mercator Projection 
#                    coordinates, and UserCoords are lat/lon.
#    getUserCoords:  Returns the point location in the user coordinates.

class pyProjPoint(DFLib.Point):
   """
   The "pyProjPoint" class derives from the DFLib::Abstract::Point interface,
   and implements its interface methods
   """
   def __init__(self,*args):
      """
      pyProjPoint constructor
      """
      super(pyProjPoint,self).__init__()
      self.myXY=DFLib.vectord(2)
      self.myMercProj=Proj('+proj=merc +datum=WGS84 +lat_ts=0')   
      self.myUserCoords=DFLib.vectord(2)
      if (len(args)!=0):
         self.setUserCoords(args[0])
      else:
         self.setXY((0,0))
   def setXY(self,aPosition):
      """
      given a vector (or list) set the location of this point in X-Y coordinates.
      """
      self.myXY[0]=aPosition[0]
      self.myXY[1]=aPosition[1]
   def getXY(self):
      """
      Return a pointer to the point's internal XY vector (an STL vector of 
      doubles).

      In C++ implementations, we return a const reference, not a pointer, to
      prevent changing the position by direct manipulation of the vector.
      I don't know how to do that in python to prevent abuse of the vector.
      It should be used only to *query* the position, never to change it.
      Do NOT use this to change the location of the point.
      """
      return self.myXY
  # This is NOT really a "clone" operation, but it serves the purpose enough
  # for DFLib
   def Clone(self):
      """
      Create a new object that is an exact copy of this one in every way.
      Typically used only by the Fix Cut Average computation routine, but
      can be used to generate new objects as needed without knowing what
      type they are.
      """
      foo=pyProjPoint(self.getUserCoords()).__disown__()
      return foo
   def getUserCoords(self):
      """
      Return a vector of coordinates of this point in a user coordinate
      system (which may be different from the XY coordinate system 
      used for computation of fixes).
      """
      self.myUserCoords[0],self.myUserCoords[1]=self.myMercProj(self.myXY[0],self.myXY[1],inverse=True)
      return self.myUserCoords
   def setUserCoords(self,uPosition):
      """
      Set the position of this point using a vector or list of coordinates in
      the user coordinate system.
      """
      self.setXY(self.myMercProj(uPosition[0],uPosition[1]))

class pyProjReport(DFLib.Report):
   """
     The pyProjReport class implements the DFLib::Abstract::Report interface.
   """
   def __init__(self,name,valid,Point,bearing,sigma):
      """
      pyProjReport constructor
      """
      super(pyProjReport,self).__init__(name,valid)
      self.myPoint=pyProjPoint(Point.getUserCoords())
      self.setBearing(bearing)
      self.setSigma(sigma)
   def setReceiverLocation(self,theLocation):
      """
      Sets the receiver location in XY coordinates.
      """
      self.myPoint.setXY(theLocation)
   def getReceiverLocation(self):
      """
      Returns the receiver location in XY coordinates.
      This is returned as a vectord (STL vector<double>) pointer.
      """
      return self.myPoint.getXY()
   def setBearing(self,theBearing):
      """
      Sets the report's bearing (in the XY space, which means it needs
      to have grid convergence and magnetic declination taken into
      account!) This is in degrees.
      """
      self.bearing=theBearing*math.pi/180.0
      while (self.bearing < 0):
         self.bearing = self.bearing + 2*math.pi
      while (self.bearing >= 2*math.pi):
         self.bearing = self.bearing - 2*math.pi
   def setSigma(self,theSigma):
      """
      Sets the standard deviation of the distribution of
      random errors expected from this type of receiver.  This is in degrees.
      """
      self.sigma = theSigma*math.pi/180.0
   def getReportBearingRadians(self):
      """
      Returns the report bearing relative to grid north (Y) in radians.
      """
      return self.bearing
   def getBearing(self):
      """
      Returns the report bearing in degrees.
      """
      return self.getReportBearingRadians()*180.0/math.pi
   def getBearingStandardDeviationRadians(self):
      """
      Returns receiver standard deviation in radians.
      """
      return self.sigma
   def getSigma(self):
      """
      Returns receiver standard deviation in degrees.
      """
      return self.getBearingStandardDeviationRadians()*180.0/math.pi

#Create 3 reports.  These correspond exactly to the "receivers3" points
#in the main DFLib directory.
# The bearings are from a run of testlsDF_proj using a transmitter
# position of -106.58274, 34.919551 and no bearing randomization


r1Point=pyProjPoint([-(106+38./60.+15.4/3600.) , (34+56./60.+48.7/3600)])
r1=pyProjReport("r1",True,r1Point,121.147,1.0)
r1p=r1.getReceiverLocation()

r2Point=pyProjPoint([-(106+18./60.+17.0/3600.) , (34+58./60.+12.5/3600)])
r2=pyProjReport("r2",True,r2Point,257.539,3.0)
r2p=r2.getReceiverLocation()

r3Point=pyProjPoint([-(106+28./60.+53.8/3600.) , (35+10./60.+33.9/3600)])
r3=pyProjReport("r3",True,r3Point,197.962,2.0)
r3p=r3.getReceiverLocation()

print('r1: (x,y)=(%f,%f)  bearing=%f sigma=%f'%(r1p[0],r1p[1],r1.getBearing(),r1.getSigma()))
print('r1: (x,y)=(%f,%f)  bearing(rad)=%f sigma(rad)=%f'%(r1p[0],r1p[1],r1.getReportBearingRadians(),r1.getBearingStandardDeviationRadians()))
print('r2: (x,y)=(%f,%f)  bearing=%f sigma=%f'%(r2p[0],r2p[1],r2.getBearing(),r2.getSigma()))
print('r2: (x,y)=(%f,%f)  bearing(rad)=%f sigma(rad)=%f'%(r2p[0],r2p[1],r2.getReportBearingRadians(),r2.getBearingStandardDeviationRadians()))

print('r3: (x,y)=(%f,%f)  bearing=%f sigma=%f'%(r3p[0],r3p[1],r3.getBearing(),r3.getSigma()))
print('r3: (x,y)=(%f,%f)  bearing(rad)=%f sigma(rad)=%f'%(r3p[0],r3p[1],r3.getReportBearingRadians(),r3.getBearingStandardDeviationRadians()))

rc=DFLib.ReportCollection()
rc.addReport(r1)
rc.addReport(r2)
rc.addReport(r3)

LSfix=pyProjPoint()
rc.computeLeastSquaresFix(LSfix)

LSv=LSfix.getXY()
LSu=LSfix.getUserCoords()
print('ls fix: x=%f y=%f'%(LSv[0],LSv[1]))
print('ls fix(LL): x=%f y=%f'%(LSu[0],LSu[1]))

FCAfix=pyProjPoint()
FCA_stddev=DFLib.vectord(2)
rc.computeFixCutAverage(FCAfix,FCA_stddev)
   
FCAv=FCAfix.getXY()
FCAu=FCAfix.getUserCoords()

print('FCA fix: x=%f y=%f'%(FCAv[0],FCAv[1]))
print('FCA fix(LL): x=%f y=%f'%(FCAu[0],FCAu[1]))
print('FCA stddevs(LL) x=%f y=%f'%(FCA_stddev[0],FCA_stddev[1]))

MLfix = LSfix.Clone()
rc.aggressiveComputeMLFix(MLfix)


MLv=MLfix.getXY()
MLu=MLfix.getUserCoords()
print('ML fix: x=%f y=%f'%(MLv[0],MLv[1]))
print('ML fix(LL): x=%f y=%f'%(MLu[0],MLu[1]))

# convert to DD MM.MMMM H
(ML_lon,ML_lat)=MLu;
ML_H_lon="E";
if (ML_lon < 0):
   ML_H_lon="W";
   ML_lon = -ML_lon;
ML_H_lat="N";
if (ML_lat < 0):
   ML_H_lat="S";
   ML_lat = -ML_lat;

(ML_lon_f,ML_lon_d) = math.modf(ML_lon);
(ML_lat_f,ML_lat_d) = math.modf(ML_lat);
print ('ML fix %d°%8.4f%s %d°%7.4f%s'%(ML_lon_d, ML_lon_f*60, ML_H_lon, ML_lat_d, ML_lat_f*60, ML_H_lat))
