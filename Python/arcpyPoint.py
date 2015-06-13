"""
   This simple Python script is meant to demonstrate how to use the DFLib
   python bindings with Point and Report classes implemented in Python.

   You do not have to build the C++ library version of DFLib to use this, only
   the Python library:

   python setup.py build_ext --inplace [--swig <path to swig>]

   This assumes you actually have installed swig, and are set up to compile
   C++ codes.  This should create the DFLib python module.


"""
import arcpy
import DFLib
import copy
import math

#  The two classes below constitute a simplistic implementation of the minimum
#  work you'd have to do to use DFLib in Python.
#
#  The "arcpyPoint" class derives from the DFLib::Abstract::Point interface,
#  and implements its interface methods:
#    __init__ : the constructor for the class
#    setXY:     Sets coordinates in "XY" space, the one that is actually used
#               by the ReportCollection class for computation.
#    getXY:     Returns an STL vector containing the XY position.
#    Clone:     Makes a new arcpyPoint and populates it with the values in the
#               current object, returning a pointer to it.  This method is ONLY
#               used by the "Fix Cut Average" computation.
#    setUserCoords:  Allows you to specify the point location in some coordinate
#                    system other than the XY coordinates used for computation.
#                    In this implementation, XY is Mercator Projection 
#                    coordinates, and UserCoords are lat/lon.
#    getUserCoords:  Returns the point location in the user coordinates.

class arcpyPoint(DFLib.Point):
   """
   The "arcpyPoint" class derives from the DFLib::Abstract::Point interface,
   and implements its interface methods
   """
   def __init__(self,sr,*args):
      """
      arcpyPoint constructor
      """
      super(arcpyPoint,self).__init__()
      self.myXY=DFLib.vectord(2)
      self.myUserCoords=DFLib.vectord(2)
      self.myUserSR=sr
      self.myMercatorSR=arcpy.SpatialReference("WGS 1984 World Mercator")

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
      foo=arcpyPoint(self.myUserSR,self.getUserCoords()).__disown__()
      return foo
   def getUserCoords(self):
      """
      Return a vector of coordinates of this point in a user coordinate
      system (which may be different from the XY coordinate system 
      used for computation of fixes).
      """
      gM=arcpy.PointGeometry(arcpy.Point(self.myXY[0],self.myXY[1]),self.myMercatorSR)
      gU=gM.projectAs(self.myUserSR)
      self.myUserCoords[0]=gU.firstPoint.X
      self.myUserCoords[1]=gU.firstPoint.Y
      return self.myUserCoords
   def setUserCoords(self,uPosition):
      """
      Set the position of this point using a vector or list of coordinates in
      the user coordinate system.
      """
      gU=arcpy.PointGeometry(arcpy.Point(uPosition[0],uPosition[1]),self.myUserSR)
      gM=gU.projectAs(self.myMercatorSR)
      self.setXY((gM.firstPoint.X,gM.firstPoint.Y))

class arcpyReport(DFLib.Report):
   """
     The arcpyReport class implements the DFLib::Abstract::Report interface.
   """
   def __init__(self,name,valid,Point,bearing,sigma):
      """
      arcpyReport constructor
      """
      super(arcpyReport,self).__init__(name,valid)
      self.myPoint=Point.Clone()
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
# The bearings are the product of a single run of testlsDFfix, randomizing
# around the true bearings to a transmitter at 106d35'W 34d55'N
# according to the sigmas associated with the transmitters
#

theSR=arcpy.SpatialReference("WGS 1984")

r1=arcpyReport("r1",True,arcpyPoint(theSR,[-(106+38./60.+15.4/3600.) , (34+56./60.+48.7/3600)]),123.866,1.0)
r1p=r1.getReceiverLocation()

r2=arcpyReport("r2",True,arcpyPoint(theSR,[-(106+18./60.+17.0/3600.) , (34+58./60.+12.5/3600)]),-104.265,3.0)
r2p=r2.getReceiverLocation()

r3=arcpyReport("r3",True,arcpyPoint(theSR,[-(106+28./60.+53.8/3600.) , (35+10./60.+33.9/3600)]),-160.354,2.0)
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

LSfix=arcpyPoint(theSR)
rc.computeLeastSquaresFix(LSfix)

LSv=LSfix.getXY()
LSu=LSfix.getUserCoords()
print 'ls fix: x=%f y=%f'%(LSv[0],LSv[1])
print 'ls fix(LL): x=%f y=%f'%(LSu[0],LSu[1])

FCAfix=arcpyPoint(theSR)
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
