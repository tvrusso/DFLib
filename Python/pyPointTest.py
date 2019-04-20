import DFLib
import copy
import math
from pyproj import Proj

class pyPoint(DFLib.Point):
   def __init__(self):
      super(pyPoint,self).__init__()
      self.myXY=DFLib.vectord(2)
      self.myXY[0]=0
      self.myXY[1]=0
   def setXY(self,aPosition):
      self.myXY[0]=aPosition[0]
      self.myXY[1]=aPosition[1]
   def getXY(self):
     return self.myXY
  # This is NOT really a "clone" operation, but it serves the purpose enough
  # for DFLib
   def Clone(self):
     foo=pyPoint().__disown__()
     foo.setXY(self.getXY())
     return foo
   def getUserCoords(self):
     return self.getXY()
   def setUserCoords(self,uPosition):
     self.setXY(uPosition)

class pyReport(DFLib.Report):
   def __init__(self,name,valid):
      super(pyReport,self).__init__(name,valid)
      self.myXY=DFLib.vectord(2)
      self.myXY[0]=0
      self.myXY[1]=0
      self.bearing=0
      self.sigma=0
   def setReceiverLocation(self,theLocation):
      self.myXY[0]=theLocation[0]
      self.myXY[1]=theLocation[1]
   def getReceiverLocation(self):
      return self.myXY
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
mercProj = Proj('+proj=merc +datum=WGS84 +lat_ts=0')   
v=DFLib.vectord(2)
v[0],v[1]=mercProj(-(106+38./60.+15.4/3600.),34+56./60.+48.7/3600)
r1=pyReport("r1",True)
r1.setReceiverLocation(v)
r1.setBearing(123.866)
r1.setSigma(1.0)
r1p=r1.getReceiverLocation()

v[0],v[1]=mercProj(-(106+18./60.+17.0/3600.),34+58./60.+12.5/3600)
r2=pyReport("r2",True)
r2.setReceiverLocation(v)
r2.setBearing(-104.265)
r2.setSigma(3.0)
r2p=r2.getReceiverLocation()

v[0],v[1]=mercProj(-(106+28./60.+53.8/3600.),35+10./60.+33.9/3600)
r3=pyReport("r3",True)
r3.setReceiverLocation(v)
r3.setBearing(-160.354)
r3.setSigma(2.0)
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

LSfix=pyPoint()
rc.computeLeastSquaresFix(LSfix)

LSv=LSfix.getXY()
print('ls fix: x=%f y=%f'%(LSv[0],LSv[1]))
print('ls fix(LL): x=%f y=%f'%mercProj(LSv[0],LSv[1],inverse=True))

FCAfix=pyPoint()
FCA_stddev=DFLib.vectord(2)
rc.computeFixCutAverage(FCAfix,FCA_stddev)
   
FCAv=FCAfix.getXY()

print('FCA fix: x=%f y=%f'%(FCAv[0],FCAv[1]))
print('FCA fix(LL): x=%f y=%f'%mercProj(FCAv[0],FCAv[1],inverse=True))
print('FCA stddevs(merc) x=%f y=%f'%(FCA_stddev[0],FCA_stddev[1]))

MLfix = LSfix.Clone()
rc.aggressiveComputeMLFix(MLfix)


MLv=MLfix.getXY()
print('ML fix: x=%f y=%f'%(MLv[0],MLv[1]))
print('ML fix(LL): x=%f y=%f'%mercProj(MLv[0],MLv[1],inverse=True))
