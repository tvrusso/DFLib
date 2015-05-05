import DFLib
import copy
import math

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
   def clone(self):
     foo=pyPoint()
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
   
v=DFLib.vectord(2)
v[0]=1.0
v[1]=2.5
p=pyPoint()
p.setXY(v)
p2=p.clone()
v[0]=3.0
v[1]=4.5
p2.setXY(v)

v2=p.getXY()
print 'v2: %f %f'%(v2[0],v2[1])
v3=p2.getXY()
print 'v3: %f %f'%(v3[0],v3[1])

r1=pyReport("r1",True)
r1.setReceiverLocation(v2)
r1.setBearing(175.0)
r1.setSigma(27.0)
r1p=r1.getReceiverLocation()

r2=pyReport("r2",True)
r2.setReceiverLocation(v3)
r2.setBearing(185.0)
r2.setSigma(18.0)
r2p=r2.getReceiverLocation()

print 'r1: (x,y)=(%f,%f)  bearing=%f sigma=%f'%(r1p[0],r1p[1],r1.getBearing(),r1.getSigma())
print 'r1: (x,y)=(%f,%f)  bearing(rad)=%f sigma(rad)=%f'%(r1p[0],r1p[1],r1.getReportBearingRadians(),r1.getBearingStandardDeviationRadians())
print 'r2: (x,y)=(%f,%f)  bearing=%f sigma=%f'%(r2p[0],r2p[1],r2.getBearing(),r2.getSigma())
print 'r2: (x,y)=(%f,%f)  bearing(rad)=%f sigma(rad)=%f'%(r2p[0],r2p[1],r2.getReportBearingRadians(),r2.getBearingStandardDeviationRadians())
