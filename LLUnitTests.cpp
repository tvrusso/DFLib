#include <iostream>
#include <vector>
#include <limits>
using namespace std;

#include "DF_LatLon_Point.hpp"

main(int argc, char **argv)
{

  double dtol=10*sqrt(numeric_limits<double>::epsilon());
  vector<double> xyVals(2);
  xyVals[0]=-106.482;
  xyVals[1]=35.0913;
  cout << " Creating a point using values x=1,y=2 " << endl;
    
  DFLib::LatLon::Point myFirstPoint(xyVals);
  cout << " Point created. " << endl;

  vector<double> xyStored = myFirstPoint.getLL();

  cout << " Retrieved XY from stored point.  Lon = " << xyStored[0]
       << " Lat = " << xyStored[1];

  if (xyStored[0] == xyVals[0] && xyStored[1] == xyVals[1])
    cout << " PASSED " << endl;
  else
    cout << " FAILED " << endl;

  xyVals = myFirstPoint.getXY();
  xyVals[0] += 200;
  xyVals[1] += 200;
  cout << " Resetting XY in stored point to X=" << xyVals[0]  << ",Y= "
       << xyVals[1] << endl;
  myFirstPoint.setXY(xyVals);
  xyStored = myFirstPoint.getXY();

  cout << " Retrieved XY from stored point.  X = " << xyStored[0]
       << " Y = " << xyStored[1];

  if (xyStored[0] == xyVals[0] && xyStored[1] == xyVals[1])
    cout << " PASSED " << endl;
  else
  {
    cout << " FAILED " << endl;
    cout << "       storedX=" << xyStored[0] << " Should be " << xyVals[0]
         << endl;
    cout << "       storedY=" << xyStored[1] << " Should be " << xyVals[1]
         << endl;
  }

  xyStored = myFirstPoint.getLL();
  cout << " Lat Lon values of stored, modified point are " 
       << " Lon=" << xyStored[0]
       << " Lat=" << xyStored[1] << endl;

  cout << " Cloning point. " << endl;
  DFLib::LatLon::Point *mySecondPointPtr = myFirstPoint.Clone();
  xyStored = mySecondPointPtr->getXY();

  cout << " Retrieved XY from cloned point.  X = " << xyStored[0]
       << " Y = " << xyStored[1];

  if (fabs(xyStored[0]-xyVals[0])< dtol && fabs(xyStored[1]-xyVals[1])<dtol)
    cout << " PASSED " << endl;
  else
  {
    cout << " FAILED " << endl;
    cout << "       storedX=" << xyStored[0] << " Should be " << xyVals[0]
         << " difference " << xyStored[0]-xyVals[0] << endl;
    cout << "       storedY=" << xyStored[1] << " Should be " << xyVals[1]
         << " difference " << xyStored[1]-xyVals[1] << endl;
    cout << "       Tolerance is " << dtol << endl;
        
  }
        
  cout << " Resetting XY in cloned point to X=0,Y=1 " << endl;
  xyVals[0]=0;
  xyVals[1]=1;
  mySecondPointPtr->setXY(xyVals);
  xyStored = mySecondPointPtr->getXY();

  cout << " Retrieved XY from cloned point.  X = " << xyStored[0]
       << " Y = " << xyStored[1];

  if (xyStored[0] == xyVals[0] && xyStored[1] == xyVals[1])
    cout << " PASSED " << endl;
  else
    cout << " FAILED " << endl;

  xyStored = mySecondPointPtr->getLL();
  cout << " Retrieved LL from modified, cloned point.  Lon = " << xyStored[0]
       << " Lat = " << xyStored[1] << endl;
    
}
