#include <iostream>
#include <vector>

using namespace std;

#include "DF_XY_Point.hpp"

main(int argc, char **argv)
{

  vector<double> xyVals(2);
  xyVals[0]=1;
  xyVals[1]=2;
  cout << " Creating a point using values x=1,y=2 " << endl;
    
  DFLib::XY::Point myFirstPoint(xyVals);
  cout << " Point created. " << endl;

  vector<double> xyStored = myFirstPoint.getXY();

  cout << " Retrieved XY from stored point.  X = " << xyStored[0]
       << " Y = " << xyStored[1];

  if (xyStored[0] == xyVals[0] && xyStored[1] == xyVals[1])
    cout << " PASSED " << endl;
  else
    cout << " FAILED " << endl;

  cout << " Resetting XY in stored point to X=2,Y=2 " << endl;
  xyVals[0]=xyVals[1]=2;
  myFirstPoint.setXY(xyVals);
  xyStored = myFirstPoint.getXY();

  cout << " Retrieved XY from stored point.  X = " << xyStored[0]
       << " Y = " << xyStored[1];

  if (xyStored[0] == xyVals[0] && xyStored[1] == xyVals[1])
    cout << " PASSED " << endl;
  else
    cout << " FAILED " << endl;

  cout << " Cloning point. " << endl;
  DFLib::XY::Point *mySecondPointPtr = myFirstPoint.Clone();
  xyStored = mySecondPointPtr->getXY();

  cout << " Retrieved XY from cloned point.  X = " << xyStored[0]
       << " Y = " << xyStored[1];

  if (xyStored[0] == xyVals[0] && xyStored[1] == xyVals[1])
    cout << " PASSED " << endl;
  else
    cout << " FAILED " << endl;
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

}
