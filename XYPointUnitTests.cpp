//    DFLib: A library of Bearings Only Target Localization algorithms
//    Copyright (C) 2009-2011  Thomas V. Russo
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile$
//
// Purpose        : Trivial set of tests to make sure the XY class methods
//                  actually do what they're expected to do.
//
// Special Notes  : 
//
// Creator        : 
//
// Creation Date  : 
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision$
//
// Revision Date  : $Date$
//
// Current Owner  : $Author$
//-------------------------------------------------------------------------
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
