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
// Purpose        : set of tests to verify that the Proj classes do what 
//                  is expected of them.
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
#include <limits>

#include "Util_Misc.hpp"
#include "DF_Proj_Point.hpp"

main(int argc, char **argv)
{

  double dtol=10*sqrt(numeric_limits<double>::epsilon());
  std::vector<double> xyVals(2);
  xyVals[0]=-106.482;
  xyVals[1]=35.0913;
  std::cout << " Creating a point using values lon="<<xyVals[0]<<",lat= " <<xyVals[1]<< std::endl;

  std::vector<std::string> projArgs;
  projArgs.push_back("proj=latlong");
  projArgs.push_back("datum=WGS84");
  std::cout << " pushed projArgs into std::vector " << std::endl;

  try
  {
    DFLib::Proj::Point myFirstPoint(xyVals,projArgs);

    std::cout << " Point created. " << std::endl;

    std::vector<double> xyStored = myFirstPoint.getUserCoords();

    std::cout << " Retrieved XY from stored point.  Lon = " << xyStored[0]
         << " Lat = " << xyStored[1];

    if (xyStored[0] == xyVals[0] && xyStored[1] == xyVals[1])
      std::cout << " PASSED " << std::endl;
    else
      std::cout << " FAILED " << std::endl;

    xyVals = myFirstPoint.getXY();
    xyVals[0] += 200;
    xyVals[1] += 200;
    std::cout << " Resetting XY in stored point to X=" << xyVals[0]  << ",Y= "
         << xyVals[1] << std::endl;
    myFirstPoint.setXY(xyVals);
    xyStored = myFirstPoint.getXY();

    std::cout << " Retrieved XY from stored point.  X = " << xyStored[0]
         << " Y = " << xyStored[1];

    if (xyStored[0] == xyVals[0] && xyStored[1] == xyVals[1])
      std::cout << " PASSED " << std::endl;
    else
    {
      std::cout << " FAILED " << std::endl;
      std::cout << "       storedX=" << xyStored[0] << " Should be " << xyVals[0]
           << std::endl;
      std::cout << "       storedY=" << xyStored[1] << " Should be " << xyVals[1]
           << std::endl;
    }

    xyStored = myFirstPoint.getUserCoords();
    std::cout << " Lat Lon values of stored, modified point are " 
         << " Lon=" << xyStored[0]
         << " Lat=" << xyStored[1] << std::endl;

    std::cout << " Cloning point. " << std::endl;
    DFLib::Proj::Point *mySecondPointPtr = myFirstPoint.Clone();
    xyStored = mySecondPointPtr->getXY();

    std::cout << " Retrieved XY from cloned point.  X = " << xyStored[0]
         << " Y = " << xyStored[1];

    if (fabs(xyStored[0]-xyVals[0])< dtol && fabs(xyStored[1]-xyVals[1])<dtol)
      std::cout << " PASSED " << std::endl;
    else
    {
      std::cout << " FAILED " << std::endl;
      std::cout << "       storedX=" << xyStored[0] << " Should be " << xyVals[0]
           << " difference " << xyStored[0]-xyVals[0] << std::endl;
      std::cout << "       storedY=" << xyStored[1] << " Should be " << xyVals[1]
           << " difference " << xyStored[1]-xyVals[1] << std::endl;
      std::cout << "       Tolerance is " << dtol << std::endl;
        
    }
        
    std::cout << " Resetting XY in cloned point to X=0,Y=1 " << std::endl;
    xyVals[0]=0;
    xyVals[1]=1;
    mySecondPointPtr->setXY(xyVals);
    xyStored = mySecondPointPtr->getXY();

    std::cout << " Retrieved XY from cloned point.  X = " << xyStored[0]
         << " Y = " << xyStored[1];

    if (xyStored[0] == xyVals[0] && xyStored[1] == xyVals[1])
      std::cout << " PASSED " << std::endl;
    else
      std::cout << " FAILED " << std::endl;

    xyStored = mySecondPointPtr->getUserCoords();
    std::cout << " Retrieved LL from modified, cloned point.  Lon = " << xyStored[0]
         << " Lat = " << xyStored[1] << std::endl;
  }
  catch (DFLib::Util::Exception x)
  {
    std::cerr << "Ooops... got exception creating myFirstPoint" 
         << x.getEmsg() << std::endl;
  }
    
}
