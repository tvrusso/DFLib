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
// Purpose        : Simplistic implementation of "Point" interface where
//                  the "user" and "XY" coordinate systems are the same.
//
// Special Notes  : Used for DFLib::XY::Report class, meaning its use is 
//                  pretty much limited to the same sort of theoretical and
//                  testing purposes of that class.
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
#include "DF_XY_Point.hpp"

#include <vector>
using namespace std;

namespace DFLib
{
  namespace XY
  {

    // class Point

    Point::Point()
    {
      myXY.resize(2,0.0);
    }

    Point::Point(const vector<double> &aPosition)
      :myXY(aPosition)
    {
    }

    Point::Point(Point &right)
      :myXY(right.myXY)
    {
    }

    void Point::setXY(const vector<double> &aPosition)
    {
      myXY = aPosition;
    }

    const vector<double> & Point::getXY()
    {
      return(myXY);
    }

    Point * Point::Clone()
    {
      Point *retPoint;
      retPoint = new Point(*this);
      return retPoint;
    }
  }
}
