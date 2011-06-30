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
#ifndef DF_XY_POINT_HPP
#define DF_XY_POINT_HPP

#include "DFLib_port.h"

#include "DF_Abstract_Point.hpp"

namespace DFLib
{
  namespace XY
  {
    class CPL_DLL Point : public DFLib::Abstract::Point
    {
    private:
      vector<double> myXY;
    public:
      /// Default
      Point();

      /// Constructor
      Point(const vector<double> &aPosition);
      /// Copy Constructor
      Point(Point &right);
      // no destructor needed
      /// Set position from X-Y
      virtual void setXY(const vector<double> &aPosition);
      /// Set X-Y position
      virtual const vector<double> &getXY();

      /// get User Coords (wrapper as required by abstract interface)
      virtual const vector<double> &getUserCoords() { return getXY();};

      /// set User Coords (wrapper as required by abstract interface)
      virtual void setUserCoords(const vector<double> &uPosition) 
      { setXY(uPosition);};

      /// Clone Self
      virtual Point *Clone();
    };
  }
}

#endif
