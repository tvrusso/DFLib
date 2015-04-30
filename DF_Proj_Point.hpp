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
// Purpose        : Implement a DFLib::Abstract::Point interface such that
//                  the "user" coordinate system is whatever the user has
//                  specified with a proj.4 spatial reference system, and the
//                  "XY" coordinate system is a mercator projection on the 
//                  WGS84 ellipsoid.
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
#ifndef DF_PROJ_POINT_HPP
#define DF_PROJ_POINT_HPP

#include "DFLib_port.h"
#include <vector>
#include <string>
#include "proj_api.h"
#include "DF_Abstract_Point.hpp"

namespace DFLib
{
  namespace Proj
  {
    class CPL_DLL Point : public DFLib::Abstract::Point
    {
    private:
      std::vector<double> theMerc;
      bool mercDirty;
      std::vector<double> theUserCoords;
      bool userDirty;
      projPJ userProj, mercProj;
    public:

      /// \brief Constructor
      ///
      /// \param uPosition coordinates <em>in user's coordinates</em>
      /// \param projArgs a vector of strings representing the Proj.4
      /// description of the user's coordinate system.
      ///
      Point(const std::vector<double> &uPosition,const std::vector<std::string> &projArgs);


      /// \brief Copy Constructor
      Point(const Point &right);

      /// \brief destructor
      ~Point();

      /// \brief assignment operator
      Point& operator=(const Point& rhs);

      /// \brief set user projection information
      ///
      /// \param projArgs vector of strings representing the new projection
      ///
      /// This method can be used to change the coordinate system of a point.
      /// When called, the Mercator representation of the point is updated,
      /// then the user projection is changed.  The next call to getUserCoords
      /// will therefore return the point's coordinates in the new coordinate
      /// system.

      void setUserProj(const std::vector<std::string> &projArgs);

      /// \brief return true if user projection is a lat/lon system
      ///
      ///  Primarily useful for deciding how to display coordinates 
      ///  (as for example whether to display them in sexigesimal or
      ///  decimal representation)
      bool isUserProjLatLong() const;

      /// \brief set mercator projection (XY) position
      /// 
      /// \param aPosition coordinates in mercator projection
      ///
      virtual void setXY(const std::vector<double> &mPosition);
      /// \brief get mercator projection (XY) position
      /// 
      /// \return vector of coordinates in mercator projection
      ///
      virtual const std::vector<double> &getXY();


      /// \brief get user position
      ///
      /// \return a vector<double> containing the coordinates of this
      /// point in the user's coordinate system.  
      ///
      /// This method works by converting the mercator coordinates of
      /// the point back to the user's coordinate system using the Proj.4
      /// cartographic projection library if the mercator coordinates 
      /// or user projection have changed since the last call. 

      virtual const std::vector<double> &getUserCoords();

      /// \brief set user position
      ///
      /// Sets the coordinates of the point in user's coordinate system,
      /// marking the mercator coordinates invalid.  Mercator coordinates
      /// will be updated by conversion with Proj.4 upon the next call to
      /// getXY().
      ///
      /// \param uPosition vector of doubles with user coordinates.

      virtual void setUserCoords(const std::vector<double> &uPosition)  ;

      virtual Point * Clone();
    private:
      void mercToUser();
      void userToMerc();
    };
  }
}
#endif // DF_PROJ_POINT_HPP
