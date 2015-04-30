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
// Purpose        : Implements the DFLib::Abstract::Report interface such that
//                  the "user" coordinate system is specified using a proj.4
//                  spatial reference system and the "XY" coordinate system
//                  is a mercator projection on the WGS84 ellipsoid.
//
// Special Notes  : Use of this class treats bearing lines to target 
//                  transmitter as loxodromes, not geodesics.  Technically,
//                  that is incorrect.  For practical VHF target localization,
//                  it is probably good enough.  To implement a Report class
//                  that implements completely correct geodesic calculations
//                  would be much harder, and I haven't done it yet.  One
//                  would require an X-Y coordinate system in which straight 
//                  lines are geodesics rather than loxodromes, and there is no
//                  conformal projection that does that.  
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

#include "DF_Proj_Point.hpp" 
#include "DF_Proj_Report.hpp"

namespace DFLib
{
  namespace Proj
  {
    

  /// \brief DF report constructor for user-specified coordinate system.
  /// \param theLocation position vector <em>in lat/lon</em> of this report.
  /// \param Bearing bearing IN DEGREES
  /// \param std_dev standard deviation in degrees
  /// \param projArgs a vector of strings to pass to pj_init in order to 
  ///        define the user coordinate system.
    Report::Report(const std::vector<double> &theLocation,
                   const double &Bearing,const double &std_dev,
                   const std::string &theName,
                   const std::vector<std::string> &projArgs)
      : DFLib::Abstract::Report(theName,true),
        bearing(Bearing*M_PI/180.0),
        sigma(std_dev*M_PI/180.0)
    {
      receiverLocation = new Point(theLocation,projArgs);
      // Make sure our bearing is *always* 0<bearing<2*pi.  If it isn't,
      // reset it:
      while (bearing < 0)
        bearing += 2*M_PI;

      while (bearing >= 2*M_PI)
        bearing -= 2*M_PI;

    }        

    /// \brief Copy constructor

     Report::Report(const DFLib::Proj::Report & right)
       : DFLib::Abstract::Report(right),
         bearing(right.bearing),
         sigma(right.sigma)
     {
       // do not copy the pointer to the point, make a copy and save a
       // pointer to it.
       receiverLocation = new Point(*(right.receiverLocation));
     }

    /// \brief assignment operator
    Report& Report::operator=(const Report& rhs)
    {
      if (this == &rhs) //  do nothing if assigning to same pointer
        return *this;

      setReportName(rhs.getReportName());
      if (rhs.isValid())
        setValid();
      else
        setInvalid();

      if (receiverLocation) delete receiverLocation;
      receiverLocation = new Point(*(rhs.receiverLocation));
      bearing=rhs.bearing;
      sigma=rhs.sigma;
      return *this;
    }

    Report::~Report()
    {
      if (receiverLocation)
        delete receiverLocation;
    }


    /// \brief return a copy (yes, a COPY) of the receiver point object
    Point Report::getReceiverPoint() const
    {
      return (*receiverLocation);
    }
  }
}
