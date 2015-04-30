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
// Purpose        : Implement the DFLib::Abstract::Report interface with a
//                  "user" coordinate system of WGS84 lat/lon and an "XY"
//                  coordinate system of WGS84 mercator.
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
#ifndef DF_LATLON_REPORT_HPP
#define DF_LATLON_REPORT_HPP
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include "DFLib_port.h"
#include <cmath>

#include "DF_Abstract_Report.hpp"

namespace DFLib
{
  namespace LatLon
  {

    class Point;

    /// \brief DF report in Lat/Lon coordinates
    class CPL_DLL Report 
      : public DFLib::Abstract::Report
    {
    private:
      Point receiverLocation;            
      double bearing,sigma;
    public:
      Report(const std::vector<double> &theLocationLL, 
                     const double &bearing,const double &std_dev,
                     const std::string &theName);
      ~Report();
      virtual  const  std::vector<double> &getReceiverLocation();
      virtual  double getReportBearingRadians() const;
      virtual  double getBearing() const;
      virtual  double getBearingStandardDeviationRadians() const;
      virtual  double getSigma() const;
      virtual  void  setReceiverLocationLL(std::vector<double> &theLocation);
      virtual  void  setReceiverLocationMercator(std::vector<double> &theLocation);
      //! set bearing in degrees
      virtual  void  setBearing(double Bearing);
      //! set standard deviation in degrees
      virtual  void  setSigma(double Sigma);
    };
  }

  /// \brief Lat/Lon DF report constructor.
  /// \param theLocation position vector <em>in lat/lon</em> of this report.
  /// \param Bearing bearing IN DEGREES
  /// \param std_dev standard deviation in degrees
  inline DFLib::LatLon::Report::Report(const std::vector<double> &theLocation,
                                       const double &Bearing,const double &std_dev,
                                       const std::string &theName)
    :  DFLib::Abstract::Report(theName,true),
       receiverLocation(theLocation),
      bearing(Bearing*M_PI/180.0),
      sigma(std_dev*M_PI/180.0)
  {

    // Make sure our bearing is *always* 0<bearing<2pi.  If it isn't,
    // reset it:
    while (bearing < 0)
      bearing += 2*M_PI;
    while (bearing >= 2*M_PI)
      bearing -= 2*M_PI;

  }        

  inline DFLib::LatLon::Report::~Report()
  {
  }

  inline double DFLib::LatLon::Report::getReportBearingRadians() const
  {
    // bearing *must* be in 0<bearing<2*pi
    return bearing;
  }
  inline double DFLib::LatLon::Report::getBearing() const
  {
    // always return bearing in range 0-360 degrees
    return (getReportBearingRadians()*180.0/M_PI);
  }

  inline double DFLib::LatLon::Report::getBearingStandardDeviationRadians() const
  {
    return sigma;
  }
  inline double DFLib::LatLon::Report::getSigma() const
  {
    return sigma*180/M_PI;
  }

  inline void DFLib::LatLon::Report::setReceiverLocationLL(std::vector<double> &theLocation)
  {
    receiverLocation.setLL(theLocation);
  }

  inline void DFLib::LatLon::Report::setReceiverLocationMercator(std::vector<double> &theLocation)
  {
    receiverLocation.setXY(theLocation);
  }

  inline const std::vector<double> & DFLib::LatLon::Report::getReceiverLocation() 
  { 
    return receiverLocation.getXY();
  }

  inline void DFLib::LatLon::Report::setBearing(double Bearing)
  {
    // bearing *must* be in 0<bearing<2*pi
    bearing=Bearing*M_PI/180.0;
    while (bearing < 0)
      bearing += 2*M_PI;
    while (bearing >= 2*M_PI)
      bearing -= 2*M_PI;
  }
  inline void DFLib::LatLon::Report::setSigma(double Sigma)
  {
    sigma=Sigma*M_PI/180;
  }
}
#endif // DF_LATLON_REPORT_HPP
