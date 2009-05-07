#ifndef DF_XY_REPORT_HPP
#define DF_XY_REPORT_HPP
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include "port.h"
#include <cmath>

#include "DF_Abstract_Report.hpp"
#include "DF_XY_Point.hpp"

namespace DFLib
{
  namespace XY
  {
    /// \brief DF report in XY coordinates
    /// This class is meant for simplistic DF codes that work in X/Y 
    ///  coordinates already.  Mostly intended for testing the interface.
    class CPL_DLL Report 
      : public DFLib::Abstract::Report
    {
    private:
      Point receiverLocation;            
      double bearing,sigma;
    public:
      Report(const vector<double> &theLocation, 
                     const double &bearing,const double &std_dev,
                     const string &theName);
      ~Report();
      virtual  const  vector<double> &getReceiverLocation();
      virtual  double getReportBearingRadians() const;
      virtual  double getBearing() const;
      virtual  double getBearingStandardDeviationRadians() const;
      virtual  double getSigma() const;
      virtual  void  setReceiverLocation(vector<double> &theLocation);
      //! set bearing in degrees
      virtual  void  setBearing(double Bearing);
      //! set standard deviation in degrees
      virtual  void  setSigma(double Sigma);
    };
  }


  /// \brief XYDF report constructor.
  /// \param theLocation position vector of this report.
  /// \param Bearing bearing IN DEGREES
  /// \param std_dev standard deviation in degrees
  inline DFLib::XY::Report::Report(const vector<double> &theLocation,
                                   const double &Bearing,const double &std_dev,
                                   const string &theName)
    : DFLib::Abstract::Report(theName,true),
      receiverLocation(theLocation),
      bearing(Bearing*M_PI/180.0),
      sigma(std_dev*M_PI/180.0)
  {

    // Make sure our bearing is *always* 0<=bearing<2pi.  If it isn't,
    // reset it:
    while (bearing < 0)
      bearing += 2*M_PI;

    while (bearing >= 2*M_PI)
      bearing -= 2*M_PI;


  }        

  inline DFLib::XY::Report::~Report()
  {
  }

  inline double DFLib::XY::Report::getReportBearingRadians() const
  {
    // bearing *must* be in 0<bearing<2*pi
    return bearing;
  }
  inline double DFLib::XY::Report::getBearing() const
  {
    // always return bearing in range 0-360 degrees
    return (getReportBearingRadians()*180.0/M_PI);
  }

  inline double DFLib::XY::Report::getBearingStandardDeviationRadians() const
  {
    return sigma;
  }
  inline double DFLib::XY::Report::getSigma() const
  {
    return sigma*180/M_PI;
  }

  inline void DFLib::XY::Report::setReceiverLocation(vector<double> &theLocation)
  {
    receiverLocation.setXY(theLocation);
  }

  inline const vector<double> & DFLib::XY::Report::getReceiverLocation() 
  { 
    return receiverLocation.getXY();
  }

  inline void DFLib::XY::Report::setBearing(double Bearing)
  {
    // bearing *must* be in 0<=bearing<2*pi
    bearing=Bearing*M_PI/180.0;
    while (bearing < 0)
      bearing += 2*M_PI;
    while (bearing >= 2*M_PI)
      bearing -= 2*M_PI;
  }
  inline void DFLib::XY::Report::setSigma(double Sigma)
  {
    sigma=Sigma*M_PI/180;
  }
}
#endif // DF_XY_REPORT_HPP
