#ifndef DF_PROJ_REPORT_HPP
#define DF_PROJ_REPORT_HPP
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include "port.h"
#include <cmath>

#include "DF_Abstract_Report.hpp"
#include "DF_Proj_Point.hpp"

namespace DFLib
{
  namespace Proj
  {
    /// \brief DF report in user-specified coordinates
    class Report 
      : public DFLib::Abstract::Report
    {
    private:
      Point *receiverLocation;            
      double bearing,sigma;
    public:
      CPL_DLL Report(const vector<double> &theLocationUser, 
                     const double &bearing,const double &std_dev,
                     const string &theName,const vector<string>&projArgs);


      CPL_DLL Report(const Report & right);
      Report & operator=(const Report& rhs);

      virtual CPL_DLL ~Report();
      virtual  const CPL_DLL  vector<double> &getReceiverLocation();
      virtual CPL_DLL  double getReportBearingRadians() const;
      virtual CPL_DLL  double getBearing() const;
      virtual CPL_DLL  double getBearingStandardDeviationRadians() const;
      virtual CPL_DLL  double getSigma() const;
      virtual CPL_DLL  void  setReceiverLocationUser(vector<double> &theLocation);
      virtual CPL_DLL  void  setReceiverLocationMercator(vector<double> &theLocation);
      //! set bearing in degrees
      virtual CPL_DLL  void  setBearing(double Bearing);
      //! set standard deviation in degrees
      virtual CPL_DLL  void  setSigma(double Sigma);
    };
  }

  inline DFLib::Proj::Report::~Report()
  {
    if (receiverLocation)
      delete receiverLocation;
  }

  inline double DFLib::Proj::Report::getReportBearingRadians() const
  {
    // bearing *must* be in 0<bearing<2*pi
    return bearing;
  }
  inline double DFLib::Proj::Report::getBearing() const
  {
    // always return bearing in range 0-360 degrees
    return (getReportBearingRadians()*180.0/M_PI);
  }

  inline double DFLib::Proj::Report::getBearingStandardDeviationRadians() const
  {
    return sigma;
  }
  inline double DFLib::Proj::Report::getSigma() const
  {
    return sigma*180/M_PI;
  }

  inline void DFLib::Proj::Report::setReceiverLocationUser(vector<double> &theLocation)
  {
    receiverLocation->setUserCoords(theLocation);
  }

  inline void DFLib::Proj::Report::setReceiverLocationMercator(vector<double> &theLocation)
  {
    receiverLocation->setXY(theLocation);
  }

  inline const vector<double> & DFLib::Proj::Report::getReceiverLocation() 
  { 
    return receiverLocation->getXY();
  }

  inline void DFLib::Proj::Report::setBearing(double Bearing)
  {
    // bearing *must* be in 0<bearing<2*pi
    bearing=Bearing*M_PI/180.0;
    while (bearing < 0)
      bearing += 2*M_PI;
    while (bearing >= 2*M_PI)
      bearing -= 2*M_PI;

  }
  inline void DFLib::Proj::Report::setSigma(double Sigma)
  {
    sigma=Sigma*M_PI/180;
  }
}
#endif // DF_PROJ_REPORT_HPP
