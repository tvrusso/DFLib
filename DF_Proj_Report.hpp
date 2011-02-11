#ifndef DF_PROJ_REPORT_HPP
#define DF_PROJ_REPORT_HPP
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include "DFLib_port.h"
#include <cmath>
using namespace std;

#include "DF_Abstract_Report.hpp"
#include "DF_Proj_Point.hpp"

namespace DFLib
{
  namespace Proj
  {
    /// \brief DF report in user-specified coordinates
    class CPL_DLL Report 
      : public DFLib::Abstract::Report
    {
    private:
      Point *receiverLocation;            
      double bearing,sigma;
    public:
      Report(const vector<double> &theLocationUser, 
                     const double &bearing,const double &std_dev,
                     const string &theName,const vector<string>&projArgs);


      Report(const Report & right);
      Report & operator=(const Report& rhs);

      virtual ~Report();
      virtual  const  vector<double> &getReceiverLocation();
      virtual Point getReceiverPoint() const;
      virtual  double getReportBearingRadians() const;
      virtual  double getBearing() const;
      virtual  double getBearingStandardDeviationRadians() const;
      virtual  double getSigma() const;
      virtual  void  setReceiverLocationUser(const vector<double> &theLocation);
      virtual  void  setReceiverLocationMercator(const vector<double> &theLocation);
      //! set bearing in degrees
      virtual  void  setBearing(double Bearing);
      //! set standard deviation in degrees
      virtual  void  setSigma(double Sigma);
      //! allow us to change the projection of the user coordinates
      virtual void setUserProj(const vector<string> &projArgs);
    };
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

  inline void DFLib::Proj::Report::setReceiverLocationUser(const vector<double> &theLocation)
  {
    receiverLocation->setUserCoords(theLocation);
  }

  inline void DFLib::Proj::Report::setReceiverLocationMercator(const vector<double> &theLocation)
  {
    receiverLocation->setXY(theLocation);
  }

  inline void DFLib::Proj::Report::setUserProj(const vector<string> &projArgs)
  {
    receiverLocation->setUserProj(projArgs);
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
