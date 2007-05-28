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
        class Report 
            : public DFLib::Abstract::Report
        {
        private:
            Point receiverLocation;            
            double bearing,sigma;
        public:
            CPL_DLL Report(const double &theLon, const double &theLat, 
                       const double &bearing,const double &std_dev);
            CPL_DLL Report(const vector<double> &theLocation, 
                       const double &bearing,const double &std_dev,
                           const string &theName);
            CPL_DLL ~Report();
            virtual  const CPL_DLL  vector<double> &getReceiverLocation();
            virtual CPL_DLL  double getReportBearingRadians() const;
            virtual CPL_DLL  double getBearing() const;
            virtual CPL_DLL  double getBearingStandardDeviationRadians() const;
            virtual CPL_DLL  double getSigma() const;
            virtual CPL_DLL  void  setReceiverLocation(vector<double> &theLocation);
            //! set bearing in degrees
            virtual CPL_DLL  void  setBearing(double Bearing);
            //! set standard deviation in degrees
            virtual CPL_DLL  void  setSigma(double Sigma);
        };
    }


    /// \brief XY DF report constructor with defaults
    inline DFLib::XY::Report::Report(const double &theX = 0,
                                     const double &theY = 0,
                                     const double &Bearing = 0,
                                     const double &std_dev = 1)
        : bearing(Bearing), sigma(std_dev*M_PI/180)
    {
        vector<double> XY(2);
        string crappe("");
        XY[0]=theX;
        XY[1]=theY;
        receiverLocation.setXY(XY);
        setReportName(crappe);
    }

    /// \brief XYDF report constructor.
    /// \param theLocation position vector of this report.
    /// \param Bearing bearing IN DEGREES
    /// \param std_dev standard deviation in degrees
    inline DFLib::XY::Report::Report(const vector<double> &theLocation,
                                  const double &Bearing,const double &std_dev,
                                     const string &theName)
        : receiverLocation(theLocation),
          bearing(Bearing*M_PI/180.0),
          sigma(std_dev*M_PI/180.0)
    {
        setReportName(theName);
        // Make sure our bearing is *always* 0<bearing<2pi.  If it isn't,
        // reset it:
        while (bearing < 0)
            bearing += 2*M_PI;
    }        

    inline DFLib::XY::Report::~Report()
    {
    }

    inline double DFLib::XY::Report::getReportBearingRadians() const
    {
        // bearing *must* be in 0<bearing<pi
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
        // bearing *must* be in 0<bearing<pi
        bearing=Bearing*M_PI/180.0;
        while (bearing < 0)
            bearing += 2*M_PI;
    }
    inline void DFLib::XY::Report::setSigma(double Sigma)
    {
        sigma=Sigma*M_PI/180;
    }
}
#endif // DF_XY_REPORT_HPP
