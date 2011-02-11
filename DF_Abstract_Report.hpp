#ifndef DF_ABSTRACT_REPORT_HPP
#define DF_ABSTRACT_REPORT_HPP
#include "DFLib_port.h"
#include <vector>
#include <cmath>
#include <string>

#include "DF_Abstract_Point.hpp"

using namespace std;

namespace DFLib 
{
  enum FixStatus {NO_DATA,GOOD_FIX,NO_FIX};


  namespace Abstract
  {
    // Forward declaration:
    /// \brief Abstract interface for DFReport.
    ///  This is the base class for DF report classes. 
    class CPL_DLL Report
    {
    private:
      string ReportName_;
      bool validReport_;
    public:
      // pure virtual functions:

      /// \brief virtual destructor because base classes need one
      virtual ~Report() {};

      /// \brief constructor for base Report class
      ///
      /// This cannot be used to create an abstract report, because this is
      /// an abstract class.
      Report(string n,bool v);

      /// \brief copy constructor for base Report class
      ///
      /// This cannot be used to create an abstract report, because this is
      /// an abstract class.
      Report(const Report & right);

      /// \brief return receiver location in double vector of XY coords
      virtual  const vector<double> &getReceiverLocation() = 0;
      /// \brief return reported bearing to target
      ///
      /// It is essential that getReportBearingRadians always return
      /// the bearing in the correct range \f$0<\theta<2\pi\f$.
      /// \return bearing in radians, always in the range \f$0<\theta<2\pi\f$.
      virtual  double getReportBearingRadians() const = 0;
      virtual  double getBearingStandardDeviationRadians() const = 0;

      ///\brief Return the name of this report
      virtual const string &getReportName() const { return ReportName_;};

      ///\brief set the name of this report
      virtual void setReportName(const string &theName) { ReportName_=theName;};

      ///\brief Set this report as valid
      virtual void setValid() { validReport_=true;};
      ///\brief Set this report as invalid
      virtual void setInvalid() { validReport_=false;};

      virtual void toggleValidity() { validReport_ = !validReport_;};

      ///\brief check this report's validity
      virtual bool isValid() const { return validReport_; };

      /// \brief compute point at which the line from this report intersects that from another
      /// \param Report2 pointer to the other report
      /// \param returnPoint reference to place to store solution point
      /// \param cutAngle on return, the angle in radians made by the two bearing lines
      /// \param fs fix status
      void computeFixCut(DFLib::Abstract::Report *Report2, 
                                 Point &returnPoint, 
                                 double &cutAngle,FixStatus &fs);
      /// \brief compute bearing from this reporting location to some other point
      ///
      /// \param aPoint point to which bearing requested.
      /// \return bearing in radians, in range \f$0<\theta<2\pi\f$
      ///
      /// The base class computes the bearing in X-Y coordinates by taking the 
      /// arctangent of "rise over run."  This is certainly the appropriate
      /// thing to do when the X-Y coordinate system is a conformal projection.
      /// When working with non-conformal projections the computation could possibly be more involved.

      double computeBearingToPoint(vector<double> &aPoint);

      /// \brief compute distance from this reporting location to some other point.
      /// \param aPoint point to which distance is requested
      /// \return distance
      ///
      /// The distance returned is the Cartesian distance in X-Y coordinates,
      /// \f$d = \sqrt{\Delta x^2 + \Delta y^2}\f$.  It is <em>not</em> the
      /// geodesic distance on the ellipsoid.
      ///  
      double computeDistanceToPoint(vector<double> &aPoint);
    };
  }

}
#endif // DF_ABSTRACT_REPORT_HPP
