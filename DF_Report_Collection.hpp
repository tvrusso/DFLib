#ifndef DF_REPORT_COLLECTION_HPP
#define DF_REPORT_COLLECTION_HPP
#include "port.h"

#include <vector>
using namespace std;

#include "Util_Abstract_Group.hpp"
#include "DF_Abstract_Report.hpp"
#include "DF_Abstract_Point.hpp"

namespace DFLib
{
  class ReportCollection : public DFLib::Abstract::Group
  {
  private: 
    vector<DFLib::Abstract::Report * > theReports;
    vector<double> evaluationPoint;
    bool f_is_valid;
    bool g_is_valid;
    bool h_is_valid;
    double function_value;
    vector<double> gradient;
    vector<vector<double> > hessian;

    // Declare the copy constructor and assignment operators, but
    // don't define them.  We should *never* copy a collection or attempt
    // to assign one to another.  This makes it illegal to do so.
    ReportCollection(ReportCollection &right);
    ReportCollection &operator=(ReportCollection &right);

  public:
    CPL_DLL ReportCollection();

    /// \brief DF Report Collection destructor
    ///
    ///  Never destroys the objects in its vector, as they might be getting
    ///  used for something else by caller
    virtual CPL_DLL ~ReportCollection();

    /// \brief destroy all reports stored in collection
    ///
    /// Provided in case our caller does NOT need the stored pointers for
    /// something else, and doesn't want to keep track of them.
    CPL_DLL void deleteReports();

    /// \brief Add a DF report to the collection
    ///
    /// \return this report's number in the collection.
    virtual CPL_DLL int addReport(DFLib::Abstract::Report * aReport);

    /// \brief return the fix cut average of this collection's reports
    ///
    /// A fix cut is the intersection of two DF reports.  The Fix Cut 
    /// Average is the average of all fix cuts from all pairs of reports 
    /// in the collection.  Fix cuts at shallow angles can be excluded by 
    /// specifying a non-zero value for minAngle (in degrees).
    ///
    /// \param FCA Returned fix cut average
    /// \param FCA_stddev standard deviation of fix cuts
    /// \param minAngle reports whose fix cut occur at less than this angle will not be included in the average.
    CPL_DLL bool computeFixCutAverage(DFLib::Abstract::Point &FCA, 
                                      vector<double> &FCA_stddev,
                                      double minAngle=0);

    /*!
      \brief Computes least squares solution of DF problem.
       
      The least squares solution is the point that has the minimum
      orthogonal distance to all bearing lines.
        
      The least squares solution is given by

      \f$ 
      P_{LS} = (A^TA)^{-1}A^Tb 
      \f$

      where  \f$A\f$ is the \f$Nx2\f$ matrix whose \f$k^{th}\f$ row
      is the unit vector orthogonal to the bearing line from receiver
      \f$k\f$.  Row \f$k\f$ of A is therefore the vector
      \f$[\cos(\theta_k), -\sin(\theta_k)]\f$ where \f$\theta_k\f$
      is the bearing from station \f$k\f$ measured clockwise from North.
      The \f$T\f$ superscript denotes matrix transpose.
        
      \f$b\f$ is the \f$N\f$ element vector whose \f$k^{th}\f$ element
      is the projection of the position vector of the \f$k^{th}\f$ 
      receiver onto the unit vector orthogonal to that receiver's 
      bearing line:  
      \f$
      b_k = [\cos(\theta_k),-sin(\theta_k)]\cdot[r_{kx},r_{ky}]
      \f$
      where \f${\bf r}_k\f$ is the position vector of the \f$k^{th}\f$
      receiver.
    */
    CPL_DLL void computeLeastSquaresFix(DFLib::Abstract::Point &LS_Fix);

    /*!
      \brief Computes Maximum Likelihood solution of DF problem

      Given a collection of DF fixes with specified standard deviation,
      computes the point of maximum likelihood (ML fix).

      The probability of a transmitter being at a particular location x,y
      given that each transmitter \f$i\f$ has heard the signal at bearing 
      \f$\theta_i\f$ is given by a multivariate gaussian probability
      distribution:

      \f$
      P(x,y) = K*\exp(-\sum_{i=0}^n (\tilde{\theta_i} - \theta_i(x,y))^2/(2\sigma_i^2))
      \f$

      where \f$\tilde{\theta_i}\f$ is the bearing measured by the
      \f$i^{th}\f$ receiver, \f$\theta_i(x,y)\f$ is the bearing from
      the \f$i^{th}\f$ receiver to the point \f$(x,y)\f$, \f$\sigma_i\f$
      is the standard deviation of the measurements by receiver \f$i\f$, and
      \f$K\f$ is a normalization coefficient.  The argument of the exponential
      is the cost function.

      The ML fix is the point that minimizes the cost function, thereby
      maximizing the probability of that point.

      Note that it is not correct to use \f$P(x,y)\f$ directly as a
      probability distribution in \f$x\f$ and \f$y\f$, because it is
      really the probability density in \f$\theta\f$ space.  To get
      the real probability density in \f$(x,y)\f$ space requires a
      change of variables factor \f$det(JJ^T)\f$ where \f$J\f$ is the
      jacobian of the transformation between \f$\theta\f$ and
      \f$(x,y)\f$ space.

      The correct transformation would be necessary if one were trying
      to compute the probability of the transmitter being in some area
      of \f$(x,y)\f$ space.  One potential future development in this 
      library would be to compute the probability on a grid around the ML fix,
      then find the convex hull containing points that sum to a particular
      confidence level.  This would be very time consuming and maybe not
      worth the trouble.

      To compute the ML fix requires an initial guess, as from the least
      squares fix.  We use the method of conjugate gradients to find the 
      minimum of the cost function.

    */

    CPL_DLL void computeMLFix(DFLib::Abstract::Point &MLFix);

    CPL_DLL double computeCostFunction(vector<double> &evaluationPoint);
    CPL_DLL void computeCostFunctionAndGradient(vector<double> &evaluationPoint,
                                                double &f,
                                                vector<double> &gradf);
    CPL_DLL void computeCostFunctionAndHessian(vector<double> &evaluationPoint,
                                               double &f,
                                               vector<double> &gradf,
                                               vector<vector<double> > &h);

    inline CPL_DLL void toggleValidity(int i)
    {
      if (i<theReports.size()&&i>=0)
        theReports[i]->toggleValidity();
    }

    inline CPL_DLL bool isValid(int i) const
    {
      return(theReports[i]->isValid());
    }

    inline CPL_DLL int size() {return theReports.size();};

    // This version returns something the caller can never use to change
    // a report from underneath us.
    inline CPL_DLL const DFLib::Abstract::Report * const getReport(int i) const
    { return (theReports[i]); };

    inline CPL_DLL const vector<double> &getReceiverLocationXY(int i) 
    { return (theReports[i]->getReceiverLocation()); };

    inline virtual CPL_DLL void setEvaluationPoint(vector<double> &ep)
    {
      evaluationPoint = ep;
      f_is_valid=false;
      g_is_valid=false;
      h_is_valid=false;
    }

    /// \return function value
    inline virtual CPL_DLL double getFunctionValue()
    {
      if (!f_is_valid)
      {
        function_value=computeCostFunction(evaluationPoint);
        f_is_valid=true;
      }
      return function_value;
    }

    /// \param g gradient returned
    /// \return function value
    inline virtual CPL_DLL double getFunctionValueAndGradient(vector<double> &g)
    {
      if (!g_is_valid)
      {
        computeCostFunctionAndGradient(evaluationPoint,
                                       function_value,gradient);
        f_is_valid=true;
        g_is_valid=true;
      }
      g=gradient;
      return function_value;
    }

    /// \param g gradient returned
    /// \param h hessian returned
    /// \return function value
    inline virtual CPL_DLL double 
    getFunctionValueAndHessian(vector<double> &g,
                               vector<vector<double> > &h)
    {
      if (!h_is_valid)
      {
        computeCostFunctionAndHessian(evaluationPoint,
                                      function_value,gradient,
                                      hessian);
        f_is_valid=true;
        g_is_valid=true;
        h_is_valid=true;
      }
      g=gradient;
      h=hessian;
      return function_value;
    }
  };
}
#endif // DF_REPORT_COLLECTION_HPP
