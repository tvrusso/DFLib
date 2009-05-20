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
  class CPL_DLL ReportCollection : public DFLib::Abstract::Group
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
    ReportCollection();

    /// \brief DF Report Collection destructor
    ///
    ///  Never destroys the objects in its vector, as they might be getting
    ///  used for something else by caller
    virtual ~ReportCollection();

    /// \brief destroy all reports stored in collection
    ///
    /// Provided in case our caller does NOT need the stored pointers for
    /// something else, and doesn't want to keep track of them.
    virtual void deleteReports();

    /// \brief Add a DF report to the collection
    ///
    /// \return this report's number in the collection.
    virtual int addReport(DFLib::Abstract::Report * aReport);

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
    bool computeFixCutAverage(DFLib::Abstract::Point &FCA, 
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
    void computeLeastSquaresFix(DFLib::Abstract::Point &LS_Fix);

    /*!
      \brief Computes Stansfield estimate of Maximum Likelihood solution of DF problem

      Given a collection of DF fixes with specified standard
      deviation, computes the point of maximum likelihood in the
      approximation of small angles.  It is based on the paper
      "Statistical Theory of D.F. Fixing" by R. G. Stansfield,
      J. I.E.E. Vol 94, Part IIA, 1947.

      Stansfield's approach approximates the bearing error
      \f$\theta_i(x,y)-\tilde{\theta}_i\f$ by
      \f$\sin(\theta_i(x,y)-\tilde{\theta}_i)\f$. 

      Fly in the ointment: Stansfield uses angles measured
      counter-clockwise from the East in those equations.  To
      translate to our usage (bearings clockwise from North) simply
      interchange cosine and sine in the final expressions below.  To keep
      consistent with Stansfield's paper, I retain his convention throughout
      this documentation.

      Starting as in the Maximum Likelihood fix,

      \f$
      P(x,y) = K*\exp(-\sum_{i=0}^n (\tilde{\theta_i} - \theta_i(x,y))^2/(2\sigma_i^2))
      \f$

      Making the substitution that \f$\tilde{\theta_i} -
      \theta_i(x,y)\approx\sin(\tilde{\theta_i} - \theta_i(x,y))\f$
      and that \f$\sin(\tilde{\theta_i} - \theta_i(x,y))=q_i/d_i\f$,
      where \f$q_i\f$ is the perpendicular distance from our bearing
      line to the point (x,y) and \f$d_i\f$ is the distance from our
      receiver site to that point, we get:

      
      \f$
      P(q) = K*\exp(-\sum_{i=0}^n (q_i)^2/(2(d_i\sigma_i)^2))
      \f$

      Unfortunately, the d_i all depend on x,y, making this a fairly ugly
      nonlinear problem to solve.  To simplify matters, Stansfield introduced
      a point which he repeatedly refers to as the actual position of the
      transmitter, making the approximation that d_i is approximately the
      distance to that point and using it instead.  But really any point
      O (for Origin) can be used instead in what follows, so long as O isn't
      too far from (x,y).

      Assume that p_i is the perpendicular distance from our bearing line i

      to O and that \f$\Delta x\f$ and \f$\Delta y\f$ are the offsets from
      point O to point Q.  Then 
      \f$
      q_i = p_i+\Delta x \sin(\tilde{\theta}_i)-\Delta y\cos(\tilde{\theta}_i)
      \f$
      by a very simple geometric argument.  Substituting this mess into the
      new cost function and solving the least squares problem, one concludes
      that the values of \f$\Delta x\f$ and \f$\Delta y\f$ that maximize the
      probability are:

      \f$
      \Delta x = \frac{1}{\lambda\mu-\nu^2}\sum_i p_i\frac{\nu\cos(\tilde{\theta}_i)-\mu\sin(\tilde{\theta}_i)}{(d_i\sigma_i)^2}
      \f$
      and
      \f$
      \Delta y = \frac{1}{\lambda\mu-\nu^2}\sum_i p_i\frac{\lambda\cos(\tilde{\theta}_i)-\nu\sin(\tilde{\theta}_i)}{(d_i\sigma_i)^2}
      \f$

      with
      \f$
      \lambda=\sum_i \frac{\sin^2(\tilde{\theta}_i)}{(d_i\sigma_i)^2}
      \f$

      \f$
      \mu=\sum_i \frac{\cos^2(\tilde{\theta}_i)}{(d_i\sigma_i)^2}
      \f$

      and 
      \f$
      \nu = \sum_i \frac{\cos(\tilde{\theta}_i)\sin(\tilde{\theta}_i)}{(d_i\sigma_i)^2}
      \f$

      Finally, it is the case that the numbers \f$\lambda\f$,
      \f$\mu\f$ and \f$\nu\f$ form the elements of the covariance
      matrix for the distribution of \f$q_i\f$ and can be used to form
      a confidence region.

      If one rotates the system of coordinates around the maximum likelihood
      fix location by an angle \f$\phi\f$ such that
      \f$
      tan(2\phi) = -\frac{2\nu}{\lambda-\mu}
      \f$
      so that 

      \f$
      \Delta x = X\cos(\phi) - Y\sin(\phi)
      \f$

      \f$
      \Delta y = X\sin\phi) + Y\cos(\phi)      
      \f$
      
      where the delta quantities are offsets from the fix location,
      then the region defined by:

      \f$
      \frac{X^2}{a^2} + \frac{Y^2}{b^2} = -2 \ln(1-P')
      \f$

      encloses the region in which there is a probability \f$P'\f$ for
      the true location to be.

      \f$a\f$ and \f$b\f$ are given by:

      \f$
      \frac{1}{a^2}=\lambda - \nu\tan(\phi)
      \f$

      and

      \f$
      \frac{1}{b^2}=\mu + \nu\tan(\phi)
      \f$

      Note that Stansfield's paper has an error in it regarding the
      two expressions above.  The error was pointed out and corrected
      in a report "Probabalistic Position-Fixing" by Steve Burnett et
      al.  This paper was a report from the Claremont Graduate School,
      Claremont McKenna College Mathematics Clinic, 1986.  The only place
      I've been able to find it is the URL:

      http://oai.dtic.mil/oai/oai?verb=getRecord&metadataPrefix=html&identifier=ADA190397      

      Remember when reading the code that Stansfield uses a different angle
      convention for bearings than DFLib does.

      And there is another subtlety that is not discussed in
      Stansfield, but is discussed in "Numerical Calculations for
      Passive Geolocation Scenarios" by Don Koks

      http://www.dsto.defence.gov.au/publications/4949/DSTO-RR-0319.pdf

      Since the Stansfield fix is supposed to be a least squares solution to
      the minimization of the cost function (and therefore the maximization of
      the probability density), the presence of the distance from receiver to
      test point in the cost function is a problem that interferes with the
      solution.  So an iterative process is employed.  Practically, this is the
      algorithm:

      0) starting with an initial guess point O, compute the distances from 
         receivers to O, call them \f$d_i\f$ and use them as an approximation to
         the distances to the test point.

      Iterate:

        1) Using the \f$d_i\f$ determined above, solve the least squares problem using
           the expressions above.  This yields an offset vector \f$(\Delta x, \Delta y)\f$ from O to the approximate fix.
         
        2) Compute the distances from receivers to \f$O+(\Delta x, \Delta y)\f$ and
           call them \f$d_i\f$ again.
           
        3) If \f$(\Delta x, \Delta y)\f$ has changed appreciably this iteration,
           return to step 1 and iterate again.  Otherwise, use \f$O+(\Delta x, \Delta y)\f$ as the Stansfield fix.

      The "changed appreciably" bit will take a little black art, I
      think.  I propose that we simply check that the norm of the
      offset vector hasn't changed more than some tolerance since the last 
      iteration.

    */
    void computeStansfieldFix(DFLib::Abstract::Point &SFix,double &a, 
                              double &b, double &phi);
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

    void computeMLFix(DFLib::Abstract::Point &MLFix);

    
    /*!
      \brief compute cost function for point x,y
      
      this returns the cost function for the transmitter being at x,y
      given the DF reports we have.  The probability density uses the
      cost function in the argument of an exponential.  Minimizing the
      cost function will therefore maximize the probability density.
      
      The cost function is the sum
      \f$
      f(x,y) = \sum_{i=0}^n (\tilde{\theta_i} - \theta_i(x,y))^2/(2\sigma_i^2)
      \f$
      
      where \f$\tilde{\theta_i}\f$ is the measured bearing from receiver
      location i and \f$\theta_i(x,y)\f$ is the bearing from receiver
      location i to point (x,y).  Care must be taken to assure that the
      bearing differences are are always kept in the range e
      \f$-\pi<\tilde{\theta_i} - \theta_i(x,y)<=\pi\f$ to avoid
      discontinuities that break the minimization operation.
    */
    
    double computeCostFunction(vector<double> &evaluationPoint);
    void computeCostFunctionAndGradient(vector<double> &evaluationPoint,
                                                double &f,
                                                vector<double> &gradf);
    void computeCostFunctionAndHessian(vector<double> &evaluationPoint,
                                               double &f,
                                               vector<double> &gradf,
                                               vector<vector<double> > &h);

    // Note, unlike size(), this one doesn't count reports that are marked
    // invalid
    int numValidReports() const;

    inline virtual void toggleValidity(int i)
    {
      if (i<theReports.size()&&i>=0)
        theReports[i]->toggleValidity();
    };

    inline bool isValid(int i) const
    {
      if (i<theReports.size() && i>=0)
        return(theReports[i]->isValid());
      return (false);
    };

    inline int size() const {return theReports.size();};

    // This version returns something the caller can never use to change
    // a report from underneath us.
    inline const DFLib::Abstract::Report * getReport(int i) const
    { 
      if (i<theReports.size() && i>=0)
        return (theReports[i]); 
      return (0);
    };

    int getReportIndex(const string &name) const;
    int getReportIndex(const DFLib::Abstract::Report *reportPtr) const;

    inline const vector<double> &getReceiverLocationXY(int i) 
    { return (theReports[i]->getReceiverLocation()); };

    inline virtual void setEvaluationPoint(vector<double> &ep)
    {
      evaluationPoint = ep;
      f_is_valid=false;
      g_is_valid=false;
      h_is_valid=false;
    };

    /// \return function value
    inline virtual double getFunctionValue()
    {
      if (!f_is_valid)
      {
        function_value=computeCostFunction(evaluationPoint);
        f_is_valid=true;
      }
      return function_value;
    };

    /// \param g gradient returned
    /// \return function value
    inline virtual double getFunctionValueAndGradient(vector<double> &g)
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
    };

    /// \param g gradient returned
    /// \param h hessian returned
    /// \return function value
    inline virtual double 
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
    };
  };
}
#endif // DF_REPORT_COLLECTION_HPP
