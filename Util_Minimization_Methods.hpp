#ifndef UTIL_MINIMIZATION_METHODS_HPP
#define UTIL_MINIMIZATION_METHODS_HPP
#include "DFLib_port.h"

#include <vector>
#include <string>

using namespace std;

namespace DFLib
{

  namespace Abstract
  {
    class Group;
  }

  namespace Util
  {
        
    class CPL_DLL Minimizer
    {
    private:
      DFLib::Abstract::Group *theGroup;
            
    public:
      inline  Minimizer(DFLib::Abstract::Group *aGroup)
        :theGroup(aGroup)
      { };
            
      /// \brief bracket minimum of function
      /// 
      /// Given scalars \f$a,b,c\f$ and an initial point \f$X0\f$
      /// and direction, find new values of \f$a, b,\f$ and \f$c\f$ 
      /// such that the minimum value of theGroup's
      ///  function in the chosen direction is bracketed 
      /// (i.e. \f$f(x0+direction*b)<f(x0+direction*a) \f$
      ///  and \f$f(x0+direction*b) < f(x0+direction*c))\f$
      /// Nearly verbatim port to C++ from Numerical Recipes in C
       void bracketMinimum(double &a, double &b, double &c, 
                                  vector<double> &X0,
                                  vector<double> &direction);
            
      /// \brief minimize by Brent's method with derivatives
      /// given a,b,c scalars that bracket the minimum in given direction
      /// from X0, search for minimum using Brent's method
      /// Returns function value, sets "xmin" to abscissa at minimum.
      /// Almost verbatim port out of Numerical Recipes in C
       double brentMinimize(double a, double b, double c,
                                   vector<double> &X0,vector<double>&direction,
                                   double &xmin);
            
      /// \brief minimize function in given direction
      /// perform a line search along the given direction for the minimum
      /// value of the function.
      /// \param X0 input: starting point  output: minimum point
      /// \param dir input: direction to search output: actual displacement
      /// \return function value at minimum.
      /// Nearly verbatim port out of Numerical Recipes in C
       double lineSearch(vector<double> &X0, vector<double> &dir); 
            
      /// \brief minimize function of vector value by method of conjugate gradients
      ///
      /// \param X0 on input starting point, on exit solution.
      /// \param ftol convergence tolerance on function.
      /// \param iter returned number of iterations taken
      /// \return value of function at minimum.
       double conjugateGradientMinimize(vector<double> &X0, double ftol,
                                               int &iter);
            

       /// \brief minimuze function of vector argument by Nelder-Mead 
       ///        simplex method.
       ///
       /// \param Simplex vector of vector of points representing simplex vertices
       /// \return index into modified simplex of vertex with lowest function 
       ///    value
       ///
       /// The Nelder-Mead simplex method is a minimization method in which
       /// a "simplex" is morphed through the landscape of the function and
       /// tends downhill.  In 2-D the simplex is a triangle, in 3-D a 
       /// tetrahedron, and in more dimensions it is a generalization of 
       /// such a shape.  In  \f$N\f$ dimensions, the simplex has \f$N+1\f$
       /// vertices.  The method keeps track of the best, worst, and "second
       /// worst" fucntion values at the vertices and at each step attempts to
       /// modify the simplex to make the collection better:
       ///
       /// Centroid point:  the centroid of the side opposite the worst point.
       ///
       /// Reflected point: The reflection of the worst point through the 
       /// centroid.
       ///
       /// Expanded point:  a point extrapolated along the line from the
       /// centroid to the reflected point at twice the distance.
       ///
       /// Contracted point:  a point halfway between the centroid and the 
       /// worst point.
       ///
       /// Reduced simplex: the simplex shrunk by half along all
       /// sides, leaving the best point unmodified.
       /// 
       /// The method improves the simplex by the following algorithm:
       ///
       /// if the reflected point is the best so far, try the expanded point
       /// If the expanded point is the best so far, replace the worst point with it.  Otherwise replace the worst point with the reflected point.
       ///
       /// If the reflected point is not the best so far, but is better than
       /// the second-worst, replace the worst point with it.
       ///
       /// If we haven't replaed the worst point yet, check the
       /// contracted point.  If it is better than the second-worst point,
       /// replace the worst point with it.
       ///
       /// If nothing has worked so far, reduce the size of the simplex by
       /// shrinking along all sides towards the current best point.
       ///
       /// The method stops when we have either exceeded the maximum number
       /// of function evaluations (here hard-coded to 5000) or the spread
       /// of values between best and worst has been reduced to a sufficiently
       /// small fraction of the average of best and worst.  
       /// The stopping criterion is (abs(diff)/average)<sqrt(machine_epsilon)

       int nelderMeadMinimize(vector<vector<double> > &Simplex);

      /// \brief Evaluate \f$F(x0+x*dir)\f$ where x0 and dir are vectors
      ///
      /// F is the group's function "computeFunctionValue"
      /// \return function value at \f$X0+x*dir\f$
       double simpleF(double &x,vector<double>&X0,vector<double>&dir);
      /// Evaluate F(x0+x*dir) and its directional derivative where x0 and 
      /// dir are vectors.  The directional derivative is the inner product
      /// of the gradient and the direction: \f$df = dir\cdot\nabla f\f$.
       double simpleFandDeriv(double &x,vector<double>&X0,vector<double>&dir,
                                     double &df);
            
            
    };
  }
}
#endif
