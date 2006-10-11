#ifndef UTIL_ABSTRACT_GROUP_HPP
#define UTIL_ABSTRACT_GROUP_HPP
#include "port.h"

#include <vector>
using namespace std;

namespace DFLib
{
    namespace Abstract
    {
        /// \class Util::Abstract::Group
        /// base class for objects to be used with minimization methods
        /// Provides an interface for returning function and derivative values
        /// at a specified point.
        class Group
        {
        public:
            /// Set the point at which to evaluate the function
            virtual CPL_DLL void setEvaluationPoint(vector<double> &ep) = 0;
            /// evaluate the function at the evaluation point
            /// \return value of function
            virtual CPL_DLL double getFunctionValue() = 0;
            /// evaluate the function and its gradient  at the evaluation point
            /// \param gradient STL vector containing gradient on return
            /// \return value of function
            virtual CPL_DLL double 
              getFunctionValueAndGradient(vector<double> &gradient) = 0;
            /// evaluate the function and its gradient  at the evaluation point
            /// \param gradient STL vector containing gradient on return
            /// \param hessian contains hessian on return
            /// \return value of function
            virtual CPL_DLL double
              getFunctionValueAndHessian(vector<double> &gradient,
                                         vector<vector<double> > &hessian)=0;
        };
    }
}
#endif
