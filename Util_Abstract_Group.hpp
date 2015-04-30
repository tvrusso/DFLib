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
// Purpose        : Abstract interface class for use with minimization methods
//                  of DFLib::Util::Minimizer
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
#ifndef UTIL_ABSTRACT_GROUP_HPP
#define UTIL_ABSTRACT_GROUP_HPP
#include "DFLib_port.h"

#include <vector>

namespace DFLib
{
  namespace Abstract
  {
    /// \class Util::Abstract::Group
    /// base class for objects to be used with minimization methods
    /// Provides an interface for returning function and derivative values
    /// at a specified point.
    class CPL_DLL Group
    {
    public:
      /// Set the point at which to evaluate the function
      virtual  void setEvaluationPoint(std::vector<double> &ep) = 0;
      /// evaluate the function at the evaluation point
      /// \return value of function
      virtual  double getFunctionValue() = 0;
      /// evaluate the function and its gradient  at the evaluation point
      /// \param gradient STL vector containing gradient on return
      /// \return value of function
      virtual  double 
      getFunctionValueAndGradient(std::vector<double> &gradient) = 0;
      /// evaluate the function and its gradient  at the evaluation point
      /// \param gradient STL vector containing gradient on return
      /// \param hessian contains hessian on return
      /// \return value of function
      virtual  double
      getFunctionValueAndHessian(std::vector<double> &gradient,
                                 std::vector<std::vector<double> > &hessian)=0;
    };
  }
}
#endif
