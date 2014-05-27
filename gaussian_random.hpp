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
// Purpose        : This class implements a random number generator that 
//                  returns a random number from a gaussian distribution
//                  with specified mean and standard deviation
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
#ifndef GAUSSIAN_RANDOM_HPP
#define GAUSSIAN_RANDOM_HPP
#include "DFLib_port.h"
#include <cmath>
#include <cstdlib>

namespace DFLib
{
  namespace Util
  {
    /// Provide a random number generator returning values from a 
    /// normal distribution of specified mean and standard deviation.
    class CPL_DLL gaussian_random_generator
    {
    private:
      double ysave;
      bool use_last;
      bool seeded;
      double mean;
      double std_dev;
      /// return double precision uniform deviates.
      double uniformRandom();
#ifdef _MSC_VER
      double y,maxran,v[98];
      int iff;
#endif
    public:
      /// Constructor with specified mean and standard deviation
       gaussian_random_generator(double mean, double std_dev);
      /// Default constructor, mean 0 and standard deviation 1
       gaussian_random_generator():mean(0),std_dev(1),use_last(false) {};
      /// Get normally distributed random deviate.
      /// \return random value from distribution
       double getRandom();
    };
  }
}
#endif
