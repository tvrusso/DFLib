//-*- mode:C++ ; c-basic-offset: 2 -*-
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
#ifdef _MSC_VER
#define _CRT_RAND_S
#include <ctime>
#endif
#include <cmath>
#include <cstdlib>
#include "Util_Misc.hpp"
#include "gaussian_random.hpp"

namespace DFLib
{
  namespace Util
  {
    gaussian_random_generator::gaussian_random_generator(double m, double s)
      :mean(m),
       std_dev(s),
#ifdef _MSC_VER
       iff(0),
#endif
       use_last(false)
    {
    }

    double gaussian_random_generator::getRandom()
    {
      double x1, x2, w, y1;

      if (use_last)      /* use value from previous call */
      {
        y1 = ysave;
        use_last = false;
      }
      else
      {
        do 
        {
          x1 = 2.0 * uniformRandom() - 1.0;
          x2 = 2.0 * uniformRandom() - 1.0;
          w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );
        
        w = sqrt( (-2.0 * log( w ) ) / w );
        y1 = x1 * w;
        ysave = x2 * w;
        use_last = true;
      }
    
      return( mean + y1 * std_dev );
    }
    /// Returns a random double from a uniform distribution
    double gaussian_random_generator::uniformRandom()
    {
#ifndef _MSC_VER
      return drand48();
#else
      double dum;
      int j;

      if (iff == 0)
      {
        iff=1;
        maxran=RAND_MAX+1.0;
        srand(time(NULL));
        for (j=0;j<98;++j) dum=rand();
        for (j=0;j<98;++j) v[j]=rand();
        y=rand();
      }
      j=98.0*y/maxran;
      if (j>97||j<0)
        throw(Exception("uniformRandom: This cannot happen."));
      y=v[j];
      v[j]=rand();
      return y/maxran;
#endif
    }
  }
}
