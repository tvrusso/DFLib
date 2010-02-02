//-*- mode:C++ ; c-basic-offset: 2 -*-
#ifdef _MSC_VER
#define _CRT_RAND_S
#include <ctime>
#endif
#include <cmath>
#include <cstdlib>
#include "Util_Misc.hpp"
#include "gaussian_random.hpp"
using namespace std;

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
        for (j=0;j<97;++j) dum=rand();
        for (j=0;j<97;++j) v[j]=rand();
        y=rand();
      }
      j=97*y/maxran;
      if (j>=97||j<0)
        throw(Exception("uniformRandom: This cannot happen."));
      y=v[j];
      v[j]=rand();
      return y/maxran;
#endif
    }
  }
}
