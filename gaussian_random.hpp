#ifndef GAUSSIAN_RANDOM_HPP
#define GAUSSIAN_RANDOM_HPP
#include "port.h"
#include <cmath>
#include <cstdlib>

namespace DFLib
{
    namespace Util
    {
        /// Provide a random number generator returning values from a 
        /// normal distribution of specified mean and standard deviation.
        class gaussian_random_generator
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
            double y,maxran,v[97];
            int iff;
#endif
        public:
            /// Constructor with specified mean and standard deviation
            CPL_DLL gaussian_random_generator(double mean, double std_dev);
            /// Default constructor, mean 0 and standard deviation 1
            CPL_DLL gaussian_random_generator():mean(0),std_dev(1),use_last(false) {};
            /// Get normally distributed random deviate.
            /// \return random value from distribution
            CPL_DLL double getRandom();
        };
    }
}
#endif
