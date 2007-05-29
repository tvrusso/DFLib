#ifndef UTIL_MISC_HPP
#define UTIL_MISC_HPP
#include "port.h"

#include <vector>
#include <string>

using namespace std;

namespace DFLib
{
    
  namespace Util
  {
        
    class Exception
    {
    private:
      string ErrorMessage;
    public:
      CPL_DLL Exception(const string &eMsg)
        :ErrorMessage(eMsg)
      {};
      const CPL_DLL string &getEmsg() const
      {
        return ErrorMessage;
      };
    };
  }
}
#endif
