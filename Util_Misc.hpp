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
        
    class CPL_DLL Exception
    {
    private:
      string ErrorMessage;
    public:
      Exception(const string &eMsg)
        :ErrorMessage(eMsg)
      {};
      const string &getEmsg() const
      {
        return ErrorMessage;
      };
    };
  }
}
#endif
