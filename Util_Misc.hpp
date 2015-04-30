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
// Purpose        : Provide a standard class of object to be thrown when
//                  DFLib methods need to throw an exception.
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
#ifndef UTIL_MISC_HPP
#define UTIL_MISC_HPP
#include "DFLib_port.h"

#include <string>


namespace DFLib
{
    
  namespace Util
  {
        
    class CPL_DLL Exception
    {
    private:
      std::string ErrorMessage;
    public:
      Exception(const std::string &eMsg)
        :ErrorMessage(eMsg)
      {};
      const std::string &getEmsg() const
      {
        return ErrorMessage;
      };
    };
  }
}
#endif
