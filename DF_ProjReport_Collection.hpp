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
// Purpose        : A subclass of DFLib::ReportCollection in which we make
//                  the explicit requirement that all reports are 
//                  DFLib::Proj::Report.
//
// Special Notes  : In some of my GUI work on qDF I found that the lax
//                  assumptions by DFLib::ReportCollection (i.e. that
//                  any implementation of DFLib::Abstract::Report could be
//                  included in the report) made for some problems.  Having
//                  a class where we could be *sure* every report in the class
//                  was of the same type was useful.  Thus, this subclass that
//                  merely wraps the other with calls that specify a concrete
//                  class as acceptable arguments rather than an abstract one.
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
#ifndef DF_PROJREPORT_COLLECTION_HPP
#define DF_PROJREPORT_COLLECTION_HPP
#include "DFLib_port.h"

#include "DF_Report_Collection.hpp"

namespace DFLib
{
  class CPL_DLL ProjReportCollection: public ReportCollection
  {
    // This class is just a report collection that refuses to take anything
    // but "DFLib::Proj::Report" objects instead of the abstract ones.

  public:
    ProjReportCollection();

    virtual ~ProjReportCollection();


    // Override the base class one...
    virtual int addReport(DFLib::Proj::Report * aReport);

  };
}
#endif
