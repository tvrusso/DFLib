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
