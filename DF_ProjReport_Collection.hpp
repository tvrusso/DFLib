#ifndef DF_PROJREPORT_COLLECTION_HPP
#define DF_PROJREPORT_COLLECTION_HPP

#include "DF_Report_Collection.hpp"

namespace DFLib
{
  class ProjReportCollection: public ReportCollection
  {
    // This class is just a report collection that refuses to take anything
    // but "DFLib::Proj::Report" objects instead of the abstract ones.

  public:
    CPL_DLL ProjReportCollection();

    virtual CPL_DLL ~ProjReportCollection();


    // Override the base class one...
    virtual CPL_DLL int addReport(DFLib::Proj::Report * aReport);

  };
}
#endif
