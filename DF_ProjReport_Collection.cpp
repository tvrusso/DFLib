//-*- mode:C++ ; c-basic-offset: 2 -*-
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <iostream>
#include <cmath>
#include "DF_Proj_Report.hpp"
#include "DF_ProjReport_Collection.hpp"

using namespace std;
namespace DFLib
{

  ProjReportCollection::ProjReportCollection()
    : ReportCollection()
  {
  }

  ProjReportCollection::~ProjReportCollection()
  {
  }
  

  int ProjReportCollection::addReport(DFLib::Proj::Report * aReport)
  {
    return (DFLib::ReportCollection::addReport(dynamic_cast<DFLib::Abstract::Report *>(aReport)));
  }

}
