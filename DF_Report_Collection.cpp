//-*- mode:C++ ; c-basic-offset: 2 -*-
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <iostream>
#include <cmath>
#include "DF_Abstract_Report.hpp"
#include "DF_Report_Collection.hpp"
#include "Util_Minimization_Methods.hpp"

using namespace std;
namespace DFLib
{
  // Class DFReportCollection

  ReportCollection::ReportCollection()
    :f_is_valid(false),
     g_is_valid(false),
     h_is_valid(false)
  {
    theReports.clear();
  }

  ReportCollection::~ReportCollection()
  {
  }

  void ReportCollection::deleteReports()
  {
    vector<DFLib::Abstract::Report *>::iterator iterReport=theReports.begin();
    vector<DFLib::Abstract::Report *>::iterator lastReport=theReports.end();
    while (iterReport != lastReport)
    {
      delete *iterReport;
      ++iterReport;
    }
    theReports.clear();
  }

  int ReportCollection::addReport(DFLib::Abstract::Report *aReport)
  {
    theReports.push_back(aReport);
    return (theReports.size()-1); // return the index to this report.
  }

  bool ReportCollection::computeFixCutAverage(DFLib::Abstract::Point &FCA,
                                              vector<double> &FCA_stddev,
                                              double minAngle)
  {
    vector<double> tempVec;
    FixStatus fs;
    bool retval;
    int numCuts=0;
    vector<double> tempFCA;
    vector<double> tempScratch;
    DFLib::Abstract::Point *tempPoint;

    // Make a point object that uses the same coordinate system that our
    // return object will use, and initialize it to 0.
    tempPoint = FCA.Clone();
    tempFCA.resize(2);
    tempFCA[0]=tempFCA[1]=0;
    tempPoint->setXY(tempFCA);

    FCA_stddev.resize(2);
    FCA_stddev[0]=FCA_stddev[1]=0;

    vector<DFLib::Abstract::Report *>::iterator iterReportI;
    vector<DFLib::Abstract::Report *>::iterator reportEnd=theReports.end();
    // loop over all reports
    for (iterReportI=theReports.begin();iterReportI!=reportEnd;
         ++iterReportI)
    {
      if ((*iterReportI)->isValid())
      {
        // loop over all reports after this one
        vector<DFLib::Abstract::Report *>::iterator iterReportJ;
        for (iterReportJ=iterReportI,++iterReportJ;
             iterReportJ != reportEnd;
             ++iterReportJ)
        {
          if ((*iterReportJ)->isValid())
          {
            double cutAngle;
            (*iterReportI)->computeFixCut(*iterReportJ,*tempPoint,cutAngle,fs);
            if (fs == GOOD_FIX && fabs(cutAngle) >= minAngle*M_PI/180.0)
            {
              numCuts++;
              tempVec = tempPoint->getUserCoords();
              tempFCA[0] += tempVec[0];
              tempFCA[1] += tempVec[1];
              FCA_stddev[0] += tempVec[0]*tempVec[0];
              FCA_stddev[1] += tempVec[1]*tempVec[1];
            }
          }
        }
      }
    }
    if (numCuts != 0) // we actually got at least one cut
    {
      tempFCA[0] /= numCuts;
      tempFCA[1] /= numCuts;

      // Do not compute fix cut average standard deviation unless there's more
      // than one cut!
      if (numCuts > 1)
      {
        FCA_stddev[0] /= numCuts;
        FCA_stddev[1] /= numCuts;
        // FCA_stddev now has <tempFCA^2>.  Now compute 
        // sqrt((<tempFCA^2>-<tempFCA>^2)), the standard deviation
        FCA_stddev[0] = sqrt((FCA_stddev[0]-tempFCA[0]*tempFCA[0]));
        FCA_stddev[1] = sqrt((FCA_stddev[1]-tempFCA[1]*tempFCA[1]));
      }
      else
      {
        FCA_stddev[0]=FCA_stddev[1]=0.0; 
      }
      retval = true;
    }
    else
    {
      retval = false;
    }

    FCA.setUserCoords(tempFCA);
    delete tempPoint;
    return retval;
  }


  /// \brief compute ML fix
  void ReportCollection::computeMLFix(DFLib::Abstract::Point &MLFix)
  {

    DFLib::Util::Minimizer bogus(this);
    vector<double> NR_fix = MLFix.getXY();
    int j;
    double tempF=bogus.conjugateGradientMinimize(NR_fix,1e-5,j);
    MLFix.setXY(NR_fix);
  }

  /// \brief compute cost function for point x,y
  ///
  /// this returns the cost function for the transmitter being at x,y
  /// given the DF reports we have.  The probability density uses the
  /// cost function in the argument of an exponential.  Minimizing the
  /// cost function will therefore maximize the probability density.

  double ReportCollection::computeCostFunction(vector<double> &evaluationPoint)
  {
    double f=0;
    vector<DFLib::Abstract::Report *>::iterator iterReport;
    vector<DFLib::Abstract::Report *>::iterator reportEnd=theReports.end();

    // Loop over all reports, sum up 
    //    (1/(2*sigma^2)*(measured_bearing-bearing_to_point)^2

    for (iterReport=theReports.begin();
         iterReport!=reportEnd;
         ++iterReport)
    {
      if ((*iterReport)->isValid())
      {
        double bearing = (*iterReport)->getReportBearingRadians();
        double bearing_to_point 
          = (*iterReport)->computeBearingToPoint(evaluationPoint);
        double sigma = (*iterReport)->getBearingStandardDeviationRadians();
        double deltatheta=bearing_to_point - bearing;

        // Make deltatheta in range -pi<deltatheta<=pi
        while (deltatheta <= -M_PI)
          deltatheta += 2*M_PI;
        while (deltatheta > M_PI)
          deltatheta -= 2*M_PI;

        f += 1/(2*sigma*sigma)*(deltatheta)*
          (deltatheta);
      }
    }
    return (f);
  }

  /// \brief compute cost function for point x,y and its gradient
  void 
  ReportCollection::computeCostFunctionAndGradient
  (
   vector<double> &evaluationPoint,
   double &f,
   vector<double> &gradient
   )
  {

    vector<DFLib::Abstract::Report *>::iterator iterReport;
    vector<DFLib::Abstract::Report *>::iterator reportEnd=theReports.end();

    f=0;
    gradient.resize(2);
    gradient[0]=gradient[1]=0;

    // Loop over all reports, sum up 
    //    (1/(2*sigma^2)*(measured_bearing-bearing_to_point)^2

    for (iterReport=theReports.begin();
         iterReport!=reportEnd;
         ++iterReport)
    {
      if ((*iterReport)->isValid())
      {
        double bearing = (*iterReport)->getReportBearingRadians();
        double bearing_to_point 
          = (*iterReport)->computeBearingToPoint(evaluationPoint);
        double sigma = (*iterReport)->getBearingStandardDeviationRadians();
        double deltatheta;
        double xr=
          (*iterReport)->getReceiverLocation()[0]-evaluationPoint[0];
        double yr=
          (*iterReport)->getReceiverLocation()[1]-evaluationPoint[1];
        double d=sqrt(xr*xr+yr*yr);
        double c=cos(bearing_to_point);
        double s=sin(bearing_to_point);
        
        deltatheta=(bearing-bearing_to_point);
        // Make deltatheta in range -pi<deltatheta<=pi
        while (deltatheta <= -M_PI)
          deltatheta += 2*M_PI;
        while (deltatheta > M_PI)
          deltatheta -= 2*M_PI;

        f += 1/(2*sigma*sigma)*(deltatheta*deltatheta);
        gradient[0] += (deltatheta)/(sigma*sigma*d)*(-c);
        gradient[1] += (deltatheta)/(sigma*sigma*d)*( s);
      }
    }
  }

  /// \brief compute cost function for point x,y its gradient, and its hessian.
  void 
  ReportCollection::computeCostFunctionAndHessian
  (
   vector<double> &evaluationPoint, 
   double &f, vector<double> &gradient, vector<vector<double> > &hessian
   )
  {

    vector<DFLib::Abstract::Report *>::iterator iterReport;
    vector<DFLib::Abstract::Report *>::iterator reportEnd=theReports.end();

    f=0;
    gradient.resize(2);
    gradient[0]=gradient[1]=0;
    hessian.resize(2);
    hessian[0].resize(2);
    hessian[1].resize(2);
    hessian[0][0]=hessian[0][1]=hessian[1][0]=hessian[1][1]=0.0;

    // Loop over all reports, sum up 
    //    (1/(2*sigma^2)*(measured_bearing-bearing_to_point)^2

    for (iterReport=theReports.begin();
         iterReport!=reportEnd;
         ++iterReport)
    {
      if ((*iterReport)->isValid())
      {
        double bearing = (*iterReport)->getReportBearingRadians();
        double bearing_to_point 
          = (*iterReport)->computeBearingToPoint(evaluationPoint);
        double sigma = (*iterReport)->getBearingStandardDeviationRadians();
        double deltatheta=(bearing-bearing_to_point);
        double xr=
          (*iterReport)->getReceiverLocation()[0]-evaluationPoint[0];
        double yr=
          (*iterReport)->getReceiverLocation()[1]-evaluationPoint[1];
        double d=sqrt(xr*xr+yr*yr);
        double c=cos(bearing_to_point);
        double s=sin(bearing_to_point);
        double coef = (1/(sigma*sigma*d*d));
        
        deltatheta=(bearing-bearing_to_point);
        // Make deltatheta in range -pi<deltatheta<=pi
        while (deltatheta <= -M_PI)
          deltatheta += 2*M_PI;
        while (deltatheta > M_PI)
          deltatheta -= 2*M_PI;
        
        f += 1/(2*sigma*sigma)*(deltatheta*deltatheta);
        gradient[0] += (deltatheta)/(sigma*sigma*d)*(-c);
        gradient[1] += (deltatheta)/(sigma*sigma*d)*( s);
        
        hessian[0][0] += coef*(c*c-s*c*deltatheta);
        hessian[0][1] += coef*(-s*c-s*s*deltatheta);
        hessian[1][0] += coef*(-s*c+c*c*deltatheta);
        hessian[1][1] += coef*(s*s-c*s*deltatheta);
      }
    }
  }

  /// \brief compute least squares solution from all df reports.
  void ReportCollection::computeLeastSquaresFix(DFLib::Abstract::Point &LS_Fix)
  {

    double atb1,atb2,a11,a12,a22;
    atb1=atb2=a11=a12=a22=0.0;
    double det;
    vector <double> LS_point;
    
    vector<DFLib::Abstract::Report *>::iterator iterReport;
    vector<DFLib::Abstract::Report *>::iterator reportEnd=theReports.end();
    
    LS_point.resize(2);
    
    for (iterReport=theReports.begin();iterReport!=reportEnd;++iterReport)
    {
      if ((*iterReport)->isValid())
      {
        double bearing=(*iterReport)->getReportBearingRadians();
        double c=cos(bearing);
        double s=sin(bearing);
        double rx=(*iterReport)->getReceiverLocation()[0];
        double ry=(*iterReport)->getReceiverLocation()[1];
        double b=rx*c-ry*s;
	
        atb1 += c*b;
        atb2 += -s*b;
        a11 += s*s;
        a12 += s*c;
        a22 += c*c;
      }
    }
    
    det = a11*a22-a12*a12;
    
    LS_point[0]=(a11*atb1+a12*atb2)/det;
    LS_point[1]=(a12*atb1+a22*atb2)/det;

    LS_Fix.setXY(LS_point);
    LS_point = LS_Fix.getXY();
  }

  int ReportCollection::numValidReports() const
  {
    int numVal=0;
    for (int i=0; i < theReports.size(); i++)
    {
      if (theReports[i]->isValid())
        numVal++;
    }
    return (numVal);
  }

  /// \brief return the index of the report that has the given name, or -1
  ///
  /// Stupid linear search, but should be OK for a realistic size of collection.
  int ReportCollection::getReportIndex(const string & name) const
  {
    int reportIndex=-1;
    for (int i=0; i< theReports.size(); i++) 
    {
      if (theReports[i]->getReportName() == name)
        reportIndex=i;
    }
    return reportIndex;
  }

  /// \brief return the index of the report that has the given pointer, or -1
  ///
  /// Stupid linear search, but should be OK for a realistic size of collection.
  int ReportCollection::getReportIndex(const DFLib::Abstract::Report *reportPtr) const
  {
    int reportIndex=-1;
    for (int i=0; i< theReports.size(); i++) 
    {
      if (theReports[i] == reportPtr)
        reportIndex=i;
    }
    return reportIndex;
  }

}
