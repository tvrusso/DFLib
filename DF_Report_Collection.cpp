//-*- mode:C++ ; c-basic-offset: 2 -*-
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <iostream>
#include <cmath>
#include <limits>
#include "DF_Abstract_Report.hpp"
#include "DF_Report_Collection.hpp"
#include "Util_Minimization_Methods.hpp"
#include "Util_Misc.hpp"

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
  void ReportCollection::aggressiveComputeMLFix(DFLib::Abstract::Point &MLFix)
  {

    DFLib::Util::Minimizer bogus(this);
    vector<double> NR_fix = MLFix.getXY();
    int j;

    // First do a quickie Nelder-Mead simplex minimize
    vector<vector<double> > Simplex(3);
    Simplex[0]=NR_fix;
    Simplex[1]=NR_fix;
    Simplex[2]=NR_fix;
    
    // get gradient of cost function at base point.
    vector<double> gradient;
    double f;
    computeCostFunctionAndGradient(Simplex[0],f,gradient);
    // normalize:
    f=sqrt(gradient[0]*gradient[0]+gradient[1]*gradient[1]);
    gradient[0]/=f;
    gradient[1]/=f;
    // perturb along the downhill direction
    Simplex[1][0] += -10*gradient[0];
    Simplex[1][1] += -10*gradient[1];

    // perturb along the direction orthogonal to gradient here
    Simplex[2][0] += -10*gradient[1];
    Simplex[2][1] +=  10*gradient[0];

    try 
    {
      int simpIndex=bogus.nelderMeadMinimize(Simplex);
      NR_fix=Simplex[simpIndex];
    }
    catch (DFLib::Util::Exception x)
    {
      cerr << " Caught exception in nelderMeadMinimize:" << endl
           << x.getEmsg() << endl;
    }
    
    double tempF=bogus.conjugateGradientMinimize(NR_fix,1e-5,j);
    MLFix.setXY(NR_fix);
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

  /// \brief Compute Stansfield fix
  void ReportCollection::computeStansfieldFix(DFLib::Abstract::Point &SFix,
                                              double &am2, double &bm2,
                                              double &phi)
  {
    vector<double> initialFix = SFix.getXY();
    vector<double> distances;
    vector<double> sines;
    vector<double> cosines;
    vector<double> sigmas;
    vector<double> p;   // Stansfield's "p_i"
    vector<double> temp(2);
    vector<double> deltas(2);

    double mu, nu, lambda;
    double lastNorm=1e100; 
    double currentNorm=1e100; // a ridiculous value to start with
    int numIters=0;
    double tol=sqrt(numeric_limits<double>::epsilon());

    vector<DFLib::Abstract::Report *>::iterator iterReport;
    vector<DFLib::Abstract::Report *>::iterator reportEnd=theReports.end();
    distances.clear();
    sines.clear();
    cosines.clear();
    distances.reserve(theReports.size());
    cosines.reserve(theReports.size());
    sines.reserve(theReports.size());
    sigmas.reserve(theReports.size());
    int i=0;

    // initialize
    for (iterReport=theReports.begin();
         iterReport!=reportEnd;
         ++iterReport)
    {
      if ((*iterReport)->isValid())
      {
        distances.push_back((*iterReport)->computeDistanceToPoint(initialFix));
        cosines.push_back(cos((*iterReport)->getReportBearingRadians()));
        sines.push_back(sin((*iterReport)->getReportBearingRadians()));
        sigmas.push_back((*iterReport)->getBearingStandardDeviationRadians());
        temp=(*iterReport)->getReceiverLocation();
        // remember difference between Stansfield and DFLib convention
        // This is the perpendicular distance between the bearing line
        // from this receiver to the point initialFix.  Cosine and sine
        // interchanged because our bearing is clockwise from north, not 
        // counterclockwise from east.
        p.push_back( cosines[i]*(initialFix[0]-temp[0])
                     -sines[i]*(initialFix[1]-temp[1]));
        ++i;
      }
    }

    // we only set these nonzero if we converge.    
    am2=bm2=0;

    // we are now ready to iterate.
    do 
    {
      mu=nu=lambda=0;
      lastNorm=currentNorm;
      // compute mu, nu, lambda
      for (int i=0; i< p.size(); ++i)
      {
        double dsigma2=(distances[i]*distances[i]*sigmas[i]*sigmas[i]);
        // again, sine and cosine opposite from Stansfield because of
        // angular convention
        mu += (sines[i]*sines[i])/dsigma2;
        nu += (cosines[i]*sines[i])/dsigma2;
        lambda += (cosines[i]*cosines[i])/dsigma2;
      }
      double denom=lambda*mu-nu*nu;
      deltas[0]=0;
      deltas[1]=0;
      for (int i=0; i< p.size(); ++i)
      {
        double dsigma2=(distances[i]*distances[i]*sigmas[i]*sigmas[i]);
        deltas[0] += p[i]*(nu*sines[i]-mu*cosines[i])/dsigma2;
        deltas[1] += p[i]*(lambda*sines[i]-nu*cosines[i])/dsigma2;
      }
      deltas[0] /= denom;
      deltas[1] /= denom;
      currentNorm=sqrt(deltas[0]*deltas[0]+deltas[1]*deltas[1]);

      temp[0]=initialFix[0]+deltas[0];
      temp[1]=initialFix[1]+deltas[1];
      int tempInt=0;
      // now we have to generate new estimates of distance:
      // initialize
      for (iterReport=theReports.begin();
           iterReport!=reportEnd;
           ++iterReport)
      {
        if ((*iterReport)->isValid())
        {
          distances[tempInt++]=((*iterReport)->computeDistanceToPoint(temp));
        }
      }
      
      ++numIters;
    } while (abs(currentNorm-lastNorm)>tol && numIters<=100);

    // we get here either because we failed to converge or because we did
    // converge.  Check.
    if (numIters > 100)
      throw(Util::Exception("Too many iterations in computeStansfieldFix"));
    else
    {
      // we converged, compute the error ellipse info and save the 
      // fix in the point we were given
      initialFix[0] += deltas[0];
      initialFix[1] += deltas[1];
      SFix.setXY(initialFix);
    
      // tan(2*phi)= -2*nu/(lambda-mu)
      phi=.5*atan2(-2*nu,lambda-mu);
      am2=(lambda-nu*tan(phi));
      bm2=(mu+nu*tan(phi));
    }
  }

  /// \brief Compute Cramer-Rao bounds
  void ReportCollection::computeCramerRaoBounds(DFLib::Abstract::Point &MLFix,
                                                double &am2, double &bm2, 
                                                double &phi)
  {
    am2=0;
    bm2=0;
    vector<double> initialFix = MLFix.getXY();
    vector<double> temp(2);
    double lambda=0;
    double mu=0;
    double nu=0;
    vector<DFLib::Abstract::Report *>::iterator iterReport;
    vector<DFLib::Abstract::Report *>::iterator reportEnd=theReports.end();
    for (iterReport=theReports.begin(); iterReport!=reportEnd; ++iterReport)
    {
      if ((*iterReport)->isValid())
      {
        temp=(*iterReport)->getReceiverLocation();
        double dx=initialFix[0]-temp[0];
        double dy=initialFix[1]-temp[1];
        double sigma=(*iterReport)->getBearingStandardDeviationRadians();
        double ds2=dx*dx+dy*dy;
        double denom=sigma*sigma*ds2*ds2;

        lambda += dy*dy/denom;
        nu += dx*dy/denom;
        mu += dx*dx/denom;
      }
    }

    phi=.5*atan2(-2*nu,lambda-mu);
    am2=(lambda-nu*tan(phi));
    bm2=(mu+nu*tan(phi));
    
  }


  /// \brief Compute Cost Function

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
