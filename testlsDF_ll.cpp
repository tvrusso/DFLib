//-*- mode:C++ ; c-basic-offset: 2 -*-
//    DFLib: A library of Bearings Only Target Localization algorithms
//    Copyright (C) 2009-2015  Thomas V. Russo
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
// Purpose        : Test harness using the DFLib "lat/lon" classes.
//
// Special Notes  : This test harness was the second in the series, and uses
//                  the DFLib::LatLon::Point and Report classes.  The LatLon
//                  class was later generalized to the "Proj" class, which is
//                  far more flexible.
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
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#include <ctime>
#endif
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <proj.h>

extern "C" {
  double proj_dmstor(const char *, char **);
}

#include "DFLib_Misc_Defs.h"
#include "Util_Misc.hpp"
#include "gaussian_random.hpp"
#include "DF_Report_Collection.hpp"
#include "DF_LatLon_Point.hpp"
#include "DF_LatLon_Report.hpp"
#include "Util_Minimization_Methods.hpp"


int main(int argc,char **argv)
{
  
  double lon,lat;
  std::vector<double> transPos(2,0.0);
  int i,j;
  char dms_string[128];
  
  std::vector<double> FCA_stddev;
  std::vector<double> NR_fix;
  DFLib::LatLon::Point LS_fix(transPos);
  DFLib::LatLon::Point FixCutAverage(transPos);
  char EW,NS;
  bool done;
  double normf,lastnormf;
  double lastf;
  
  DFLib::ReportCollection rColl;
  DFLib::Util::Minimizer bogus(&rColl);
  
  // newton-raphson temporaries
  double f;
  std::vector<double> gradf;
  std::vector<std::vector<double> > jac;

  std::ofstream gnuplotFile("testlsDFfix.gnuplot");
  std::ofstream gridFile("function.grid");
  gnuplotFile << "set angles degrees" << std::endl;
  gnuplotFile << "set parametric" << std::endl;
  gnuplotFile.precision(16); gnuplotFile.width(20);

  std::string progName(argv[0]);
  argv++;
  argc--;

  time_t seed;
  if (argc > 0)
  {
    std::string testArg(argv[0]);

    if (testArg == "--seed")
    {
      argv++;
      argc--;
      seed=atoi(argv[0]);
      argv++;
      argc--;
      std::cerr << " using seed " << seed << std::endl;
    }
    else
    {
      seed=time(NULL);
      std::cerr << " using time " << seed << " as random number seed." << std::endl;
    }
  }
  else
  {
    seed=time(NULL);
    std::cerr << " using time " << seed << " as random number seed." << std::endl;
  }
#ifdef _MSC_VER
  srand(seed);
#else
  srand48(seed);
#endif

  if (argc < 2)
  {
    std::cerr << "Usage: " << progName << " <trans lon> <trans lat> " << std::endl;
    std::cerr << " Remember to pipe list of receiver lon/lats into stdin!" << std::endl;
    exit(1);
  }

  lon=proj_dmstor(argv[0],NULL);
  lat=proj_dmstor(argv[1],NULL);

  std::cout << "Transmitter location in decimal degrees: Lon: " << lon*RAD_TO_DEG
       << " Lat: " << lat*RAD_TO_DEG << std::endl;

  transPos[0]=lon*RAD_TO_DEG;
  transPos[1]=lat*RAD_TO_DEG;
  DFLib::LatLon::Point transPoint(transPos);
  transPos=transPoint.getXY();

  std::cout << " Transmitter " << " at X=" << transPos[0]
       << " Y= " << transPos[1] << std::endl;


  // Now read receiver lon/lats from stdin.  These are in dms format per
  // proj.4 standard, space delimited.
  while (!std::cin.eof())
  {
    double temp_sigma;
    char junk_space;
    DFLib::LatLon::Report *reportPtr;
    std::vector<double> tempVector(2);
    std::string reportName;
    std::ostringstream ost;

    std::cin.get(dms_string,sizeof(dms_string),' ');
    if (std::cin.eof())
      break;
    lon=proj_dmstor(dms_string,NULL);
    // get the space: 
    std::cin.get(junk_space);
    std::cin.get(dms_string,sizeof(dms_string),' ');
    if (std::cin.eof())
      break;
    lat=proj_dmstor(dms_string,NULL);
    std::cin.get(junk_space);

    std::cin >>  temp_sigma;

    if (std::cin.eof())
      break;
    DFLib::Util::gaussian_random_generator rand_gen(0,temp_sigma);

    std::cout << " Got receiver number " << rColl.size()
         << " Position = " << lon*RAD_TO_DEG << " " << lat*RAD_TO_DEG 
         << " With standard deviation " << temp_sigma
         << std::endl;
    tempVector[0]=lon*RAD_TO_DEG;
    tempVector[1]=lat*RAD_TO_DEG;
    double bearing=0; // temporary

    ost << "report " << rColl.size();
    reportPtr = new DFLib::LatLon::Report(tempVector,bearing,temp_sigma,
                                          ost.str());



    bearing=reportPtr->computeBearingToPoint(transPos)*RAD_TO_DEG;
    std::cout << " True bearing to transmitter is " << bearing << std::endl;

    bearing += rand_gen.getRandom();
    std::cout << " after randomizing, bearing to transmitter is " << bearing << std::endl;
    reportPtr->setBearing(bearing);

    rColl.addReport(reportPtr);


  }
  gnuplotFile << "plot [t=0:40000] ";

  std::cout << "Receiver locations in mercator: " << std::endl;
  for (i=0;i<rColl.size();++i)
  {
    std::vector<double> receiverLoc = 
      rColl.getReceiverLocationXY(i);
    double rb=dynamic_cast<DFLib::LatLon::Report const *>(rColl.getReport(i))->getBearing();
    double rbr=rColl.getReport(i)->getReportBearingRadians();
    std::cout << " Receiver " << i << " at X=" << receiverLoc[0]
         << " Y= " << receiverLoc[1] << std::endl;
    std::cout << "  bearing from this receiver to transmitter is "
         <<  rb
         << " degrees (" << rbr << " radians)"<< std::endl;

    if (i != 0)
      gnuplotFile << ",";
    gnuplotFile << receiverLoc[0] << "+sin("<<rb<<")*t,"
                << receiverLoc[1] << "+cos("<<rb<<")*t with lines title \"station " << i << "\" ";

  }
  gnuplotFile << std::endl;

  rColl.computeLeastSquaresFix(LS_fix);
  rColl.computeFixCutAverage(FixCutAverage,FCA_stddev);
  std::vector <double> FCA_point_ll=FixCutAverage.getLL();
  std::cerr<< "Fix Cut Average at lon " << FCA_point_ll[0] << " lat " << FCA_point_ll[1] << std::endl;

  std::vector <double> LS_point=LS_fix.getXY();
  std::vector <double> FCA_point=FixCutAverage.getXY();

  gnuplotFile << "replot " << LS_point[0] << "," << LS_point[1] << " with points title \"LS Fix\"" << std::endl;
  gnuplotFile << "replot " << transPos[0] << "," << transPos[1] << " with points title \"Actual Location\"" << std::endl;
  gnuplotFile << "replot " << FCA_point[0] << "," << FCA_point[1] << " with points title \"Fix Cut Average\"" << std::endl;

  for(i = 0; i<rColl.size() ; ++i)
  {
    const std::vector<double> &receiverLoc = 
      rColl.getReceiverLocationXY(i);
    gnuplotFile << "replot " << receiverLoc[0] << "," << receiverLoc[1] 
                << " with points title \"Station "<< i << "\"" << std::endl;
  }

  std::cout << " Mercator coordinates of LS fix: " 
       << "X = " << LS_point[0] << " Y = " << LS_point[1] << std::endl;
  std::cout << " Mercator coordinates of Fix Cut Average: " 
       << "X = " << FCA_point[0] << " Y = " << FCA_point[1] << std::endl;
  
  std::vector<double> latlon=LS_fix.getUserCoords();

  EW='E';
  NS='N';

  if (lon < 0)
  {
    latlon[0] *= -1;
    EW = 'W';
  }

  if (lat < 0)
  {
    latlon[1] *= -1;
    NS = 'S';
  }

  std::cout << "  Longitude of LS fix: " << (int) latlon[0] << "d" 
       << (latlon[0]-(int)latlon[0])*60 << "\"" << EW << std::endl;

  std::cout << "  Latitude of LS fix: " << (int) latlon[1] << "d" 
       << (latlon[1]-(int)latlon[1])*60 << "\"" << NS << std::endl;

  //  for (double minCutAngle=0; minCutAngle < 50; minCutAngle += 5.0)
  for (double minCutAngle=0; minCutAngle < 5; minCutAngle += 5.0)
  {
    if (rColl.computeFixCutAverage(FixCutAverage,FCA_stddev,minCutAngle))
    {
      std::vector<double> latlon=FixCutAverage.getUserCoords();
	
	
      EW='E';
      NS='N';
	
      if (latlon[0] < 0)
      {
        latlon[0] *= -1;
        EW = 'W';
      }
	
      if (latlon[1] < 0)
      {
        latlon[1] *= -1;
        NS = 'S';
      }
	
      std::cout << "  Longitude of Fix Cut Average (min cut angle="<< minCutAngle <<"): " << (int) latlon[0] << "d" 
           << (latlon[0]-(int)latlon[0])*60 << "\"" << EW << std::endl;
	
      std::cout << "  Latitude of Fix Cut Average (min cut angle="<< minCutAngle <<"): " << (int) latlon[1] << "d" 
           << (latlon[1]-(int)latlon[1])*60 << "\"" << NS << std::endl;
	
      std::cout << "   Std Dev of FCA = (" << FCA_stddev[0] << " , " 
           << " , " << FCA_stddev[1] << ")" << std::endl;
    }
  }

  //
  //  NR_fix.resize(2);
  //  // write out a grid of 10 meter "pixels" showing function values
  //  for (i=-100;i<=100;++i)
  //  {
  //    for (j=-100;j<=100;++j)
  //    {
  //      NR_fix[0] = LS_fix[0]+10.0*i;
  //      NR_fix[1] = LS_fix[1]+10.0*j;
  //      gridFile << rColl.computeCostFunction(NR_fix) << " ";
  //    }
  //    gridFile << std::endl;
  //  }

  // Now try Conjugate Gradients on Jml, always starting from OV fix.
  NR_fix=LS_fix.getXY();
  done = false;
  j=0;

  //  rColl.computeCostFunctionAndHessian(NR_fix,f,gradf,jac);

  {
    double tempF;

    try
    {
      tempF=bogus.conjugateGradientMinimize(NR_fix,1e-5,j);
      std::cout << " CG says minimum at  " << NR_fix[0] << "," << NR_fix[1] 
           << " where the function is " << tempF << std::endl;
    }
    catch (DFLib::Util::Exception x)
    {
      std::cerr << " Ooops... got exception " << x.getEmsg() << std::endl;
    }
  }

  std::cout << "Final C-G ML result took " << j  << " iterations: X=" << NR_fix[0] << " Y=" << NR_fix[1] << std::endl;
  gnuplotFile << "replot " << NR_fix[0] << "," << NR_fix[1] << " with points title \"ML Fix\"" << std::endl;
  DFLib::LatLon::Point NRPoint(NR_fix); // bogus point for now
  NRPoint.setXY(NR_fix);

  latlon=NRPoint.getUserCoords();

  EW='E';
  NS='N';

  if (latlon[0] < 0)
  {
    latlon[0] *= -1;
    EW = 'W';
  }

  if (latlon[1] < 0)
  {
    latlon[1] *= -1;
    NS = 'S';
  }

  std::cout << "  Longitude of ML fix: " << (int) latlon[0] << "d" 
       << (latlon[0]-(int)latlon[0])*60 << "\"" << EW << std::endl;

  std::cout << "  Latitude of ML fix: " << (int) latlon[1] << "d" 
       << (latlon[1]-(int)latlon[1])*60 << "\"" << NS << std::endl;

  gnuplotFile << "pause -1" << std::endl;
  gnuplotFile.close();
}
