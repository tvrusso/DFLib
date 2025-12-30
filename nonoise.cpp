//-*- mode:C++ ; c-basic-offset: 2 -*-
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
// Filename       : nonoise.cpp
//
// Purpose        :  Another test harness.
//
// Special Notes  : This is essentially testlsDFfix with the addition
//                  of random bearing errors commented out.  I don't
//                  even remember what I used this for.
//-------------------------------------------------------------------------
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#include <ctime>
#endif
#include <cmath>
#include <iostream>
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
#include "DF_XY_Point.hpp"
#include "DF_XY_Report.hpp"
#include "Util_Minimization_Methods.hpp"

PJ *convertPJ;

void convertMercToLatLon(std::vector<double> &merc, double &lon, double &lat)
{
  PJ_COORD data;
  PJ_COORD newCoord;
  data.xy.x = merc[0];
  data.xy.y = merc[1];

  newCoord=proj_trans(convertPJ,PJ_INV,data);

  if (std::isnan(newCoord.lp.lam) || std::isnan(newCoord.lp.phi))
  {
    cerr << "Converting " << merc[0] << ", " << merc[1] <<
      " to lat/lon failed" << std::endl;
    exit(1);
  }

  // Now output lat/lon of fix in degrees+decimal minutes
  lon=newCoord.lp.lam*RAD_TO_DEG;
  lat=newCoord.lp.phi*RAD_TO_DEG;
}

void convertLatLonToMerc(std::vector<double> &merc, double &lon, double &lat)
{
  PJ_COORD data;
  PJ_COORD newCoord;

  data.lp.lam = lon;
  data.lp.phi = lat;

  newCoord=proj_trans(convertPJ,PJ_FWD,data);

  if (std::isnan(newCoord.xy.x) || std::isnan(newCoord.xy.y))
  {
    cerr << "Converting " << merc[0] << ", " << merc[1] <<
      " to lat/lon failed" << std::endl;
    exit(1);
  }

  // Now output lat/lon of fix in degrees+decimal minutes
  merc[0]=newCoord.xy.x;
  merc[1]=newCoord.xy.y;
}

int main(int argc,char **argv)
{

  double lon,lat;
  std::vector<double> transPos(2);
  int i,j;
  char dms_string[128];

  std::vector<double> FixCutAverage,FCA_stddev;
  std::vector<double> LS_fix;
  std::vector<double> NR_fix;
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

  ofstream gnuplotFile("testlsDFfix.gnuplot");
  gnuplotFile << "set angles degrees" << std::endl;
  gnuplotFile << "set parametric" << std::endl;
  gnuplotFile.precision(16); gnuplotFile.width(20);


#ifdef _MSC_VER
  srand(time(NULL));
#else
  srand48(time(NULL));
#endif
  if (argc < 3)
  {
    cerr << "Usage: " << argv[0] << " <trans lon> <trans lat> " << std::endl;
    cerr << " Remember to pipe list of receiver lon/lats into stdin!" << std::endl;
    exit(1);
  }

  std::string latlon_args="+proj=latlong +datum=WGS84";
  std::string mercator_args="+proj=merc +ellps=WGS84 +lat_ts=0";

  char *latlon_argv= new char [latlon_args.size()+1];
  strcpy(latlon_argv,latlon_args.c_str());
  char *mercator_argv= new char [mercator_args.size()+1];
  strcpy(mercator_argv,mercator_args.c_str());

  convertPJ = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
                                     latlon_argv,
                                     mercator_argv,
                                     0);
  delete [] latlon_argv;
  delete [] mercator_argv;

  lon=proj_dmstor(argv[1],NULL);
  lat=proj_dmstor(argv[2],NULL);

  std::cout << "Transmitter location in decimal degrees: Lon: " << lon*RAD_TO_DEG
       << " Lat: " << lat*RAD_TO_DEG << std::endl;

  convertLatLonToMerc(transPos,lon,lat);

  std::cout << " Transmitter " << " at X=" << transPos[0]
       << " Y= " << transPos[1] << std::endl;


  // Now read receiver lon/lats from stdin.  These are in dms format per
  // proj.4 standard, space delimited.
  while (!std::cin.eof())
  {
    double temp_sigma;
    char junk_space;
    DFLib::XY::Report *reportPtr;
    std::vector<double> tempVector(2);

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
    convertLatLonToMerc(tempVector,lon,lat);

    double bearing=atan2((transPos[0]-tempVector[0]),(transPos[1]-tempVector[1]))*RAD_TO_DEG;
    std::cout << " True bearing to transmitter is " << bearing << std::endl;

    //    bearing += rand_gen.getRandom();
    std::cout << "  no randomizing done "  << std::endl;


    reportPtr = new DFLib::XY::Report(tempVector,bearing,temp_sigma);
    rColl.addReport(reportPtr);

  }
  gnuplotFile << "plot [t=0:40000] ";

  std::cout << "Receiver locations in mercator: " << std::endl;
  for (i=0;i<rColl.size();++i)
  {
    const std::vector<double> &receiverLoc =
      rColl.getReport(i)->getReceiverLocation();
    double rb=dynamic_cast<DFLib::XY::Report const *>(rColl.getReport(i))->getBearing();
    double rbr=rColl.getReport(i)->getReportBearingRadians();
    std::cout << " Receiver " << i << " at X=" << receiverLoc[0]
         << " Y= " << receiverLoc[1] << std::endl;
    std::cout << "  bearing from this receiver to transmitter is "
         <<  rb
         << " degrees (" << rbr << " radians)"<< std::endl;

    if (i != 0)
      gnuplotFile << ",";
    gnuplotFile << receiverLoc[0] << "+sin("<<rb<<")*t,"
                << receiverLoc[1] << "+cos("<<rb<<")*t with lines title \"station " << i << "\' ";

  }
  gnuplotFile << std::endl;

  rColl.computeLeastSquaresFix(LS_fix);
  rColl.computeFixCutAverage(FixCutAverage,FCA_stddev);


  gnuplotFile << "replot " << LS_fix[0] << "," << LS_fix[1] << " with points title \"LS Fix\"" << std::endl;
  gnuplotFile << "replot " << transPos[0] << "," << transPos[1] << " with points title \"Actual Location\"" << std::endl;
  gnuplotFile << "replot " << FixCutAverage[0] << "," << FixCutAverage[1] << " with points title \"Fix Cut Average\"" << std::endl;

  for(i = 0; i<rColl.size() ; ++i)
  {
    const std::vector<double> &receiverLoc =
      rColl.getReport(i)->getReceiverLocation();
    gnuplotFile << "replot " << receiverLoc[0] << "," << receiverLoc[1]
                << " with points title \"Station "<< i << "\'" << std::endl;
  }

  std::cout << " Mercator coordinates of LS fix: "
       << "X = " << LS_fix[0] << " Y = " << LS_fix[1] << std::endl;
  std::cout << " Mercator coordinates of Fix Cut Average: "
       << "X = " << FixCutAverage[0] << " Y = " << FixCutAverage[1] << std::endl;
  std::cout << " Fix Cut Average Standard Deviations: "
       << "X = " << FCA_stddev[0] << " Y = " << FCA_stddev[1] << std::endl;

  convertMercToLatLon(LS_fix,lon,lat);

  EW='E';
  NS='N';

  if (lon < 0)
  {
    lon *= -1;
    EW = 'W';
  }

  if (lat < 0)
  {
    lat *= -1;
    NS = 'S';
  }

  std::cout << "  Longitude of LS fix: " << (int) lon << "d"
       << (lon-(int)lon)*60 << "\'" << EW << std::endl;

  std::cout << "  Latitude of LS fix: " << (int) lat << "d"
       << (lat-(int)lat)*60 << "\'" << NS << std::endl;

  for (double minCutAngle=0; minCutAngle < 50; minCutAngle += 5.0)
  {
    rColl.computeFixCutAverage(FixCutAverage,FCA_stddev,minCutAngle);
    convertMercToLatLon(FixCutAverage,lon,lat);

    EW='E';
    NS='N';

    if (lon < 0)
    {
      lon *= -1;
      EW = 'W';
    }

    if (lat < 0)
    {
      lat *= -1;
      NS = 'S';
    }

    std::cout << "  Longitude of Fix Cut Average (min cut angle="<< minCutAngle <<"): " << (int) lon << "d"
         << (lon-(int)lon)*60 << "\'" << EW << std::endl;

    std::cout << "  Latitude of Fix Cut Average (min cut angle="<< minCutAngle <<"): " << (int) lat << "d"
         << (lat-(int)lat)*60 << "\'" << NS << std::endl;
  }
  // Now try Conjugate Gradients on Jml, always starting from OV fix.
  NR_fix=LS_fix;
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
      cerr << " Ooops... got exception " << x.getEmsg() << std::endl;
    }
  }

  std::cout << "Final C-G ML result took " << j  << " iterations: X=" << NR_fix[0] << " Y=" << NR_fix[1] << std::endl;
  gnuplotFile << "replot " << NR_fix[0] << "," << NR_fix[1] << " with points title \"ML Fix\"" << std::endl;
  convertMercToLatLon(NR_fix,lon,lat);

  EW='E';
  NS='N';

  if (lon < 0)
  {
    lon *= -1;
    EW = 'W';
  }

  if (lat < 0)
  {
    lat *= -1;
    NS = 'S';
  }

  std::cout << "  Longitude of ML fix: " << (int) lon << "d"
       << (lon-(int)lon)*60 << "\'" << EW << std::endl;

  std::cout << "  Latitude of ML fix: " << (int) lat << "d"
       << (lat-(int)lat)*60 << "\'" << NS << std::endl;

#if 0
  while (!done && j < 100)
  {
    double determJ;
    double deltax,deltay;
    bool sd_step = false;

    determJ=jac[0][0]*jac[1][1]-jac[0][1]*jac[1][0];
    normf = sqrt(gradf[0]*gradf[0]+gradf[1]*gradf[1]);

    std::cout << " Before step: f = " << f << std::endl;

    if (j!= 0)
    {
      std::cout << "   normf = " << normf
           << "   lastnormf=" << lastnormf
           << "   difference = " << normf-lastnormf << std::endl;
      std::cout << "   Function value = " << f << std::endl;
      std::cout << "   Gradf=("<<gradf[0]<<","<<gradf[1]<<")"<<std::endl;
    }
    if (j!=0 && fabs(determJ) <1e-10)
    {
      std::cout << " Oops --- determinant getting small: " << determJ << std::endl;
    }

    // Newton direction
    // compute deltax, deltay from Jdelta=-f:
    // using analytic inverse of 2x2 matrix:
    deltax = (jac[1][1]*(-gradf[0])-jac[0][1]*(-gradf[1]))/determJ;
    deltay = (-jac[1][0]*(-gradf[0])+jac[0][0]*(-gradf[1]))/determJ;

    lastnormf=normf;
    lastf=f;

    NR_fix[0] += deltax;
    NR_fix[1] += deltay;

    rColl.computeCostFunctionAndHessian(NR_fix,f,gradf,jac);


    j++;
    std::cout << " N-R Iteration " << j << "X " << NR_fix[0] << " Y "
         << NR_fix[1]
         << " dx " << deltax << " dy " << deltay << std::endl;
    std::cout << "   f = " << f << std::endl;

    if ( deltax*deltax + deltay*deltay < 1e-4) // norm less than 1e-8
      done=true;
  }

  std::cout << "Final N-R result: X=" << NR_fix[0] << " Y=" << NR_fix[1] << std::endl;
  gnuplotFile << "replot " << NR_fix[0] << "," << NR_fix[1] << " with points title \"N-R ML Fix\"" << std::endl;
  convertMercToLatLon(NR_fix,lon,lat);

  EW='E';
  NS='N';

  if (lon < 0)
  {
    lon *= -1;
    EW = 'W';
  }

  if (lat < 0)
  {
    lat *= -1;
    NS = 'S';
  }

  std::cout << "  Longitude of NR fix: " << (int) lon << "d"
       << (lon-(int)lon)*60 << "\'" << EW << std::endl;

  std::cout << "  Latitude of NR fix: " << (int) lat << "d"
       << (lat-(int)lat)*60 << "\'" << NS << std::endl;

#endif

  gnuplotFile << "pause -1" << std::endl;
  gnuplotFile.close();
}
