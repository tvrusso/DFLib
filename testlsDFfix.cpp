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
// Filename       : $RCSfile$
//
// Purpose        : Test fixture for DFLib methods
//
// Special Notes  : This was the first test harness, and used the trivial 
//                  XY style of reports and points.  All conversion from
//                  lat/lon to XY was done by hand in this early test.  Later
//                  test harnesses used more elaborate implementations of the
//                  Abstract::Point and Abstract::Report interfaces and 
//                  are simpler.
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
#include <fstream>
#include <sstream>
#include <vector>
#include <proj_api.h>

extern "C" {
  double dmstor(const char *, char **);
}

#include "Util_Misc.hpp"
#include "gaussian_random.hpp"
#include "DF_Report_Collection.hpp"
#include "DF_XY_Point.hpp"
#include "DF_XY_Report.hpp"
#include "Util_Minimization_Methods.hpp"

using namespace std;

char *latlon_argv[2]={"proj=latlong",
                      "datum=WGS84"};
char *mercator_argv[3]={"proj=merc",
                        "ellps=WGS84",
                        "lat_ts=0"};
projPJ latlonProj, mercProj;

void convertMercToLatLon(vector<double> &merc, double &lon, double &lat)
{
  projUV data;
  double z;
  data.u = merc[0];
  data.v = merc[1];
  z=0;
  if (pj_transform(mercProj,latlonProj,1,0,&(data.u),&(data.v),&z) != 0)
  {
    cerr << "Converting " << merc[0] << ", " << merc[1] << 
      " to lat/lon failed" << endl;
    exit(1);
  }
  
  // Now output lat/lon of fix in degrees+decimal minutes
  lon=data.u*RAD_TO_DEG;
  lat=data.v*RAD_TO_DEG;
}

void convertLatLonToMerc(vector<double> &merc, double &lon, double &lat)
{
  projUV data;
  double z;
  data.u = lon;
  data.v = lat;
  z=0;
  if (pj_transform(latlonProj,mercProj,1,0,&(data.u),&(data.v),&z) != 0)
  {
    cerr << "Converting " << merc[0] << ", " << merc[1] << 
      " to lat/lon failed" << endl;
    exit(1);
  }
  
  // Now output lat/lon of fix in degrees+decimal minutes
  merc[0]=data.u;
  merc[1]=data.v;
}

int main(int argc,char **argv)
{
  
  double lon,lat;
  vector<double> transPos(2,0.0);
  int i,j;
  char dms_string[128];
  
  vector<double> FCA_stddev;
  vector<double> NR_fix;

  DFLib::XY::Point LS_fix(transPos);
  DFLib::XY::Point FixCutAverage(transPos);

  char EW,NS;
  bool done;
  double normf,lastnormf;
  double lastf;
  
  DFLib::ReportCollection rColl;
  DFLib::Util::Minimizer bogus(&rColl);
  
  // newton-raphson temporaries
  double f;
  vector<double> gradf;
  vector<vector<double> > jac;

  ofstream gnuplotFile("testlsDFfix.gnuplot");
  ofstream gridFile("function.grid");
  gnuplotFile << "set angles degrees" << endl;
  gnuplotFile << "set parametric" << endl;
  gnuplotFile.precision(16); gnuplotFile.width(20);

#ifdef _MSC_VER
  srand(time(NULL));
#else
  srand48(time(NULL));
#endif
  if (argc < 3)
  {
    cerr << "Usage: " << argv[0] << " <trans lon> <trans lat> " << endl;
    cerr << " Remember to pipe list of receiver lon/lats into stdin!" << endl;
    exit(1);
  }

  if (!(latlonProj = pj_init(2,latlon_argv)))
  {
    printf("Using from definition: ");
    for( i = 0; i < 2; i++ )
      printf( "%s ", latlon_argv[i] );
    printf( "\n" );
    
    printf("Projection initialization error\n");
    exit(1);
  }

  if (!(mercProj = pj_init(3,mercator_argv)))
  {
    printf("Using from definition: ");
    for( i = 0; i < 3; i++ )
      printf( "%s ", mercator_argv[i] );
    printf( "\n" );
    
    printf("Projection initialization error\n");
    exit(1);
  }

  lon=dmstor(argv[1],NULL);
  lat=dmstor(argv[2],NULL);

  cout << "Transmitter location in decimal degrees: Lon: " << lon*RAD_TO_DEG
       << " Lat: " << lat*RAD_TO_DEG << endl;

  convertLatLonToMerc(transPos,lon,lat);

  cout << " Transmitter " << " at X=" << transPos[0]
       << " Y= " << transPos[1] << endl;


  // Now read receiver lon/lats from stdin.  These are in dms format per
  // proj.4 standard, space delimited.
  while (!cin.eof())
  {
    double temp_sigma;
    char junk_space;
    DFLib::XY::Report *reportPtr;
    vector<double> tempVector(2);
    ostringstream ostr;

    cin.get(dms_string,sizeof(dms_string),' ');
    if (cin.eof())
      break;
    lon=dmstor(dms_string,NULL);
    // get the space: 
    cin.get(junk_space);
    cin.get(dms_string,sizeof(dms_string),' ');
    if (cin.eof())
      break;
    lat=dmstor(dms_string,NULL);
    cin.get(junk_space);

    cin >>  temp_sigma;

    if (cin.eof())
      break;
    DFLib::Util::gaussian_random_generator rand_gen(0,temp_sigma);

    cout << " Got receiver number " << rColl.size()
         << " Position = " << lon*RAD_TO_DEG << " " << lat*RAD_TO_DEG 
         << " With standard deviation " << temp_sigma
         << endl;
    ostr << "report " << rColl.size();
    convertLatLonToMerc(tempVector,lon,lat);

    double bearing=atan2((transPos[0]-tempVector[0]),(transPos[1]-tempVector[1]))*RAD_TO_DEG;
    cout << " True bearing to transmitter is " << bearing << endl;

    bearing += rand_gen.getRandom();
    cout << " after randomizing, bearing to transmitter is " << bearing << endl;


    reportPtr = new DFLib::XY::Report(tempVector,bearing,temp_sigma,ostr.str());
    rColl.addReport(reportPtr);

  }
  gnuplotFile << "plot [t=0:40000] ";

  cout << "Receiver locations in mercator: " << endl;
  for (i=0;i<rColl.size();++i)
  {
    vector<double> receiverLoc = 
      rColl.getReceiverLocationXY(i);
    double rb=dynamic_cast<DFLib::XY::Report const *>(rColl.getReport(i))->getBearing();
    double rbr=rColl.getReport(i)->getReportBearingRadians();
    cout << " Receiver " << i << " at X=" << receiverLoc[0]
         << " Y= " << receiverLoc[1] << endl;
    cout << "  bearing from this receiver to transmitter is "
         <<  rb
         << " degrees (" << rbr << " radians)"<< endl;

    if (i != 0)
      gnuplotFile << ",";
    gnuplotFile << receiverLoc[0] << "+sin("<<rb<<")*t,"
                << receiverLoc[1] << "+cos("<<rb<<")*t with lines title \"station " << i << "\" ";

  }
  gnuplotFile << endl;

  rColl.computeLeastSquaresFix(LS_fix);
  rColl.computeFixCutAverage(FixCutAverage,FCA_stddev);
  vector <double>LS_point = LS_fix.getXY();
  vector <double>FCA_point = FixCutAverage.getXY();


  gnuplotFile << "replot " << LS_point[0] << "," << LS_point[1] << " with points title \"LS Fix\"" << endl;
  gnuplotFile << "replot " << transPos[0] << "," << transPos[1] << " with points title \"Actual Location\"" << endl;
  gnuplotFile << "replot " << FCA_point[0] << "," << FCA_point[1] << " with points title \"Fix Cut Average\"" << endl;

  for(i = 0; i<rColl.size() ; ++i)
  {
    const vector<double> &receiverLoc = 
      rColl.getReceiverLocationXY(i);
    gnuplotFile << "replot " << receiverLoc[0] << "," << receiverLoc[1] 
                << " with points title \"Station "<< i << "\"" << endl;
  }

  cout << " Mercator coordinates of LS fix: " 
       << "X = " << LS_point[0] << " Y = " << LS_point[1] << endl;
  cout << " Mercator coordinates of Fix Cut Average: " 
       << "X = " << FCA_point[0] << " Y = " << FCA_point[1] << endl;
  cout << " Fix Cut Average Standard Deviations: " 
       << "X = " << FCA_stddev[0] << " Y = " << FCA_stddev[1] << endl;

  convertMercToLatLon(LS_point,lon,lat);

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

  cout << "  Longitude of LS fix: " << (int) lon << "d" 
       << (lon-(int)lon)*60 << "\"" << EW << endl;

  cout << "  Latitude of LS fix: " << (int) lat << "d" 
       << (lat-(int)lat)*60 << "\"" << NS << endl;

  for (double minCutAngle=0; minCutAngle < 50; minCutAngle += 5.0)
  {
    rColl.computeFixCutAverage(FixCutAverage,FCA_stddev,minCutAngle);
    FCA_point = FixCutAverage.getXY();
    convertMercToLatLon(FCA_point,lon,lat);

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

    cout << "  Longitude of Fix Cut Average (min cut angle="<< minCutAngle <<"): " << (int) lon << "d" 
         << (lon-(int)lon)*60 << "\"" << EW << endl;

    cout << "  Latitude of Fix Cut Average (min cut angle="<< minCutAngle <<"): " << (int) lat << "d" 
         << (lat-(int)lat)*60 << "\"" << NS << endl;
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
  //    gridFile << endl;
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
      cout << " CG says minimum at  " << NR_fix[0] << "," << NR_fix[1] 
           << " where the function is " << tempF << endl;
    }
    catch (DFLib::Util::Exception x)
    {
      cerr << " Ooops... got exception " << x.getEmsg() << endl;
    }
  }

  cout << "Final C-G ML result took " << j  << " iterations: X=" << NR_fix[0] << " Y=" << NR_fix[1] << endl;
  gnuplotFile << "replot " << NR_fix[0] << "," << NR_fix[1] << " with points title \"ML Fix\"" << endl;
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

  cout << "  Longitude of ML fix: " << (int) lon << "d" 
       << (lon-(int)lon)*60 << "\"" << EW << endl;

  cout << "  Latitude of ML fix: " << (int) lat << "d" 
       << (lat-(int)lat)*60 << "\"" << NS << endl;

#if 0
  while (!done && j < 100)
  {
    double determJ;
    double deltax,deltay;
    bool sd_step = false;

    determJ=jac[0][0]*jac[1][1]-jac[0][1]*jac[1][0];
    normf = sqrt(gradf[0]*gradf[0]+gradf[1]*gradf[1]);

    cout << " Before step: f = " << f << endl;

    if (j!= 0) 
    {
      cout << "   normf = " << normf
           << "   lastnormf=" << lastnormf
           << "   difference = " << normf-lastnormf << endl;
      cout << "   Function value = " << f << endl;
      cout << "   Gradf=("<<gradf[0]<<","<<gradf[1]<<")"<<endl;
    }
    if (j!=0 && fabs(determJ) <1e-10)
    {
      cout << " Oops --- determinant getting small: " << determJ << endl;
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
    cout << " N-R Iteration " << j << "X " << NR_fix[0] << " Y " 
         << NR_fix[1] 
         << " dx " << deltax << " dy " << deltay << endl;
    cout << "   f = " << f << endl;

    if ( deltax*deltax + deltay*deltay < 1e-4) // norm less than 1e-8
      done=true;
  }

  cout << "Final N-R result: X=" << NR_fix[0] << " Y=" << NR_fix[1] << endl;
  gnuplotFile << "replot " << NR_fix[0] << "," << NR_fix[1] << " with points title \"N-R ML Fix\"" << endl;
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

  cout << "  Longitude of NR fix: " << (int) lon << "d" 
       << (lon-(int)lon)*60 << "\"" << EW << endl;

  cout << "  Latitude of NR fix: " << (int) lat << "d" 
       << (lat-(int)lat)*60 << "\"" << NS << endl;
      
#endif

  gnuplotFile << "pause -1" << endl;
  gnuplotFile.close();
}
