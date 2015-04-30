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
// Purpose        : Test harness for DFLib methods using the "Proj" classes.
//
// Special Notes  : This is the last test harness I wrote.  It uses the
//                  DFLib::Proj implementations of DFLib::Abstract::Point
//                  and DFLib::Abstract::Report.
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
// we need this for dmstor and various conversion factors
#include <proj_api.h>
extern "C" {
  double dmstor(const char *, char **);
}

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef WIN32
#include <cfloat>
inline bool isnan(double v) {return _isnan(v)!=0;}
inline bool isinf(double v) {return !_finite(v);}
#endif

#include "Util_Misc.hpp"
#include "gaussian_random.hpp"
#include "DF_Proj_Point.hpp"
#include "DF_Report_Collection.hpp"
#include "DF_Proj_Report.hpp"

#if 0
#ifdef HAVE_GDAL_H
#include <gdal_priv.h>
#endif
#endif


int main(int argc,char **argv)
{
  /*!
    \brief A simple demo program to test out DFLib fixes using the DFLIB::Proj::Report class

    testlsDFfix_proj takes a longitude/latitude pair in PROJ.4 format
    as command line arguments.  This coordinate pair represents the
    actual location of a transmitter.  It also takes as standard in a
    text file of receiver locations with their standard deviations.
    Several example receivers files are distributed with DFLib.

    The program then computes the actual bearings from the given
    receiver locations to the known transmitter location, and adds
    random errors to the bearings.  The random errors are selected
    from a Gaussian distribution with the standard deviations given in
    the receivers file.

    From the randomized bearings, the code computes DF fixes using all
    the methods available in DFLib: the Fix Cut Average, Least Squares
    (orthogonal vectors), Maximum Likelihood and Stansfield fixes.  It 
    outputs a file "testlsDFfix.gnuplot" of Gnuplot commands to plot the
    DF problem and the various fixes.  Error ellipses for the Stansfield
    solution are also plotted.

    The program also outputs a file "testlsDFfix.grasspoints" with the
    station locations and fix locations in a format that can be read by
    GRASS GIS's v.in.ascii. 

    Code to dump out a GeoTIFF file of the Maximum Likelihood method's 
    cost function is commented out with preprocessor "if 0" statements, and
    requires GDAL.  Only the Autoconf build method is currently set up to
    detect and use GDAL.  This stuff is commented out because in at least one
    version of GDAL it lead to reports of memory access issues from Valgrind.
    But it can be useful to uncomment it when trying to understand why some
    geometries of DF problem lead to the Maximum Likelihood and Stansfield
    fixes failing to converge (one finds that the surface is very flat near
    the minimum, and the various minimization methods can't do anything but
    shoot off to infinity).

  */
  double lon,lat;
  std::vector<double> transPos(2,0.0);
  int i,j;
  char dms_string[128];
  
  std::vector<double> FCA_stddev;
  std::vector<double> NR_fix;
  std::vector<std::string> projArgs;
  projArgs.push_back("proj=latlong");
  projArgs.push_back("datum=WGS84");
  DFLib::Proj::Point LS_fix(transPos,projArgs);
  DFLib::Proj::Point FixCutAverage(transPos,projArgs);
  char EW,NS;
  bool done;
  double normf,lastnormf;
  double lastf;
  
  DFLib::ReportCollection rColl;
  
  // newton-raphson temporaries
  double f;
  std::vector<double> gradf;
  std::vector<std::vector<double> > jac;

  std::ofstream gnuplotFile("testlsDFfix.gnuplot");
  std::ofstream gridFile("function.grid");
  std::ofstream pointsFile("testlsDFfix.grasspoints");
  gnuplotFile << "set angles degrees" << std::endl;
  gnuplotFile << "set size square" << std::endl;
  gnuplotFile << "set parametric" << std::endl;
  gnuplotFile.precision(16); gnuplotFile.width(20);
  pointsFile.precision(16); pointsFile.width(20);
  std::cout.precision(16); std::cout.width(20);

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

#if 0
#ifdef HAVE_LIBGDAL  
  testArg=argv[0];
  std::string geotiffName;
  bool geotiffRequested=false;

  if (testArg == "--geotiff")
  {
    std::cout<< "You asked for a geotiff ... not yet." << std::endl;
    argv++;
    argc--;
    geotiffRequested=true;
    geotiffName=argv[0];
    std::cout << "Geotiff name: " << geotiffName << std::endl;
    argv++;
    argc--;
  }
#endif
#endif

  if (argc < 2)
  {
    std::cerr << "Usage: " << progName;
#if 0
#ifdef HAVE_LIBGDAL    
    std::cerr << " [--geotiff <geotiffname>]";
#endif
#endif
    std::cerr << " <trans lon> <trans lat> " << std::endl;
    std::cerr << " Remember to pipe list of receiver lon/lats into stdin!" << std::endl;
    exit(1);
  }

  lon=dmstor(argv[0],NULL);
  lat=dmstor(argv[1],NULL);

  std::cout << "Transmitter location in decimal degrees: Lon: " << lon*RAD_TO_DEG
       << " Lat: " << lat*RAD_TO_DEG << std::endl;

  transPos[0]=lon*RAD_TO_DEG;
  transPos[1]=lat*RAD_TO_DEG;
  std::cout << " making point from transPos " << std::endl;
  std::cout << " projArgs is currently: " << std::endl;
  for (int junk=0; junk<projArgs.size();++junk)
  {
    std::cout << "   " << projArgs[junk];
  }

  DFLib::Proj::Point transPoint(transPos,projArgs);
  transPos=transPoint.getXY();

  std::cout << " Transmitter " << " at X=" << transPos[0]
       << " Y= " << transPos[1] << std::endl;


  // Now read receiver lon/lats from stdin.  These are in dms format per
  // proj.4 standard, space delimited.
  while (!std::cin.eof())
  {
    double temp_sigma;
    char junk_space;
    DFLib::Proj::Report *reportPtr;
    std::vector<double> tempVector(2);
    std::string reportName;
    std::ostringstream ost;

    std::cin.get(dms_string,sizeof(dms_string),' ');
    if (std::cin.eof())
      break;
    lon=dmstor(dms_string,NULL);
    // get the space: 
    std::cin.get(junk_space);
    std::cin.get(dms_string,sizeof(dms_string),' ');
    if (std::cin.eof())
      break;
    lat=dmstor(dms_string,NULL);
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
    reportPtr = new DFLib::Proj::Report(tempVector,bearing,temp_sigma,
                                          ost.str(),projArgs);



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
    double rb=dynamic_cast<DFLib::Proj::Report const *>(rColl.getReport(i))->getBearing();
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

    pointsFile << i << "|"<<receiverLoc[0]<<"|"<<receiverLoc[1]<<"|Receiver "<<i<<std::endl;
  }
  gnuplotFile << std::endl;

  rColl.computeLeastSquaresFix(LS_fix);
  rColl.computeFixCutAverage(FixCutAverage,FCA_stddev);
  std::vector <double> LS_point=LS_fix.getXY();
  std::vector <double> FCA_point=FixCutAverage.getXY();

  gnuplotFile << "replot " << LS_point[0] << "," << LS_point[1] << " with points title \"LS Fix\"" << std::endl;
  gnuplotFile << "replot " << transPos[0] << "," << transPos[1] << " with points title \"Actual Location\"" << std::endl;
  gnuplotFile << "replot " << FCA_point[0] << "," << FCA_point[1] << " with points title \"Fix Cut Average\"" << std::endl;

  pointsFile << 100 << "|"<<LS_point[0]<<"|"<<LS_point[1]<<"|LSFix" <<std::endl;
  pointsFile << 101 << "|"<<FCA_point[0]<<"|"<<FCA_point[1]<<"|FCA" <<std::endl;
  pointsFile << 666 << "|"<<transPos[0]<<"|"<<transPos[1]<<"|Transmitter" <<std::endl;

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

  /*
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
  */

  NR_fix.resize(2);

//   LS_point=LS_fix.getXY();
//   // write out a grid of 5 meter "pixels" showing function values
//   for (i=-500;i<=500;++i)
//   {
//     for (j=-500;j<=500;++j)
//     {
//       NR_fix[0] = LS_point[0]+10.0*i;
//       NR_fix[1] = LS_point[1]+10.0*j;
//       gridFile << rColl.computeCostFunction(NR_fix) << " ";
//     }
//     gridFile << std::endl;
//   }

#if 0
#ifdef HAVE_LIBGDAL
#define RASTSIZ (2049)
#define RWID ((RASTSIZ-1)/2)
#define PIXSIZ (80.0)
  if (geotiffRequested)
  {
    // create a RASTSIZxRASTSIZ buffer centered on LS point 
    double rasterBuff[RASTSIZ*RASTSIZ];
    for (i=RWID; i>=-RWID; --i)
    {
      NR_fix[1]=LS_point[1]+PIXSIZ*i;
      for (j=-RWID; j<=RWID; ++j)
      {
        NR_fix[0]=LS_point[0]+PIXSIZ*j;
        rasterBuff[RASTSIZ*(RWID-i)+RWID+j]=rColl.computeCostFunction(NR_fix);
      }
    }
    // Now create a geotiff and write out our buffer

    const char *pszFormat = "GTiff";
    GDALDriver *poDriver;

    GDALAllRegister();
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

    if( poDriver != NULL )
    {
      GDALDataset *poDstDS;       
      char **papszOptions = NULL;
      
      poDstDS = poDriver->Create( geotiffName.c_str(), RASTSIZ, RASTSIZ, 1, 
                                  GDT_Float64, papszOptions );

      
      double adfGeoTransform[6] = { LS_point[0]-RWID*PIXSIZ, PIXSIZ, 0, LS_point[1]+RWID*PIXSIZ, 0, -PIXSIZ };
      char *pszSRS_WKT="PROJCS[\"unnamed\",GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433],AUTHORITY[\"EPSG\",\"4326\"]],PROJECTION[\"Mercator_1SP\"],PARAMETER[\"central_meridian\",0],PARAMETER[\"scale_factor\",1],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]]]";

      GDALRasterBand *poBand;

      poDstDS->SetGeoTransform( adfGeoTransform );
      poDstDS->SetProjection( pszSRS_WKT );
      poBand = poDstDS->GetRasterBand(1);
      poBand->RasterIO( GF_Write, 0, 0, RASTSIZ, RASTSIZ, 
                        rasterBuff, RASTSIZ, RASTSIZ, GDT_Float64, 0, 0 );    

      /* Once we're done, close properly the dataset */
      GDALClose( (GDALDatasetH) poDstDS );

    }
    else
    {
      std::cerr << "Could not open GeoTiff driver." << std::endl;
    }
  }
#endif
#endif
  // Now try Conjugate Gradients on Jml, always starting from OV fix.
  j=0;

  DFLib::Proj::Point NRPoint=LS_fix;

  rColl.computeMLFix(NRPoint);

  // Now let's check how far we are from a selected receiver
  const DFLib::Proj::Report *r0=dynamic_cast<const DFLib::Proj::Report *>(rColl.getReport(0));
  DFLib::Proj::Point rPoint0=r0->getReceiverPoint();

  // assure we're in WGS84 lat/lon like our ML point
  rPoint0.setUserProj(projArgs);
  std::vector<double> r0_coords=rPoint0.getUserCoords();

  NR_fix = NRPoint.getXY();
  latlon=NRPoint.getUserCoords();

  bool retryFix=false;
  bool fixFailed=false;
  double haversin_d;
  if (!(isinf(latlon[0]) || isinf(latlon[1]) || isnan(latlon[0]) || isnan(latlon[1])
        ||isinf(NR_fix[0]) || isinf(NR_fix[1]) || isnan(NR_fix[0]) || isnan(NR_fix[1])))
  {
    // compute very rough distance on sphere with haversine formula:
    double dlon=(latlon[0]-r0_coords[0])/RAD_TO_DEG;
    double dlat=(latlon[1]-r0_coords[1])/RAD_TO_DEG;
    double haversin_a=sin(dlat/2.0)*sin(dlat/2.0)+cos(latlon[1]/RAD_TO_DEG)*cos(r0_coords[1]/RAD_TO_DEG)*sin(dlon/2)*sin(dlon/2);
    double haversin_c=2*atan2(sqrt(haversin_a),sqrt(1-haversin_a));
    haversin_d=3596*haversin_c;   // miles, give or take
    std::cout << " latlon[0]-r0_coords[0]= " << latlon[0]-r0_coords[0];
    std::cout << " latlon[1]-r0_coords[1]= " << latlon[1]-r0_coords[1];
    std::cout << " dlon="<<dlon <<std::endl;
    std::cout << " dlat="<<dlat <<std::endl;
    std::cout << " haversin_a="<<haversin_a <<std::endl;
    std::cout << " haversin_c="<<haversin_c <<std::endl;
    std::cout << " haversin_d="<<haversin_d <<std::endl;
    if (haversin_d>100)    // don't freakin' trust it
    {
      retryFix=true;
    }
  }
  else
  {
    std::cout << " Simple attempt returned bogus numbers. " << std::endl;
    std::cout << "latlon[0]=" << latlon[0] << std::endl;
    std::cout << "latlon[1]=" << latlon[1] << std::endl;
    std::cout << "NR_fix[0]=" << NR_fix[0] << std::endl;
    std::cout << "NR_fix[1]=" << NR_fix[1] << std::endl;
    retryFix=true;
  }
  if (retryFix)
  {
    NRPoint=LS_fix;
    rColl.aggressiveComputeMLFix(NRPoint);

    NR_fix = NRPoint.getXY();
    latlon=NRPoint.getUserCoords();
    
    if (!(isinf(latlon[0]) || isinf(latlon[1]) || isnan(latlon[0]) || isnan(latlon[1])
          ||isinf(NR_fix[0]) || isinf(NR_fix[1]) || isnan(NR_fix[0]) || isnan(NR_fix[1])))
    {
      double dlon=(latlon[0]-r0_coords[0])/RAD_TO_DEG;
      double dlat=(latlon[1]-r0_coords[1])/RAD_TO_DEG;
      double haversin_a=sin(dlat/2.0)*sin(dlat/2.0)+cos(latlon[1]/RAD_TO_DEG)*cos(r0_coords[1]/RAD_TO_DEG)*sin(dlon/2)*sin(dlon/2);
      double haversin_c=2*atan2(sqrt(haversin_a),sqrt(1-haversin_a));
      haversin_d=3596*haversin_c;   // miles, give or take
      std::cout << " latlon[0]-r0_coords[0]= " << latlon[0]-r0_coords[0];
      std::cout << " latlon[1]-r0_coords[1]= " << latlon[1]-r0_coords[1];
      std::cout << " dlon="<<dlon <<std::endl;
      std::cout << " dlat="<<dlat <<std::endl;
      std::cout << " haversin_a="<<haversin_a <<std::endl;
      std::cout << " haversin_c="<<haversin_c <<std::endl;
      std::cout << " haversin_d="<<haversin_d <<std::endl;
      if (haversin_d>100)
      {
        fixFailed=true;
      }
      else
      {
        fixFailed=false;
      }
    }
    else
    {
      fixFailed=true;
      std::cout << " Second attempt returned bogus numbers. " << std::endl;
      std::cout << "latlon[0]=" << latlon[0] << std::endl;
      std::cout << "latlon[1]=" << latlon[1] << std::endl;
      std::cout << "NR_fix[0]=" << NR_fix[0] << std::endl;
      std::cout << "NR_fix[1]=" << NR_fix[1] << std::endl;
    }
    if (fixFailed)
    {
      std::cout << " more aggressive attempt still failed to get a reasonable fix."
           << std::endl;
    }
  }
  if (!fixFailed)
  {
    std::cout << " ML Fix is about " << haversin_d
         << " miles from a receiver." << std::endl;
    gnuplotFile << "replot " << NR_fix[0] << "," << NR_fix[1] << " with points title \"ML Fix\"" << std::endl;
    pointsFile << 102 << "|"<<NR_fix[0]<<"|"<<NR_fix[1]<<"|ML" <<std::endl;
    
    std::cout << " getting user coordinates " << std::endl;
    
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
    
    latlon=NRPoint.getXY();
    std::cout << " Mercator coords of ML fix" << latlon[0] << " , " << latlon[1] 
         << std::endl;
    std::cout << "  Offset from LS by " << latlon[0]-LS_point[0] << " , " 
         << latlon[1]-LS_point[1] << std::endl; 
    
  }
  // Lastly, try to compute the stansfield solution using LS fix as staring
  // guess:
  DFLib::Proj::Point StansfieldPoint=LS_fix;
  double am2,bm2,phi;
  try
  {
    rColl.computeStansfieldFix(StansfieldPoint,am2,bm2,phi);

    NR_fix = StansfieldPoint.getXY();
    gnuplotFile << "replot " << NR_fix[0] << "," << NR_fix[1] << " with points title \"Stansfield Fix\"" << std::endl;
    pointsFile << 103 << "|"<<NR_fix[0]<<"|"<<NR_fix[1]<<"|Stansfield" <<std::endl;
    
    std::cout << " getting user coordinates " << std::endl;
    
    latlon=StansfieldPoint.getUserCoords();
    
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
    
    std::cout << "  Longitude of Stansfield fix: " << (int) latlon[0] << "d" 
         << (latlon[0]-(int)latlon[0])*60 << "\"" << EW << std::endl;
    
    std::cout << "  Latitude of Stansfield fix: " << (int) latlon[1] << "d" 
         << (latlon[1]-(int)latlon[1])*60 << "\"" << NS << std::endl;
    
    latlon=StansfieldPoint.getXY();
    std::cout << " Mercator coords of Stansfield fix" << latlon[0] << " , " << latlon[1] 
         << std::endl;
    std::cout << "  Offset from LS by " << latlon[0]-LS_point[0] << " , " 
         << latlon[1]-LS_point[1] << std::endl; 

    if (am2>0 && bm2>0)
    {
      double a=sqrt(1/am2);
      double b=sqrt(1/bm2);
      std::cout << " Stansfield ellipse: a="<<a << "  b=" << b << " phi=" << phi
           << " rotation in degrees="<<phi*RAD_TO_DEG 
           << std::endl;
      
      double rho=sqrt(-2*log(.5));   // 50% confidence interval
      double cosphi=cos(phi);
      double sinphi=sin(phi);
      gnuplotFile << "replot " << NR_fix[0] << "+"<<a<<"*"<<rho<<"*"<<cosphi
                  <<"*cos(360.0/40000.0*t)-"<<b<<"*"<<rho<<"*"<<sinphi
                  <<"*sin(360.0/40000.0*t),"
                  <<NR_fix[1] << "+"<<a<<"*"<<rho<<"*"<<sinphi
                  <<"*cos(360.0/40000.0*t)+"<<b<<"*"<<rho<<"*"<<cosphi
                  <<"*sin(360.0/40000.0*t) w l title \"50% Stansfield confidence\"" << std::endl;
      rho=sqrt(-2*log(.25));   // 75% confidence interval
      gnuplotFile << "replot " << NR_fix[0] << "+"<<a<<"*"<<rho<<"*"<<cosphi
                  <<"*cos(360.0/40000.0*t)-"<<b<<"*"<<rho<<"*"<<sinphi
                  <<"*sin(360.0/40000.0*t),"
                  <<NR_fix[1] << "+"<<a<<"*"<<rho<<"*"<<sinphi
                  <<"*cos(360.0/40000.0*t)+"<<b<<"*"<<rho<<"*"<<cosphi
                  <<"*sin(360.0/40000.0*t) w l title \"75% Stansfield confidence\"" << std::endl;
      
      rho=sqrt(-2*log(.05));   // 95% confidence interval
      gnuplotFile << "replot " << NR_fix[0] << "+"<<a<<"*"<<rho<<"*"<<cosphi
                  <<"*cos(360.0/40000.0*t)-"<<b<<"*"<<rho<<"*"<<sinphi
                  <<"*sin(360.0/40000.0*t),"
                  <<NR_fix[1] << "+"<<a<<"*"<<rho<<"*"<<sinphi
                  <<"*cos(360.0/40000.0*t)+"<<b<<"*"<<rho<<"*"<<cosphi
                  <<"*sin(360.0/40000.0*t) w l title \"95% Stansfield confidence\"" << std::endl;
    }
    else
    {
      std::cerr << "Stansfield ellipse parameters came back bogus: 1/a^2=" << am2
           << " and 1/b^2="<< bm2 << std::endl;
    }
  }
  catch(DFLib::Util::Exception x)
  {
    std::cerr << "Ooops, got exception trying to compute Stansfield fix:" << std::endl
         << x.getEmsg()
         << std::endl;
  }


  gnuplotFile << "pause -1" << std::endl;
  gnuplotFile.close();
}
