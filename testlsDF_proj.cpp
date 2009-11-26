//-*- mode:C++ ; c-basic-offset: 2 -*-
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
#include <projects.h>
#undef XY

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
#include "DF_Report_Collection.hpp"
#include "DF_Proj_Report.hpp"

#if 0
#ifdef HAVE_GDAL_H
#include <gdal_priv.h>
#endif
#endif

using namespace std;


int main(int argc,char **argv)
{
  
  double lon,lat;
  vector<double> transPos(2,0.0);
  int i,j;
  char dms_string[128];
  
  vector<double> FCA_stddev;
  vector<double> NR_fix;
  vector<string> projArgs;
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
  vector<double> gradf;
  vector<vector<double> > jac;

  ofstream gnuplotFile("testlsDFfix.gnuplot");
  ofstream gridFile("function.grid");
  ofstream pointsFile("testlsDFfix.grasspoints");
  gnuplotFile << "set angles degrees" << endl;
  gnuplotFile << "set size square" << endl;
  gnuplotFile << "set parametric" << endl;
  gnuplotFile.precision(16); gnuplotFile.width(20);
  pointsFile.precision(16); pointsFile.width(20);
  cout.precision(16); cout.width(20);

  string progName(argv[0]);
  argv++;
  argc--;

  time_t seed;
  if (argc > 0)
  {
    string testArg(argv[0]);

    if (testArg == "--seed")
    {
      argv++;
      argc--;
      seed=atoi(argv[0]);
      argv++;
      argc--;
      cerr << " using seed " << seed << endl;
    }
  }
  else
  {
    seed=time(NULL);
    cerr << " using time " << seed << " as random number seed." << endl;
  }
#ifdef _MSC_VER
  srand(seed);
#else
  srand48(seed);
#endif

#if 0
#ifdef HAVE_LIBGDAL  
  testArg=argv[0];
  string geotiffName;
  bool geotiffRequested=false;

  if (testArg == "--geotiff")
  {
    cout<< "You asked for a geotiff ... not yet." << endl;
    argv++;
    argc--;
    geotiffRequested=true;
    geotiffName=argv[0];
    cout << "Geotiff name: " << geotiffName << endl;
    argv++;
    argc--;
  }
#endif
#endif

  if (argc < 2)
  {
    cerr << "Usage: " << progName;
#if 0
#ifdef HAVE_LIBGDAL    
    cerr << " [--geotiff <geotiffname>]";
#endif
#endif
    cerr << " <trans lon> <trans lat> " << endl;
    cerr << " Remember to pipe list of receiver lon/lats into stdin!" << endl;
    exit(1);
  }

  lon=dmstor(argv[0],NULL);
  lat=dmstor(argv[1],NULL);

  cout << "Transmitter location in decimal degrees: Lon: " << lon*RAD_TO_DEG
       << " Lat: " << lat*RAD_TO_DEG << endl;

  transPos[0]=lon*RAD_TO_DEG;
  transPos[1]=lat*RAD_TO_DEG;
  cout << " making point from transPos " << endl;
  cout << " projArgs is currently: " << endl;
  for (int junk=0; junk<projArgs.size();++junk)
  {
    cout << "   " << projArgs[junk];
  }

  DFLib::Proj::Point transPoint(transPos,projArgs);
  transPos=transPoint.getXY();

  cout << " Transmitter " << " at X=" << transPos[0]
       << " Y= " << transPos[1] << endl;


  // Now read receiver lon/lats from stdin.  These are in dms format per
  // proj.4 standard, space delimited.
  while (!cin.eof())
  {
    double temp_sigma;
    char junk_space;
    DFLib::Proj::Report *reportPtr;
    vector<double> tempVector(2);
    string reportName;
    ostringstream ost;

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
    tempVector[0]=lon*RAD_TO_DEG;
    tempVector[1]=lat*RAD_TO_DEG;
    double bearing=0; // temporary

    ost << "report " << rColl.size();
    reportPtr = new DFLib::Proj::Report(tempVector,bearing,temp_sigma,
                                          ost.str(),projArgs);



    bearing=reportPtr->computeBearingToPoint(transPos)*RAD_TO_DEG;
    cout << " True bearing to transmitter is " << bearing << endl;

    bearing += rand_gen.getRandom();
    cout << " after randomizing, bearing to transmitter is " << bearing << endl;
    reportPtr->setBearing(bearing);

    rColl.addReport(reportPtr);


  }
  gnuplotFile << "plot [t=0:40000] ";

  cout << "Receiver locations in mercator: " << endl;
  for (i=0;i<rColl.size();++i)
  {
    vector<double> receiverLoc = 
      rColl.getReceiverLocationXY(i);
    double rb=dynamic_cast<DFLib::Proj::Report const *>(rColl.getReport(i))->getBearing();
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

    pointsFile << i << "|"<<receiverLoc[0]<<"|"<<receiverLoc[1]<<"|Receiver "<<i<<endl;
  }
  gnuplotFile << endl;

  rColl.computeLeastSquaresFix(LS_fix);
  rColl.computeFixCutAverage(FixCutAverage,FCA_stddev);
  vector <double> LS_point=LS_fix.getXY();
  vector <double> FCA_point=FixCutAverage.getXY();

  gnuplotFile << "replot " << LS_point[0] << "," << LS_point[1] << " with points title \"LS Fix\"" << endl;
  gnuplotFile << "replot " << transPos[0] << "," << transPos[1] << " with points title \"Actual Location\"" << endl;
  gnuplotFile << "replot " << FCA_point[0] << "," << FCA_point[1] << " with points title \"Fix Cut Average\"" << endl;

  pointsFile << 100 << "|"<<LS_point[0]<<"|"<<LS_point[1]<<"|LSFix" <<endl;
  pointsFile << 101 << "|"<<FCA_point[0]<<"|"<<FCA_point[1]<<"|FCA" <<endl;
  pointsFile << 666 << "|"<<transPos[0]<<"|"<<transPos[1]<<"|Transmitter" <<endl;

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
  
  vector<double> latlon=LS_fix.getUserCoords();

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

  cout << "  Longitude of LS fix: " << (int) latlon[0] << "d" 
       << (latlon[0]-(int)latlon[0])*60 << "\"" << EW << endl;

  cout << "  Latitude of LS fix: " << (int) latlon[1] << "d" 
       << (latlon[1]-(int)latlon[1])*60 << "\"" << NS << endl;

  /*
  //  for (double minCutAngle=0; minCutAngle < 50; minCutAngle += 5.0)
  for (double minCutAngle=0; minCutAngle < 5; minCutAngle += 5.0)
  {
    if (rColl.computeFixCutAverage(FixCutAverage,FCA_stddev,minCutAngle))
    {
      vector<double> latlon=FixCutAverage.getUserCoords();
	
	
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
	
      cout << "  Longitude of Fix Cut Average (min cut angle="<< minCutAngle <<"): " << (int) latlon[0] << "d" 
           << (latlon[0]-(int)latlon[0])*60 << "\"" << EW << endl;
	
      cout << "  Latitude of Fix Cut Average (min cut angle="<< minCutAngle <<"): " << (int) latlon[1] << "d" 
           << (latlon[1]-(int)latlon[1])*60 << "\"" << NS << endl;
	
      cout << "   Std Dev of FCA = (" << FCA_stddev[0] << " , " 
           << " , " << FCA_stddev[1] << ")" << endl;
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
//     gridFile << endl;
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
      cerr << "Could not open GeoTiff driver." << endl;
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
  vector<double> r0_coords=rPoint0.getUserCoords();

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
    cout << " latlon[0]-r0_coords[0]= " << latlon[0]-r0_coords[0];
    cout << " latlon[1]-r0_coords[1]= " << latlon[1]-r0_coords[1];
    cout << " dlon="<<dlon <<endl;
    cout << " dlat="<<dlat <<endl;
    cout << " haversin_a="<<haversin_a <<endl;
    cout << " haversin_c="<<haversin_c <<endl;
    cout << " haversin_d="<<haversin_d <<endl;
    if (haversin_d>100)    // don't freakin' trust it
    {
      retryFix=true;
    }
  }
  else
  {
    cout << " Simple attempt returned bogus numbers. " << endl;
    cout << "latlon[0]=" << latlon[0] << endl;
    cout << "latlon[1]=" << latlon[1] << endl;
    cout << "NR_fix[0]=" << NR_fix[0] << endl;
    cout << "NR_fix[1]=" << NR_fix[1] << endl;
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
      cout << " latlon[0]-r0_coords[0]= " << latlon[0]-r0_coords[0];
      cout << " latlon[1]-r0_coords[1]= " << latlon[1]-r0_coords[1];
      cout << " dlon="<<dlon <<endl;
      cout << " dlat="<<dlat <<endl;
      cout << " haversin_a="<<haversin_a <<endl;
      cout << " haversin_c="<<haversin_c <<endl;
      cout << " haversin_d="<<haversin_d <<endl;
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
      cout << " Second attempt returned bogus numbers. " << endl;
      cout << "latlon[0]=" << latlon[0] << endl;
      cout << "latlon[1]=" << latlon[1] << endl;
      cout << "NR_fix[0]=" << NR_fix[0] << endl;
      cout << "NR_fix[1]=" << NR_fix[1] << endl;
    }
    if (fixFailed)
    {
      cout << " more aggressive attempt still failed to get a reasonable fix."
           << endl;
    }
  }
  if (!fixFailed)
  {
    cout << " ML Fix is about " << haversin_d
         << " miles from a receiver." << endl;
    gnuplotFile << "replot " << NR_fix[0] << "," << NR_fix[1] << " with points title \"ML Fix\"" << endl;
    pointsFile << 102 << "|"<<NR_fix[0]<<"|"<<NR_fix[1]<<"|ML" <<endl;
    
    cout << " getting user coordinates " << endl;
    
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
    
    cout << "  Longitude of ML fix: " << (int) latlon[0] << "d" 
         << (latlon[0]-(int)latlon[0])*60 << "\"" << EW << endl;
    
    cout << "  Latitude of ML fix: " << (int) latlon[1] << "d" 
         << (latlon[1]-(int)latlon[1])*60 << "\"" << NS << endl;
    
    latlon=NRPoint.getXY();
    cout << " Mercator coords of ML fix" << latlon[0] << " , " << latlon[1] 
         << endl;
    cout << "  Offset from LS by " << latlon[0]-LS_point[0] << " , " 
         << latlon[1]-LS_point[1] << endl; 
    
  }
  // Lastly, try to compute the stansfield solution using LS fix as staring
  // guess:
  DFLib::Proj::Point StansfieldPoint=LS_fix;
  double am2,bm2,phi;
  try
  {
    rColl.computeStansfieldFix(StansfieldPoint,am2,bm2,phi);

    NR_fix = StansfieldPoint.getXY();
    gnuplotFile << "replot " << NR_fix[0] << "," << NR_fix[1] << " with points title \"Stansfield Fix\"" << endl;
    pointsFile << 103 << "|"<<NR_fix[0]<<"|"<<NR_fix[1]<<"|Stansfield" <<endl;
    
    cout << " getting user coordinates " << endl;
    
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
    
    cout << "  Longitude of Stansfield fix: " << (int) latlon[0] << "d" 
         << (latlon[0]-(int)latlon[0])*60 << "\"" << EW << endl;
    
    cout << "  Latitude of Stansfield fix: " << (int) latlon[1] << "d" 
         << (latlon[1]-(int)latlon[1])*60 << "\"" << NS << endl;
    
    latlon=StansfieldPoint.getXY();
    cout << " Mercator coords of Stansfield fix" << latlon[0] << " , " << latlon[1] 
         << endl;
    cout << "  Offset from LS by " << latlon[0]-LS_point[0] << " , " 
         << latlon[1]-LS_point[1] << endl; 

    if (am2>0 && bm2>0)
    {
      double a=sqrt(1/am2);
      double b=sqrt(1/bm2);
      cout << " Stansfield ellipse: a="<<a << "  b=" << b << " phi=" << phi
           << " rotation in degrees="<<phi*RAD_TO_DEG 
           << endl;
      
      double rho=sqrt(-2*log(.5));   // 50% confidence interval
      double cosphi=cos(phi);
      double sinphi=sin(phi);
      gnuplotFile << "replot " << NR_fix[0] << "+"<<a<<"*"<<rho<<"*"<<cosphi
                  <<"*cos(360.0/40000.0*t)-"<<b<<"*"<<rho<<"*"<<sinphi
                  <<"*sin(360.0/40000.0*t),"
                  <<NR_fix[1] << "+"<<a<<"*"<<rho<<"*"<<sinphi
                  <<"*cos(360.0/40000.0*t)+"<<b<<"*"<<rho<<"*"<<cosphi
                  <<"*sin(360.0/40000.0*t) w l title \"50% Stansfield confidence\"" << endl;
      rho=sqrt(-2*log(.25));   // 75% confidence interval
      gnuplotFile << "replot " << NR_fix[0] << "+"<<a<<"*"<<rho<<"*"<<cosphi
                  <<"*cos(360.0/40000.0*t)-"<<b<<"*"<<rho<<"*"<<sinphi
                  <<"*sin(360.0/40000.0*t),"
                  <<NR_fix[1] << "+"<<a<<"*"<<rho<<"*"<<sinphi
                  <<"*cos(360.0/40000.0*t)+"<<b<<"*"<<rho<<"*"<<cosphi
                  <<"*sin(360.0/40000.0*t) w l title \"75% Stansfield confidence\"" << endl;
      
      rho=sqrt(-2*log(.05));   // 95% confidence interval
      gnuplotFile << "replot " << NR_fix[0] << "+"<<a<<"*"<<rho<<"*"<<cosphi
                  <<"*cos(360.0/40000.0*t)-"<<b<<"*"<<rho<<"*"<<sinphi
                  <<"*sin(360.0/40000.0*t),"
                  <<NR_fix[1] << "+"<<a<<"*"<<rho<<"*"<<sinphi
                  <<"*cos(360.0/40000.0*t)+"<<b<<"*"<<rho<<"*"<<cosphi
                  <<"*sin(360.0/40000.0*t) w l title \"95% Stansfield confidence\"" << endl;
    }
    else
    {
      cerr << "Stansfield ellipse parameters came back bogus: 1/a^2=" << am2
           << " and 1/b^2="<< bm2 << endl;
    }
  }
  catch(DFLib::Util::Exception x)
  {
    cerr << "Ooops, got exception trying to compute Stansfield fix:" << endl
         << x.getEmsg()
         << endl;
  }


  gnuplotFile << "pause -1" << endl;
  gnuplotFile.close();
}
