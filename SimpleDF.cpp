#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <projects.h>
#undef XY

#include "DF_Report_Collection.hpp"
#include "DF_Proj_Report.hpp"
#include "Util_Misc.hpp"

using namespace std;
void printCoords(const vector<double> &latlon, const string &text);

int main(int argc, char **argv)
{

  /*!
    \brief A very simple text-based DF fixing program demonstrating how to use DFLib classes.

    SimpleDF takes as command line argument the name of a text file
    with lines in the following format:

    NAME LONGITUDE LATITUDE BEARING Datum STDDEV VALID

    where longitude and latitude are in PROJ.4 format, Bearing is magnetic
    bearing IN CENTRAL NEW MEXICO IN 2009, Datum is 1 for WGS84/NAD83 and 0 for
    NAD27, stddev is the standard deviation of bearing error expected for the 
    equipment and user, and VALID is 0 for invalid and 1 for valid reports.

    The magnetic declination is hard-coded here as 9.8 degrees, which is the
    correct declination in Albuquerque, NM in mid-2009.  The user must change 
    the appropriate line of code to fix this for the region and date.

    The DF reports in the input file are input into a DFLib::ReportCollection,
    and the various DF fixes computed.  The fixes are printed to standard 
    output.

    This primitive program shows how to use the DFLib report classes.
    It does very little error checking of the fixes after they're
    computed.  testlsDF_proj.cpp actually does more error checking.
    This program is not intended as a finished product, but as a demo.

    The program was not updated to compute Stansfield fix after the Stansfield
    fix was implemented in DFLib.

    Neither was the program updated to compute the Maximum Likelihood
    fix more aggressively than by the conjugate gradient method.  See
    the program "testlsDF_proj.cpp for how to call the more aggressive
    ML computation when the simple conjugate gradients method fails
    due to flatness of the cost surface.

  */

  string stationName;
  string tempstr;
  vector<double> stationPos(2,0);
  int datum; // 0=NAD27, 1=NAD83/WGS84
  double bearing; // magnetic
  double sd; // standard deviation
  vector<string> NAD27_args;
  vector<string> WGS84_args;
  DFLib::ReportCollection rColl;
  DFLib::Proj::Report *reportPtr;
  int validity;

  NAD27_args.push_back("proj=latlong");
  WGS84_args.push_back("proj=latlong");

  NAD27_args.push_back("datum=NAD27");
  WGS84_args.push_back("datum=WGS84");

  if (argc < 2)
  {
    cerr << "Usage: " << argv[0] << " <inputfile>" << endl;
    exit(1);
  }

  ifstream infile(argv[1],ifstream::in);
  
  if (!infile.good())
  {
    cout << " Failed to open file " << argv[1] << endl;
  }
  else
  {
    while (!infile.eof())
    {
      infile >> stationName;
      if (!infile.eof())
      {
        infile >> tempstr; // lon
        stationPos[0]=dmstor(tempstr.c_str(),NULL)*RAD_TO_DEG;
        infile >> tempstr; // lat
        stationPos[1]=dmstor(tempstr.c_str(),NULL)*RAD_TO_DEG;
        infile >> bearing;
        bearing += 9.8;  // hard coded magnetic declination
        infile >> datum;
        infile >> sd;
        infile >> validity;

        cout << "Station " << stationName << " at (" 
             << stationPos[0] << "," << stationPos[1] << ")" 
             << " with datum "; 
        if (datum==0)
          cout << "NAD27";
        else
          cout << "WGS84/NAD83" ;
        
        cout << " bearing " << bearing << " sd " << sd << endl;
        cout << string((validity==1)?"VALID":"IGNORE") << endl;

        reportPtr = new DFLib::Proj::Report(stationPos,bearing,sd,stationName,
                                            (datum==0)?NAD27_args:WGS84_args);
        if (validity == 0)
          reportPtr->setInvalid();

        rColl.addReport(reportPtr);
      }
    }
    cout << " Got " << rColl.size() << " reports " << endl;
    // Give it a bogus initial point
    DFLib::Proj::Point LSFix(stationPos,WGS84_args);
    try {
      rColl.computeLeastSquaresFix(LSFix);
    }
    catch (DFLib::Util::Exception x)
    {
      cerr << " Ooops .... got exception trying to compute LS fix: " 
           << x.getEmsg()
           << endl;
    }

    DFLib::Proj::Point FCA=LSFix;
    vector<double> FCA_stddev(2);
    try 
    {
      rColl.computeFixCutAverage(FCA,FCA_stddev);
    }
    catch (DFLib::Util::Exception x)
    {
      cerr << " Ooops .... got exception trying to compute FCA: " 
           << x.getEmsg()
           << endl;
    }

    DFLib::Proj::Point MLFix=LSFix;
    try 
    {
      rColl.computeMLFix(MLFix);
    }
    catch (DFLib::Util::Exception x)
    {
      cerr << " Ooops .... got exception trying to compute ML Fix: " 
           << x.getEmsg()
           << endl;
    }

    vector<double> LS_point=LSFix.getUserCoords();
    printCoords(LS_point,string("Least Squares Fix"));

    vector<double> FCA_point=FCA.getUserCoords();
    printCoords(FCA_point,string("Fix Cut Average"));
    cout << " Standard deviation of FCA is " << FCA_stddev[0] << " longitude"
         << " and " << FCA_stddev[1] << " latitude." << endl;
    vector<double> ML_point=MLFix.getUserCoords();
    printCoords(ML_point,string("Maximum Likelihood Fix"));
  }
}

void printCoords(const vector<double> &latlon,const string &text)
{
  int latfac=1;  int lonfac=1;
  char EW,NS;

  EW='E';
  NS='N';

  if (latlon[0] < 0)
  {
    lonfac=-1;
    EW='W';
  }

  if (latlon[1] < 0)
  {
    latfac=-1;
    NS='S';
  }

  cout << " Longitude of " << text << ": " 
       << static_cast<int>(latlon[0]*lonfac)
       << "d" 
       << (latlon[0]*lonfac-static_cast<int>(latlon[0]*lonfac))*60 
       << "'" << EW 
       << endl;

  cout << " Latitude of " << text << ": " 
       << static_cast<int>(latlon[1]*latfac)
       << "d" 
       << (latlon[1]*latfac-static_cast<int>(latlon[1]*latfac))*60 
       << "'" << NS 
       << endl;

}
