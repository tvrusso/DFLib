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
// Purpose        : Provide a trivial DF fixing program that uses DFLib
//                  to produce fixes from an input file of DF reports.
//
// Special Notes  : To use, one must provide an input text file in a rigid 
//                  format.  The format is described in comments below.
//                  This was a quick hack to show that  code could actually
//                  use DFLib to produce useful fixes from real data, but its 
//                  interface is so crude it is not really useful in real work.
//                  It's just a precursor to a real DF fixing program using 
//                  DFLib.
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
#endif

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <proj_api.h>

extern "C" {
  double dmstor(const char *, char **);
}

#include "DF_Proj_Point.hpp"
#include "DF_Report_Collection.hpp"
#include "DF_Proj_Report.hpp"
#include "Util_Misc.hpp"

using namespace std;
void printCoords(const vector<double> &latlon, const string &text);

// The magnetic declination.  This is correct in Albuquerque, NM in mid-2009.
// CHANGE THIS for your own environment!
// In this usage, east declination is positive, west negative.
#define MAG_DEC (9.8)

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
      // Read in a single line of the input file.  We do next to no error
      // checking.  We assume that if we can read a report name, we have the
      // whole line.  This needs to be much more robust in a production code.
      infile >> stationName;
      if (!infile.eof())
      {
        infile >> tempstr; // lon
        stationPos[0]=dmstor(tempstr.c_str(),NULL)*RAD_TO_DEG;
        infile >> tempstr; // lat
        stationPos[1]=dmstor(tempstr.c_str(),NULL)*RAD_TO_DEG;
        infile >> bearing;
        // The bearing as input is magnetic.  Convert to true by adding in
        // the declination (assumes East declination is positive, here)
        bearing += MAG_DEC;
        infile >> datum;
        infile >> sd;
        infile >> validity;

        // We now have the data for this line, so print it all out and
        // create a new report object.
        cout << "Station " << stationName << " at (" 
             << stationPos[0] << "," << stationPos[1] << ")" 
             << " with datum "; 
        if (datum==0)
          cout << "NAD27";
        else
          cout << "WGS84/NAD83" ;
        
        cout << " bearing " << bearing << " sd " << sd << endl;
        cout << string((validity==1)?"VALID":"IGNORE") << endl;

        // Here's the guts.  We create a report object with our characteristics
        // and with the appropriate datum.
        reportPtr = new DFLib::Proj::Report(stationPos,bearing,sd,stationName,
                                            (datum==0)?NAD27_args:WGS84_args);
        if (validity == 0)
          reportPtr->setInvalid();

        // Now add the new report to the collection
        rColl.addReport(reportPtr);
      }
    }
    cout << " Got " << rColl.size() << " reports " << endl;

    // computeLeastSquaresFix doesn't use the actual value of the point
    // we pass to it, it merely returns the answer.  We give it a bogus
    // position just to initialize the point properly.
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

    // Fix cut average also doesn't use the input point except as a place to
    // store the answer.  We initialize it to something valid by merely copying
    // the LSFix point into FCA.
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

    // computeMLFix *DOES* use the input point as a starting guess for its
    // minimization search.  Initialize to the LS fix.  It can also FAIL if
    // the geometry of bearing measurements leads to a very flat cost surface,
    // and we *should* be testing the result before reporting it.  But we're
    // not.  See testlsDF_proj.cpp for a clumsy example of testing the 
    // ML fix before proceeding.
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

    // That's it.  We've computed all the fixes, now we just display them.
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
