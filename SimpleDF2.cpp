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
// Special Notes  : This is basically the same program as SimpleDF, but using
//                  the DFLib::ProjReportCollection instead of 
//                  ReportCollection.
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
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <proj_api.h>

extern "C" {
  double dmstor(const char *, char **);
}


#include "DF_Proj_Point.hpp"
#include "DF_Proj_Report.hpp"
#include "DF_ProjReport_Collection.hpp"
#include "Util_Misc.hpp"

using namespace std;
void printCoords(const vector<double> &latlon, const string &text);

int main(int argc, char **argv)
{

  string stationName;
  string tempstr;
  vector<double> stationPos(2,0);
  int datum; // 0=NAD27, 1=NAD83/WGS84
  double bearing; // magnetic
  double sd; // standard deviation
  vector<string> NAD27_args;
  vector<string> WGS84_args;
  DFLib::ProjReportCollection rColl;
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
