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
// Purpose        : Abstract interface class for the Report classes.
//
// Special Notes  : 
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
// This is needed for MSVC++:
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include <iostream>
#include "DF_Abstract_Report.hpp"
#include "DF_Abstract_Point.hpp"
#include <cmath>


namespace DFLib
{
  DFLib::Abstract::Report::Report(std::string n, bool v)
    : ReportName_(n),
      validReport_(v)
  { }

  DFLib::Abstract::Report::Report(const DFLib::Abstract::Report & right)
    : ReportName_(right.ReportName_),
      validReport_(right.validReport_)
  { }

  void DFLib::Abstract::Report::computeFixCut(DFLib::Abstract::Report *Report2, 
                                              Point &returnPoint,
                                              double &cutAngle,
                                              FixStatus &fs)
  {
    std::vector<double> p2,p2p;
    std::vector<double> rp;
    double theta2;
    double thetaprime;
    double phi;
    double DEG_TO_RAD=M_PI/180.0;


    // translate myself to origin, drag Report2 along for ride:
    p2p=Report2->getReceiverLocation();

    p2p[0] -= getReceiverLocation()[0];
    p2p[1] -= getReceiverLocation()[1];

    // Thetas are always in 0<theta<2PI for the arithmetic to work:
    thetaprime  = getReportBearingRadians();
    theta2  = Report2->getReportBearingRadians();

    // Now rotate counter clockwise about the origin by thetaprime
    p2.resize(2);
    p2[0] = p2p[0]*cos(thetaprime)-p2p[1]*sin(thetaprime);
    p2[1] = p2p[1]*cos(thetaprime)+p2p[0]*sin(thetaprime);
    theta2 -= thetaprime;
    // convert theta2 to -PI<theta2<PI so our tests below work
    while (theta2 > M_PI)
      theta2 -= 2*M_PI;

    while (theta2 < -M_PI)
      theta2 += 2*M_PI;

    cutAngle = fabs(theta2);

    // if his point is in left half-plane, reflect around y axis  
    if (p2[0]<0)
    {
      p2[0] *= -1;
      theta2 *= -1;
    }

    rp.resize(2);
    rp[0]=rp[1]=0.0;

    // Special cases.  We treat degenerate fixes as "no fix" for now.

    if ((p2[0] == 0 && p2[1] == 0)    // he and I are in same spot!
        || p2[0] ==0             // he is on my beam or back beam (degenerate)
        || fabs(theta2) < 1e-6          // our beams are parallel, no fix
        || fabs(fabs(theta2)-M_PI)<1e-6) // parallel, opp direction
    {
      fs=DFLib::NO_FIX;
    }
    else
    {
      // Compute angle that line from me to him makes with X axis:
      phi = atan2(p2[1],p2[0]);
            
      // now we have my bearing as y axis, me at origin, 
      // his bearing and his position in appropriate half-plane.
      // The fix, if there is one, will be on positive Y axis.
            
      if ( theta2 > 0 || (theta2 <= -(phi+M_PI/2.0)))
        // his bearing points away from Y axis, or points at negative Y
        // (treat degenerate case where he points right
        // along bearing to me as no fix, too)
      {
        fs=DFLib::NO_FIX;
      } 
      else // his beam intersects positive y axis, we have a fix!
      {
        double xfp,yfp;
        // theta2 is negative, M_PI/2+theta2 is positive angle between 
        // horizontal and beam
        rp[1]=p2[1]+p2[0]*tan(M_PI/2+theta2);
        fs=DFLib::GOOD_FIX;
        // now rotate clockwise by my theta:
        xfp = rp[0]*cos(thetaprime)+rp[1]*sin(thetaprime);
        yfp = rp[1]*cos(thetaprime)-rp[0]*sin(thetaprime);


        rp[0]=xfp+getReceiverLocation()[0]; // translate back.
        rp[1]=yfp+getReceiverLocation()[1];
      }
    }
    returnPoint.setXY(rp);
  }

  double DFLib::Abstract::Report::computeBearingToPoint(std::vector<double> &aPoint)
  {
    double dx=aPoint[0]-getReceiverLocation()[0];
    double dy=aPoint[1]-getReceiverLocation()[1];
    double bearingToPoint=atan2(dx , dy);
    while (bearingToPoint < 0)
      bearingToPoint += 2*M_PI;
    while (bearingToPoint > 2*M_PI)
      bearingToPoint -= 2*M_PI;
    return (bearingToPoint);
  }

  double DFLib::Abstract::Report::computeDistanceToPoint(std::vector<double> &aPoint)
  {
    double dx=aPoint[0]-getReceiverLocation()[0];
    double dy=aPoint[1]-getReceiverLocation()[1];
    return (sqrt(dx*dx+dy*dy));
  }

}
