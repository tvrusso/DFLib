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
// Purpose        : Implement a DFLib::Abstract::Point interface such that
//                  the "user" coordinate system is WGS84 lat/lon
//                  "XY" coordinate system is a mercator projection on the 
//                  WGS84 ellipsoid.
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
#include "DF_LatLon_Point.hpp"
#include "Util_Misc.hpp"

#include <vector>
#include <iostream>
using namespace std;

namespace DFLib
{
  namespace LatLon
  {

    // class Point
    Point::Point()
      : llDirty(true),
        mercDirty(false)
    {
      char *latlon_argv[2]={"proj=latlong",
                            "datum=WGS84"};
      char *mercator_argv[3]={"proj=merc",
                              "ellps=WGS84",
                              "lat_ts=0"};
      if (!(latlonProj = pj_init(2,latlon_argv)))
      {
        throw(Util::Exception("Failed to initialize lat/lon projection"));
      }

      if (!(mercProj = pj_init(3,mercator_argv)))
      {
        throw(Util::Exception("Failed to initialize mercator projection"));
      }

      theLatLon.resize(2,0.0);
      theMerc.resize(2,0.0);
    }

    Point::Point(const vector<double> &aPosition)
      : theLatLon(aPosition),
        llDirty(true),
        mercDirty(false)

    {
      // Create the proj.4 stuff.  This is highly inefficient, as it should
      // live in the base class, not each object.  Fix that.

      char *latlon_argv[2]={"proj=latlong",
                            "datum=WGS84"};
      char *mercator_argv[3]={"proj=merc",
                              "ellps=WGS84",
                              "lat_ts=0"};
      if (!(latlonProj = pj_init(2,latlon_argv)))
      {
        throw(Util::Exception("Failed to initialize lat/lon projection"));
      }

      if (!(mercProj = pj_init(3,mercator_argv)))
      {
        throw(Util::Exception("Failed to initialize mercator projection"));
      }

      // Don't bother trying to force the mercator --- we'll do that if
      // we query, because we're setting llDirty to true, just initialize
      // to junk.
      theMerc.resize(2,0.0);
    }

    Point::Point(const Point &right)
      : latlonProj(right.latlonProj),
        mercProj(right.mercProj)
    {
      if (right.mercDirty)
      {
        // Right's mercator's been changed without its LL being 
        // updated, so copy its mercator values and say we're dirty.
        theMerc=right.theMerc;
        mercDirty=true;
      } 
      else
      {
        // Just copy its LL and mark dirty.
        theLatLon=right.theLatLon;
        llDirty=true;
      }
    }

    void Point::setXY(const vector<double> &aPosition)
    {
      theMerc = aPosition;
      // Essential to set llDirty=false, otherwise getXY will try to 
      // convert ll and we'll be wrong.
      mercDirty=true;  
      llDirty=false;
    }

    const vector<double> & Point::getXY()
    {
      if (llDirty)
      {
        // our LL has been changed, so update
        llToMerc();
      }
      return(theMerc);
    }

    void Point::setLL(const vector<double> &llPosition)
    {
      theLatLon = llPosition;
      // If we set ll, must now trump mercator.
      llDirty=true; mercDirty=false;
    }

    const vector<double> & Point::getLL()
    {
      if (mercDirty)
      {
        // our LL has been changed, so update
        mercToLL();
      }
      return(theLatLon);
    }
    
    Point * Point::Clone()
    {
      Point *retPoint;
      retPoint = new Point(*this);
      return retPoint;
    }

    void Point::llToMerc()
    {
      projUV data;
      double z;
      data.u = theLatLon[0]*DEG_TO_RAD;
      data.v = theLatLon[1]*DEG_TO_RAD;
      z=0;
      if (pj_transform(latlonProj,mercProj,1,0,&(data.u),&(data.v),&z) != 0)
      {
        throw(Util::Exception("Failure converting LL to Mercator"));
      }
      theMerc.resize(2);
      theMerc[0]=data.u;
      theMerc[1]=data.v;

      llDirty=false;
    }      


    void Point::mercToLL()
    {
      projUV data;
      double z;
      data.u = theMerc[0];
      data.v = theMerc[1];
      z=0;
      if (pj_transform(mercProj,latlonProj,1,0,&(data.u),&(data.v),&z) != 0)
      {
        throw(Util::Exception("Failure converting Mercator to LL"));
      }
      theLatLon.resize(2);
      theLatLon[0]=data.u*RAD_TO_DEG;
      theLatLon[1]=data.v*RAD_TO_DEG;

      mercDirty=false;
    }      
  }
}
