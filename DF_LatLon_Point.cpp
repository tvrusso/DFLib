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
// Filename       : DF_LatLon_Point.cpp
//
// Purpose        : Implement a DFLib::Abstract::Point interface such that
//                  the "user" coordinate system is WGS84 lat/lon
//                  "XY" coordinate system is a mercator projection on the
//                  WGS84 ellipsoid.
//
//-------------------------------------------------------------------------
#include "DF_LatLon_Point.hpp"
#include "Util_Misc.hpp"
#include "DFLib_Misc_Defs.h"

#include <vector>
#include <iostream>
#include <cmath>

namespace DFLib
{
  namespace LatLon
  {

    // class Point
    Point::Point()
      : llDirty(true),
        mercDirty(false)
    {
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

      if (convertPJ == 0)
      {
        throw(Util::Exception("Failed to initialize crs_to_crs"));
      }

      theLatLon.resize(2,0.0);
      theMerc.resize(2,0.0);

      delete [] latlon_argv;
      delete [] mercator_argv;
    }

    Point::Point(const std::vector<double> &aPosition)
      : theLatLon(aPosition),
        llDirty(true),
        mercDirty(false)

    {
      // Create the proj.4 stuff.  This is highly inefficient, as it should
      // live in the base class, not each object.  Fix that.


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

      if (convertPJ == 0)
      {
        throw(Util::Exception("Failed to initialize crs_to_crs"));
      }

      // Don't bother trying to force the mercator --- we'll do that if
      // we query, because we're setting llDirty to true, just initialize
      // to junk.
      theMerc.resize(2,0.0);
      delete [] latlon_argv;
      delete [] mercator_argv;
    }

    Point::Point(const Point &right)
    {
      PJ_PROJ_INFO theProjDef = proj_pj_info(right.convertPJ);
      convertPJ = proj_create(PJ_DEFAULT_CTX,theProjDef.definition);

      if (convertPJ == 0)
      {
        throw(Util::Exception("Failed to initialize transform in copy constructor"));
      }

      if (right.mercDirty)
      {
        // Right's mercator's been changed without its LL being
        // updated, so copy its mercator values and say we're dirty.
        theMerc=right.theMerc;
        mercDirty=true;
        llDirty=false;
      }
      else
      {
        // Just copy its LL and mark dirty.
        theLatLon=right.theLatLon;
        llDirty=true;
        mercDirty=false;
      }
    }

    Point::~Point()
    {
      proj_destroy(convertPJ);
    }

    void Point::setXY(const std::vector<double> &aPosition)
    {
      theMerc = aPosition;
      // Essential to set llDirty=false, otherwise getXY will try to
      // convert ll and we'll be wrong.
      mercDirty=true;
      llDirty=false;
    }

    bool Point::isUserProjRadians() const
    {
      return (proj_angular_input(convertPJ,PJ_FWD));
    }

    const std::vector<double> & Point::getXY()
    {
      if (llDirty)
      {
        // our LL has been changed, so update
        llToMerc();
      }
      return(theMerc);
    }

    void Point::setLL(const std::vector<double> &llPosition)
    {
      theLatLon = llPosition;
      // If we set ll, must now trump mercator.
      llDirty=true; mercDirty=false;
    }

    const std::vector<double> & Point::getLL()
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
      PJ_COORD data;
      PJ_COORD newCoord;

      if (isUserProjRadians())
      {
        data.lp.lam = theLatLon[0]*DEG_TO_RAD;
        data.lp.phi = theLatLon[1]*DEG_TO_RAD;
      }
      else
      {
        data.xy.x = theLatLon[0];
        data.xy.y = theLatLon[1];
      }
      newCoord=proj_trans(convertPJ,PJ_FWD,data);
      if (std::isnan(newCoord.xy.x) || std::isnan(newCoord.xy.y))
      {
        throw(Util::Exception("Failure converting LL to Mercator"));
      }
      theMerc.resize(2);
      theMerc[0]=newCoord.xy.x;
      theMerc[1]=newCoord.xy.y;

      llDirty=false;
    }


    void Point::mercToLL()
    {
      PJ_COORD data;
      PJ_COORD newCoord;

      data.xy.x = theMerc[0];
      data.xy.y = theMerc[1];
      newCoord=proj_trans(convertPJ,PJ_INV,data);
      if (std::isnan(newCoord.lp.lam) || std::isnan(newCoord.lp.phi))
      {
        throw(Util::Exception("Failure converting Mercator to LL"));
      }
      theLatLon.resize(2);
      if (isUserProjRadians())
      {
        theLatLon[0]=data.lp.lam*RAD_TO_DEG;
        theLatLon[1]=data.lp.phi*RAD_TO_DEG;
      }
      else
      {
        theLatLon[0]=newCoord.xy.x;
        theLatLon[1]=newCoord.xy.y;
      }
      mercDirty=false;
    }
  }
}
