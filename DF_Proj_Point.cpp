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
// Filename       : DF_Proj_Point.cpp
//
// Purpose        : Implement a DFLib::Abstract::Point interface such that
//                  the "user" coordinate system is whatever the user has
//                  specified with a proj.4 spatial reference system, and the
//                  "XY" coordinate system is a mercator projection on the
//                  WGS84 ellipsoid.
//
//-------------------------------------------------------------------------
#include "DF_Proj_Point.hpp"
#include "Util_Misc.hpp"
#include "DFLib_Misc_Defs.h"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>

namespace DFLib
{
  namespace Proj
  {

    // class Point
    Point::Point(const std::vector<double> &uPosition,const std::vector<std::string> &projArgs)
      : theUserCoords(uPosition),
        userDirty(true),
        mercDirty(false)
    {
      // Create the proj.4 stuff.  This is highly inefficient, as it should
      // live in the base class, not each object.  Fix that.

      std::string mercator_args="+proj=merc +ellps=WGS84 +lat_ts=0";
      char *mercator_argv= new char [mercator_args.size()+1];
      strcpy(mercator_argv,mercator_args.c_str());

      int numUserArgs = projArgs.size();
      std::string user_args;

      if (projArgs[0][0] != '+')
        user_args = "+"+projArgs[0];
      else
        user_args = projArgs[0];

      for (int i=1; i<numUserArgs; ++i)
      {
        if (projArgs[i][0] != '+')
        {
          user_args += " +";
        }
        else
        {
          user_args += " ";
        }
        user_args += projArgs[i];
      }

      char *user_argv = new char [user_args.size()+1];
      strcpy(user_argv,user_args.c_str());

      convertPJ = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
                                                 user_argv,
                                                 mercator_argv,
                                                 NULL);

      if (convertPJ == 0)
      {
        throw(Util::Exception("Failed to initialize crs_to_crs"));
      }

      // Don't bother trying to force the mercator --- we'll do that if
      // we query, because we're setting userDirty to true, just initialize
      // to junk.
      theMerc.resize(2,0.0);

      delete [] user_argv;
      delete [] mercator_argv;

    }

    Point::Point(const Point &right)
    {

      // Must *COPY* the definition, not the pointer to it!
      // Get the text version of the definition
      PJ_PROJ_INFO theProjDef = proj_pj_info(right.convertPJ);
      // generate a new projUJ pointer
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
        userDirty=false;
      }
      else
      {
        // Just copy its LL and mark dirty.
        theUserCoords=right.theUserCoords;
        userDirty=true;
        mercDirty=false;
      }
    }

    Point::~Point()
    {
      proj_destroy(convertPJ);
    }

    Point& Point::operator=(const Point& rhs)
    {
      if (this == &rhs) return *this;

      if (convertPJ)
      {
        convertPJ=proj_destroy(convertPJ);
        convertPJ=0;
      }

      PJ_PROJ_INFO theProjDef = proj_pj_info(rhs.convertPJ);
      convertPJ = proj_create(PJ_DEFAULT_CTX,theProjDef.definition);


      mercDirty=rhs.mercDirty;
      userDirty=rhs.userDirty;

      // This is OK, because we're copying the dirty bools, too.
      // we have to copy both of these, otherwise we risk copying only the
      // ones that haven't been updated for consistency!
      theMerc=rhs.theMerc;
      theUserCoords=rhs.theUserCoords;

      return (*this);
    }

    void Point::setXY(const std::vector<double> &mPosition)
    {
      theMerc = mPosition;
      // Essential to set userDirty=false, otherwise getXY will try to
      // convert ll and we'll be wrong.
      mercDirty=true;
      userDirty=false;
    }

    void Point::setUserProj(const std::vector<std::string> &projArgs)
    {
      int numUserArgs = projArgs.size();
      std::string user_args;
      std::string mercator_args="+proj=merc +ellps=WGS84 +lat_ts=0";
      char *mercator_argv= new char [mercator_args.size()+1];
      strcpy(mercator_argv,mercator_args.c_str());


      if (projArgs[0][0] != '+')
        user_args = "+"+projArgs[0];
      else
        user_args = projArgs[0];

      for (int i=1; i<numUserArgs; ++i)
      {
        if (projArgs[i][0] != '+')
        {
          user_args += " +";
        }
        else
        {
          user_args += " ";
        }
        user_args += projArgs[i];
      }

      char *user_argv = new char [user_args.size()+1];
      strcpy(user_argv,user_args.c_str());

      // before we clobber userProj, if we have already have valid userProj
      //  already then we might also have
      // valid data, and need to make sure our coordinates are correct in the
      // new system:
      if (convertPJ)
      {
        if (userDirty)
        {
          userToMerc(); // make sure our mercator coordinates are up-to-date
        }
        proj_destroy(convertPJ); // free the existing projection
        convertPJ = NULL;
      }

      convertPJ = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
                                                 user_argv,
                                                 mercator_argv,
                                                 NULL);
      PJ_PROJ_INFO theProjDef = proj_pj_info(convertPJ);

      if (convertPJ == 0)
      {
        throw(Util::Exception("Failed to initialize crs_to_crs"));
      }

      // Now, we have just changed the projection, so mark mercDirty as if
      // we had changed the mercator coordinates ourselves.
      mercDirty=true;

      delete [] user_argv;
      delete [] mercator_argv;

    }

    bool Point::isUserProjRadians() const
    {
      return (proj_angular_input(convertPJ,PJ_FWD));
    }

    const std::vector<double> & Point::getXY()
    {
      if (userDirty)
      {
        // our LL has been changed, so update
        userToMerc();
      }
      return(theMerc);
    }

    void Point::setUserCoords(const std::vector<double> &llPosition)
    {
      theUserCoords = llPosition;
      // If we set ll, must now trump mercator.
      userDirty=true; mercDirty=false;
    }

    const std::vector<double> & Point::getUserCoords()
    {
      if (mercDirty)
      {
        // our XY has been changed, so update
        mercToUser();
      }
      return(theUserCoords);
    }

    Point * Point::Clone()
    {
      Point *retPoint;
      retPoint = new Point(*this);
      return retPoint;
    }

    void Point::userToMerc()
    {
      PJ_COORD data;
      PJ_COORD newCoord;

      if (isUserProjRadians())
      {
        data.lp.lam = theUserCoords[0]*DEG_TO_RAD;
        data.lp.phi = theUserCoords[1]*DEG_TO_RAD;
      }
      else
      {
        data.xy.x = theUserCoords[0];
        data.xy.y = theUserCoords[1];
      }

      newCoord=proj_trans(convertPJ,PJ_FWD,data);

      if (std::isnan(newCoord.xy.x) || std::isnan(newCoord.xy.y))
      {
        throw(Util::Exception("Failure converting user coords to Mercator"));
      }

      theMerc.resize(2);
      theMerc[0]=newCoord.xy.x;
      theMerc[1]=newCoord.xy.y;

      userDirty=false;
    }


    void Point::mercToUser()
    {
      PJ_COORD data;
      PJ_COORD newCoord;

      data.xy.x = theMerc[0];
      data.xy.y = theMerc[1];
      newCoord=proj_trans(convertPJ,PJ_INV,data);
      if (std::isnan(newCoord.lp.lam) || std::isnan(newCoord.lp.phi))
      {
        throw(Util::Exception("Failure converting Mercator to user coords"));
      }

      theUserCoords.resize(2);
      if (isUserProjRadians())
      {
        theUserCoords[0]=newCoord.lp.lam*RAD_TO_DEG;
        theUserCoords[1]=newCoord.lp.phi*RAD_TO_DEG;
      }
      else
      {
        theUserCoords[0]=newCoord.xy.x;
        theUserCoords[1]=newCoord.xy.y;
      }
      mercDirty=false;
    }
  }
}
