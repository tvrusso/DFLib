#ifndef DF_XY_POINT_HPP
#define DF_XY_POINT_HPP

#include "DFLib_port.h"

#include "DF_Abstract_Point.hpp"

namespace DFLib
{
  namespace XY
  {
    class CPL_DLL Point : public DFLib::Abstract::Point
    {
    private:
      vector<double> myXY;
    public:
      /// Default
      Point();

      /// Constructor
      Point(const vector<double> &aPosition);
      /// Copy Constructor
      Point(Point &right);
      // no destructor needed
      /// Set position from X-Y
      virtual void setXY(const vector<double> &aPosition);
      /// Set X-Y position
      virtual const vector<double> &getXY();

      /// get User Coords (wrapper as required by abstract interface)
      virtual const vector<double> &getUserCoords() { return getXY();};

      /// set User Coords (wrapper as required by abstract interface)
      virtual void setUserCoords(const vector<double> &uPosition) 
      { setXY(uPosition);};

      /// Clone Self
      virtual Point *Clone();
    };
  }
}

#endif
