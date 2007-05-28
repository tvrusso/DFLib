#ifndef DF_XY_POINT_HPP
#define DF_XY_POINT_HPP

#include "port.h"

#include "DF_Abstract_Point.hpp"

namespace DFLib
{
    namespace XY
    {
        class Point : public DFLib::Abstract::Point
        {
        private:
            vector<double> myXY;
        public:
            /// Default
            CPL_DLL Point();

            /// Constructor
            CPL_DLL Point(const vector<double> &aPosition);
            /// Copy Constructor
            CPL_DLL Point(Point &right);
            // no destructor needed
            /// Set position from X-Y
            virtual CPL_DLL void setXY(const vector<double> &aPosition);
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
