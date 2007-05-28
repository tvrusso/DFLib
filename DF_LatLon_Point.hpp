#ifndef DF_LATLON_POINT_HPP
#define DF_LATLON_POINT_HPP

#include <vector>
#include "projects.h"
// projects.h rudely defines this, and we don't want it:
#undef XY
#include "DF_Abstract_Point.hpp"
using namespace std;

namespace DFLib
{
    namespace LatLon
    {
        class Point : public DFLib::Abstract::Point
        {
        private:
            vector<double> theMerc;
            bool mercDirty;
            vector<double> theLatLon;
            bool llDirty;
            projPJ latlonProj, mercProj;
        public:
            /// \brief Default Constructor
            ///
            Point();

            /// \brief Constructor
            ///
            /// \param llPosition coordinates <em>in Lat Lon</em>
            Point(const vector<double> &llPosition);


            /// \brief Copy Constructor
            Point(const Point &right);

            /// \brief set mercator projection (XY) position
            /// 
            /// \param aPosition coordinates in mercator projection
            ///
            virtual void setXY(const vector<double> &aPosition);
            /// \brief get mercator projection (XY) position
            /// 
            /// \return vector of coordinates in mercator projection
            ///
            virtual const vector<double> &getXY();


	  /// \brief get user position
	  ///
	  /// This is just a wrapper for getLL as needed by the abstract
	  /// interface

	  virtual const vector<double> &getUserCoords() { return getLL(); };

	  /// \brief set user position
	  ///
	  /// This is just a wrapper for setLL as needed by the abstract
	  /// interface

	  virtual void setUserCoords(const vector<double> &uPosition)  
	  { setLL(uPosition); };
	  
            /// \brief set lat/lon position
            /// 
            /// \param llPosition coordinates in Lat/Lon
            ///
            void setLL(const vector<double> &llPosition);
            /// \brief get Lat/Lon position
            /// 
            /// \return vector of coordinates in lat/lon
            ///
            const vector<double> &getLL();


            Point * Clone();
        private:
            void mercToLL();
            void llToMerc();
        };
    }
}
#endif // DF_LATLON_POINT_HPP
