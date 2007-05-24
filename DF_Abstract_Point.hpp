#ifndef DF_ABSTRACT_POINT_HPP
#define DF_ABSTRACT_POINT_HPP

#include <vector>
using namespace std;

namespace DFLib
{
    namespace Abstract
    {
        /// \brief Abstract base class for generic "point" objects
        ///
        /// The purpose of this extraordinarily simple class is to
        /// provide a generic interface for low-level DF classes to
        /// obtain the X-Y coordinates of a DF receiver location (or
        /// for returning a DF fix point) without regard to how the user
        /// actually provides the data or expects to view it.
        ///
        /// The abstract interface works only with the X-Y representations
        /// of the coordinates.  It is expected that in most cases for
        /// DF work this will be Mercator projection coordinates, but 
        /// it could be any planar representation that is suitable for the
        /// user's work.
        ///

        class Point
        {
        public:
            /// \brief Set X-Y coordinates of the point
            ///
            /// \param aPosition a vector containing the X and Y coordinates
            ///

            virtual void setXY(const vector<double> &aPosition)=0;
            
            /// \brief Get X-Y coordinates of the point
            /// 
            /// Note that this is NOT declared as a const member function,
            /// because implementations \e could, for efficiency, store 
            /// both X-Y coordinates and another system of coordinates.
            /// In that case, getXY might trigger a coordinate system 
            /// transformation, and might need to tinker with internals.
            ///
            /// The reference returned, however, \e is const.  One is not
            /// allowed to tinker with the internal representation through
            /// it.
            ///
            /// \return reference to an STL vector containing the coordinates 
            ///
            virtual const vector<double> &getXY() = 0;

	  /// \brief Get coordinates in the user's coordinate system
	  ///
	  
	  virtual const vector<double> &getUserCoords() = 0;

	  /// \brief Set coordinates in the user's coordinate system
	  ///
	  
	  virtual void setUserCoords(const vector<double> &uPosition) = 0;

            /// \brief Clone Self
            ///
            /// Make a copy of yourself, and return a pointer to the copy.
            /// This is useful as a prototype: pass a pointer to any concrete 
            /// implementation to a function, and that function has a way of
            /// making other objects of precisely the same type by cloning.
            /// 
            /// \return pointer to
            virtual Point *Clone() = 0; 
        };
    }
}
#endif // DF_ABSTRACT_POINT_HPP
