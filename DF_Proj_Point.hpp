#ifndef DF_PROJ_POINT_HPP
#define DF_PROJ_POINT_HPP

#include "port.h"
#include <vector>
#include "projects.h"
// projects.h rudely defines this, and we don't want it:
#undef XY
#include "DF_Abstract_Point.hpp"
using namespace std;

namespace DFLib
{
  namespace Proj
  {
    class CPL_DLL Point : public DFLib::Abstract::Point
    {
    private:
      vector<double> theMerc;
      bool mercDirty;
      vector<double> theUserCoords;
      bool userDirty;
      projPJ userProj, mercProj;
    public:

      /// \brief Constructor
      ///
      /// \param uPosition coordinates <em>in user's coordinates</em>
      /// \param projArgs a vector of strings representing the Proj.4
      /// description of the user's coordinate system.
      ///
      Point(const vector<double> &uPosition,const vector<string> &projArgs);


      /// \brief Copy Constructor
      Point(const Point &right);

      /// \brief destructor
      ~Point();

      /// \brief assignment operator
      Point& operator=(const Point& rhs);

      /// \brief set user projection information
      ///
      /// \param projArgs vector of strings representing the new projection
      ///
      /// This method can be used to change the coordinate system of a point.
      /// When called, the Mercator representation of the point is updated,
      /// then the user projection is changed.  The next call to getUserCoords
      /// will therefore return the point's coordinates in the new coordinate
      /// system.

      void setUserProj(const vector<string> &projArgs);

      /// \brief return true if user projection is a lat/lon system
      ///
      ///  Primarily useful for deciding how to display coordinates 
      ///  (as for example whether to display them in sexigesimal or
      ///  decimal representation)
      bool isUserProjLatLong() const;

      /// \brief set mercator projection (XY) position
      /// 
      /// \param aPosition coordinates in mercator projection
      ///
      virtual void setXY(const vector<double> &mPosition);
      /// \brief get mercator projection (XY) position
      /// 
      /// \return vector of coordinates in mercator projection
      ///
      virtual const vector<double> &getXY();


      /// \brief get user position
      ///
      /// \return a vector<double> containing the coordinates of this
      /// point in the user's coordinate system.  
      ///
      /// This method works by converting the mercator coordinates of
      /// the point back to the user's coordinate system using the Proj.4
      /// cartographic projection library if the mercator coordinates 
      /// or user projection have changed since the last call. 

      virtual const vector<double> &getUserCoords();

      /// \brief set user position
      ///
      /// Sets the coordinates of the point in user's coordinate system,
      /// marking the mercator coordinates invalid.  Mercator coordinates
      /// will be updated by conversion with Proj.4 upon the next call to
      /// getXY().
      ///
      /// \param uPosition vector of doubles with user coordinates.

      virtual void setUserCoords(const vector<double> &uPosition)  ;

      virtual Point * Clone();
    private:
      void mercToUser();
      void userToMerc();
    };
  }
}
#endif // DF_PROJ_POINT_HPP
