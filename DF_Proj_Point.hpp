#ifndef DF_PROJ_POINT_HPP
#define DF_PROJ_POINT_HPP

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
    class Point : public DFLib::Abstract::Point
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
      /// \param llPosition coordinates <em>in user's coordinates</em>
      Point(const vector<double> &uPosition,const vector<string> &projArgs);


      /// \brief Copy Constructor
      Point(const Point &right);

      /// \brief destructor
      ~Point();

      /// \brief assignment operator
      Point& operator=(const Point& rhs);

      /// \brief set user projection information
      void setUserProj(vector<string> &projArgs);

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
      /// This is just a wrapper for getLL as needed by the abstract
      /// interface

      virtual const vector<double> &getUserCoords();

      /// \brief set user position
      ///
      /// This is just a wrapper for setLL as needed by the abstract
      /// interface

      virtual void setUserCoords(const vector<double> &uPosition)  ;

      virtual Point * Clone();
    private:
      void mercToUser();
      void userToMerc();
    };
  }
}
#endif // DF_PROJ_POINT_HPP
