#include "DF_XY_Point.hpp"

#include <vector>
using namespace std;

namespace DFLib
{
    namespace XY
    {

        // class Point

        Point::Point(const vector<double> &aPosition)
            :myXY(aPosition)
        {
        }

        Point::Point(Point &right)
            :myXY(right.myXY)
        {
        }

        void Point::setXY(const vector<double> &aPosition)
        {
            myXY = aPosition;
        }

        const vector<double> & Point::getXY()
        {
            return(myXY);
        }

        Point * Point::Clone()
        {
            Point *retPoint;
            retPoint = new Point(*this);
            return retPoint;
        }
    }
}
