//
// Created by t-idkess on 09-Apr-18.
//

#ifndef ROD_PATH_H
#define ROD_PATH_H


#include "CGAL_defines.h"
#include <vector>
#include <utility>

using namespace std;

class Path {
public:
    /*!
    * This class should be used to represent a translation and rotation inside the path.
    */
    class PathMovement {
    public:
        /*!
         * The location of the rod at the end of the movement.
         */
        Point_2 location;
        /*!
         * The rotation of the rod at the end of the movement. (should be [0,2PI] where 0 is right and PI/2 is up)
         */
        double rotation;
        /*!
         * The rotation direction of the robot ( clockwise or counterclockwise).
         */
        CGAL::Orientation orientation;

        PathMovement() : location(0, 0), rotation(0), orientation(CGAL::CLOCKWISE) {}

        PathMovement(Point_2 location, double rotation, CGAL::Orientation orientation1) :
                location(std::move(location)), rotation(rotation), orientation(orientation1) {}

        friend ostream &operator<<(ostream &os, const PathMovement &dt);
    };

    Path(vector<PathMovement> path);

    bool
    verify(FT rodLength, Point_2 rodStartPoint, double rodStartRotation, Point_2 rodEndPoint, double rodEndRotation,
           vector<Polygon_2> obstacles, int pathVerifierNumberOfPoints);

    friend ostream &operator<<(ostream &os, const Path &pt);

private:
    vector<PathMovement> _path;
};


#endif //ROD_PATH_H
