//
// Created by t-idkess on 09-Apr-18.
//

#ifndef ROD_MYRODPATHFINDER_H
#define ROD_MYRODPATHFINDER_H

#include <vector>
#include "CGAL_defines.h"
#include "Path.h"

using namespace std;

class MyRodPathFinder {
public:
    vector<Path::PathMovement>
    getPath(FT rodLength, Point_2 rodStartPoint, double rodStartRotation, Point_2 rodEndPoint, double rodEndRotation,
            vector<Polygon_2> obstacles);
};


#endif //ROD_MYRODPATHFINDER_H
