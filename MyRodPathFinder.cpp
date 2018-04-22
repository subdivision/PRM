//
// Created by t-idkess on 09-Apr-18.
//

#include "MyRodPathFinder.h"

vector<Path::PathMovement>
MyRodPathFinder::getPath(FT rodLength, Point_2 rodStartPoint, double rodStartRotation, Point_2 rodEndPoint,
                         double rodEndRotation, vector<Polygon_2> obstacles) {
    //todo:
    vector<Path::PathMovement> res;
    res.emplace_back(rodStartPoint, rodStartRotation, CGAL::CLOCKWISE);
    res.emplace_back(rodStartPoint, ((3 * CGAL_PI) / 2), CGAL::CLOCKWISE);
    res.emplace_back(Point_2(0, 1), (CGAL_PI / 2), CGAL::CLOCKWISE);
    res.emplace_back(Point_2(1, 1), 0, CGAL::CLOCKWISE);
    res.emplace_back(rodEndPoint, rodEndRotation, CGAL::COUNTERCLOCKWISE);
    return res;
}
