//
// Created by t-idkess on 09-Apr-18.
//

#ifndef ROD_MYRODPATHFINDER_H
#define ROD_MYRODPATHFINDER_H

#include <vector>
#include <random>
#include <math.h>

#include "NaiveQueryHandler.h"

#include "CGAL_defines.h"
#include "Path.h"
#include "MyQueryHandler.h"

#define NUM_OF_POINTS 100000
#define RADIUS 0.5

using namespace std;

class cPoint{
public:
    Point_2 point;
    double rotation;
    bool visited = false;
    cPoint* last;
    vector<cPoint*> neighbors;
};

class MyRodPathFinder {

    Point_2 rodStartPoint, rodEndPoint;
    double rodStartRotation, rodEndRotation;
    vector<Polygon_2> obstacles;
    FT rodLength;
    cPoint start, end;
    vector<cPoint> cPoints;
    uniform_real_distribution<double> xUnif, yUnif, rUnif;
    std::default_random_engine re;

    void setDistributions(FT rodLength, vector<Polygon_2>& obstacles);
    cPoint getRandomPoint();
    FT cPointDistance(cPoint a, cPoint b);
    bool checkConnectCPoint(cPoint a, cPoint b, MyQueryHandler& queryHandler);

    void setRandomPoints(unsigned long n, MyQueryHandler& queryHandler);
    void addStartAndEndPoints();
    void connectPoints(FT r, MyQueryHandler queryHandler);


    bool findPath();
    vector<Path::PathMovement> fetchPath();


public:
    vector<Path::PathMovement> getPath(FT rodLength, Point_2 rodStartPoint, double rodStartRotation,
                                       Point_2 rodEndPoint, double rodEndRotation, vector<Polygon_2>& obstacles);
};


#endif //ROD_MYRODPATHFINDER_H
