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

#define NUM_OF_POINTS 1000
#define RADIUS 0.5
#define STEP_QUERIES 50

using namespace std;

class cPoint{
public:
    Point_2 point;
    double rotation;
    bool visited = false;
    bool inQueue = false;
    cPoint* last = nullptr;
};

class MyRodPathFinder {

    Point_2 rodStartPoint, rodEndPoint;
    double rodStartRotation, rodEndRotation;
    vector<Polygon_2> obstacles;
    FT rodLength;
    cPoint startCPoint, endCPoint;
    vector<cPoint> cPoints;
    uniform_real_distribution<double> xUnif, yUnif, rUnif;
    std::default_random_engine re;

    void setDistributions(FT rodLength, vector<Polygon_2>& obstacles);
    FT cPointDistance(cPoint *a, cPoint *b);
    bool checkConnectCPoint(cPoint *a, cPoint *b, IQueryHandler& queryHandler);

    void setRandomPoints(unsigned long n, IQueryHandler& queryHandler);

    bool findPath(IQueryHandler& queryHandler);
    vector<Path::PathMovement> fetchPath();

    FT pointsDistance(Point_2 a, Point_2 b);


public:
    vector<Path::PathMovement> getPath(FT rodLength, Point_2 rodStartPoint, double rodStartRotation,
                                       Point_2 rodEndPoint, double rodEndRotation, vector<Polygon_2>& obstacles);
};


#endif //ROD_MYRODPATHFINDER_H
