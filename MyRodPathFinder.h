//
// Created by t-idkess on 09-Apr-18.
//

#ifndef ROD_MYRODPATHFINDER_H
#define ROD_MYRODPATHFINDER_H

#include <vector>
#include <random>
#include <math.h>
#include <map>

#include "NaiveQueryHandler.h"

#include "CGAL_defines.h"
#include "Path.h"
#include "MyQueryHandler.h"

#define NUM_OF_POINTS 30000
#define RADIUS 0.5
//#define STEP_QUERIES 100

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
    const double STEP_QUERIES = 1000;
    vector<Polygon_2> obstacles;
    FT rodLength;
    cPoint startCPoint, endCPoint;
    uniform_real_distribution<double> xUnif, yUnif, rUnif;
    std::default_random_engine re;
    map<Point_2, cPoint> cMap;
    Point_set PSet;

    void setDistributions(FT rodLength, vector<Polygon_2>& obstacles);
    FT cPointDistance(cPoint *a, cPoint *b);
    bool checkConnectCPoint(cPoint *a, cPoint *b, IQueryHandler& queryHandler);
    Point_2 endRodPoint(cPoint *a);
    Point_2 endRodPoint(Point_2 a, double dir);

    void setRandomPoints(unsigned long n, IQueryHandler& queryHandler);

    bool findPath(IQueryHandler& queryHandler);
    vector<Path::PathMovement> fetchPath();

    FT pointsDistance(Point_2 a, Point_2 b);


public:
    vector<Path::PathMovement> getPath(FT rodLength, Point_2 rodStartPoint, double rodStartRotation,
                                       Point_2 rodEndPoint, double rodEndRotation, vector<Polygon_2>& obstacles);
};


#endif //ROD_MYRODPATHFINDER_H
