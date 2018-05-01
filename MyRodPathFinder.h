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

#define NUM_OF_POINTS 5000
//#define RADIUS 2
#define STEP_QUERIES 5
#define VERIFY_QUERIES 1000

using namespace std;

struct cPoint;

struct cPoint{
public:
    cPoint();
    cPoint(Point_2 p, Point_2 e, double r);
    Point_2 point;
    Point_2 endPoint;
    double rotation;
    int inQueue = 0;
    int cost = 0;
    cPoint* last = nullptr;
};



class MyRodPathFinder {
    double RADIUS = 2;
    typedef pair<cPoint*, cPoint*> Edge;
    int legalCounter = 0;
    int runIndex = 1; //first run
    vector<Polygon_2> obstacles;
    FT rodLength;
    cPoint startCPoint, endCPoint;
    uniform_real_distribution<double> xUnif, yUnif, rUnif;
    std::default_random_engine re;
    multimap<Point_2, cPoint> cMap;
    Point_set PSet;
    map<Edge, int> edges;

    void setDistributions(FT rodLength, vector<Polygon_2>& obstacles);
    double cPointDistance(cPoint *a, cPoint *b);
    bool checkConnectCPoint(cPoint *a, cPoint *b, IQueryHandler& queryHandler, int queries);
    bool checkConnectCPointWrapper(cPoint *a, cPoint *b, IQueryHandler &queryHandler);
    Point_2 endRodPoint(Point_2 a, double dir);

    void setRandomPoints(unsigned long n, IQueryHandler& queryHandler);


    bool findPath(IQueryHandler& queryHandler);
    bool checkPath(IQueryHandler& queryHandler);
    vector<Path::PathMovement> fetchPath();

    FT pointsDistance(Point_2 a, Point_2 b);


public:
    vector<Path::PathMovement> getPath(FT rodLength, Point_2 rodStartPoint, double rodStartRotation,
                                       Point_2 rodEndPoint, double rodEndRotation, vector<Polygon_2>& obstacles);
};


#endif //ROD_MYRODPATHFINDER_H
