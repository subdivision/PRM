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

#define NUM_OF_POINTS 10000
#define RADIUS 1
#define STEP_QUERIES 1000
#define END_RADIUS 5

using namespace std;

struct cPoint;

typedef pair<cPoint*, double> Edge;

struct CmpEdges
{
    bool operator()(const Edge lhs, const Edge rhs) const;
};

struct cPoint{
public:
    cPoint();
    cPoint(Point_2 p, Point_2 e, double r);
    Point_2 point;
    Point_3 point3;
    Point_2 endPoint;
    double rotation;
    double distanceToEnd;
    double distance = 0;
    bool visited = false;
    cPoint* last = nullptr;
    set<Edge, CmpEdges> edges;
};


struct CmpCPointsPtrs
{
    bool operator()(const cPoint* lhs, const cPoint* rhs) const;
};

class MyRodPathFinder {
    int checks=0;
    int edgesNum = 0;
    int legalCounter = 0;
    vector<Polygon_2> obstacles;
    FT rodLength;
    cPoint startCPoint, endCPoint;

    uniform_real_distribution<double> xUnif, yUnif, rUnif;
    std::default_random_engine re;

    map<Point_3, cPoint> cMap;
    set<cPoint*, CmpCPointsPtrs> queue;
    Tree tree;

    void setDistributions(FT rodLength, vector<Polygon_2>& obstacles);
    double cPointDistance(cPoint *a, cPoint *b);
    bool checkConnectCPoint(cPoint *a, cPoint *b, IQueryHandler& queryHandler);

    Point_2 endRodPoint(Point_2 a, double dir);
    void setRandomPoints(unsigned long n, IQueryHandler& queryHandler);


    bool findPath(IQueryHandler& queryHandler);
    void addEdge(cPoint *current, cPoint *temp);
    void addNeighbors(cPoint *current);
    bool connectCPoint(cPoint *current, IQueryHandler& queryHandler);

    vector<Path::PathMovement> fetchPath();

public:
    vector<Path::PathMovement> getPath(FT rodLength, Point_2 rodStartPoint, double rodStartRotation,
                                       Point_2 rodEndPoint, double rodEndRotation, vector<Polygon_2>& obstacles);
};


#endif //ROD_MYRODPATHFINDER_H
