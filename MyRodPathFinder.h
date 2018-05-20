//
// Created by t-idkess on 09-Apr-18.
//

#ifndef ROD_MYRODPATHFINDER_H
#define ROD_MYRODPATHFINDER_H

#include <vector>
#include <random>
#include <math.h>
#include <map>
#include <queue>

#include "NaiveQueryHandler.h"

#include "CGAL_defines.h"
#include "Path.h"
#include "MyQueryHandler.h"

#define NUM_OF_POINTS_PER_SQUARE 100
#define RADIUS 1.5
#define STEP_QUERIES 100
#define END_RADIUS 10

using namespace std;

struct cPoint;
struct Edge;

struct CmpEdges
{
    bool operator()(const Edge lhs, const Edge rhs) const;
};

struct cPoint{
public:
    cPoint();
    cPoint(Point_2 p, double r);
    Point_2 point;
    Point_3 point3;
    double startX, startY;
    double rotation;
    double heuristic=-1;
    double distance = 0;
    bool visited = false;
    int state = 0;
    cPoint* last = nullptr;
};

struct Edge{
        cPoint *from, *to;
        double distance;
};

class MyRodPathFinder {
    //statistics:
    int numberOfRandomConfiguration = 0;
    int legalConfiguration = 0;
    int discoveredConfigurations = 0;
    int processedConfigurations = 0;

    int checkedEdges = 0;
    int forbiddenEdges = 0;
    int numOfEdges = 0;
    int queueMaxSize = 0;
    int endEdgesChecked = 0;

    int pathLength = 0;


    vector<Polygon_2> obstacles;
    FT rodLength;
    cPoint startCPoint, endCPoint;

    uniform_real_distribution<double> xUnif, yUnif, rUnif;
    std::default_random_engine re;

    map<Point_3, cPoint> cMap;
    priority_queue<Edge, vector<Edge>, CmpEdges> queue;
    Tree tree;

    void setDistributions(FT rodLength, vector<Polygon_2>& obstacles);
    double cPointDistance(cPoint *a, cPoint *b);
    double heuristic(cPoint *cp);
    bool checkConnectCPoint(cPoint *a, cPoint *b, IQueryHandler& queryHandler);

    Point_2 endRodPoint(Point_2 a, double dir);
    void setRandomPoints(IQueryHandler& queryHandler);


    bool findPath(IQueryHandler& queryHandler);
    void addEdge(cPoint *current, cPoint *temp);
    void addNeighbors(cPoint *current);
    bool checkEdge(Edge& edge, IQueryHandler& queryHandler);

    vector<Path::PathMovement> fetchPath();


public:
    vector<Path::PathMovement> getPath(FT rodLength, Point_2 rodStartPoint, double rodStartRotation,
                                       Point_2 rodEndPoint, double rodEndRotation, vector<Polygon_2>& obstacles);
    void printStatistics();
    void exportCPoints();

};


#endif //ROD_MYRODPATHFINDER_H
