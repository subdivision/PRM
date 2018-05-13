//
// Created by t-idkess on 09-Apr-18.
//

#include "MyRodPathFinder.h"
#include <queue>

cPoint::cPoint(Point_2 p, Point_2 e, double r) :
        point(p), point3(p[0],p[1],r), endPoint(e), rotation(r) {}

cPoint::cPoint() {
}

bool CmpCPointsPtrs::operator()(const cPoint *lhs, const cPoint *rhs) const {

    if(lhs->distanceToEnd != rhs->distanceToEnd)
        return lhs->distanceToEnd < rhs->distanceToEnd;

    if(lhs->distance != rhs->distance)
        return lhs->distance < rhs->distance;

    return lhs->point3 < rhs->point3;
}

bool CmpEdges::operator()(const Edge lhs, const Edge rhs) const {
    if(lhs.second != rhs.second)
        return lhs.second < rhs.second;

    return lhs.first->point3 < rhs.first->point3;
}

vector<Path::PathMovement> MyRodPathFinder::getPath(FT rodLength, Point_2 rodStartPoint,
                                                    double rodStartRotation, Point_2 rodEndPoint,
                                                    double rodEndRotation, vector<Polygon_2> &obstacles) {
    this->rodLength = rodLength;
    //NaiveQueryHandler queryHandler(rodLength, obstacles);
    MyQueryHandler queryHandler(rodLength, obstacles);

    if (!queryHandler.isLegalConfiguration(rodStartPoint, rodStartRotation))
        throw "start point is not legal";
    if (!queryHandler.isLegalConfiguration(rodEndPoint, rodEndRotation))
        throw "end point is not legal";

    re.seed(std::chrono::system_clock::now().time_since_epoch().count());

    startCPoint = {rodStartPoint, endRodPoint(rodStartPoint, rodStartRotation), rodStartRotation};
    endCPoint = {rodEndPoint, endRodPoint(rodEndPoint, rodEndRotation), rodEndRotation};

    setDistributions(rodLength, obstacles);

    setRandomPoints(NUM_OF_POINTS, queryHandler);
    if(findPath(queryHandler))
        return fetchPath();

    throw "no path found";
}

void MyRodPathFinder::setDistributions(FT rodLength, vector<Polygon_2> &obstacles) {
    Polygon_set_2 freeSpace;
    freeSpace.join(obstacles.begin(), obstacles.end());
    Arrangement_2 arr = freeSpace.arrangement();
    FT d = rodLength;

    Arr_VrtxCIter iVrtx = arr.vertices_begin();
    FT mostLeft = startCPoint.point.x() < endCPoint.point.x() ? startCPoint.point.x() : endCPoint.point.x();
    FT mostRight = startCPoint.point.x() < endCPoint.point.x() ? endCPoint.point.x() : startCPoint.point.x();
    FT mostUp = startCPoint.point.y() < endCPoint.point.y() ? endCPoint.point.y() : startCPoint.point.y();
    FT mostDown = startCPoint.point.y() < endCPoint.point.y() ? startCPoint.point.y() : endCPoint.point.y();

    for (; iVrtx != arr.vertices_end(); ++iVrtx) {
        if (iVrtx->point().x() < mostLeft)
            mostLeft = iVrtx->point().x();
        if (iVrtx->point().x() > mostRight)
            mostRight = iVrtx->point().x();
        if (iVrtx->point().y() < mostDown)
            mostDown = iVrtx->point().y();
        if (iVrtx->point().y() > mostUp)
            mostUp = iVrtx->point().y();
    }

    xUnif = uniform_real_distribution<double>(CGAL::to_double(mostLeft - d/10), CGAL::to_double(mostRight + d/10));
    yUnif = uniform_real_distribution<double>(CGAL::to_double(mostDown - d/10), CGAL::to_double(mostUp + d/10));
    rUnif = uniform_real_distribution<double>(0, 2 * M_PI);
}

void MyRodPathFinder::setRandomPoints(unsigned long n, IQueryHandler &queryHandler) {
    for (int i = 1; i <= n; i++) {
        Point_2 p = {xUnif(re), yUnif(re)};
        double d = rUnif(re);
        if (queryHandler.isLegalConfiguration(p, d)) {
            Point_3 p3(p[0],p[1],d);
            tree.insert(p3);
            cPoint cp = cPoint(p, endRodPoint(p, d), d);
            cp.distanceToEnd = cPointDistance(&cp, &this->endCPoint);
            cMap.insert(std::pair<Point_3, cPoint>(p3, cp));
            legalCounter++;
        }
    }
    cout << "number of legal positions " << legalCounter << endl;
}

void MyRodPathFinder::addEdge(cPoint *current, cPoint *temp) {
    if(temp->visited)
        return;

    this->queue.insert(temp);

    double newDistance = current->distance + cPointDistance(current, temp);
    temp->edges.insert(Edge(current, newDistance));
}

bool MyRodPathFinder::connectCPoint(cPoint *current, IQueryHandler& queryHandler) {
    //cout << "\t\tedges " << current->edges.size() << endl;
    for(auto edgeIt = current->edges.begin(); edgeIt != current->edges.end();)
    {
        checks++;
        //cout << "\t\t check edge from " << edgeIt->first->point << " dis " << edgeIt->second;
        if(checkConnectCPoint(edgeIt->first,current,queryHandler))
        {
            current->last = edgeIt->first;
            current->distance = edgeIt->second;
            current->visited = true;
            //cout << " passed!\n";
            return true;
        }

        edgeIt = current->edges.erase(edgeIt);
        //cout << "\n";
    }
    //cout << " failed!\n";
    return false;
}

void MyRodPathFinder::addNeighbors(cPoint *current) {
    list<Point_3> L;
    list<Point_3>::const_iterator it;
    Fuzzy_sphere rc(current->point3, RADIUS);


    tree.search(std::back_inserter(L), rc);

    edgesNum += L.size();
    //cout << "adding " << L.size() << " edges\n";
    for (it = L.begin(); it != L.end(); it++) {
        cPoint* temp = &(cMap[*it]);
        addEdge(current, temp);
    }
}

bool MyRodPathFinder::findPath(IQueryHandler &queryHandler) {
    cPoint* refEndPt = &(this->endCPoint);

    addNeighbors(&(this->startCPoint));

    while (!queue.empty()) {
        cPoint *current = *(queue.begin());
        //cout << "current " << current->point << " " << current->distanceToEnd;
        if(!connectCPoint(current,queryHandler))
        {
            queue.erase(current);
            continue;
        }
        if(current->distanceToEnd < END_RADIUS) {
            checks++;
            if (checkConnectCPoint(current, &this->endCPoint, queryHandler)) {
                this->endCPoint.last = current;
                cout << "found, checks " << checks << " edgse " << edgesNum << endl;
                return true;
            }
        }

        addNeighbors(current);

        queue.erase(current);
    }
    cout << "not found, checks " << checks << " edgse " << edgesNum << endl;
    return false;
}

bool MyRodPathFinder::checkConnectCPoint(cPoint *a, cPoint *b, IQueryHandler &queryHandler) {
    Vector_2 pointsVector(a->point, b->point);
    double d = b->rotation - a->rotation;
    if (d > M_PI)
        d -= 2 * M_PI;
    else if (d < -M_PI)
        d += 2 * M_PI;

    double currentDir = a->rotation;
    double stepDir = d / STEP_QUERIES;

    FT currentX = a->point.x();
    FT stepX = pointsVector.x() / STEP_QUERIES;

    FT currentY = a->point.y();
    FT stepY = pointsVector.y() / STEP_QUERIES;
    for (int i = 1; i < STEP_QUERIES; i++) {
        currentDir += stepDir;
        currentX += stepX;
        currentY += stepY;

        if (!queryHandler.isLegalConfiguration({currentX, currentY}, currentDir))
            return false;
    }
    return true;
}

vector<Path::PathMovement> MyRodPathFinder::fetchPath() {
    vector<Path::PathMovement> tempVector;
    cPoint *temp = &this->endCPoint;
    cPoint *start = &this->startCPoint;
    while (temp != start) {
        CGAL::Orientation orientation = CGAL::CLOCKWISE;
        double d = temp->last->rotation - temp->rotation;
        if ((d > 0 && d < M_PI) || d < -M_PI)
            tempVector.push_back({temp->point, temp->rotation, CGAL::CLOCKWISE});
        else
            tempVector.push_back({temp->point, temp->rotation, CGAL::COUNTERCLOCKWISE});

        temp = temp->last;
    }
    vector<Path::PathMovement> path;
    path.push_back({this->startCPoint.point, this->startCPoint.rotation, CGAL::CLOCKWISE});
    for (int i = static_cast<int>(tempVector.size() - 1); i >= 0; i--)
        path.push_back(tempVector[i]);

    cout << "path lenght " << path.size() << endl;
    return path;
}

double MyRodPathFinder::cPointDistance(cPoint *a, cPoint *b) {
    double xdiff = CGAL::to_double(a->point.x()) - CGAL::to_double(b->point.x());
    double ydiff = CGAL::to_double(a->point.y()) - CGAL::to_double(b->point.y());
    double rdiff = a->rotation - b->rotation;

    return sqrt(xdiff*xdiff + ydiff*ydiff + rdiff*rdiff);
}

Point_2 MyRodPathFinder::endRodPoint(Point_2 a, double dir) {
    Vector_2 a_dir = {cos(dir), sin(dir)};
    return a + (a_dir * rodLength);
}







