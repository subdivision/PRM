//
// Created by t-idkess on 09-Apr-18.
//

#include "MyRodPathFinder.h"
#include <queue>

cPoint::cPoint(Point_2 p, Point_2 e, double r) :
        point(p), point3(p[0],p[1],r), endPoint(e), rotation(r) {}

cPoint::cPoint() {
}

bool CmpEdges::operator()(const Edge lhs, const Edge rhs) const {
    if(lhs.distance != rhs.distance)
        return lhs.distance > rhs.distance;

    if(lhs.from->point3 != lhs.from->point3)
        return lhs.from->point3 > rhs.from->point3;

    return lhs.to->point3 > rhs.to->point3;
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

    setRandomPoints(queryHandler);
    if(findPath(queryHandler))
        return fetchPath();

    throw "no path found";
}

void MyRodPathFinder::exportCPoints()
{
    ofstream file;
    file.open("points", ios_base::out | ios_base::trunc);
    for (auto it=cMap.begin(); it!=cMap.end(); ++it)
    {
        file << it->second.state << " ";
        file << CGAL::to_double(it->second.point.x()) << " " << CGAL::to_double(it->second.point.y()) <<
             " " << it->second.rotation << endl;
    }
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
    double addition = max(1.0, CGAL::to_double(d));
    mostLeft -= addition;
    mostRight += addition;
    mostDown -= addition;
    mostUp += addition;

    int boxSize = int(CGAL::to_double((mostRight-mostLeft) * (mostUp-mostDown)));
    numberOfRandomConfiguration = NUM_OF_POINTS_PER_SQUARE*boxSize;

    xUnif = uniform_real_distribution<double>(CGAL::to_double(mostLeft), CGAL::to_double(mostRight));
    yUnif = uniform_real_distribution<double>(CGAL::to_double(mostDown), CGAL::to_double(mostUp));
    rUnif = uniform_real_distribution<double>(0, 2 * M_PI);
}

void MyRodPathFinder::setRandomPoints(IQueryHandler &queryHandler) {
    for (int i = 1; i <= numberOfRandomConfiguration; i++) {
        Point_2 p = {xUnif(re), yUnif(re)};
        double d = rUnif(re);
        if (queryHandler.isLegalConfiguration(p, d)) {
            Point_3 p3(p[0],p[1],d);
            tree.insert(p3);
            cPoint cp = cPoint(p, endRodPoint(p, d), d);
            cp.distanceToEnd = cPointDistance(&cp, &this->endCPoint);
            cMap.insert(std::pair<Point_3, cPoint>(p3, cp));
            legalConfiguration++;
        }
    }
}

void MyRodPathFinder::addEdge(cPoint *current, cPoint *temp) {
    if(temp->visited)
        return;
    if(temp->state == 0)
    {
        temp->state = 1;
        discoveredConfigurations++;
    }

    double newDistance = current->distance + cPointDistance(current, temp) + temp->distanceToEnd;

    this->queue.push({current, temp, newDistance});

}

bool MyRodPathFinder::checkEdge(Edge& edge, IQueryHandler& queryHandler) {
    if(edge.to->visited)
        return false;

    checkedEdges++;
    if(checkConnectCPoint(edge.from,edge.to,queryHandler))
    {
        edge.to->last = edge.from;
        edge.to->distance = edge.distance - edge.from->distanceToEnd;
        edge.to->visited = true;
        edge.to->state = 2;
        processedConfigurations++;
        return true;
    }
    forbiddenEdges++;

    return false;
}

void MyRodPathFinder::addNeighbors(cPoint *current) {
    list<Point_3> L;
    list<Point_3>::const_iterator it;
    Fuzzy_sphere rc(current->point3, RADIUS);


    tree.search(std::back_inserter(L), rc);

    numOfEdges += L.size();
    for (it = L.begin(); it != L.end(); it++) {
        cPoint* temp = &(cMap[*it]);
        addEdge(current, temp);
    }
}

bool MyRodPathFinder::findPath(IQueryHandler &queryHandler) {
    cPoint* refEndPt = &(this->endCPoint);

    addNeighbors(&(this->startCPoint));

    while (!queue.empty()) {
        if(queue.size() > queueMaxSize)
            queueMaxSize = queue.size();
        Edge currentEdge = queue.top();
        queue.pop();
        if(!checkEdge(currentEdge,queryHandler))
            continue;

        processedConfigurations++;
        cPoint* currentCpoint = currentEdge.to;
        if(currentCpoint->distanceToEnd < END_RADIUS) {
            endEdgesChecked++;
            if (checkConnectCPoint(currentCpoint, refEndPt, queryHandler)) {
                this->endCPoint.last = currentCpoint;
                return true;
            }
        }

        addNeighbors(currentCpoint);

    }
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
        temp->state=3;
        temp = temp->last;
    }
    vector<Path::PathMovement> path;
    path.push_back({this->startCPoint.point, this->startCPoint.rotation, CGAL::CLOCKWISE});
    for (int i = static_cast<int>(tempVector.size() - 1); i >= 0; i--)
        path.push_back(tempVector[i]);

    pathLength = path.size();
    return path;
}

double MyRodPathFinder::cPointDistance(cPoint *a, cPoint *b) {
    double xdiff = CGAL::to_double(a->point.x()) - CGAL::to_double(b->point.x());
    double ydiff = CGAL::to_double(a->point.y()) - CGAL::to_double(b->point.y());

    double xdiff2 = CGAL::to_double(a->endPoint.x()) - CGAL::to_double(b->endPoint.x());
    double ydiff2 = CGAL::to_double(a->endPoint.y()) - CGAL::to_double(b->endPoint.y());

 //   double rdiff = a->rotation - b->rotation;

    return xdiff*xdiff + ydiff*ydiff + xdiff2*xdiff2 + ydiff2*ydiff2;
}

Point_2 MyRodPathFinder::endRodPoint(Point_2 a, double dir) {
    Vector_2 a_dir = {cos(dir), sin(dir)};
    return a + (a_dir * rodLength);
}

void MyRodPathFinder::printStatistics() {
    cout << "\n***** statistics *****\n";
    cout << "Configurations:\n";
    cout << "Number of random configurations - " << numberOfRandomConfiguration << endl;
    cout << "Number of legal configuration - " << legalConfiguration << endl;
    cout << "Number of discovered configurations - " << discoveredConfigurations << endl;
    cout << "Number of processed configuration - " << processedConfigurations << endl;

    cout << "Edges:\n";
    cout << "Number of edges discovered - " << numOfEdges << endl;
    cout << "Number of edges checked - " << checkedEdges << endl;
    cout << "Nubmer of forbidden edges - " << forbiddenEdges << endl;
    cout << "Nubmer of edges to the end checked - " << endEdgesChecked << endl;
    cout << "\nA* queue max size " << queueMaxSize << endl;

    if(pathLength>0)
        cout << "path length " << pathLength << endl;

}







