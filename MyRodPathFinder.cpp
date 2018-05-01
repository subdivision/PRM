//
// Created by t-idkess on 09-Apr-18.
//

#include "MyRodPathFinder.h"
#include <queue>

cPoint::cPoint(Point_2 p, Point_2 e, double r) : point(p), endPoint(e), rotation(r) {}

cPoint::cPoint() {
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

    int runs[] = {1, 2, 5, 10, 20};
    double raduises[] = {8, 5, 3, 2, 1};
    for (int i; i<8;i++) {
        //run index denote whice points is already in queue.
        setRandomPoints(NUM_OF_POINTS * runs[i], queryHandler);
        this->RADIUS = raduises[i];
        while (findPath(queryHandler)) {
            if (checkPath(queryHandler))
                return fetchPath();

            cout << "verifing failed\n";
            runIndex++;
        }
        runIndex++;
    }

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
            PSet.insert(p);
            cMap.insert(std::pair<Point_2, cPoint>(p, cPoint(p, endRodPoint(p, d), d)));
            legalCounter++;
        }
    }
    cout << "number of legal positions " << legalCounter << endl;
}

bool MyRodPathFinder::findPath(IQueryHandler &queryHandler) {
    cPoint* refEndPt = &(this->endCPoint);

    priority_queue<cPoint *, vector<cPoint *>, function<bool(cPoint *, cPoint *)> >
            queue([&](cPoint *p1, cPoint *p2) -> bool {
        return (cPointDistance(p1, refEndPt) > cPointDistance(p2, refEndPt));
                });

    queue.push(&(this->startCPoint));
    int edges = 0;
    while (!queue.empty()) {
        cPoint *current = queue.top();
        queue.pop();
        if (checkConnectCPointWrapper(current, &this->endCPoint, queryHandler)) {
            this->endCPoint.last = current;
            cout << "path found! edges " << edges << endl;
            return true;
        }

        list<Point_set_Vertex_handle> L;
        list<Point_set_Vertex_handle>::const_iterator it;
        Circle_2 rc(current->point, RADIUS);

        PSet.range_search(rc, std::back_inserter(L));
        for (it = L.begin(); it != L.end(); it++) {
            auto range = cMap.equal_range((*it)->point());
            for (auto i = range.first; i != range.second; ++i) {
                cPoint *temp = &(i->second);
                if (checkConnectCPointWrapper(current, temp, queryHandler)) {
                    queue.push(temp);
                    temp->inQueue = runIndex;
                    temp->last = current;
                    edges++;
                }
            }
        }
    }
    cout << "path not found! edges " << edges << endl;
    return false;
}

bool MyRodPathFinder::checkConnectCPointWrapper(cPoint *a, cPoint *b, IQueryHandler &queryHandler) {
    if (b->inQueue == runIndex)
        return false;
    Edge edge(a, b);
    auto it = edges.find(edge);
    if (it != edges.end())
        return it->second > 0;

    bool ans = false;
    if (cPointDistance(a, b) < RADIUS)
        ans = checkConnectCPoint(a, b, queryHandler, STEP_QUERIES+a->cost+b->cost);
    int res = ans? STEP_QUERIES : -1;
    edges[edge] = res;
    edges[make_pair(b, a)] = res;

    return ans;
}

bool MyRodPathFinder::checkConnectCPoint(cPoint *a, cPoint *b, IQueryHandler &queryHandler, int queries) {
    Vector_2 pointsVector(a->point, b->point);
    double d = b->rotation - a->rotation;
    if (d > M_PI)
        d -= 2 * M_PI;
    else if (d < -M_PI)
        d += 2 * M_PI;

    double currentDir = a->rotation;
    double stepDir = d / queries;

    FT currentX = a->point.x();
    FT stepX = pointsVector.x() / queries;

    FT currentY = a->point.y();
    FT stepY = pointsVector.y() / queries;
    for (int i = 1; i < queries; i++) {
        currentDir += stepDir;
        currentX += stepX;
        currentY += stepY;

        if (!queryHandler.isLegalConfiguration({currentX, currentY}, currentDir))
            return false;
    }
    return true;
}

bool MyRodPathFinder::checkPath(IQueryHandler &queryHandler) {
    cPoint *temp = &this->endCPoint;
    cPoint *start = &this->startCPoint;
    bool res = true;
    while (temp != start) {
        Edge edge(temp, temp->last);
        auto it = edges.find(edge);
        if (it != edges.end() && it->second >= VERIFY_QUERIES)
        {
            temp = temp->last;
            continue;
        }

        if (!checkConnectCPoint(temp, temp->last, queryHandler, VERIFY_QUERIES)) {
            cout << "verfing failed between " << temp->point << " to " << temp->last->point << endl;

            edges[make_pair(temp, temp->last)] = -1;
            edges[make_pair(temp->last, temp)] = -1;
            temp->cost += 50;
            temp->last->cost += 50;

            res = false;
        } else
        {
            edges[make_pair(temp, temp->last)] = VERIFY_QUERIES;
            edges[make_pair(temp->last, temp)] = VERIFY_QUERIES;
        }
        temp = temp->last;

    }
    return res;
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

    for (Path::PathMovement &t:path)
        cout << t << endl;
    return path;
}

double MyRodPathFinder::cPointDistance(cPoint *a, cPoint *b) {
    double xdiff = CGAL::to_double(a->point.x()) - CGAL::to_double(b->point.x());
    double ydiff = CGAL::to_double(a->point.y()) - CGAL::to_double(b->point.y());
    double rdiff = (a->rotation-b->rotation);

    return xdiff*xdiff + ydiff*ydiff + rdiff*rdiff;
}

FT MyRodPathFinder::pointsDistance(Point_2 a_point, Point_2 b_point) {
    FT distance = (a_point.x() - b_point.x()) * (a_point.x() - b_point.x()) +
                  (a_point.y() - b_point.y()) * (a_point.y() - b_point.y());

    return sqrt(CGAL::to_double(distance));
}

Point_2 MyRodPathFinder::endRodPoint(Point_2 a, double dir) {
    Vector_2 a_dir = {cos(dir), sin(dir)};
    return a + (a_dir * rodLength);
}

