//
// Created by t-idkess on 09-Apr-18.
//

#include "MyRodPathFinder.h"

vector<Path::PathMovement> MyRodPathFinder::getPath(FT rodLength, Point_2 rodStartPoint,
                                                    double rodStartRotation, Point_2 rodEndPoint,
                                                    double rodEndRotation, vector<Polygon_2>& obstacles)
{
    NaiveQueryHandler queryHandler(rodLength, obstacles);
    if(!queryHandler.isLegalConfiguration(rodStartPoint, rodStartRotation))
        throw "start point is not legal";
    if(!queryHandler.isLegalConfiguration(rodEndPoint, rodEndRotation))
        throw "end point is not legal";

    startCPoint = {rodStartPoint,rodStartRotation};
    endCPoint = {rodEndPoint, rodEndRotation};

    setDistributions(rodLength, obstacles);
    setRandomPoints(NUM_OF_POINTS, queryHandler);

    if(findPath(queryHandler))
        return fetchPath();

    throw "no path found";
}

void MyRodPathFinder::setDistributions(FT rodLength, vector<Polygon_2>& obstacles) {
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

    xUnif = uniform_real_distribution<double>(CGAL::to_double(mostLeft-d), CGAL::to_double(mostRight+d));
    yUnif = uniform_real_distribution<double>(CGAL::to_double(mostDown-d), CGAL::to_double(mostUp+d));
    rUnif = uniform_real_distribution<double>(0, 2*M_PI);
}

void MyRodPathFinder::setRandomPoints(unsigned long n, IQueryHandler& queryHandler) {
    re.seed(std::chrono::system_clock::now().time_since_epoch().count());
    list<Point_2> Lr;
    int counter=0;
    for(int i=1; i<=n; i++)
    {
        cPoint temp = {{xUnif(re),yUnif(re)},rUnif(re)};
        if(queryHandler.isLegalConfiguration(temp.point, temp.rotation)) {
            Lr.push_back(temp.point);
            cMap[temp.point] = temp;
            counter++;
        }
    }
    PSet.insert(Lr.begin(),Lr.end());
    cout << "number of legal positions " << counter << endl;
}

bool MyRodPathFinder::findPath(IQueryHandler& queryHandler) {
    list<cPoint*> queue;
    queue.push_back(&(this->startCPoint));
    int edges = 0;
    while(!queue.empty()) {
        cPoint* current = queue.front();
        current->visited = true;
        queue.pop_front();
        if(checkConnectCPoint(current,&this->endCPoint,queryHandler))
        {
            this->endCPoint.last = current;
            cout << "path found! edges " << edges <<endl;
            return true;
        }
        list<Point_set_Vertex_handle> L;
        list<Point_set_Vertex_handle>::const_iterator it;
        Circle_2 rc(current->point,RADIUS);

        PSet.range_search(rc, std::back_inserter(L));

        for (it=L.begin();it != L.end(); it++)
        {
            cPoint* temp = &(cMap[(*it)->point()]);
            if (checkConnectCPoint(current, temp, queryHandler))
            {
                queue.push_back(temp);
                edges++;
            }
        }
    }
    cout << "path not found! edges " << edges <<endl;
    return false;
}

vector<Path::PathMovement> MyRodPathFinder::fetchPath() {
    vector<Path::PathMovement> tempVector;
    cPoint *temp = &this->endCPoint;
    cPoint *start = &this->startCPoint;
    while(temp != start)
    {
        CGAL::Orientation orientation = CGAL::CLOCKWISE;
        double d =  temp->last->rotation - temp->rotation;
        if((d > 0 && d < M_PI) || d < -M_PI)
            tempVector.push_back({temp->point,temp->rotation, CGAL::CLOCKWISE});
        else
            tempVector.push_back({temp->point,temp->rotation, CGAL::COUNTERCLOCKWISE});

        temp = temp->last;
    }
    vector<Path::PathMovement> path;
    path.push_back({this->startCPoint.point,this->startCPoint.rotation,CGAL::CLOCKWISE});
    for(int i= static_cast<int>(tempVector.size() - 1); i >= 0; i--)
        path.push_back(tempVector[i]);

    for(Path::PathMovement& t:path)
        cout << t << endl;
    return path;
}

bool MyRodPathFinder::checkConnectCPoint(cPoint *a, cPoint *b, IQueryHandler &queryHandler) {
    if(b->inQueue)
        return false;
    if(cPointDistance(a,b) > RADIUS)
        return false;


    Vector_2 pointsVector(a->point,b->point);
    double d =  a->rotation - b->rotation;
    if((d > 0 && d < M_PI) || d < -M_PI)
    {
        d = b->rotation - a->rotation;
        if(d > M_PI)
            d -= 2*M_PI;
    } else {
        d = b->rotation - a->rotation;
        if(d<0)
            d += 2*M_PI;
    }

    for(int i=1; i<=STEP_QUERIES; i++)
    {
        double dir = a->rotation + d*i/STEP_QUERIES;
        if(dir >= 2*M_PI)
            dir -= 2*M_PI;
        else if(dir < 0)
            dir += 2*M_PI;

        Point_2 p(a->point.x() + pointsVector.x() * (i/STEP_QUERIES),
                  a->point.y() + pointsVector.y() * (i/STEP_QUERIES));

        if(!queryHandler.isLegalConfiguration(p, dir))
            return false;
    }

    b->inQueue = true;
    b->last = a;
    return true;


}

FT MyRodPathFinder::cPointDistance(cPoint *a, cPoint *b) {
    return pointsDistance(a->point,b->point) + pointsDistance(endRodPoint(a),endRodPoint(b));
}

FT MyRodPathFinder::pointsDistance(Point_2 a_point, Point_2 b_point)
{
    FT distance = (a_point.x() - b_point.x()) * (a_point.x() - b_point.x()) +
                  (a_point.y() - b_point.y()) * (a_point.y() - b_point.y());

    return sqrt(CGAL::to_double(distance));;
}

Point_2 MyRodPathFinder::endRodPoint(cPoint *a)
{
    return endRodPoint(a->point, a->rotation);
}

Point_2 MyRodPathFinder::endRodPoint(Point_2 a, double dir)
{
    Vector_2 a_dir = {cos(dir), sin(dir)};
    return a + (a_dir * rodLength);
}
