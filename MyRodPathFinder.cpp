//
// Created by t-idkess on 09-Apr-18.
//

#include "MyRodPathFinder.h"

vector<Path::PathMovement> MyRodPathFinder::getPath(FT rodLength, Point_2 rodStartPoint,
                                                    double rodStartRotation, Point_2 rodEndPoint,
                                                    double rodEndRotation, vector<Polygon_2>& obstacles)
{
    cout << "starting\n";
    MyQueryHandler queryHandler(rodLength, obstacles);
    cout << "setDistributions\n";
    setDistributions(rodLength, obstacles);
    cout << "setRandomPoints\n";
    setRandomPoints(NUM_OF_POINTS, queryHandler);

    startCPoint = {rodStartPoint,rodStartRotation};
    endCPoint = {rodEndPoint, rodEndRotation};

    if(findPath(queryHandler))
        return fetchPath();

    throw "no path found";
}

void MyRodPathFinder::setDistributions(FT rodLength, vector<Polygon_2>& obstacles) {
    Polygon_set_2 freeSpace;
    freeSpace.join(obstacles.begin(), obstacles.end());
    Arrangement_2 arr = freeSpace.arrangement();
    FT d = rodLength * 1.5;

    Arr_VrtxCIter iVrtx = arr.vertices_begin();
    FT mostLeft = iVrtx->point().x();
    FT mostRight = iVrtx->point().x();
    FT mostUp = iVrtx->point().y();
    FT mostDown = iVrtx->point().y();

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
    for(int i=1; i<=n; i++)
    {
        cPoint temp = {{xUnif(re),yUnif(re)},rUnif(re)};
        cout << temp.point << " " << temp.rotation << "\n";
        if(queryHandler.isLegalConfiguration(temp.point, temp.rotation))
            cPoints.push_back(temp);
    }
    cout << "number of legal positions " << cPoints.size() << endl;
}

bool MyRodPathFinder::findPath(IQueryHandler& queryHandler) {
    list<cPoint*> queue;
    queue.push_back(&(this->startCPoint));
    int edges = 0;
    while(!queue.empty()) {
        cPoint* temp = queue.front();
        temp->visited = true;
        queue.pop_front();
        if(checkConnectCPoint(temp,&this->endCPoint,queryHandler))
        {
            this->endCPoint.last = temp;
            cout << "path found! edges " << edges <<endl;
            return true;
        }
        for(int i=0; i<cPoints.size(); i++) {
            if (checkConnectCPoint(temp, &(cPoints[i]), queryHandler))
            {
                queue.push_back(&(cPoints[i]));
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
    if(b->visited || b->inQueue)
        return false;
    if(cPointDistance(a,b) > RADIUS)
        return false;

    /*if(abs(a->rotation - b->rotation) > M_PI)
    {

    }*/

    b->inQueue = true;
    b->last = a;
    return true;


}

FT MyRodPathFinder::cPointDistance(cPoint *a, cPoint *b) {
    Vector_2 a_dir = {cos(a->rotation), sin(a->rotation)};
    Vector_2 b_dir = {cos(b->rotation), sin(b->rotation)};
    Point_2 a_end = a->point + (a_dir * rodLength);
    Point_2 b_end = b->point + (b_dir * rodLength);
    return pointsDistance(a->point,b->point) + pointsDistance(a_end,b_end);
}

FT MyRodPathFinder::pointsDistance(Point_2 a_point, Point_2 b_point)
{
    FT distance = (a_point.x() - b_point.x()) * (a_point.x() - b_point.x()) +
                  (a_point.y() - b_point.y()) * (a_point.y() - b_point.y());

    return sqrt(CGAL::to_double(distance));;
}
