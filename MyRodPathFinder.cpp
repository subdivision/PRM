//
// Created by t-idkess on 09-Apr-18.
//

#include "MyRodPathFinder.h"

vector<Path::PathMovement> MyRodPathFinder::getPath(FT rodLength, Point_2 rodStartPoint,
                                                    double rodStartRotation, Point_2 rodEndPoint,
                                                    double rodEndRotation, vector<Polygon_2>& obstacles)
{
    this->rodLength = rodLength;
    //NaiveQueryHandler queryHandler(rodLength, obstacles);
    MyQueryHandler queryHandler(rodLength, obstacles);

    if(!queryHandler.isLegalConfiguration(rodStartPoint, rodStartRotation))
        throw "start point is not legal";
    if(!queryHandler.isLegalConfiguration(rodEndPoint, rodEndRotation))
        throw "end point is not legal";

    startCPoint = {rodStartPoint, endRodPoint(rodStartPoint, rodStartRotation),rodStartRotation};
    endCPoint = {rodEndPoint, endRodPoint(rodEndPoint, rodEndRotation), rodEndRotation};

    setDistributions(rodLength, obstacles);

    int runs[] = {1, 10, 20, 50, 100, 200, 500, 1000};
    for(int run:runs)
    {
        runIndex++; //run index denote whice points is already in queue.
        setRandomPoints(NUM_OF_POINTS*run, queryHandler);

        if(findPath(queryHandler))
            return fetchPath();
    }

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
    for(int i=1; i<=n; i++)
    {
        Point_2 p = {xUnif(re),yUnif(re)};
        double d = rUnif(re);
        cPoint temp = {p,endRodPoint(p,d),d};
        if(queryHandler.isLegalConfiguration(temp.point, temp.rotation)) {
            PSet.insert(p);
            cMap[temp.point] = temp;
            legalCounter++;
        }
    }
    cout << "number of legal positions " << legalCounter << endl;
}

bool MyRodPathFinder::findPath(IQueryHandler& queryHandler) {
    list<cPoint*> queue;
    queue.push_back(&(this->startCPoint));
    int edges = 0;
    while(!queue.empty()) {
        cPoint* current = queue.front();
        queue.pop_front();
        if(checkConnectCPointWrapper(current,&this->endCPoint,queryHandler))
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
            if (checkConnectCPointWrapper(current, temp, queryHandler))
            {
                queue.push_back(temp);
                temp->inTree = true;
                temp->inQueue = runIndex;
                temp->last = current;
                edges++;
            }
        }
    }
    cout << "path not found! edges " << edges <<endl;
    return false;
}
bool MyRodPathFinder::checkConnectCPointWrapper(cPoint *a, cPoint *b, IQueryHandler &queryHandler) {
    if(b->inQueue == runIndex)
        return false;
    auto it = a->edges.find(b);
    if (it != a->edges.end())
        return it->second;
    if(b->inTree)
        return false;
    bool ans = false;
    if(cPointDistance(a,b) < RADIUS)
        ans = checkConnectCPoint(a, b, queryHandler);

    a->edges.insert ( std::pair<cPoint*,bool>(b,ans) );
    b->edges.insert ( std::pair<cPoint*,bool>(a,ans) );

    return ans;
}

bool MyRodPathFinder::checkConnectCPoint(cPoint *a, cPoint *b, IQueryHandler &queryHandler) {
    Vector_2 pointsVector(a->point,b->point);
    double d =  b->rotation - a->rotation;
    if( d > M_PI)
        d -= 2*M_PI;
    else if(d<-M_PI)
        d += 2*M_PI;

    double currentDir = a->rotation;
    double stepDir = d/STEP_QUERIES;

    FT currentX = a->point.x();
    FT stepX = pointsVector.x() /STEP_QUERIES;

    FT currentY = a->point.y();
    FT stepY = pointsVector.y() / STEP_QUERIES;
    for(int i=1; i<STEP_QUERIES; i++)
    {
        currentDir += stepDir;
        currentX += stepX;
        currentY += stepY;

        if(!queryHandler.isLegalConfiguration({currentX,currentY}, currentDir))
            return false;
    }
    return true;
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

FT MyRodPathFinder::cPointDistance(cPoint *a, cPoint *b) {
    return pointsDistance(a->point,b->point) + pointsDistance(a->endPoint,b->endPoint);
}

FT MyRodPathFinder::pointsDistance(Point_2 a_point, Point_2 b_point)
{
    FT distance = (a_point.x() - b_point.x()) * (a_point.x() - b_point.x()) +
                  (a_point.y() - b_point.y()) * (a_point.y() - b_point.y());

    return sqrt(CGAL::to_double(distance));;
}

Point_2 MyRodPathFinder::endRodPoint(Point_2 a, double dir)
{
    Vector_2 a_dir = {cos(dir), sin(dir)};
    return a + (a_dir * rodLength);
}

bool CmpCPointPtrs::operator()(const cPoint *lhs, const cPoint *rhs) const {
    if(lhs->point != rhs->point)
        return lhs->point < rhs->point;

    lhs->endPoint < rhs->endPoint;
}
