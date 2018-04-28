//
// Created by t-idkess on 09-Apr-18.
//

#include "MyRodPathFinder.h"

vector<Path::PathMovement> MyRodPathFinder::getPath(FT rodLength, Point_2 rodStartPoint,
                                                    double rodStartRotation, Point_2 rodEndPoint,
                                                    double rodEndRotation, vector<Polygon_2>& obstacles)
{
    cout << "starting\n";
    NaiveQueryHandler queryHandler(rodLength, obstacles);
    cout << "setDistributions\n";
    setDistributions(rodLength, obstacles);
    cout << "setRandomPoints\n";
    setRandomPoints(NUM_OF_POINTS, queryHandler);

    startCPoint = {rodStartPoint,rodStartRotation};
    endCPoint = {rodEndPoint, rodEndRotation};

    if(findPath())
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
        cout << "temp point " << temp.point << " direction " << temp.rotation << endl;
        if(queryHandler.isLegalConfiguration(temp.point, temp.rotation))
            cPoints.push_back(temp);
    }
    cout << "number of legal positions " << cPoints.size() << endl;
}
