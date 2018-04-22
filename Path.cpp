//
// Created by t-idkess on 09-Apr-18.
//

#include "Path.h"
#include "NaiveQueryHandler.h"


Path::Path(vector<Path::PathMovement> path) : _path(std::move(path)) {}

bool
Path::verify(FT rodLength, Point_2 rodStartPoint, double rodStartRotation, Point_2 rodEndPoint, double rodEndRotation,
             vector<Polygon_2> obstacles, int pathVerifierNumberOfPoints) {
    if (_path.front().location != rodStartPoint || _path.front().rotation != rodStartRotation)
        return cout << "Failure: path start is incorrect" << endl, false;
    if (_path.back().location != rodEndPoint || _path.back().rotation != rodEndRotation)
        return cout << "Failure: path end is incorrect" << endl, false;

    NaiveQueryHandler tester(rodLength, obstacles);

    srand(time(NULL));
    vector<pair<int, double>> checks;
    while (pathVerifierNumberOfPoints--)
        checks.emplace_back(rand() % (_path.size() - 1), (double) rand() / RAND_MAX);
    sort(checks.begin(), checks.end());
    for (auto &pr: checks) {
        int index = pr.first;
        double partOf = pr.second;
        Point_2 p = _path[index].location + (_path[index + 1].location - _path[index].location) * partOf;
        double MOD = ((double) 2) * CGAL_PI;
        double r1 = fmod(_path[index].rotation, MOD), r2 = fmod(_path[index + 1].rotation, MOD), r3;
        if (r1 < 0) r1 += MOD;
        if (r2 < 0) r2 += MOD;
        if (_path[index + 1].orientation == CGAL::CLOCKWISE && r2 > r1) {
            r1 += MOD;
        } else if (_path[index + 1].orientation == CGAL::COUNTERCLOCKWISE && r1 > r2) {
            r2 += MOD;
        }
        r3 = fmod((r1 + ((r2 - r1) * partOf)), MOD);
        if (r3 < 0) r2 += MOD;
        if (r1 == r2) r3 = r1;
        if (!tester.isLegalConfiguration(p, r3))
            return
                    cout << "Failure: while moving from " << _path[index] << " to " << _path[index + 1]
                         << " (indexes in path: " << index << "->" << (index + 1) << "), at (point, rotation): (("
                         << p.x().to_double() << ", " << p.y().to_double() << "), "
                         << r3
                         << ")." << endl, false;
    }
    return true;

}

ostream &operator<<(ostream &os, const Path::PathMovement &dt) {
    os << "((" << dt.location.x().to_double() << ", " << dt.location.y().to_double() << "), " << dt.rotation << ", "
       << (dt.orientation == CGAL::CLOCKWISE ? "CLOCKWISE" : "COUNTERCLOCKWISE") << ")";
    return os;
}

ostream &operator<<(ostream &os, const Path &pt) {
    os << pt._path.size() << endl;
    for (auto &p: pt._path) {
        os << p.location.x().to_double() << " " << p.location.y().to_double() << " " << p.rotation << " "
           << (p.orientation == CGAL::CLOCKWISE ? "-1" : "1") << endl;
    }
    return os;
}
