//
// Created by t-idkess on 09-Apr-18.
//

#include <vector>
#include <boost/timer.hpp>
#include "CGAL_defines.h"
#include "MyQueryHandler.h"
#include "MyRodPathFinder.h"

using namespace std;

void getRod(FT &rodLength, Point_2 &rodStartPoint, double &rodStartRotation, Point_2 &rodEndPoint, double &rodEndRotation) {
    FT x, y;
    cin >> rodLength >> x >> y >> rodStartRotation;
    rodStartPoint = {x, y};
    cin >> x >> y >> rodEndRotation;
    rodEndPoint = {x, y};
}

void getObstacles(vector<Polygon_2> &obstacles) {
    int numberOfObstacles, numberOfVertices;
    FT x, y;
    cin >> numberOfObstacles;
    while (numberOfObstacles--) {
        Polygon_2 p;
        cin >> numberOfVertices;
        while (numberOfVertices--) {
            cin >> x >> y;
            p.push_back({x, y});
        }
        obstacles.push_back(move(p));
    }
}

/*!
 * Notice that you can provide a path to a output file
 */
int main(int argc, char *argv[]) {
    /*!
     * This value describes how many checks should be made on the path before determining it as correct.
     * (The checks are randomize and might fail and not fail depending on the point along the path that were choosen).
     * If the checks take up to much time, reduce the value.
     */
    //freopen("tests/test1","r",stdin);
    const int pathVerifierNumberOfPoints = 1000;
    FT rodLength;
    Point_2 rodStartPoint, rodEndPoint;
    double rodStartRotation, rodEndRotation;
    vector<Polygon_2> obstacles;

    getRod(rodLength, rodStartPoint, rodStartRotation, rodEndPoint, rodEndRotation);
    getObstacles(obstacles);

    MyRodPathFinder myRodPathFinder;
    try {
        boost::timer timer;
        Path path = myRodPathFinder.getPath(rodLength, rodStartPoint, rodStartRotation, rodEndPoint, rodEndRotation,
                                            obstacles);
        cout << "Path creation time: " << timer.elapsed() << endl;
        bool verifierResults = path.verify(rodLength, rodStartPoint, rodStartRotation, rodEndPoint, rodEndRotation,
                                           obstacles, pathVerifierNumberOfPoints);
        cout << "Path verifying: " << (verifierResults ? "SUCCESS!" : "FAILED!") << endl;
        if (argc != 1) { // print out the resulting path.
            ofstream file;
            file.open(argv[1], ios_base::out | ios_base::trunc);
            if (!file) {
                cerr << "Couldn't open output file: " << argv[1];
            } else {
                file << obstacles.size() << endl;
                for (auto &obs:obstacles) {
                    file << obs.size();
                    for (auto it = obs.vertices_begin(); it != obs.vertices_end(); ++it) {
                        file << " " << CGAL::to_double(it->x()) << " " << CGAL::to_double(it->y());
                    }
                    file << endl;
                }
                file << CGAL::to_double(rodLength) << endl;
                file << path << endl;
            }
        }
    }
    catch (const char* c)
    {
        cout << "ERROR: " << c << endl;
        return 0;
    }
}

