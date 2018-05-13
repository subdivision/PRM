#include <iostream>
#include <boost/timer.hpp>
#include "CGAL_defines.h"
#include "MyQueryHandler.h"
#include "NaiveQueryHandler.h"


using namespace std;

void getQueryCount(int &numberOfQueries) {
    cin >> numberOfQueries;
}

void getRodLength(FT &rodLength) {
    cin >> rodLength;
}

void getObstacles(vector<Polygon_2> &obstacles) {
    FT x, y;
    int numberOfPolygons, numberOfVertices;
    cin >> numberOfPolygons;
    while (numberOfPolygons--) {
        cin >> numberOfVertices;
        Polygon_2 p;
        while (numberOfVertices--) {
            cin >> x >> y;
            p.push_back({x, y});
        }
        if (p.is_clockwise_oriented()) p.reverse_orientation();
        obstacles.emplace_back(move(p));
    }
}

void handleQuery(int queryNum, FT& timerSum, FT& timerSum2, int &errors, IQueryHandler &queryHandler, IQueryHandler &tester) {
    double rotation;
    FT x, y;
    cin >> x >> y >> rotation;
    boost::timer timer;
    bool result = queryHandler.isLegalConfiguration({x, y}, rotation);
    timerSum += timer.elapsed();
    boost::timer timer2;
    bool verifier = tester.isLegalConfiguration({x, y}, rotation);
    timerSum2 += timer2.elapsed();

    if (result != verifier) {
        cout << "Failure in query: #" << queryNum << endl;
        errors++;
    }
}

int main() {
    ios::sync_with_stdio(false);
    FT rodLength = 7;
    vector<Polygon_2> obstacles;

    getRodLength(rodLength);
    getObstacles(obstacles);

    boost::timer timer;
    MyQueryHandler queryHandler(rodLength, obstacles);
    cout << "queryHandler creation took " << timer.elapsed() << " seconds" << endl;
    NaiveQueryHandler naiveTester(rodLength, obstacles);

    int numberOfQueries;
    getQueryCount(numberOfQueries);
    FT timerSum = 0;
    FT timerSum2 = 0;
    int errors = 0;
    cout << "start queries" << endl;
    for (auto q = 1; q <= numberOfQueries; ++q)
        handleQuery(q, timerSum, timerSum2, errors, queryHandler, naiveTester);

    cout << "Error rate: " << errors << "/" << numberOfQueries << endl;
    cout << "Mean time: " << CGAL::to_double(timerSum / numberOfQueries) << " secs" << endl;
    cout << "Naive time: " << CGAL::to_double(timerSum2 / numberOfQueries) << " secs" << endl;
    return 0;
}


