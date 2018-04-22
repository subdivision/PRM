//
// Created by t-idkess on 08-Apr-18.
//

#include "NaiveQueryHandler.h"


NaiveQueryHandler::NaiveQueryHandler(const FT &rodLength, const vector<Polygon_2> &obstacles) : _rodLength(rodLength),
                                                                                                _obstacles(obstacles) {
}

bool NaiveQueryHandler::_isLegalConfiguration(const Point_2 &point, const Vector_2 &direction, const double rotation) {
    Segment_2 s = {point, point + (direction * _rodLength)};
    vector<Point_2> points;
    for (Polygon_2 &p:_obstacles) {
        if (CGAL::bounded_side_2(p.vertices_begin(), p.vertices_end(), s.source()) == CGAL::ON_BOUNDED_SIDE ||
            CGAL::bounded_side_2(p.vertices_begin(), p.vertices_end(), s.target()) == CGAL::ON_BOUNDED_SIDE)
            return false;
        points.clear();
        points.push_back(s.source());
        points.push_back(s.target());
        for (auto it = p.edges_begin(); it != p.edges_end(); ++it) {
            auto intersection = CGAL::intersection(*it, s);
            if (intersection) {
                if (const Point_2 *pp = boost::get<Point_2>(&*intersection)) {
                    points.push_back(*pp);
                }
            }
        }
        sort(points.begin(), points.end());
        for (auto i = 0; i < points.size() - 1; ++i) {
            if (points.at(i) == points.at(i + 1))continue;
            Point_2 pp = points.at(i) + (points.at(i + 1) - points.at(i)) / 2;
            if (CGAL::bounded_side_2(p.vertices_begin(), p.vertices_end(), pp) == CGAL::ON_BOUNDED_SIDE) return false;
        }
    }

    return true;
}
