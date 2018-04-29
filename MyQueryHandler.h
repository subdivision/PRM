#ifndef RODQUERY_MYQUERYHANDLER_H
#define RODQUERY_MYQUERYHANDLER_H

#include "IQueryHandler.h"

class polygon_split_observer : public CGAL::Arr_observer<Arrangement_2> {
    void after_split_face(Face_handle f1, Face_handle f2, bool) override;
};


class MyQueryHandler : public IQueryHandler {
public:
    MyQueryHandler(const FT &rodLength, const vector<Polygon_2> &obstacles);

protected:
    bool _isLegalConfiguration(const Point_2 &point, const Vector_2 &direction, const double rotation) override;

private:
    FT _rodLength;
    Arrangement_2 _arr;
    Landmarks_pl pl;

    void verticalDecomposition(Kernel &ker);

    void addVerticalSegment(Vertex_handle v, CGAL::Object obj, Kernel &ker);

    void addFrame();
};

#endif //RODQUERY_MYQUERYHANDLER_H
