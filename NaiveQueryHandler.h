//
// Created by t-idkess on 08-Apr-18.
//

#ifndef INC_3_2_NAIVETESTER_H
#define INC_3_2_NAIVETESTER_H


#include "IQueryHandler.h"
#include "CGAL_defines.h"
#include <algorithm>
#include <vector>

class NaiveQueryHandler : public IQueryHandler {
private:
    vector<Polygon_2> _obstacles;
    FT _rodLength;
public:
    NaiveQueryHandler(const FT &rodLength, const vector<Polygon_2> &obstacles);

protected:
    bool _isLegalConfiguration(const Point_2 &point, const Vector_2 &direction, const double rotation) override;
};


#endif //INC_3_2_NAIVETESTER_H
