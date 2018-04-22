//
// Created by t-idkess on 08-Apr-18.
//

#ifndef RODQUERY_RODQUERYINTERFACE_H
#define RODQUERY_RODQUERYINTERFACE_H

#include <vector>
#include "CGAL_defines.h"

using namespace std;

class IQueryHandler {
public:
    /*!
   * @param point: the current point of the robot.
   * @param rotation: [0,2PI] where 0 is right and PI/2 is up. (counterClockwise)
   * @return
   * true if the configuration does not violate any obstacles. false otherwise.
   */
    bool isLegalConfiguration(const Point_2 &point, const double rotation) {
        return _isLegalConfiguration(point, {cos(rotation), sin(rotation)}, rotation);
    }

protected:
    /*!
     * @param point: the current point of the robot.
     * @param direction: direction of the robot (unit vector of size 1!).
     * @param rotation: the rotation value.
     * @return
     * true if the configuration does not violate any obstacles. false otherwise.
     */
    virtual bool _isLegalConfiguration(const Point_2 &point, const Vector_2 &direction, const double rotation)=0;
};

#endif //RODQUERY_RODQUERYINTERFACE_H
