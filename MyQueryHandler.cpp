//
// Created by t-idkess on 08-Apr-18.
//

#include "MyQueryHandler.h"

using namespace std;

//-----------------------------------------------------------------------------
MyQueryHandler::MyQueryHandler( const FT&                rodLength,
                                const vector<Polygon_2>& obstacles):
_rodLength( rodLength ),
_verifier( rodLength, obstacles )
{}

//-----------------------------------------------------------------------------
bool
MyQueryHandler::_isLegalConfiguration( const Point_2&  point,
                                       const Vector_2& direction,
                                       const double    rotation )
{
  Point_2 endPt = point + (direction * _rodLength);
  return _verifier.isFree( point, endPt );
}
//============================= END OF FILE ===================================
