#ifndef SEGEMENT_CHECK_H
#define SEGEMENT_CHECK_H

#include <vector>
#include <list>
#include <map>
#include <set>
#include <math.h>

#include "CGAL_defines.h"

using namespace std;

//=============================================================================
class polygon_split_observer : public CGAL::Arr_observer<Arrangement_2>
{
    void after_split_face(Face_handle f1, Face_handle f2, bool) override;
};

//=============================================================================
class SegmentCheck
{
private:

  FT                                _rodLength;
  Polygon_set_2                     _freeSpace;
  Arrangement_2                     _arr;
  void verticalDecomposition( Kernel&        ker );
  void addVerticalSegment( Vertex_handle     v,
                           CGAL::Object      obj,
                           Kernel&           ker );

  //Face_const_handle getFace( const Landmarks_pl& pl,
  //                           const Point_2&      p ) const;

  FT pointsDistance(Point_2 a, Point_2 b) const;
  Segment_2 getSegment(Halfedge_const_handle edge) const;
  Segment_2 getSegment(Point_2 a, Point_2 b) const;

  Halfedge_const_handle
  findNextHEdge( Face_const_handle&           hCurrFace,
                 const Point_2&               endPt,
                 const Segment_2&             querySegment,
                 bool&                        bDone,
                 bool&                        bResult,
                 const Halfedge_const_handle* pAnchorHEdge );

  Face_const_handle chooseFace( Halfedge_const_handle hHEdge,
                                const Point_2&        endPt );
  Face_const_handle chooseFace( Vertex_const_handle hVrtx,
                                const Vector_2&     direction );

  void addFrame( const FT& rodLength );

public:
  SegmentCheck( const FT& rodLength, const vector<Polygon_2>& obstacles );
  bool isFree( const Point_2& startPt, const Vector_2& direction ) ;
};

#endif 
