#ifndef PATHFINDER_PATHPLANNER_H
#define PATHFINDER_PATHPLANNER_H

#include <vector>
#include <list>
#include <map>
#include <set>
#include <math.h>

#include "CGAL_defines.h"

using namespace std;

//=============================================================================
class PointNode
{
public:
    explicit PointNode( FT                    distance,
                        Point_2               point,
                        PointNode*            prev,
                        Halfedge_const_handle edge );
    PointNode() = default;

    FT distance;
    Point_2 point;
    PointNode* prev;
    Halfedge_const_handle hedge;
    bool processed = false;
    vector<PointNode*> crossedSegments;
};

//=============================================================================
struct CmpfaceNodePtrs
{
    bool operator()(const PointNode* lhs, const PointNode* rhs) const;
};

//=============================================================================
class polygon_split_observer : public CGAL::Arr_observer<Arrangement_2>
{
    void after_split_face(Face_handle f1, Face_handle f2, bool) override;
};

//=============================================================================
class SegmentCheck
{
private:
  static const int EDGE_DIVIDING_PARAM;

  Polygon_set_2           _freeSpace;
  Arrangement_2           _arr;

  void addFrame( const FT& rodLength );

  void verticalDecomposition( Kernel&        ker );
  void addVerticalSegment( Vertex_handle     v,
                           CGAL::Object      obj,
                           Kernel&           ker );

  Face_const_handle getFace( const Landmarks_pl& pl,
                             const Point_2&      p ) const;

  void setFacesPath(const Point_2& startPt,
                      const Point_2& endPt,
                      set<PointNode*, CmpfaceNodePtrs>& queue,
                      map<Point_2, PointNode>&          pointsMap,
                      map<Halfedge_const_handle, vector<Point_2>>& edgesMap ) const;

  vector<Point_2>
  getEdgePoints( Halfedge_const_handle hHedge,
                 map<Halfedge_const_handle, vector<Point_2>>& edgesMap ) const;

  void addFacesToQueue( PointNode*                        pPointNode,
                        const Point_2&                    endPt,
                        Face_const_handle                 hEndFace,
                        set<PointNode*, CmpfaceNodePtrs>& queue,
                        map<Point_2, PointNode>&          pointsMap,
                        map<Halfedge_const_handle, vector<Point_2>>& edgesMap) const;

  void addFaceToQueue( PointNode*                        pPointNode,
                       Face_const_handle                 hFace,
                       set<PointNode*, CmpfaceNodePtrs>& queue,
                       map<Point_2, PointNode>&          pointsMap,
                       map<Halfedge_const_handle, vector<Point_2>>& edgesMap ) const;

  void addPointToQueue( PointNode*                        pointNode,
                        Point_2                           tempPoint,
                        Halfedge_const_handle             tempEdge,
                        set<PointNode*, CmpfaceNodePtrs>& queue,
                        map<Point_2, PointNode>&          pointsMap ) const;

  void tryToImprove( PointNode*                         pointNode,
                     Point_2                            tempPoint,
                     set<PointNode*, CmpfaceNodePtrs>&  queue,
                     map<Point_2, PointNode>&           pointsMap )const;

  void addStartPathToQueue( const Point_2&                    startPt,
                            const Point_2&                    endPt,
                            Face_const_handle&                hStartFace,
                            Face_const_handle&                hEndFace,
                            set<PointNode*, CmpfaceNodePtrs>& queue,
                            map<Point_2, PointNode>&          pointsMap,
                            map<Halfedge_const_handle, vector<Point_2>>& edgesMap ) const;

  FT pointsDistance(Point_2 a, Point_2 b) const;
  Segment_2 getSegment(Halfedge_const_handle edge) const;
  Segment_2 getSegment(Point_2 a, Point_2 b) const;

  //vector<Point_2> reversedPath(Arrangement_2& arr, Kernel& ker);

public:
    SegmentCheck( const FT& rodLength, const vector<Polygon_2>& obstacles );
    bool isFree( const Point_2& startPt, const Point_2& endPt ) const;
};

#endif //PATHFINDER_PATHPLANNER_H
