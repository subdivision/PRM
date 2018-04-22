#include "SegmentCheck.h"

const int SegmentCheck::EDGE_DIVIDING_PARAM = 5;

//=============================================================================
PointNode::PointNode( FT                     distance,
                      Point_2                point,
                      PointNode*             prev,
                      Halfedge_const_handle  hHedge ):
distance(distance), point(point), prev(prev), hedge(hHedge)
{}

//=============================================================================
bool
CmpfaceNodePtrs::operator()( const PointNode* lhs,
                             const PointNode* rhs) const
{
    if(lhs->distance != rhs->distance)
        return lhs->distance < rhs->distance;

    return lhs->point < rhs->point;
}

//=============================================================================
void
polygon_split_observer::after_split_face( Face_handle f1, Face_handle f2, bool )
{
  f2->set_contained(f1->contained());
}

//=============================================================================
SegmentCheck::SegmentCheck( const FT& rodLength,
                            const vector<Polygon_2>& obstacles ):
_traits()
{
  _freeSpace.join( obstacles.begin(), obstacles.end() );

  //complement to get the free space
  _freeSpace.complement();

  //set the free space as arrangment
  _arr = _freeSpace.arrangement();
  addFrame( rodLength );

  //ensure that when face split two side safe their property (inside/outside)
  polygon_split_observer observer;
  observer.attach( _arr );
  _pKer = &_traits;
  this->verticalDecomposition( _arr, *_pKer);
  observer.detach();

}

//-----------------------------------------------------------------------------
void
SegmentCheck::addVerticalSegment( Arrangement_2& arr,
                                  Vertex_handle  v,
                                  CGAL::Object   obj,
                                  Kernel&        ker )
{
  X_monotone_curve_2 seg;
  Vertex_const_handle vh;
  Halfedge_const_handle hh;
  Face_const_handle fh;
  Vertex_handle v2;

  if( CGAL::assign(vh, obj) )
  {
    // The given feature is a vertex.
    seg = X_monotone_curve_2(v->point(), vh->point());
    v2 = arr.non_const_handle(vh);
  }
  else if( CGAL::assign(hh, obj) )
  {
    // The given feature is a halfedge.
    if (hh->is_fictitious()) //We ignore fictitious halfedges.
        return;

    // Check whether v lies in the interior of the x-range of the edge (in
    // which case this edge should be split).
    const typename Kernel::Compare_x_2 cmp_x = ker.compare_x_2_object();
    if (cmp_x(v->point(), hh->target()->point()) == CGAL::EQUAL)
    {
        // In case the target of the edge already has the same x-coordinate as
        // the vertex v, just connect these two vertices.
        seg = X_monotone_curve_2(v->point(), hh->target()->point());
        v2 = arr.non_const_handle(hh->target());
    }
    else
    {
      // Compute the vertical projection of v onto the segment associated
      // with the halfedge. Split the edge and connect v with the split point.
      Line_2 Line;
      Line_2 supp_line(hh->source()->point(), hh->target()->point());
      Line_2 vert_line(v->point(), Point_2(v->point().x(), v->point().y() + 1));
      Point_2  point;
      CGAL::assign(point, ker.intersect_2_object()(supp_line, vert_line));
      seg = X_monotone_curve_2(v->point(), point);
      arr.split_edge(arr.non_const_handle(hh),
                     X_monotone_curve_2(hh->source()->point(), point),
                     X_monotone_curve_2(point, hh->target()->point()));
      v2 = arr.non_const_handle(hh->target());
    }
  }
  else
  {
    // Ignore faces and empty objects.
    return;
  }

    // Add the vertical segment to the arrangement using its two end vertices.
    arr.insert_at_vertices(seg, v, v2);
}

//-----------------------------------------------------------------------------
void
SegmentCheck::verticalDecomposition( Arrangement_2& arr,
                                     Kernel&        ker )
{
  typedef pair<Vertex_const_handle, pair<CGAL::Object, CGAL::Object>> Vd_entry;

  // For each vertex in the arrangment, locate the feature that lies
  // directly below it and the feature that lies directly above it.
  list<Vd_entry>   vd_list;
  CGAL::decompose(arr, back_inserter(vd_list));

  // Go over the vertices (given in ascending lexicographical xy-order),
  // and add segements to the feautres below and above it.
  const typename Kernel::Equal_2 equal = ker.equal_2_object();
  typename list<Vd_entry>::iterator  it, prev = vd_list.end();
  for (it = vd_list.begin(); it != vd_list.end(); ++it)
  {
    // If the feature above the previous vertex is not the current vertex,
    // add a vertical segment to the feature below the vertex.
    Vertex_const_handle v;
    if(  ( prev == vd_list.end() )               ||
        !CGAL::assign( v, prev->second.second  ) ||
        !equal( v->point(), it->first->point() )   )
    {
      addVerticalSegment(arr, arr.non_const_handle(it->first),
                         it->second.first, ker);
    }
    // Add a vertical segment to the feature above the vertex.
    addVerticalSegment( arr, arr.non_const_handle(it->first),
                        it->second.second, ker);
    prev = it;
  }
}

//-----------------------------------------------------------------------------
Face_const_handle
SegmentCheck::getFace( const Landmarks_pl &pl, const Point_2 &p ) const
{
  CGAL::Object obj = pl.locate(p); //find p in pl

  Vertex_const_handle vertex;
  if( CGAL::assign( vertex, obj ) )
  {
    Arr_FaceCIter iFace = _arr.faces_begin();
    for( ; iFace != _arr.faces_end(); ++iFace )
    {
      if( iFace == _arr.unbounded_face() )
        continue;

      ccb_haledge_circulator first = iFace->outer_ccb();
      ccb_haledge_circulator circ = first;
      do
      {
        Halfedge_const_handle temp = circ;
        if( temp->source()->point() == vertex->point() )
        {
          if( iFace->contained() )
            return iFace;
          else
             break;
        }
      }
      while( ++circ != first );
    }
    throw "point is not in a legal position - "
          "on a vertex between obstacles faces";
  }

  // Check it's a halfedge
  Halfedge_const_handle  hHEdge;
  if( CGAL::assign( hHEdge, obj ) )
  {
    if( hHEdge->face()->contained() )
      return hHEdge->face();
    else if( hHEdge->twin()->face()->contained() )
      return hHEdge->twin()->face();
    throw "point is not in a legal position - "
          "on an edge between two obstacles faces";
  }

  // Check whether the point is contained inside a free bounded face.
  Face_const_handle hFace;
  if( CGAL::assign( hFace, obj ) ) //if obj is face
  {
    if( hFace->contained() )
      return hFace;
  }
  throw "point is not in a legal position - "
        "inside an obstacle face";
}
//-----------------------------------------------------------------------------
FT
SegmentCheck::pointsDistance( Point_2 a_point, Point_2 b_point ) const
{
    FT distance = (a_point.x() - b_point.x()) * (a_point.x() - b_point.x()) +
               (a_point.y() - b_point.y()) * (a_point.y() - b_point.y());

    return sqrt(CGAL::to_double(distance));;
}

//-----------------------------------------------------------------------------
vector<Point_2>
SegmentCheck::getEdgePoints( Halfedge_const_handle hHEdge,
                             map<Halfedge_const_handle, vector<Point_2>>& edgesMap )
                                                                          const
{
  auto search = edgesMap.find(hHEdge);

  //if edge in database take from it
  if(search != edgesMap.end())
      return edgesMap[hHEdge];

  Point_2 source = hHEdge->source()->point();
  Point_2 target = hHEdge->target()->point();
  FT distance = pointsDistance(source, target);
  FT xDistance = target.x() - source.x();
  FT yDistance = target.y() - source.y();
  vector<Point_2> points;
  for(int i=0; i<=SegmentCheck::EDGE_DIVIDING_PARAM; i++)
      points.push_back({ source.x()
                         + i*xDistance/SegmentCheck::EDGE_DIVIDING_PARAM,
                         source.y()
                         + i*yDistance/SegmentCheck::EDGE_DIVIDING_PARAM });

  edgesMap[hHEdge] = points;

  return points;
}

//-----------------------------------------------------------------------------
void
SegmentCheck::addPointToQueue( PointNode*                        pointNode,
                               Point_2                           tempPoint,
                               Halfedge_const_handle             tempEdge,
                               set<PointNode*, CmpfaceNodePtrs>& queue,
                               map<Point_2, PointNode>&          pointsMap )
                                                                         const
{
    if(tempPoint == pointNode->point)
        return;
    auto search = pointsMap.find(tempPoint);
    if(search != pointsMap.end()) //if face already exist try to improve
      tryToImprove( pointNode, tempPoint, queue, pointsMap );
    else
    {
      FT tempDistance = pointNode->distance
                        + pointsDistance(pointNode->point, tempPoint);
      pointsMap[tempPoint] = PointNode(tempDistance, tempPoint,
                                       pointNode, tempEdge);
      queue.insert(&(pointsMap[tempPoint]));
    }
}

//-----------------------------------------------------------------------------
void
SegmentCheck::addFaceToQueue( PointNode*                        pPointNode,
                              Face_const_handle                 hFace,
                              set<PointNode*, CmpfaceNodePtrs>& queue,
                              map<Point_2, PointNode>&          pointsMap,
                              map<Halfedge_const_handle, vector<Point_2>>& edgesMap )
                                                                          const
{
  ccb_haledge_circulator ciFirst = hFace->outer_ccb();
  ccb_haledge_circulator ciCurr  = ciFirst;
  do
  {
    Halfedge_const_handle hTempHEdge = ciCurr;
    if( !hTempHEdge->twin()->face()->contained() )
      continue;
    for( Point_2 point: getEdgePoints( hTempHEdge, edgesMap ) )
      addPointToQueue( pPointNode, point, hTempHEdge, queue, pointsMap);
  }
  while( ++ciCurr != ciFirst );
}

//-----------------------------------------------------------------------------
void
SegmentCheck::tryToImprove( PointNode*                         pointNode,
                            Point_2                            tempPoint,
                            set<PointNode*, CmpfaceNodePtrs>&  queue,
                            map<Point_2, PointNode>&           pointsMap )const
{
    PointNode* temp = &(pointsMap[tempPoint]);
    if(temp->processed)
        return;
    FT tempDistance = pointNode->distance
                      + pointsDistance(pointNode->point, tempPoint);
    if(tempDistance < temp->distance)
    {
        queue.erase(temp); //remove from set beacuse it's in wrong position
        temp = &(pointsMap[tempPoint]);
        temp->distance = tempDistance;
        temp->prev = pointNode;
        queue.insert(temp); //insert in the right position
    }
}

//-----------------------------------------------------------------------------
void
SegmentCheck::addFacesToQueue( PointNode*                        pPointNode,
                               const Point_2&                    endPt,
                               Face_const_handle                 hEndFace,
                               set<PointNode*, CmpfaceNodePtrs>& queue,
                               map<Point_2, PointNode>&          pointsMap,
                               map<Halfedge_const_handle, vector<Point_2>>& edgesMap )
                                                                          const
{
  Face_const_handle hFirstFace  = pPointNode->hedge->face();
  Face_const_handle hSecondFace = pPointNode->hedge->twin()->face();
  if( hFirstFace == hEndFace || hSecondFace == hEndFace )
  {
    auto search = pointsMap.find( endPt );
    if( search != pointsMap.end() )
      tryToImprove( pPointNode, endPt, queue, pointsMap );
    else{
      FT tempDistance = pPointNode->distance
                        + pointsDistance(pPointNode->point, endPt);
      pointsMap[endPt] = PointNode( tempDistance, endPt, pPointNode,
                                    Halfedge_handle());
      queue.insert(&(pointsMap[endPt]));
    }
    return;
  }
  addFaceToQueue( pPointNode, hFirstFace, queue, pointsMap, edgesMap  );
  addFaceToQueue( pPointNode, hSecondFace, queue, pointsMap, edgesMap );
}

//-----------------------------------------------------------------------------
void
SegmentCheck::addStartPathToQueue( const Point_2&                    startPt,
                                   const Point_2&                    endPt,
                                   Face_const_handle&                hStartFace,
                                   Face_const_handle&                hEndFace,
                                   set<PointNode*, CmpfaceNodePtrs>& queue,
                                   map<Point_2, PointNode>&          pointsMap,
                                   map<Halfedge_const_handle, vector<Point_2>>& edgesMap )
                                                                          const
{
  pointsMap[startPt] = PointNode(0, startPt, nullptr, Halfedge_handle());
  PointNode* pointNode = &( pointsMap[startPt] );
  if( hStartFace == hEndFace)
  {
      pointsMap[endPt] = PointNode(0, endPt, pointNode, Halfedge_handle());
      queue.insert(&(pointsMap[endPt]));
      return;
  }

  ccb_haledge_circulator ciFirst = hStartFace->outer_ccb();
  ccb_haledge_circulator ciCurr  = ciFirst;
  do
  {
    Halfedge_const_handle hTempHEdge = ciCurr;
    for( Point_2 point: getEdgePoints( hTempHEdge, edgesMap ) )
    {
      if(point == pointNode->point && hTempHEdge->twin()->face()->contained())
      {
          pointNode->hedge = hTempHEdge;
          queue.insert( &(pointsMap[startPt]) );
          return;
      }
    }
  }
  while( ++ciCurr != ciFirst );

  do
  {
    Halfedge_const_handle hTempHEdge = ciCurr;
    if( !hTempHEdge->twin()->face()->contained() )
      continue;
    for( Point_2 point: getEdgePoints( hTempHEdge, edgesMap ) )
      addPointToQueue( pointNode, point, hTempHEdge, queue, pointsMap );
  }
  while( ++ciCurr != ciFirst );
}

//-----------------------------------------------------------------------------
// run BFS from start_face to end_face
void
SegmentCheck::setFacesPath( const Point_2&                         startPt,
                            const Point_2&                         endPt,
                            set<PointNode*, CmpfaceNodePtrs>&      queue,
                            map<Point_2, PointNode>&               pointsMap,
                            map<Halfedge_const_handle, vector<Point_2>>& edgesMap )
                                                                          const
{
  Landmarks_pl pl( _arr );

  //start and end faces:
  Face_const_handle hStartFace = getFace( pl, startPt );
  Face_const_handle hEndFace   = getFace( pl, endPt   );
  if( hStartFace == hEndFace)
    return;

  addStartPathToQueue( startPt, endPt, hStartFace, hEndFace,
                       queue, pointsMap, edgesMap );
  while(!queue.empty())
  {
    PointNode* pointNode = *(queue.begin());
    if( pointNode->point == endPt )
      return;
    pointNode->processed = true;
    this->addFacesToQueue( pointNode, endPt, hEndFace,
                           queue, pointsMap, edgesMap );
    queue.erase( pointNode );
  }
  throw "no path found!";
}

//-----------------------------------------------------------------------------
Segment_2
SegmentCheck::getSegment( Halfedge_const_handle hHEdge ) const
{
  return getSegment( hHEdge->source()->point(), hHEdge->target()->point() );
}

//-----------------------------------------------------------------------------
Segment_2
SegmentCheck::getSegment( Point_2 a, Point_2 b ) const
{
  return {a,b};
}

//-----------------------------------------------------------------------------
/*
vector<Point_2>
SegmentCheck::reversedPath()
{ //create path from BFS results
    vector<Point_2> path;
    PointNode *node = &(pointsMap[this->endPoint]);
    while (node != nullptr)
    {
        path.push_back(node->point);
        PointNode *prev = node->prev;
        if(prev->prev == nullptr)
        {
            path.push_back(prev->point);
            break;
        }
        PointNode *prevprev = prev->prev;
        while(prevprev != nullptr)
        {
            if(CGAL::intersection(getSegment(prev->edge), getSegment(node->point, prevprev->point)))
            {
                bool segmentFine = true;
                for(PointNode *temp:node->crossedSegments)
                {
                    if(!CGAL::intersection(getSegment(temp->edge), getSegment(node->point, prevprev->point)))
                    {
                        segmentFine = false;
                        break;
                    }
                }
                if(segmentFine)
                    node->crossedSegments.push_back(prev);
                else
                    break;
            } else
                break;
            prev = prevprev;
            if(prev->prev == nullptr)
            {
                path.push_back(prev->point);
                prev = prev->prev;
                break;
            }
            prevprev = prev->prev;
        }
        node = prev;
    }
    return path;
}
*/
//-----------------------------------------------------------------------------
bool
SegmentCheck::isFree( const Point_2& startPt, const Point_2& endPt ) const
{
  //bfs maps:
  //use set because need to delete efficiently
  set<PointNode*, CmpfaceNodePtrs>      queue;
  map<Point_2, PointNode>               pointsMap;
  //use to improve finding all interesting point on edge
  map<Halfedge_const_handle, vector<Point_2>> edgesMap;

  setFacesPath( startPt, endPt, queue, pointsMap, edgesMap );
//  vector<Point_2> reversedPath = this->reversedPath();
//  vector<Point_2> path;
//  for(int i = static_cast<int>(reversedPath.size())-1; i >= 0; i--)
//      path.push_back(Point_2(reversedPath[i]));

  //return path;
  return true;
}
//-----------------------------------------------------------------------------
void SegmentCheck::addFrame( const FT& rodLength )
{
  FT d = rodLength * 1.5;

  Arr_VrtxCIter iVrtx = _arr.vertices_begin();
  FT mostLeft  = iVrtx->point().x();
  FT mostRight = iVrtx->point().x();
  FT mostUp    = iVrtx->point().y();
  FT mostDown  = iVrtx->point().y();

  for(; iVrtx != _arr.vertices_end(); ++iVrtx )
  {
    if( iVrtx->point().x() < mostLeft )
      mostLeft = iVrtx->point().x();
    if( iVrtx->point().x() > mostRight )
      mostRight = iVrtx->point().x();
    if( iVrtx->point().y() < mostDown )
      mostDown = iVrtx->point().y();
    if( iVrtx->point().y() > mostUp )
      mostUp = iVrtx->point().y();
  }
  Point_2     upperLeft(  mostLeft  - d, mostUp   + d ),
              upperRight( mostRight + d, mostUp   + d ),
              lowerRight( mostRight + d, mostDown - d ),
              lowerLeft(  mostLeft  - d, mostDown - d );

  Segment_2   upperBound(upperLeft, upperRight),
              rightBound(upperRight, lowerRight),
              lowerBound(lowerRight, lowerLeft),
              leftBound(lowerLeft, upperLeft);

  Halfedge_handle tempEdge = _arr.insert_in_face_interior(
                                           upperBound, _arr.unbounded_face());
  Vertex_handle startVertex = tempEdge->source();
  tempEdge = _arr.insert_from_left_vertex(rightBound, tempEdge->target());
  tempEdge = _arr.insert_from_right_vertex(lowerBound, tempEdge->target());
  tempEdge = _arr.insert_at_vertices(leftBound, tempEdge->target(), startVertex);

  tempEdge->twin()->face()->set_contained(true);
  _arr.unbounded_face()->set_contained(false);
}
//============================== END OF FILE ==================================