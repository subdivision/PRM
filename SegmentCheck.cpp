#include "SegmentCheck.h"

const int SegmentCheck::EDGE_DIVIDING_PARAM = 2;
void print_arrangement( const Arrangement_2& arr );
void print_segment( const Point_2& p1, const Point_2& p2 );

//=============================================================================
PointNode::PointNode( FT                     distance,
                      Point_2                point,
                      PointNode*             prev,
                      Halfedge_const_handle  hHedge ):
distance(distance),
point(point),
prev(prev),
hedge(hHedge),
processed(false)
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
SegmentCheck::~SegmentCheck()
{
  delete[] _pNeigHEdges;
  delete[] _pNeigMtrx;
}
//-----------------------------------------------------------------------------
SegmentCheck::SegmentCheck( const FT& rodLength,
                            const vector<Polygon_2>& obstacles ):
_freeSpace()
{
  _freeSpace.join( obstacles.begin(), obstacles.end() );

  //complement to get the free space
  _freeSpace.complement();

  //set the free space as arrangment
  _arr = _freeSpace.arrangement();

  addFrame( rodLength );

  //ensure that when face split two side safe their property (inside/outside)
  Polygon_set_2::Traits_2 traits;
  polygon_split_observer observer;
  observer.attach(_arr);
  Kernel* ker = &traits;
  this->verticalDecomposition( *ker );
  observer.detach();
  buildNeigMatrices();
}

//-----------------------------------------------------------------------------
void
SegmentCheck::initFreeFacesVec()
{
  Arr_FaceCIter iFace = _arr.faces_begin();
  int nIdx = 0;
  for( ; iFace != _arr.faces_end(); ++iFace )
  {
    if( (*iFace).is_unbounded() )
      continue;
    if( !(*iFace).contained() )
      continue;
    Face_const_handle hFace = iFace;
    _freeFaces.push_back( hFace );
    _faceToIdx[hFace] = nIdx;
    ++nIdx;
  }
}

//-----------------------------------------------------------------------------
bool
SegmentCheck::twoFacesHaveCommonEdge( const Face_const_handle& hFace1,
                                      const Face_const_handle& hFace2,
                                      Halfedge_const_handle&   hRes )
{
  ccb_haledge_circulator circ1 = hFace1->outer_ccb();
  ccb_haledge_circulator circ2 = hFace2->outer_ccb();
  ccb_haledge_circulator curr1 = circ1;
  do
  {
    Halfedge_const_handle hHe1 = curr1;
    ccb_haledge_circulator curr2 = circ2;
    do
    {
      Halfedge_const_handle hHe2 = curr2;
      Point_2 p1s = hHe1->source()->point();
      Point_2 p2t = hHe2->target()->point();
      Point_2 p1t = hHe1->target()->point();
      Point_2 p2s = hHe2->source()->point();

      if( p1s == p2t && p1t == p2s && p1s[0] == p1t[0] )
      {
        hRes = hHe1;
        return true;
      }
    }while( ++curr2 != circ2 );
  }while( ++curr1 != circ1 );
  return false;
}
//-----------------------------------------------------------------------------
void
SegmentCheck::buildNeigMatrices()
{
  initFreeFacesVec();
  int i = 0, j = 0;
  int nMtrxSize = _freeFaces.size();
  _pNeigMtrx = new char[nMtrxSize * nMtrxSize];
  _pNeigHEdges = new Halfedge_const_handle[nMtrxSize * nMtrxSize];
  vector<Face_const_handle>::const_iterator iFace = _freeFaces.begin();
  for( ; iFace != _freeFaces.end(); ++iFace, ++i )
  {
    vector<Face_const_handle>::const_iterator iOtherFace = iFace;
    for( ++iOtherFace, j = i+1; iOtherFace != _freeFaces.end(); ++iOtherFace, ++j )
    {
      Halfedge_const_handle hCmnHEdge;
      bool bHaveCommon = twoFacesHaveCommonEdge( *iFace, *iOtherFace, hCmnHEdge);
      _pNeigMtrx[i*nMtrxSize + j] = bHaveCommon;
      _pNeigMtrx[j*nMtrxSize + i] = bHaveCommon;
      _pNeigHEdges[i*nMtrxSize + j] = hCmnHEdge;
      _pNeigHEdges[j*nMtrxSize + i] = hCmnHEdge;
    }
  }
}

//-----------------------------------------------------------------------------

extern int dijkstra_compute_path(char*, int, int, int, vector<int>& );

void
SegmentCheck::retrieveHEdges(Face_const_handle              hStartFace,
                             Face_const_handle              hEndFace,
                             vector<Halfedge_const_handle>& crossedHEdges)
{
  int nMtrxSize = _freeFaces.size();
  int i = _faceToIdx[hStartFace];
  int j = _faceToIdx[hEndFace];
  int s = i < j ? i : j;
  int d = i > j ? i : j;

  if( _pNeigMtrx[ i * nMtrxSize + j ] )
  {
    // Direct neighbors
    crossedHEdges.push_back( _pNeigHEdges[i * nMtrxSize + j] );
    return;
  }

  auto iIdxPath = _pathsBetweenFaces.find(pair<int, int>(s,d));
  vector<int> CurrIdxPath;
  if( iIdxPath == _pathsBetweenFaces.end() )
  {
    int nCurrLen = dijkstra_compute_path( _pNeigMtrx, nMtrxSize,
                                          s, d, CurrIdxPath );
    _pathsBetweenFaces[pair<int, int>(s,d)] = CurrIdxPath;
  }
  else
  {
    CurrIdxPath = iIdxPath->second;
  }

  vector<int>::const_iterator iPrevIdx = CurrIdxPath.begin();
  vector<int>::const_iterator iCurrIdx = CurrIdxPath.begin();
  for( ++iCurrIdx; iCurrIdx != CurrIdxPath.end(); ++iCurrIdx, ++iPrevIdx )
  {
    Halfedge_const_handle hHEdge = _pNeigHEdges[*iPrevIdx * nMtrxSize + *iCurrIdx];
    crossedHEdges.push_back(hHEdge);
  }
}

//-----------------------------------------------------------------------------
void
SegmentCheck::addVerticalSegment( Vertex_handle  v,
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
    v2 = _arr.non_const_handle(vh);
  }
  else if( CGAL::assign(hh, obj) )
  {
    // The given feature is a halfedge.
    if (hh->is_fictitious()) //We ignore fictitious halfedges.
    {
      return;
    }

    // Check whether v lies in the interior of the x-range of the edge (in
    // which case this edge should be split).
    const typename Kernel::Compare_x_2 cmp_x = ker.compare_x_2_object();
    if (cmp_x(v->point(), hh->target()->point()) == CGAL::EQUAL)
    {
        // In case the target of the edge already has the same x-coordinate as
        // the vertex v, just connect these two vertices.
        seg = X_monotone_curve_2(v->point(), hh->target()->point());
        v2 = _arr.non_const_handle(hh->target());
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
      _arr.split_edge(_arr.non_const_handle(hh),
                     X_monotone_curve_2(hh->source()->point(), point),
                     X_monotone_curve_2(point, hh->target()->point()));
      v2 = _arr.non_const_handle(hh->target());
    }
  }
  else
  {
    // Ignore faces and empty objects.
    return;
  }

    // Add the vertical segment to the arrangement using its two end vertices.
  _arr.insert_at_vertices(seg, v, v2);
}

//-----------------------------------------------------------------------------
void
SegmentCheck::verticalDecomposition( Kernel& ker )
{
  typedef pair<Vertex_const_handle, pair<CGAL::Object, CGAL::Object>> Vd_entry;

  // For each vertex in the arrangment, locate the feature that lies
  // directly below it and the feature that lies directly above it.
  list<Vd_entry>   vd_list;
  CGAL::decompose(_arr, back_inserter(vd_list));

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
      addVerticalSegment( _arr.non_const_handle(it->first),
                         it->second.first, ker);
    }
    // Add a vertical segment to the feature above the vertex.
    addVerticalSegment( _arr.non_const_handle(it->first),
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
      cout << "addPointToQueue " << tempPoint << "-> ("
           << tempEdge->source()->point() << ") - ("
           << tempEdge->target()->point() << ")" << endl;
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
    cout << "Improving to" << tempPoint << endl;
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
bool
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
    return true;

  addStartPathToQueue( startPt, endPt, hStartFace, hEndFace,
                       queue, pointsMap, edgesMap );
  while(!queue.empty())
  {
    PointNode* pointNode = *(queue.begin());
    if( pointNode->point == endPt )
      return true;
    pointNode->processed = true;
    this->addFacesToQueue( pointNode, endPt, hEndFace,
                           queue, pointsMap, edgesMap );
    queue.erase( pointNode );
  }
  return false;
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
void
SegmentCheck::reversedPath( const Point_2&                 startPt,
                            const Point_2&                 endPt,
                            map<Point_2, PointNode>&       pointsMap,
                            vector<Halfedge_const_handle>& crossedHEdges ) const
{
  //create path from BFS results
  PointNode* pCurr = &(pointsMap[endPt]);
  while( pCurr != nullptr )
  {
    //cout << "Hedge at 0x" << static_cast<long> (pCurr->hedge) << endl;
    //crossedHEdges.push_back(pCurr->hedge);
    PointNode* pPrev = pCurr->prev;
    crossedHEdges.push_back(pPrev->hedge);

    if(pPrev->prev == nullptr)
    {
      //path.push_back(prev->point);
      break;
    }
    PointNode* pPrevPrev = pPrev->prev;
    while(pPrevPrev != nullptr)
    {
      if( CGAL::intersection( getSegment( pPrev->hedge ),
                              getSegment( pCurr->point, pPrevPrev->point ) ) )
      {
        bool segmentFine = true;
        for( PointNode* temp : pCurr->crossedSegments )
        {
          if( !CGAL::intersection( getSegment( temp->hedge ),
                                   getSegment( pCurr->point, pPrevPrev->point )))
          {
            segmentFine = false;
            break;
          }
        }
        if(segmentFine)
          pCurr->crossedSegments.push_back(pPrev);
        else
          break;
      }
      else
        break;
      pPrev = pPrevPrev;
      if( pPrev->prev == nullptr )
      {
        //path.push_back(pPrev->point);
        pPrev = pPrev->prev;
        break;
      }
      pPrevPrev = pPrev->prev;
    } // end of (pPrevPrev != nullptr)
    pCurr = pPrev;
  }
}

//-----------------------------------------------------------------------------
bool
SegmentCheck::isFree( const Point_2& startPt, const Point_2& endPt )
{
try{
  cout << "--------------------------" << endl;
  print_segment(startPt, endPt);
  //bfs maps:
  //use set because need to delete efficiently
//  set<PointNode*, CmpfaceNodePtrs>      queue;
//  map<Point_2, PointNode>               pointsMap;
//  //use to improve finding all interesting point on edge
//  map<Halfedge_const_handle, vector<Point_2>> edgesMap;
//
//  if( !setFacesPath( startPt, endPt, queue, pointsMap, edgesMap ) )
//    return false;
//
//
//  vector<Halfedge_const_handle> CrossedHEdges;
//  reversedPath( startPt, endPt, pointsMap, CrossedHEdges);

  Landmarks_pl pl( _arr );

  //start and end faces:
  Face_const_handle hStartFace = getFace( pl, startPt );
  Face_const_handle hEndFace   = getFace( pl, endPt   );
  if( hStartFace == hEndFace)
  {
    cout << "Same face. True." << endl;
    return true;
  }

  vector<Halfedge_const_handle> CrossedHEdges;
  retrieveHEdges(hStartFace, hEndFace, CrossedHEdges);

  if( 0 == CrossedHEdges.size() )
  {
    cout << "No Hedges crossed. But not same face? False" << endl;
    return false;
  }

  bool bVertQuery = startPt.x() == endPt.x();
  FT k = 0;
  if( !bVertQuery )
    k = ( startPt.y() - endPt.y() ) / ( startPt.x() - endPt.x() );
  bool bResult = true;
  cout << "CrossedHEdges.size() = " <<  CrossedHEdges.size() << endl;
  for( Halfedge_const_handle iCurrHEdge : CrossedHEdges)
  {
    cout << "Testing HEdge ";
    print_segment( iCurrHEdge->source()->point(), iCurrHEdge->target()->point() );
    //if( !bResult )
    //  break;
    FT x0 = iCurrHEdge->source()->point().x();
    FT y0 = iCurrHEdge->source()->point().y();
    FT x1 = iCurrHEdge->target()->point().x();
    FT y1 = iCurrHEdge->target()->point().y();
    FT maxy = y0 > y1 ? y0 : y1;
    FT miny = y0 < y1 ? y0 : y1;
    if( x0 != x1 )
    {
      cout << "Not vertical hedge" << endl;
    }
    else
    {
      if( bVertQuery )
      {
        bResult = x0 == startPt.x();
      }
      else
      {
        FT yInx = k * (x0 - endPt.x()) + endPt.y();
        bResult = (miny <= yInx && yInx <= maxy);
      }
    }
  }  
  cout << "bResult = " << ( bResult ? "True" : "False" ) << endl;
  return bResult;
}
catch(const char* pExcText)
{
  cout << pExcText << " False." << endl;
  return false;
}
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
  tempEdge = _arr.insert_from_right_vertex(rightBound, tempEdge->target());
  tempEdge = _arr.insert_from_right_vertex(lowerBound, tempEdge->target());
  tempEdge = _arr.insert_at_vertices(leftBound, tempEdge->target(), startVertex);
  tempEdge->twin()->face()->set_contained(true);
  _arr.unbounded_face()->set_contained(false);
}
//============================== END OF FILE ==================================
// Debuf stuff. 
//----------------------------------------------------------------------------
//
// Iterating through a DCEL face boundary
void print_ccb( ccb_haledge_circulator circ )
{
  ccb_haledge_circulator curr = circ;
  int k = 1;
  do{
    ++k;
  }while(++curr != circ);
  cout << k << " ";
  const Point_2& pt = curr->source()->point();
  cout << CGAL::to_double(pt[0]) << " " 
       << CGAL::to_double(pt[1]) << " ";
  do {
    const Point_2& pt = curr->target()->point();
    cout << CGAL::to_double(pt[0]) << " " 
         << CGAL::to_double(pt[1]) << " ";
  } while (++curr != circ);
  std::cout << std::endl;
}


void print_arr_face(Face_const_handle f)
{
  // Print the outer boundary.
  if (f->is_unbounded())
    cout << "Unbounded face. " << std::endl;
  else {
    std::cout << "Outer boundary: ";
    print_ccb (f->outer_ccb());
  }
  // Print the boundary of each of the holes.
  Arrangement_2::Hole_const_iterator hi;
  int index = 1;
  for (hi = f->holes_begin(); hi != f->holes_end(); ++hi, ++index) {
    std::cout << " Hole #" << index << ": ";
    print_ccb (*hi);
  }
}

void print_arrangement( const Arrangement_2& arr )
{
  Arrangement_2::Face_const_iterator iFace = arr.faces_begin();
  for( ; iFace != arr.faces_end(); ++iFace )
  {
    if ((*iFace).is_unbounded())
      continue;
    print_arr_face( iFace );
  }
}

void print_segment( const Segment_2& s )
{
  const Point_2& p1 = s.source();
  const Point_2& p2 = s.target();
  double x1 = CGAL::to_double( p1[0] );
  double y1 = CGAL::to_double( p1[1] );
  double x2 = CGAL::to_double( p2[0] );
  double y2 = CGAL::to_double( p2[1] );
  cout << "[" << x1 << ", " << y1 << " - " << x2 << ", " << y2 << "]" << endl;
}

void print_segment( const Point_2& p1, const Point_2& p2 )
{
  double x1 = CGAL::to_double( p1[0] );
  double y1 = CGAL::to_double( p1[1] );
  double x2 = CGAL::to_double( p2[0] );
  double y2 = CGAL::to_double( p2[1] );
  cout << "[" << x1 << ", " << y1 << " - " << x2 << ", " << y2 << "]" << endl;
}

void print_point( const Point_2& p)
{
  double x1 = CGAL::to_double( p[0] );
  double y1 = CGAL::to_double( p[1] );
  cout << "(" << x1 << ", " << y1 << ")" << endl;
}


