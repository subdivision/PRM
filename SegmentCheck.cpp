#include "SegmentCheck.h"

void print_arr_face(Face_const_handle f);
void print_arrangement( const Arrangement_2& arr );
void print_segment( const Point_2& p1, const Point_2& p2 );
void print_segment( const Segment_2& s );
void print_halfedge(Halfedge_const_handle hHEdge );
void print_point( const Point_2& p);
void print_direction( const Direction_2& d );

//=============================================================================
void
polygon_split_observer::after_split_face( Face_handle f1, Face_handle f2, bool )
{
  f2->set_contained(f1->contained());
}

//=============================================================================
FaceData::FaceData():
nInx(-1), 
bStartProperInside(false), 
bStartOnBndry(false), 
bEndProperInside(false),
bEndOnBndry(false){}

//----------------------------------------------------------------------------
void 
FaceData::reset() 
{ 
  nInx = -1; 
  vecStrHEdges.clear(); 
  vecEndHEdges.clear(); 
  vecInxHEdges.clear(); 
  vecInxPts.clear(); 
}

//----------------------------------------------------------------------------
bool 
FaceData::hasData() 
{ 
  return -1 != nInx; 
}

//-----------------------------------------------------------------------------
#define b2s(b) (b ? "True":"False")
#define prv(v , vname, prfn) \
if( 0 == v.size() ) \
{\
  cout << vname << " empty" << endl;\
}\
else\
{\
  cout << vname << endl;\
  for( auto o : v )\
  {\
      prfn( o );\
  }\
}

void
FaceData::print()
{
  cout << "~~~~~~~~~~~~~~~~~" << endl;
  print_arr_face(hFace);
  cout << "nInx = " << nInx << endl;
  cout << "bStartProperInside = " << b2s( bStartProperInside ) << endl;
  cout << "bStartOnBndry      = " << b2s( bStartOnBndry ) << endl;
  cout << "bEndProperInside   = " << b2s(bEndProperInside) << endl; 
  cout << "bEndOnBndry        = " << b2s(bEndOnBndry) << endl;
  prv(vecStrHEdges, "Start hedges", print_halfedge )
  prv(vecEndHEdges, "End   hedges", print_halfedge )
  prv(vecInxHEdges, "Intersection hedges", print_halfedge )
  prv(vecInxPts, "Intersection points", print_point)
}

//=============================================================================
SegmentCheck::SegmentCheck( const FT& rodLength,
                            const vector<Polygon_2>& obstacles ):
_rodLength(rodLength),
_freeSpace()
{
  _freeSpace.join( obstacles.begin(), obstacles.end() );

  //complement to get the free space
  _freeSpace.complement();

  //set the free space as arrangment
  _arr = _freeSpace.arrangement();

  addFrame();

  //ensure that when face split two side safe their property (inside/outside)
  Polygon_set_2::Traits_2 traits;
  polygon_split_observer observer;
  observer.attach(_arr);
  Kernel* ker = &traits;
  verticalDecomposition( *ker );
  observer.detach();
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
FT
SegmentCheck::pointsDistance( Point_2 a_point, Point_2 b_point ) const
{
    FT distance = (a_point.x() - b_point.x()) * (a_point.x() - b_point.x()) +
               (a_point.y() - b_point.y()) * (a_point.y() - b_point.y());

    return sqrt(CGAL::to_double(distance));;
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
bool
SegmentCheck::pointProperlyOut( const Point_2& pt )
{
  FT x = pt.x();
  FT y = pt.y();
  return x < _minx || x > _maxx || y < _miny || y > _maxy;
}

//----------------------------------------------------------------------------
bool
SegmentCheck::pointOnTheOutBndry( const Point_2& pt )
{
  FT x = pt.x();
  FT y = pt.y();
  bool bUp = y == _maxy && _minx <= x && x <= _maxx;
  bool bDw = y == _miny && _minx <= x && x <= _maxx;
  bool bLf = x == _minx && _miny <= y && y <= _maxy;
  bool bRg = x == _maxx && _miny <= y && y <= _maxy;
  return bUp || bDw || bLf || bRg;
}
//----------------------------------------------------------------------------
FaceData
SegmentCheck::getUnboundedFaceData( Face_const_handle& hFace,
                                    const Segment_2& querySegment )
{
  FaceData res;
  res.hFace = hFace;
  res.nInx = 0;
  Point_2 startPt = querySegment.source();
  Point_2 endPt   = querySegment.target(); 
  res.bStartProperInside = pointProperlyOut( startPt );
  res.bEndProperInside = pointProperlyOut( endPt );
  res.bStartOnBndry = pointOnTheOutBndry( startPt );
  res.bEndOnBndry = pointOnTheOutBndry( endPt );
  ccb_haledge_circulator ciFirstHEdge = *(hFace->holes_begin()); 
  ccb_haledge_circulator ciCurrHedge = ciFirstHEdge;
  do
  {
    Halfedge_const_handle hHEdgeToCheck = ciCurrHedge;
    Segment_2 currSeg = getSegment(hHEdgeToCheck);
    CGAL::cpp11::result_of<Intersect_2(Segment_2, Segment_2)>::type
      result = intersection( currSeg, querySegment );
    const Point_2* p = nullptr;
    if( result && ( p = boost::get<Point_2>(&*result) ) )
    {
      ++res.nInx;
      res.vecInxHEdges.push_back( hHEdgeToCheck );
      res.vecInxPts.push_back( *p );
    }
  }
  while( ++ciCurrHedge != ciFirstHEdge );
 
  return res;
}

//----------------------------------------------------------------------------
FaceData 
SegmentCheck::getFaceData( Face_const_handle&  hFace,
                           const Segment_2&    querySegment )
{
  ccb_haledge_circulator ciFirstHEdge = hFace->outer_ccb();
  if( hFace->is_unbounded() )
    return getUnboundedFaceData( hFace, querySegment );
  FaceData res;
  res.hFace = hFace;
  res.nInx = 0;
  ccb_haledge_circulator ciCurrHedge = ciFirstHEdge;
  Point_2 startPt = querySegment.source();
  Point_2 endPt   = querySegment.target(); 
  int nStartPtOnLeft = 0;
  int nEndPtOnLeft   = 0;
  int nOfEdges = 0; 
  do
  {
    Halfedge_const_handle hHEdgeToCheck = ciCurrHedge;
    Segment_2 currSeg = getSegment(hHEdgeToCheck);

    if( currSeg.has_on( endPt ) )
    {
      ++res.nInx;
      res.bEndOnBndry = true;
      res.vecEndHEdges.push_back(hHEdgeToCheck);
    }
    if( currSeg.has_on( startPt ) )
    {
      ++res.nInx;
      res.bStartOnBndry = true;
      res.vecStrHEdges.push_back(hHEdgeToCheck);
    }
    CGAL::cpp11::result_of<Intersect_2(Segment_2, Segment_2)>::type
      result = intersection( currSeg, querySegment );
    const Point_2* p = nullptr;
    if( result && (p = boost::get<Point_2>(&*result)) )
    {
      if( *p != startPt && *p != endPt )
      {
        ++res.nInx;
        res.vecInxHEdges.push_back( hHEdgeToCheck );
        res.vecInxPts.push_back(*p);
      }
    }
    Point_2 srcPt = ciCurrHedge->source()->point();
    Point_2 dstPt = ciCurrHedge->target()->point();
    nStartPtOnLeft += CGAL::left_turn( srcPt, dstPt, startPt);
    nEndPtOnLeft   += CGAL::left_turn( srcPt, dstPt, endPt ); 
    ++nOfEdges;
  }
  while( ++ciCurrHedge != ciFirstHEdge );
  res.bStartProperInside = nStartPtOnLeft == nOfEdges;
  res.bEndProperInside   = nEndPtOnLeft   == nOfEdges;
  return res;  
}

//----------------------------------------------------------------------------
Halfedge_const_handle
SegmentCheck::getVerticalHEdge( const vector<Halfedge_const_handle>& v )
{
  for( Halfedge_const_handle h : v )
  {
    FT x0 = h->source()->point().x();
    FT x1 = h->target()->point().x();
    if( x0 == x1 )
      return h;
  }
  return v[0];
}

//----------------------------------------------------------------------------
Halfedge_const_handle
SegmentCheck::findNextHEdge( const FaceData&   f,
                             const Segment_2&  querySegment,
                             bool&             bDone,
                             bool&             bResult,
                             FT&               dInxEnd )
{
  cout << "\tInside findNextHEdge"<< endl;
  Halfedge_const_handle hRes; 
  if( !f.hFace->contained() )
  {
    cout << "\tForbidden face. Done"<< endl;
    bDone = true;
    bResult = false;
  }
  else if( f.vecInxPts.size() > 0 && ( f.bEndProperInside || f.bEndOnBndry ) )
  {
    cout << "\tFace with endPt found"<< endl;  
    bDone = true;
    bResult = f.hFace->contained();
  }
  else if(    ( f.bStartProperInside || f.bStartOnBndry ) 
           && ( f.bEndProperInside   || f.bEndOnBndry   ) )
  {
    cout << "\tSegment inclusion"<< endl;  
    bDone = true;
    bResult = f.hFace->contained();  
  }
  else if( 0 != f.vecInxPts.size() )
  {
    cout << "\t Checking proper intersections" << endl;
    Point_2 endPt = querySegment.target(); 
    FT minDist = pointsDistance(f.vecInxPts[0], endPt );
    Point_2 inxPt = f.vecInxPts[0];
    hRes = f.vecInxHEdges[0];
    for( int i = 0; i < f.vecInxPts.size(); ++i )
    {
      FT dPt =  pointsDistance(f.vecInxPts[i], endPt );  
      if( minDist > dPt )
      {
        minDist = dPt;
        hRes    = f.vecInxHEdges[i];
        inxPt   = f.vecInxPts[i];
      }
    }
    dInxEnd = minDist;
    cout << "\t Best intersection ";
    print_point( inxPt );
    cout << "\t on the hedge ";
    print_halfedge( hRes );
    cout << "\t Distance to endPt is " << CGAL::to_double(minDist);
  }
  else if( f.bStartOnBndry )
  {
    dInxEnd = _rodLength;
    hRes = getVerticalHEdge( f.vecStrHEdges );
    cout << "\t Start on boundary case." << endl;
    cout << "\t Result hedge ";
    print_halfedge(hRes);
  }
  else
  {
    cout << "\tUnrecognized situation" << endl;
  }
  return hRes;
}

//----------------------------------------------------------------------------
bool 
SegmentCheck::getCommonFace( Halfedge_const_handle hH1,
                             Halfedge_const_handle hH2,
                             Face_const_handle* pRes )
{
  bool bRes = true;
  if( hH1->face() == hH2->face() )
    *pRes = hH1->face();
  else if( hH1->twin()->face() == hH2->face() )
    *pRes = hH2->face();
  else if( hH1->face() == hH2->twin()->face() )
    *pRes = hH1->face();
  else if( hH1->twin()->face() == hH2->twin()->face() )
    *pRes = hH1->twin()->face();
  else
    bRes = false;
  return bRes;
}

//----------------------------------------------------------------------------
Face_const_handle 
SegmentCheck::chooseFace( Vertex_const_handle hVrtx,
                          const Vector_2&     direction )
{
  Direction_2 theDir(direction);
  cout << "Inside chooseFace with direction";
  print_direction( theDir );

  Arrangement_2::Halfedge_around_vertex_const_circulator ciStart = 
                                                 hVrtx->incident_halfedges();
  Arrangement_2::Halfedge_around_vertex_const_circulator ciCurr = ciStart;
  Arrangement_2::Halfedge_around_vertex_const_circulator ciPrev = ciStart;
  Direction_2 currDir( getSegment( hVrtx->point(), 
                                   ciCurr->source()->point()));
  cout << "\tFirst direction to test";
  print_direction(currDir);
  if( currDir == theDir )
  {
    cout << "Direction overlaps 1st edge" << endl;
    if( ciCurr->face()->contained() )
      return ciCurr->face();
    if( ciCurr->twin()->face()->contained() )
      return ciCurr->twin()->face(); 
  }
  Direction_2 prevDir = currDir;
  ++ciCurr; 
  do
  {
    Direction_2 currDir(getSegment(hVrtx->point(), 
                                   ciCurr->source()->point()));
    cout << "\tTesting between" << endl;
    cout << "\t"; 
    print_direction(prevDir);
    cout << "\t";
    print_direction(currDir);
    if( currDir == theDir )
    {
      cout << "\t Direction overlaps an edge ";
      if( ciCurr->face()->contained() )
        return ciCurr->face();
      if( ciCurr->twin()->face()->contained() )
        return ciCurr->twin()->face();
      else
        return ciCurr->face(); 
    }
    else if( theDir.counterclockwise_in_between( currDir, prevDir ) )
    {
      Face_const_handle hCmnFace;
      if( !getCommonFace( ciCurr, ciPrev, &hCmnFace ) )
        cout << "\tNo common face found" << endl;
      cout << "\t Common face found" << endl;
      return hCmnFace;
    }
    prevDir = currDir;
    ++ciCurr;
    ++ciPrev;
  }
  while( ciPrev != ciStart );
  cout << "\t Search failed" << endl;
  return ciCurr->face();
}


//----------------------------------------------------------------------------
bool
SegmentCheck::isFree( const Point_2& startPt, const Vector_2& direction )
{
  cout << "--------------------------" << endl;
  Point_2 endPt = startPt + (direction * _rodLength);
  Segment_2 querySegment = getSegment(startPt, endPt);
  print_segment(startPt, endPt);

  /* ========================= PLAN B =======================
  vector<CGAL::Object>  vecZoneElems;
  Face_handle     hFace;
  Halfedge_handle hHEdge;
  Landmarks_pl pl( _arr );
  bool bResult = true;
  CGAL::zone( _arr, querySegment, std::back_inserter(vecZoneElems), pl );
  for ( int i = 0; i < vecZoneElems.size(); ++i )
  {
    cout << "\t ----" << endl;  
    if( CGAL::assign( hFace, vecZoneElems[i] ) )
    {
      cout << "\t Face ";
      print_arr_face( hFace );
      if( !hFace->contained() )
      {
        cout << "\t Face isn't contained. Result = false" << endl;
        bResult = false;
        break;
      }
    }
    else if( CGAL::assign( hHEdge, vecZoneElems[i]) )
    {
      cout << "\t Edge ";
      print_halfedge(hHEdge);
      if(    !hHEdge->face()->contained() 
          && !hHEdge->twin()->face()->contained() )
      {
        cout << "\t Both faces aren't contained. Result = false." << endl;
        bResult = false;
        break;
      }
    }
    Vertex_handle hVrtx;
    if( CGAL::assign( hVrtx, vecZoneElems[i] ) )
    {
      cout << "\t Vertex at ";
      print_point( hVrtx->point() ); 
    }
  }
  cout << "Result = " << (bResult ? "true"  : "false" )<< endl;
  return bResult;
  ============================ PLAN B END HERE============================*/

  bool bDone = false;
  bool bResult = true;
  Landmarks_pl pl( _arr );
  CGAL::Object obj = pl.locate(startPt);
  Face_const_handle     hFace;
  Halfedge_const_handle hHEdge;
  Vertex_const_handle   hVrtx;
  while( !bDone )
  {
    if (CGAL::assign(hHEdge, obj))
    {
      cout << "Processing halfedge";
      print_halfedge( hHEdge );
      Face_const_handle hLFace = hHEdge->face();
      Face_const_handle hRFace = hHEdge->twin()->face();
      FaceData lftData = getFaceData( hLFace, querySegment );
      lftData.print();
      FaceData rghData = getFaceData( hRFace, querySegment );
      rghData.print();
      FT dL, dR;
      Halfedge_const_handle hLftHEdge = findNextHEdge( lftData, 
                                                       querySegment, 
                                                       bDone, bResult, 
                                                       dL );
      if( bDone )
        break;
      Halfedge_const_handle hRghHEdge = findNextHEdge( rghData, 
                                                       querySegment, 
                                                       bDone, bResult, 
                                                       dR );
      if( bDone )
        break;

      obj = CGAL::make_object( dL < dR ? hLftHEdge : hRghHEdge );
    }

    //-----------------------------------------------------
    if (CGAL::assign(hFace, obj)) 
    {
      cout << "Handling a face ";
      print_arr_face( hFace );
      FaceData faceData = getFaceData( hFace, querySegment );
      faceData.print();
      FT d;
      hHEdge = findNextHEdge( faceData, querySegment, 
                              bDone, bResult, d );
      if( bDone )
        break;
      obj = CGAL::make_object( hHEdge );
    }
    //------------------------------------------------------
    if( CGAL::assign( hVrtx, obj ) )
    {
      cout << "Handling a vertex" << endl;
      print_point(hVrtx->point());
      hFace = chooseFace( hVrtx, direction );
      obj = CGAL::make_object( hFace );
    }
    cout << "** Here" << endl;
  }//end of while(bDone)

  cout << "bResult = " << ( bResult ? "True" : "False" ) << endl;
  return bResult;
}
//-----------------------------------------------------------------------------
void SegmentCheck::addFrame()
{
  FT d = _rodLength;// * 1.5;

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
  _arr.unbounded_face()->set_contained(true);
  _maxx = mostRight + d;
  _minx = mostLeft  - d;
  _maxy = mostUp    + d;
  _miny = mostDown  - d;
}
//============================== END OF FILE ==================================
// Debug stuff. 
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

void print_halfedge(Halfedge_const_handle hHEdge )
{
  print_segment( hHEdge->source()->point(), hHEdge->target()->point() );
}

void print_point( const Point_2& p)
{
  double x1 = CGAL::to_double( p[0] );
  double y1 = CGAL::to_double( p[1] );
  cout << "(" << x1 << ", " << y1 << ")" << endl;
}

void print_direction( const Direction_2& d )
{
  double x1 = CGAL::to_double( d.dx() );
  double y1 = CGAL::to_double( d.dy() );
  cout << "(" << x1 << ", " << y1 << ")" << endl;
}
