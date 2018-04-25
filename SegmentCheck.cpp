#include "SegmentCheck.h"

void print_arr_face(Face_const_handle f);
void print_arrangement( const Arrangement_2& arr );
void print_segment( const Point_2& p1, const Point_2& p2 );
void print_segment( const Segment_2& s );
void print_halfedge(Halfedge_const_handle hHEdge );
void print_point( const Point_2& p);

//=============================================================================
void
polygon_split_observer::after_split_face( Face_handle f1, Face_handle f2, bool )
{
  f2->set_contained(f1->contained());
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

  addFrame( rodLength );

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
/*
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
*/
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
Halfedge_const_handle
SegmentCheck::findNextHEdge( Face_const_handle&           hCurrFace,
                             const Point_2&               endPt,
                             const Segment_2&             querySegment,
                             bool&                        bDone,
                             bool&                        bResult,
                             const Halfedge_const_handle* pAnchorHEdge )
{
  ccb_haledge_circulator ciFirstHEdge = hCurrFace->outer_ccb();
  ccb_haledge_circulator ciCurrHedge = ciFirstHEdge;
  Halfedge_const_handle hResult;
  cout << "Inside findNextHEdge" << endl;
  do
  {
    Halfedge_const_handle hHEdgeToCheck = ciCurrHedge;
    if( nullptr != pAnchorHEdge && *pAnchorHEdge == hHEdgeToCheck )
    {
      cout << "Skipping anchor hedge " << endl;
      continue;
    }

    Segment_2 currSeg = getSegment(hHEdgeToCheck);
    cout << "\t Current segment ";
    print_segment( currSeg );
    if (currSeg.has_on(endPt))
    {
      // The end point is on the hHEdgeToCheck
      cout << "\t End point is on the current segment ";
      bResult = true;
      bDone = true;
      hResult = hHEdgeToCheck;
      break;
    }
    CGAL::cpp11::result_of<Intersect_2(Segment_2, Segment_2)>::type
        result = intersection( currSeg, querySegment );
    if (result)
    {
      cout << " Intersection found ";
      if (const Segment_2 *s = boost::get<Segment_2>(&*result))
      {
        // Two segments overlap. But the end point isn't there.
        print_segment( *s );
        continue;
      } else
      {
        // This is the next hedge to jump to
        cout << " with a point" << endl;
        hResult = hHEdgeToCheck;
        break;
      }
    }
  }
  while( ++ciCurrHedge != ciFirstHEdge );
  return hResult;
}

//----------------------------------------------------------------------------
Face_const_handle 
SegmentCheck::chooseFace( Halfedge_const_handle hHEdge,
                          const Point_2&        endPt )
{
  Point_2 srcPt = hHEdge->source()->point();
  Point_2 dstPt = hHEdge->target()->point();
  Point_2 leftPt = hHEdge->next()->target()->point();
  Point_2 rightPt = hHEdge->twin()->prev()->source()->point();
  cout << "Choosing face" << endl;
  cout << "\t srcPt";
  print_point( srcPt );
  cout << "\t dstPt";
  print_point( dstPt );
  cout << "\tleftPt";
  print_point(leftPt);
  cout << "\trightPt";
  print_point(rightPt);
  Face_const_handle hLeftFace = hHEdge->face();
  Face_const_handle hRightFace = hHEdge->twin()->face();
  bool bLft = CGAL::left_turn(    srcPt, dstPt, leftPt  );
  bool bLftColl = CGAL::collinear(srcPt, dstPt, leftPt  );
  bool bRgh = CGAL::left_turn(    srcPt, dstPt, rightPt );
  bool bRghColl = CGAL::collinear(srcPt, dstPt, rightPt );
  bool bEnd = CGAL::left_turn(    srcPt, dstPt, endPt   );
  bool bEndColl = CGAL::collinear(srcPt, dstPt, endPt   );
  cout << "bLft     = " << bLft     << " bRgh     = " << bRgh     << " bEnd     = " << bEnd << endl;
  cout << "bLftColl = " << bLftColl << " bRghColl = " << bRghColl << " bEndColl = " << bEndColl << endl;
  if( bLftColl && bEndColl && !bRghColl )
    return hLeftFace;
  if( bRghColl && bEndColl && !bLftColl )
    return hRightFace;
  if( bLft == bEnd && bRgh != bEnd )
    return hLeftFace;
  if( bRgh == bEnd && bLft != bEnd )
    return hRightFace;
}
  
//----------------------------------------------------------------------------
Face_const_handle 
SegmentCheck::chooseFace( Vertex_const_handle hVrtx,
                          const Vector_2&     direction )
{
  Direction_2 theDir(direction);
  Arrangement_2::Halfedge_around_vertex_const_circulator ciStart = 
                                                 hVrtx->incident_halfedges();
  Arrangement_2::Halfedge_around_vertex_const_circulator ciPrev = ciStart; 
  Arrangement_2::Halfedge_around_vertex_const_circulator ciCurr = ciStart;
  Direction_2 currDir( getSegment( hVrtx->point(), 
                                   ciCurr->source()->point()));
  if( currDir == theDir )
  {
    cout << "Direction overlaps 1st edge" << endl;
    if( ciCurr->face()->contained() )
      return ciCurr->face();
    if( ciCurr->twin()->face()->contained() )
      return ciCurr->twin()->face(); 
  }
  ++ciCurr;
  Direction_2 prevDir( getSegment( hVrtx->point(), 
                                   ciPrev->source()->point()) );
  do
  {
    Direction_2 currDir(getSegment(hVrtx->point(), 
                                   ciCurr->source()->point()));
    if( currDir == theDir || 
        theDir.counterclockwise_in_between( currDir, prevDir ) )
    {
      if( ciCurr->face()->contained() )
        return ciCurr->face();
      if( ciCurr->twin()->face()->contained() )
        return ciCurr->twin()->face(); 
    }
    prevDir = currDir;
  }
  while( ++ciCurr != ciStart );
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
  Face_const_handle hCurrFace;
  Halfedge_const_handle hHEdge;
  Vertex_const_handle hVrtx;
  while( !bDone )
  {
    // Check if it's a halfedge
    if (CGAL::assign(hHEdge, obj))
    {
      cout << "Processing halfedge";
      print_halfedge( hHEdge );
      hCurrFace = chooseFace( hHEdge, endPt );
      if( !hCurrFace->contained() )
      {
        cout << "Got to forbidden face"<< endl;
        bResult = false;
        bDone = true;
        break;
      }
      if (hCurrFace->is_unbounded())
      {
        cout << "Got to unbounded face" << endl;
        bResult = true;
        bDone = true;
        break;
      }
      Halfedge_const_handle hNextHEdge = findNextHEdge(hCurrFace, endPt, querySegment, bDone, bResult, &hHEdge);
      if( hNextHEdge == hHEdge )
      {
        cout << "No next intesection found."<< endl;
        bResult = true;
        bDone = true;
        break;
      }
      obj = CGAL::make_object( hNextHEdge );
    } // end of if "obj is an Hedge?"

    if (CGAL::assign(hCurrFace, obj)) //if obj is face
    {
      cout << "Handling a face" << endl;
      if (!hCurrFace->contained())
      {
        cout << "Got to forbidden face" << endl;
        bDone = true;
        bResult = false;
        break;
      }
      hHEdge = findNextHEdge( hCurrFace, endPt, querySegment, bDone, bResult, nullptr );
      obj = CGAL::make_object( hHEdge );
      continue;
    }

    if( CGAL::assign( hVrtx, obj ) )
    {
      cout << "Handling a vertex" << endl;
      print_point(hVrtx->point());
      hCurrFace = chooseFace( hVrtx, direction );
      obj = CGAL::make_object( hCurrFace );
    }
    cout << "Here" << endl;
  }//end of while(bDone)

  cout << "bResult = " << ( bResult ? "True" : "False" ) << endl;
  return bResult;
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
  _arr.unbounded_face()->set_contained(true);
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


