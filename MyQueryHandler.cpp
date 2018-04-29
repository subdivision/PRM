#include "MyQueryHandler.h"

using namespace std;

void polygon_split_observer::after_split_face(Face_handle f1, Face_handle f2, bool) {
    f2->set_contained(f1->contained());
}

MyQueryHandler::MyQueryHandler(const FT &rodLength, const vector<Polygon_2> &obstacles) :
        _rodLength(rodLength) {
    Polygon_set_2 freeSpace;
    freeSpace.join(obstacles.begin(), obstacles.end());

    //complement to get the free space
    freeSpace.complement();

    //set the free space as arrangment
    _arr = freeSpace.arrangement();

    addFrame();

    //ensure that when face split two side safe their property (inside/outside)
    Polygon_set_2::Traits_2 traits;
    polygon_split_observer observer;
    observer.attach(_arr);
    Kernel *ker = &traits;
    verticalDecomposition(*ker);
    observer.detach();
    _arr.unbounded_face()->set_contained(true);

    pl.attach(_arr);
}

bool MyQueryHandler::_isLegalConfiguration(const Point_2 &point, const Vector_2 &direction, const double rotation) {
    Point_2 endPt = point + (direction * _rodLength);
    Segment_2 querySegment(point, endPt);

    vector<CGAL::Object> vecZoneElems;
    Face_handle hFace;
    Halfedge_handle hHEdge;

    CGAL::zone(_arr, querySegment, std::back_inserter(vecZoneElems), pl);
    for (int i = 0; i < vecZoneElems.size(); ++i) {
        if (CGAL::assign(hFace, vecZoneElems[i]) && !hFace->contained())
            return false;
    }
    return true;
}



void MyQueryHandler::addVerticalSegment(Vertex_handle v, CGAL::Object obj, Kernel &ker) {
    X_monotone_curve_2 seg;
    Vertex_const_handle vh;
    Halfedge_const_handle hh;
    Face_const_handle fh;
    Vertex_handle v2;

    if (CGAL::assign(vh, obj)) {
        // The given feature is a vertex.
        seg = X_monotone_curve_2(v->point(), vh->point());
        v2 = _arr.non_const_handle(vh);
    } else if (CGAL::assign(hh, obj)) {
        // The given feature is a halfedge.
        if (hh->is_fictitious()) //We ignore fictitious halfedges.
        {
            return;
        }

        // Check whether v lies in the interior of the x-range of the edge (in
        // which case this edge should be split).
        const typename Kernel::Compare_x_2 cmp_x = ker.compare_x_2_object();
        if (cmp_x(v->point(), hh->target()->point()) == CGAL::EQUAL) {
            // In case the target of the edge already has the same x-coordinate as
            // the vertex v, just connect these two vertices.
            seg = X_monotone_curve_2(v->point(), hh->target()->point());
            v2 = _arr.non_const_handle(hh->target());
        } else {
            // Compute the vertical projection of v onto the segment associated
            // with the halfedge. Split the edge and connect v with the split point.
            Line_2 Line;
            Line_2 supp_line(hh->source()->point(), hh->target()->point());
            Line_2 vert_line(v->point(), Point_2(v->point().x(), v->point().y() + 1));
            Point_2 point;
            CGAL::assign(point, ker.intersect_2_object()(supp_line, vert_line));
            seg = X_monotone_curve_2(v->point(), point);
            _arr.split_edge(_arr.non_const_handle(hh),
                            X_monotone_curve_2(hh->source()->point(), point),
                            X_monotone_curve_2(point, hh->target()->point()));
            v2 = _arr.non_const_handle(hh->target());
        }
    } else {
        // Ignore faces and empty objects.
        return;
    }

    // Add the vertical segment to the arrangement using its two end vertices.
    _arr.insert_at_vertices(seg, v, v2);
}

void MyQueryHandler::verticalDecomposition(Kernel &ker) {
    typedef pair<Vertex_const_handle, pair<CGAL::Object, CGAL::Object>> Vd_entry;

    // For each vertex in the arrangment, locate the feature that lies
    // directly below it and the feature that lies directly above it.
    list<Vd_entry> vd_list;
    CGAL::decompose(_arr, back_inserter(vd_list));

    // Go over the vertices (given in ascending lexicographical xy-order),
    // and add segements to the feautres below and above it.
    const typename Kernel::Equal_2 equal = ker.equal_2_object();
    typename list<Vd_entry>::iterator it, prev = vd_list.end();
    for (it = vd_list.begin(); it != vd_list.end(); ++it) {
        // If the feature above the previous vertex is not the current vertex,
        // add a vertical segment to the feature below the vertex.
        Vertex_const_handle v;
        if ((prev == vd_list.end()) ||
            !CGAL::assign(v, prev->second.second) ||
            !equal(v->point(), it->first->point())) {
            addVerticalSegment(_arr.non_const_handle(it->first),
                               it->second.first, ker);
        }
        // Add a vertical segment to the feature above the vertex.
        addVerticalSegment(_arr.non_const_handle(it->first),
                           it->second.second, ker);
        prev = it;
    }
}

void MyQueryHandler::addFrame() {
    FT d = _rodLength;// * 1.5;

    Arr_VrtxCIter iVrtx = _arr.vertices_begin();
    FT mostLeft = iVrtx->point().x();
    FT mostRight = iVrtx->point().x();
    FT mostUp = iVrtx->point().y();
    FT mostDown = iVrtx->point().y();

    for (; iVrtx != _arr.vertices_end(); ++iVrtx) {
        if (iVrtx->point().x() < mostLeft)
            mostLeft = iVrtx->point().x();
        if (iVrtx->point().x() > mostRight)
            mostRight = iVrtx->point().x();
        if (iVrtx->point().y() < mostDown)
            mostDown = iVrtx->point().y();
        if (iVrtx->point().y() > mostUp)
            mostUp = iVrtx->point().y();
    }

    Point_2 upperLeft(mostLeft - d, mostUp + d),
            upperRight(mostRight + d, mostUp + d),
            lowerRight(mostRight + d, mostDown - d),
            lowerLeft(mostLeft - d, mostDown - d);
    Segment_2 upperBound(upperLeft, upperRight),
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
