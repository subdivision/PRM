//
// Created by t-idkess on 08-Apr-18.
//

#ifndef RODQUERY_CGAL_DEFINES_H
#define RODQUERY_CGAL_DEFINES_H

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_set_2.h>
#include <CGAL/Point_set_2.h>
#include <CGAL/basic.h>
#include <CGAL/intersections.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/graph_traits_Dual_Arrangement_2.h>
#include <CGAL/Arr_face_index_map.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Polygon_with_holes_2.h>

typedef typename CGAL::Gmpq Number_type;
typedef typename CGAL::Cartesian<Number_type> Kernel;
typedef typename Kernel::FT FT;
typedef typename Kernel::Point_2 Point_2;
typedef typename Kernel::Vector_2 Vector_2;
typedef typename Kernel::Segment_2 Segment_2;
typedef typename Kernel::Direction_2 Direction_2;
typedef typename Kernel::Line_2 Line_2;
typedef typename CGAL::Circle_2<Kernel> Circle_2;
typedef typename CGAL::Point_set_2<Kernel> Point_set;
typedef typename CGAL::Polygon_2<Kernel> Polygon_2;
typedef CGAL::Point_set_2<Kernel>::Vertex_handle Point_set_Vertex_handle;
typedef typename CGAL::Polygon_set_2<Kernel> Polygon_set_2;
typedef Polygon_set_2::Arrangement_2 Arrangement_2;

typedef typename Arrangement_2::Face_handle                     Face_handle;
typedef typename Arrangement_2::Vertex_handle                   Vertex_handle;
typedef typename Arrangement_2::Vertex_const_handle             Vertex_const_handle;
typedef typename Arrangement_2::Face_const_handle               Face_const_handle;
typedef typename Arrangement_2::Face_handle                     Face_handle;
typedef typename Arrangement_2::Halfedge_handle                 Halfedge_handle;
typedef typename Arrangement_2::Halfedge_const_handle           Halfedge_const_handle;
typedef typename Arrangement_2::X_monotone_curve_2              X_monotone_curve_2;
typedef typename Arrangement_2::Vertex_const_iterator           Arr_VrtxCIter;
typedef typename Arrangement_2::Ccb_halfedge_const_circulator   ccb_haledge_circulator;
typedef typename CGAL::Arr_landmarks_point_location<Arrangement_2> Landmarks_pl;

typedef Kernel::Intersect_2 Intersect_2;

#endif //RODQUERY_CGAL_DEFINES_H
