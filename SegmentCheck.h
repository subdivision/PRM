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
class polygon_split_observer : public CGAL::Arr_observer<Arrangement_2> {
    void after_split_face(Face_handle f1, Face_handle f2, bool) override;
};

//=============================================================================
class FaceData {
public:
    FaceData();

    void reset();

    bool hasData();

    void print();

    int nInx;
    bool bStartProperInside;
    bool bStartOnBndry;
    bool bEndProperInside;
    bool bEndOnBndry;
    Face_const_handle hFace;
    vector<Halfedge_const_handle> vecStrHEdges;
    vector<Halfedge_const_handle> vecEndHEdges;
    vector<Halfedge_const_handle> vecInxHEdges;
    vector<Point_2> vecInxPts;
};

//=============================================================================
class SegmentCheck {
private:

    FT _rodLength;
    FT _maxx, _minx, _maxy, _miny;
    Polygon_set_2 _freeSpace;
    Arrangement_2 _arr;
    Landmarks_pl pl;


    void verticalDecomposition(Kernel &ker);

    void addVerticalSegment(Vertex_handle v,
                            CGAL::Object obj,
                            Kernel &ker);

    bool pointProperlyOut(const Point_2 &pt);

    bool pointOnTheOutBndry(const Point_2 &pt);

    FaceData getUnboundedFaceData(Face_const_handle &hFace,
                                  const Segment_2 &querySegment);

    FaceData getFaceData(Face_const_handle &hFace,
                         const Segment_2 &querySegment);

    FT pointsDistance(Point_2 a, Point_2 b) const;

    Segment_2 getSegment(Halfedge_const_handle edge) const;

    Segment_2 getSegment(Point_2 a, Point_2 b) const;


    Halfedge_const_handle
    getVerticalHEdge(const vector<Halfedge_const_handle> &v);

    Halfedge_const_handle
    findNextHEdge(const FaceData &faceData,
                  const Segment_2 &querySegment,
                  bool &bDone,
                  bool &bResult,
                  FT &dInxEnd);


    Face_const_handle chooseFace(Vertex_const_handle hVrtx,
                                 const Vector_2 &direction);

    bool getCommonFace(Halfedge_const_handle hH1,
                       Halfedge_const_handle hH2,
                       Face_const_handle *pRes);

    void addFrame();

public:
    SegmentCheck(const FT &rodLength, const vector<Polygon_2> &obstacles);

    bool isFree(const Point_2 &startPt, const Vector_2 &direction);
};

#endif 
