#ifndef H_SEGMENTS_ALGO_H
#define H_SEGMENTS_ALGO_H

#include <vector>
#include "math_lib.h"


typedef enum {Inside, Outside, onAB, onAC, onBC, isA, isB, isC } point_triangle_placement;
inline bool IsOnBoundary(point_triangle_placement in) { return in > Outside; }
inline bool IsOnVertex(point_triangle_placement in) { return in > onBC; }
inline bool IsOnEdges(point_triangle_placement in) { return in > Outside && in < isA; }

point_triangle_placement IsPointInTriangle(const SVec2f& pnt, const SVec2f& A,
                                           const SVec2f& B, const SVec2f& C);


struct exit_type
{
 point_triangle_placement ans;
 point_triangle_placement where_end;
 double cedge_res;
 double ncedge_res;
 exit_type(point_triangle_placement in_ans = Outside, 
           double ce = -1.0, double nce = -1.0)
   : ans(in_ans), cedge_res(ce), ncedge_res(nce)  { }
 
 exit_type(const exit_type& src) 
   : ans(src.ans), where_end(src.where_end), 
     cedge_res(src.cedge_res), ncedge_res(src.ncedge_res) { }
 
 const exit_type& operator=(const exit_type& src)
 {
   ans = src.ans;
   where_end = src.where_end;
   cedge_res = src.cedge_res;
   ncedge_res = src.ncedge_res;
   return * this;
 }
};

exit_type
  TriangleToSegmentRelationShip(const SVec2f& beg, const SVec2f& end, 
                                point_triangle_placement enter_point,
                                const SVec2f& A, const SVec2f& B, const SVec2f& C,
                                double alpha_beg = 0.0, double alpha_end = 1.0);

struct S2Dsegment
{
  SVec2f begin;
  SVec2f end;
  int node_num;
  int edge_num;
  bool IsCodirectional(const S2Dsegment& seg) const 
  {
    return segment_vector() * seg.segment_vector() > 0.0;
  }

  S2Dsegment() : node_num(-1), edge_num(-1) { }
  SVec2f segment_vector() const 
  {
    return end - begin;
  }
  S2Dsegment(const SVec2f& b, const SVec2f & e, int Inode_num = -1, int Iedge_num=-1) : 
      begin(b), end(e), node_num(Inode_num), edge_num(Iedge_num) { }

  S2DLine GetLine() const
  {
    // computing equation parameters, for line involves segment
    double l_a(begin.y - end.y);
    double l_b(-begin.x + end.x);
    return S2DLine(l_a, l_b, -begin.x * l_a - begin.y * l_b);
  }
};


struct IntersectedPair
{
  int n1, n2;
  double u1, u2;
  IntersectedPair() { n1 = n2 = -1;  u1 = u2 = -1.0; }
  IntersectedPair(int in1, double iu1, int in2, double iu2) 
  { 
    n1 = in1;
    n2 = in2;
    u1 = iu1;
    u2 = iu2;
  }
};


int AreSegmentsIntersecting(const SVec2f& s1_beg, const SVec2f& s1_end,
                            const SVec2f& s2_beg, const SVec2f& s2_end,
                            double& u1, double& u2);

double distance_point2segment(const SVec2f& s_beg, const SVec2f& s_end, const SVec2f& pt);

inline int AreSegmentsIntersecting(const SVec2f& s1_beg, const SVec2f& s1_end,
                            const SVec2f& s2_beg, const SVec2f& s2_end,
                            double& u1, double& u2,
                            const double& a_beg, const double& a_end)
{
  int ans(AreSegmentsIntersecting(s1_beg, s1_end, s2_beg, s2_end, u1, u2));
  if(ans == 1)
    return a_beg <= u1 && a_end >= u1 ? 1 : 0;
  else
    return ans;
}

inline int AreSegmentsIntersecting(const S2Dsegment& s1, const S2Dsegment& s2,
                             double& u1, double& u2)
{
  return AreSegmentsIntersecting(s1.begin, s1.end, s2.begin, s2.end, u1, u2);
}



int GetIntersectedPairs(std :: vector<S2Dsegment>& segments1, 
                        std :: vector<S2Dsegment>& segments2,
                        std :: vector<IntersectedPair>& pairs);

#endif
