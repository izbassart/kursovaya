#include "segment_algo.h"
#include "math.h"


S2DLine :: S2DLine(const S2Dsegment& seg)
{
  // computing equation parameters, for line involves segment
  a_ = seg.begin.y - seg.end.y;
  b_ = -seg.begin.x + seg.end.x;
  c_ = -seg.begin.x * a_ - seg.begin.y * b_;
}

S2DLine :: S2DLine(const SVec2f& beg, const SVec2f& end)
{
  a_ = beg.y - end.y;
  b_ = -beg.x + end.x;
  c_ = -beg.x * a_ - beg.y * b_;
}

S2DLine :: S2DLine(const SVec3f& beg, const SVec3f& end)
{
  a_ = beg.y - end.y;
  b_ = -beg.x + end.x;
  c_ = -beg.x * a_ - beg.y * b_;
}



// В этой функции считаем, что ребро 100% входит в треугольник
// Выводы считаются для ребер AB, AC и BC.
// относительно заданных A,B и С.
exit_type
  TriangleToSegmentRelationShip(const SVec2f& beg, const SVec2f& end, 
                                point_triangle_placement enter_point, 
                                const SVec2f& A, const SVec2f& B, const SVec2f& C,
                                double alpha_beg, double alpha_end)
{
  S2DLine line(beg, end);
  exit_type answer;
   //Начала и концы отрезка ребра
  SVec2f beg_pnt((1.0 - alpha_beg) * beg + alpha_beg * end);
  SVec2f end_pnt((1.0 - alpha_end) * beg + alpha_end * end);

  answer.where_end = IsPointInTriangle(end_pnt, A, B, C);
  if(answer.where_end == Inside)
  {
    answer.ans = Inside;
  }
  else if(answer.where_end == Outside)
  {
    // Анализ на пересесение и касание 
    if(IsOnVertex(enter_point))
    {
      if(enter_point == isA)
      {
        AreSegmentsIntersecting(beg, end, B, C, 
                                answer.cedge_res, answer.ncedge_res);
        if(answer.ncedge_res == 0.0)
        {
          answer.ans = isB;
        }
        else if(answer.ncedge_res == 1.0)
        {
          answer.ans = isC;
        }
        else
        {
          answer.ans = onBC;
        }
      }
      else if(enter_point == isB)
      {
        AreSegmentsIntersecting(beg, end, A, C, 
                                answer.cedge_res, answer.ncedge_res);
        if(answer.ncedge_res == 0.0)
        {
          answer.ans = isA;
        }
        else if(answer.ncedge_res == 1.0)
        {
          answer.ans = isC;
        }
        else
        {
          answer.ans = onAC;
        }
      }
      else if(enter_point == isC)
      {
        AreSegmentsIntersecting(beg, end, A, B, 
                                answer.cedge_res, answer.ncedge_res);
        if(answer.ncedge_res == 0.0)
        {
          answer.ans = isA;
        }
        else if(answer.ncedge_res == 1.0)
        {
          answer.ans = isB;
        }
        else
        {
          answer.ans = onAB;
        }
      }
    }
    else
    {
      double out_fAB(-1.0), out_fAC(-1.0), out_fBC(-1.0);
      double self_AB(-1.0), self_AC(-1.0), self_BC(-1.0);
      bool ansAB(AreSegmentsIntersecting(beg, end, A, B, out_fAB, self_AB,
                                         alpha_beg, alpha_end) == 1);
      bool ansAC(AreSegmentsIntersecting(beg, end, A, C, out_fAC, self_AC,
                                         alpha_beg, alpha_end) == 1);
      bool ansBC(AreSegmentsIntersecting(beg, end, B, C, out_fBC, self_BC,
                                         alpha_beg, alpha_end) == 1);
     
      if(enter_point == Inside)
      {
        if(ansAB && ansAC)
        {
          answer.ans = isA;
          answer.cedge_res = out_fAC;
        }
        else if(ansBC && ansAC)
        {
          answer.ans = isC;
          answer.cedge_res = out_fAC;
        }
        else if(ansAB && ansBC)
        {
          answer.ans = isB;
          answer.cedge_res = out_fAB;
        }
        else if(ansAB)
        {
          answer.ans = onAB;
          answer.cedge_res = out_fAB;
          answer.ncedge_res = self_AB;
        }
        else if(ansAC)
        {
          answer.ans = onAC;
          answer.cedge_res = out_fAC;
          answer.ncedge_res = self_AC;
        } 
        else if(ansBC)
        {
          answer.ans = onBC;
          answer.cedge_res = out_fBC;
          answer.ncedge_res = self_BC;
        }
      }
      else if(IsOnEdges(enter_point))
      {
        S2DLine line(beg, end);
        bool accrossA(line.sign(A) == 0);
        bool accrossB(line.sign(B) == 0);
        bool accrossC(line.sign(C) == 0);

        if(enter_point == onAB)
        {
          if(accrossA)
          {
            if(ansBC)
            {
              answer.ans = isB;
              answer.cedge_res = out_fBC;
            }
            else if(ansAC)
            {
              answer.ans = isA;
              answer.cedge_res = out_fAC;
            } 
          } 
          else if(accrossC)
          {
            answer.ans = isC;
            answer.cedge_res = out_fBC;
          }
          else
          {
            if(ansBC)
            {
              answer.ans = onBC;
              answer.cedge_res = out_fBC;
              answer.ncedge_res = self_BC;
            }
            else if(ansAC)
            {
              answer.ans = onAC;
              answer.cedge_res = out_fAC;
              answer.ncedge_res = self_AC;
            }
          }
        }
        else if(enter_point == onAC)
        {
          if(accrossA)
          {
            if(ansBC)
            {
              answer.ans = isC;
              answer.cedge_res = out_fBC;
            }
            else if(ansAB)
            {
              answer.ans = isA;
              answer.cedge_res = out_fAB;
            } 
          } 
          else if(accrossB)
          {
            answer.ans = isB;
            answer.cedge_res = out_fBC;
          }
          else 
          {
            if(ansBC)
            {
              answer.ans = onBC;
              answer.cedge_res = out_fBC;
              answer.ncedge_res = self_BC;
            }
            else if(ansAB)
            {
              answer.ans = onAB;
              answer.cedge_res = out_fAB;
              answer.ncedge_res = self_AB;
            }
          }
        }
        else if(enter_point == onBC)
        {
          if(accrossB)
          {
            if(ansAC)
            {
              answer.ans = isC;
              answer.cedge_res = out_fAC;
            }
            else if(ansAB)
            {
              answer.ans = isB;
              answer.cedge_res = out_fAB;
            } 
          } 
          else if(accrossA)
          {
            answer.ans = isA;
            answer.cedge_res = out_fAB;
          }
          else 
          {
            if(ansAC)
            {
              answer.ans = onAC;
              answer.cedge_res = out_fAC;
              answer.ncedge_res = self_AC;
            }
            else if(ansAB)
            {
              answer.ans = onAB;
              answer.cedge_res = out_fAB;
              answer.ncedge_res = self_AB;
            }
          }
        }
      }
      else if(enter_point == Outside)
      {
        answer.ans = Outside;
        if(ansAB && self_AB > 0.0)
        {
          answer.ans = onAB;
          answer.ncedge_res = self_AB;
          answer.cedge_res = out_fAB;
        }

        if(ansBC && self_BC > 0.0 && self_BC > answer.cedge_res)
        {
          answer.ans = onBC;
          answer.ncedge_res = self_BC;
          answer.cedge_res = out_fBC;
        }

        if(ansAC && self_AC > 0.0 && self_AC > answer.cedge_res)
        {
          answer.ans = onAC;
          answer.ncedge_res = self_AC;
          answer.cedge_res = out_fAC;
        }
      }
    }
  }
  else if(IsOnVertex(answer.where_end))
  {
    answer.ans = answer.where_end;
    answer.cedge_res = 1.0;
  }
  else if(IsOnEdges(answer.where_end))
  {
    answer.cedge_res = 1.0;
    answer.ans = answer.where_end;
    if(answer.where_end == onAB)
      answer.ncedge_res = (end_pnt - A).length() / (B - A).length(); 
    else if(answer.where_end == onAC)
      answer.ncedge_res = (end_pnt - A).length() / (C - A).length(); 
    else if(answer.where_end == onBC)
      answer.ncedge_res = (end_pnt - B).length() / (C - B).length(); 
  }
  return answer;
}

point_triangle_placement IsPointInTriangle(const SVec2f& pnt, const SVec2f& A,
                                           const SVec2f& B, const SVec2f& C)
{
  S2DLine AB(A,B), AC(A,C), BC(B,C);
  int sAB(AB.sign(pnt) * AB.sign(C));
  int sAC(AC.sign(pnt) * AC.sign(B));
  int sBC(BC.sign(pnt) * BC.sign(A));
  
  if(sAB == 1 && sAC == 1 && sBC== 1)
    return Inside;
   if(sAB == 0 && sAC == 1 && sBC== 1)
     return onAB;
   if(sAB == 1 && sAC == 0 && sBC== 1)
     return onAC;
   if(sAB == 1 && sAC == 1 && sBC== 0)
     return onBC;
  else if(pnt == A)
    return isA;
  else if(pnt == B)
    return isB;
  else if(pnt == C)
    return isC;
  else
    return Outside;
}

// return_type 0 - non-intersected
//             1 - intersected
//             2 - on the same line

int AreSegmentsIntersecting(const SVec2f& s1_beg, const SVec2f& s1_end,
                            const SVec2f& s2_beg, const SVec2f& s2_end,
                            double& u1, double& u2)
{

  double znam((s2_beg.y - s2_end.y) * (s1_end.x - s1_beg.x) - 
              (s2_beg.x - s2_end.x) * (s1_end.y - s1_beg.y));

  double chisl((s2_beg.y - s2_end.y) * (s2_beg.x - s1_beg.x) -
               (s2_beg.x - s2_end.x) * (s2_beg.y - s1_beg.y));

  if(fabs(znam) < EPS_PRECISION)
  {
    if(fabs(chisl) < EPS_PRECISION)
      return 2;
    else
      return 0;
  }

  u1 = chisl / znam;

  u2 = ((s1_end.x - s1_beg.x) * (s2_beg.y - s1_beg.y) -
            (s1_end.y - s1_beg.y) * (s2_beg.x - s1_beg.x)) / znam;

  return (u1 >= 0.0) && (u1 <= 1.0) && (u2 >= 0.0) && (u2 <= 1.0) ? 1 : 0;
}


int GetIntersectedPairs(std :: vector<S2Dsegment>& segments1, 
                        std :: vector<S2Dsegment>& segments2,
                        std :: vector<IntersectedPair>& pairs)
{
  pairs.clear();
  int n, q;
  double u1, u2;
  for(n = 0; n < int(segments1.size()); n++)
    for(q = 0; q < int(segments2.size()); q++)
    {
      if(AreSegmentsIntersecting(segments1[n], segments2[q], u1, u2))
        pairs.push_back(IntersectedPair(n, u1, q, u2));
    }
                            
  return int(pairs.size());
}

double distance_point2segment(const SVec2f& s_beg, const SVec2f& s_end, const SVec2f& pt)
{
  SVec2f AB(s_end - s_beg);
  SVec2f AC(pt - s_beg);
  SVec2f BC(pt - s_end);
  double lab = AB.length();
  // проверка: тупоугольный ли  трегульник
  if( (AB * AC < 0) || (AB * BC > 0) || (lab == 0.0))
    return std :: min(AC.length(), BC.length());
  else
  {
    return fabs(AC.x * BC.y - AC.y * BC.x) / lab;
  }
}

