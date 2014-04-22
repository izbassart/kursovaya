#include "./math_lib.h"
#include "segment_algo.h"
#include <set>
#include <stack>
#include <algorithm>

double SVec2f :: org_x = 0.0;
double SVec2f :: org_y = 0.0;

Complex_SVec3f operator*(const double& f, const Complex_SVec3f& vec)
{
  return Complex_SVec3f(f * vec.re, f * vec.im);
}

Complex_SVec3f operator^(const SVec3f& f, const Complex_SVec3f& vec)
{
  return Complex_SVec3f(f ^ vec.re, f ^ vec.im);
}

complex_double operator*(const SVec3f& f, const Complex_SVec3f& vec)
{
  return complex_double(f * vec.re, f * vec.im);
}

SVec3f :: SVec3f(const SVec4f & v)
{
  x = v.x / v.t;
  y = v.y / v.t;
  z = v.z / v.t;
}

std :: string SVec3f :: get_string() const
{
  std :: stringstream str;
  str << x << ' ' << y << ' ' << z;
  return str.str();
}

void SVec3f :: set_from_string(const std :: string& str)
{
  std :: stringstream in_str(str);
  in_str >> x >> y >> z;
}

double mixed_product(const SVec3f& a, const SVec3f& b, const SVec3f& c) 
{
	return a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) + 
		a.z * (b.x * c.y - b.y * c.x);
}

std :: ostream& operator<<(std :: ostream& st, const SVec3f& vec)
{
    st << vec.x << ' ' << vec.y << ' ' << vec.z;
	return st;
}

std :: istream& operator>>(std :: istream& st, SVec3f& vec)
{
  st >> vec.x >> vec.y >> vec.z;
	return st;
}

SVec3f operator*(const double& f, const SVec3f& r) 
{
  return r * f;
}

std :: ostream& operator<<(std :: ostream& st, const SMatrix3f& matr)
{
  for(size_t i(0); i < 3; i++)
	{
		st << matr(i,0) << ' ' << matr(i,1) << ' ' << matr(i,2) << std :: endl;
	}
	return st;
}

SVec3f operator*(const SVec3f& vec, const SMatrix3f& matr)
{
   SVec3f ret_vec;
	for(size_t i = 0; i < 3; i++)
	{
	double sum(0.0);
	for(int q(0); q < 3; q++)
		sum += vec.elem[q] * matr(i,q);

	ret_vec.elem[i] = sum;
	}
	return ret_vec;
}

std :: ostream& operator<<(std :: ostream& st, const SVec4f& vec)
{
    st << vec.x << ' ' << vec.y << ' ' << vec.z << ' ' << vec.t;
	return st;
}

SVec4f operator*(const double& f, const SVec4f& r) 
{
  return r * f;
}

std :: ostream& operator<<(std :: ostream& st, const SMatrix4f& matr)
{
  for(size_t i(0); i < 4; i++)
	{
		st << matr(i,0) << ' ' << matr(i,1) << ' ' << matr(i,2) << ' ' << matr(i,3) << std :: endl;
	}
	return st;
}

SVec4f operator*(const SVec4f& vec, const SMatrix4f& matr)
{
   SVec4f ret_vec;
	for(size_t i = 0; i < 4; i++)
	{
	double sum(0.0);
	for(int q(0); q < 4; q++)
		sum += vec.elem[q] * matr(i,q);

	ret_vec.elem[i] = sum;
	}
	return ret_vec;
}

bool SMatrix3f :: Solve(const SVec3f& in, SVec3f& out) const
{
 double d(Determinant());
 if(fabs(d) > EPS_PRECISION)
 {
   double d1(in[0] * (m_data[4] * m_data[8] - m_data[7] * m_data[5]) 
           + in[1] * (m_data[7] * m_data[2] - m_data[1] * m_data[8])
           + in[2] * (m_data[1] * m_data[5] - m_data[4] * m_data[2]));

   double d2(in[0] * (m_data[6] * m_data[5] - m_data[3] * m_data[8]) 
           + in[1] * (m_data[0] * m_data[8] - m_data[2] * m_data[6])
           + in[2] * (m_data[2] * m_data[3] - m_data[0] * m_data[5]));

   double d3(in[0] * (m_data[3] * m_data[7] - m_data[4] * m_data[6]) 
           + in[1] * (m_data[1] * m_data[6] - m_data[0] * m_data[7])
           + in[2] * (m_data[0] * m_data[4] - m_data[1] * m_data[3]));


   out[0] = d1 / d;
   out[1] = d2 / d;
   out[2] = d3 / d;
   return true;
 }
 return false;
}


OrthoProjector :: OrthoProjector(const SVec3f& ld, const SVec3f& up)
{
  e1 = (ld ^ up);
  if(e1.length() < EPS_PRECISION)
    e1 = ld ^ SVec3f(up.y, -up.x, 0.0);
  e2 = (e1 ^ ld);
  e1.Normalize(); 
  e2.Normalize();
}

SVec3f OrthoProjector :: GetProjection(const SVec3f& vec) const
{
  SVec3f ret = (vec * e1) * e1 + (vec * e2) * e2;
  return ret;
}

SMatrix3f OrthoProjection(const SVec3f& vec)
{
  SMatrix3f matrix;
  int i, j1, j2;
  
  if(fabs(vec.x) > fabs(vec.y))
  {
    if(fabs(vec.x) > fabs(vec.z))
    {
      i = 0;
      j1 = 1;
      j2 = 2;
    }
    else
    {
      i = 2;
      j1 = 0;
      j2 = 1;
    }
  }
  else
  {
    if(fabs(vec.y) > fabs(vec.z))
    {
      i = 1;
      j1 = 0;
      j2 = 2;
    }
    else
    {
      i = 2;
      j1 = 0;
      j2 = 1;
    }
  }

  matrix.SetZero();
  if(fabs(vec[i]) > EPS_PRECISION)
  {
    SVec3f f1(0.0,0.0,0.0), f2(0.0,0.0,0.0);
    f1[i] = -vec[j1] / vec[i]; f1[j1] = 1.0; 
    f2[i] = -vec[j2] / vec[i]; f2[j2] = 1.0; 
    f1.Normalize(); 
    f2 = f2 - (f2 * f1) * f1; 
    f2.Normalize();
    for(int q(0); q < 3; q++)
    {
      matrix(0, q) = f1[q];
      matrix(1, q) = f2[q];
    }
  }

  return matrix;
}

 SMatrix3f OrthoProjection(const SVec3f& vec, const SVec3f& up)
 {
  SMatrix3f matrix;
  matrix.SetZero();
  SVec3f f1(vec ^ up);
  if(f1.length() < EPS_PRECISION)
    return OrthoProjection(vec);
  SVec3f f2(f1 ^ vec);
  f1.Normalize(); 
  f2.Normalize();
  for(int q(0); q < 3; q++)
  { 
    matrix(0, q) = f2[q];
    matrix(1, q) = f1[q];
  }
  return matrix;
 }


 
 SPlate :: SPlate(const SVec3f& a, const SVec3f& b, const SVec3f& c) 
   : v_a(a), v_b(b), v_c(c)
 {
   normal.x = (a.y - b.y) * (b.z - c.z) - (a.z - b.z) * (b.y - c.y);
   normal.y = -(a.x - b.x) * (b.z - c.z) + (a.z - b.z) * (b.x - c.x);
   normal.z = (a.x - b.x) * (b.y - c.y) - (a.y - b.y) * (b.x - c.x);
   D = - normal * a;
 }

 bool SPlate :: InfinitePlateIntersection(const SVec3f& beg, const SVec3f& end)
 {
  return sign(beg) * sign(end) < 0;
 }

 bool SPlate :: GetPlateIntersection(const SVec3f& beg, const SVec3f& end, SVec3f& res) const
 {
   SVec3f vec(end - beg);
   double znam(vec * normal);
   if(!same_val(znam, 0))
   {
     double tau(-((beg * normal) + D) / znam);
     if((tau >= 0) && (tau <= 1.0))
     {
       res = beg + tau * vec;
       return true;
     }
   }
   return false;
 }

 // in case of intersection returs the ray parameter tau
 double SPlate :: GetRayPlateIntersection(const SVec3f& beg, const SVec3f& dir, SVec3f& res) const
 {
   double znam(dir * normal);
   if(!same_val(znam, 0))
   {
     double tau(-((beg * normal) + D) / znam);
     if(tau >= 0)
     {
       res = beg + tau * dir;
       return tau;
     }
   }
   return -1.0;
 }



 int SPlate :: PlateSegmentIntersections(const SVec3f& beg, const SVec3f& end)
 {
  int s1(sign(beg)), s2(sign(end));
  SVec3f isk(end - beg);
  SVec3f e1(v_a - beg);
  SVec3f e2(v_b - beg);
  SVec3f e3(v_c - beg);
    
 /* double d1(sign(v_a));
  double d2(sign(v_b));
  double d3(sign(v_c));*/

  SMatrix3f m(e1, e2, e3);
  SVec3f pet;
  bool intersect(m.Solve(isk, pet) && pet.x > 0.0 && pet.y > 0.0 && pet.z > 0.0);

  if(intersect)
  {
    if(s1 * s2 > 0)
    {
      return 2;
    }
    else
    if(s1 == 0.0 || s2 == 0.0) // не очень аккуратно
    {
      return 3;//true;
    }
    else
    {
      return 1;
    }
  }
  return 0;
}

double SPlate :: DistanceToPointByVector(const SVec3f& point, 
                                         const SVec3f& vect)const 
{
  double zn(normal * vect);
  return fabs(zn) > EPS_PRECISION ? - vect.length() * (D + normal * point) / zn : 0;
}

 // Угол поворота против часовой стрелке. Если < 0, то ошибка
 // Код оттестирован
 double SVec2f :: operator%(const SVec2f& oth) const
 {
   double d(length());
   if(d > EPS_PRECISION)
   {
     double m(1.0 / d);
     double sin_fi((-oth.y * x + oth.x * y) * m);
     double cos_fi((y * oth.y + x * oth.x) * m);
     return sin_fi >= 0.0 ? acos(cos_fi) : 2 * M_PI - acos(cos_fi);
   }
   return -1.0;
 }

 bool Solve2x2(const SVec2f& e1, const SVec2f& e2, 
              const SVec2f& right, SVec2f& result)
 {
   double coef(e1.x * e2.y - e2.x * e1.y);
   if(coef != 0.0)
   {
     double inverse_coef(1.0 / coef);
     result.x = inverse_coef * (right.x * e2.y - right.y * e2.x);
     result.y = inverse_coef * (e1.x * right.y - e1.y * right.x);
     return true;
   }
   return false;
 }

 struct __edge
 {
   int beg;
   int end;
   __edge() { beg = end = -1; }
   __edge(int b, int e) :  beg(b), end(e) { if(beg > end) std :: swap(beg, end); }
   bool operator==(const __edge& s) const
   {
     return beg == s.beg && end == s.end;
   }
   bool operator<(const __edge& e) const
   {
     return (beg < e.beg) || (beg == e.beg && end < e.end);
   }
 };

 bool MakeTriangulation(const std :: vector<SVec3f>& points, 
                         std :: vector<IntTriade>& triangles)
 {
   static std :: vector<SVec2f> points2D;
   points2D.resize(points.size());
   if(points.size() >= 3) 
   {
     SVec3f norm((points[2] - points[0]) ^ (points[1] - points[0]));
     std :: vector<SVec3f> :: const_iterator p_it(points.begin());
     std :: vector<SVec2f> :: iterator p2d_it(points2D.begin());
     if(Is_Equal(norm.z, 0.0))
     {
       for(; p_it != points.end(); p_it++, p2d_it++)
       {
         p2d_it->x = p_it->x;
         p2d_it->y = p_it->z;
       }
     }
     else
     {
       for(; p_it != points.end(); p_it++, p2d_it++)
       {
         p2d_it->x = p_it->x;
         p2d_it->y = p_it->y;
       }
     }
     return MakeTriangulation(points2D, triangles);
   }
   return false;
 }

 bool MakeTriangulation(const std :: vector<SVec3f>& points, 
                             std :: vector<IntTriade>& triangles,
                             const SMatrix3f& project)
 {
   static std :: vector<SVec2f> points2D;
   points2D.resize(points.size());
   if(points.size() >= 3) 
   {
     std :: vector<SVec3f> :: const_iterator p_it(points.begin());
     std :: vector<SVec2f> :: iterator p2d_it(points2D.begin());
     for(; p_it != points.end(); p_it++, p2d_it++)
        *p2d_it = project * (*p_it);
     return MakeTriangulation(points2D, triangles);
   }
   return false;
 }

enum { inside, outside, boundary };         // положение точки
//     ВНУТРИ, ВНЕ,     НА ГРАНИЦЕ
enum { touching, crossing, inessential };   // положение ребра} ;

int edgeType (const SVec2f &a, const SVec2f &b, const SVec2f &e)
{
  switch (a.classify(b, e)) {
    case Left:
      return ((b.y<a.y)&&(a.y<=e.y)) ? crossing : inessential; 
    case Right:
      return ((e.y<a.y)&&(a.y<=b.y)) ? crossing : inessential; 
    case Between:
    case Begin:
    case End:
      return touching;
    default:
      return inessential;
  }
}

int IsInside(const SVec2f& probe_point, const std :: vector<SVec2f>& points)
{
  bool parity(false);
  size_t sz(points.size());
  for (size_t i = 0;  i < sz; i++) 
  { 
    switch (edgeType(probe_point, points[i], points[(i + 1) % sz])) 
    {
      case touching:
        return boundary;
      case crossing:
        parity = !parity;
    }
  }
  return (parity ? inside : outside);
}

struct __variant
{
  int var_num;
  double _cs;
  __variant() : var_num(-1), _cs(-2.0) { }
  __variant(int n, double c) : var_num(n), _cs(c) { }
  bool operator<(const __variant& v) const
  {
    return _cs < v._cs;
  }
};

bool cycle_check(int a, int b, int num)
{
  return (labs(a % num - b % num) == 1);    
}

bool IsEdgeInside(const std :: vector<SVec2f>& points, int b, int e)
{
  if(cycle_check(b, e, points.size()))
    return true;
  double d1, d2;
  for(size_t q = 0; q < points.size(); q++)
  {
    size_t qn((q + 1) % points.size());
    if(q != b && qn != b && e != q && qn != e)
    {
      if(AreSegmentsIntersecting(points[q], points[qn],
                                 points[b], points[e], d1, d2) != 0)
        return false;
    }
  }
  return IsInside((points[b] + points[e]) / 2.0, points) != outside; 
  
}


bool MakeTriangulation(const std :: vector<SVec2f>& points, 
                        std :: vector<IntTriade>& triangles)
{
  triangles.clear();
  std :: set<__edge >   edges;
  std :: set<int> vertecies;
  if(points.size() < 3) 
    return false;
  for(int q(0); q < int(points.size() - 1); q++)
  {
    edges.insert(__edge(q, q + 1));
    vertecies.insert(q); 
  }

  vertecies.erase(0);
  std :: stack<std :: pair<__edge, int> > now_edges;
  //now_edges.push(__edge(0, points.size() - 1));
  // Первичный выбор точки
  static std :: vector<__variant> vars;
  vars.resize(vertecies.size()); 
  std :: set<int> :: iterator vb_it(vertecies.begin());
  std :: vector<__variant> :: iterator var_it(vars.begin());
  for(; vb_it != vertecies.end(); vb_it++, var_it++)
  {
    int t_var(*vb_it);   
    var_it->var_num = t_var;
    SVec2f v1(points.front() - points[t_var]); 
    v1.Normalize();
    SVec2f v2(points.back() - points[t_var]);
    v2.Normalize();
    var_it->_cs = v1 ^ v2;
  }
  std :: sort(vars.begin(), vars.end());
  for(size_t q(0); q < vars.size(); q++)
  {
    int var_num(vars[q].var_num);
    if(/*IsInside((points.front() + points.back() + points[var_num]) / 3.0, points) == inside*/IsEdgeInside(points, 0, var_num) && IsEdgeInside(points, points.size() - 1, var_num))
    {
      vertecies.erase(var_num);
      __edge be(var_num, 0), ee(var_num, int(points.size() - 1));
      std :: set<__edge> :: iterator s_it(edges.find(be));
      if(s_it != edges.end())
        edges.erase(s_it);
      else
        now_edges.push(std :: pair<__edge, int>(be, int(points.size() - 1)));

      s_it = edges.find(ee);
      if(s_it != edges.end())
        edges.erase(s_it);
      else
        now_edges.push(std :: pair<__edge, int>(ee, 0));
      triangles.push_back(IntTriade(0, int(points.size() - 1), var_num));
      break;
    }
  }

  while(!now_edges.empty())
  {
    std :: set<int> :: iterator v_it(vertecies.begin());
    int best_var(-1);
    double best_cos(2.0);
    std :: pair<__edge, int> t_p(now_edges.top());
    const __edge& t_edge(t_p.first);
    S2DLine line(points[t_edge.beg], points[t_edge.end]);
    int l_sign(line.sign(points[t_p.second]));
    now_edges.pop();
    for(; v_it != vertecies.end(); v_it++)
    {
      int t_var(*v_it);
      const SVec2f& pnt(points[t_var]);
      if(l_sign * line.sign(pnt) < 0 && IsEdgeInside(points, t_var, t_edge.beg) && 
         IsEdgeInside(points, t_var, t_edge.end))
      {
        double _cs((points[t_edge.beg] - pnt) ^ (points[t_edge.end] - pnt));
        if(_cs < best_cos)
        {
          best_cos = _cs;
          best_var = t_var;
        }
      }
    }
    if(best_var > 0)
    {
      vertecies.erase(best_var);
      __edge be(best_var, t_edge.beg), ee(best_var, t_edge.end);
      std :: set<__edge> :: iterator s_it(edges.find(be));
      if(s_it != edges.end())
        edges.erase(s_it);
      else
        now_edges.push(std :: pair<__edge, int>(be, t_edge.end));

      s_it = edges.find(ee);
      if(s_it != edges.end())
        edges.erase(s_it);
      else
        now_edges.push(std :: pair<__edge, int>(ee, t_edge.beg));
      triangles.push_back(IntTriade(t_edge.beg, t_edge.end, best_var));
    }
    else
    {
#ifdef _MSC_VER
     FILE * pf;
     fopen_s(&pf, "triangle_error", "a");
#else
      FILE * pf = fopen("triangle_error", "a");
#endif
      fprintf(pf, "----------------------------------------------------------------\n");
      for(size_t q = 0; q < points.size(); q++)
        fprintf(pf, "%f %f \n", points[q].x, points[q].y);
      fclose(pf);
      return false; // ОШИБКА
    }
  }

  return true;
}

void SMatrix3f :: SetByColumns(const SVec3f& col1, const SVec3f& col2, 
                          const SVec3f& col3)
{
  m_data[0] = col1.x; m_data[3] = col1.y; m_data[6] = col1.z;
  m_data[1] = col2.x; m_data[4] = col2.y; m_data[7] = col2.z;
  m_data[2] = col3.x; m_data[5] = col3.y; m_data[8] = col3.z;
}

void SMatrix4f :: SetByColumns(const SVec4f& col1, const SVec4f& col2, 
                          const SVec4f& col3, const SVec4f& col4)
{
  m_data[0] = col1.x; m_data[4] = col1.y; m_data[8] = col1.z; m_data[12] = col1.t;
  m_data[1] = col2.x; m_data[5] = col2.y; m_data[9] = col2.z; m_data[13] = col2.t;
  m_data[2] = col3.x; m_data[6] = col3.y; m_data[10] = col3.z; m_data[14] = col3.t;
  m_data[3] = col4.x; m_data[7] = col4.y; m_data[11] = col4.z; m_data[15] = col4.t;
}

SMatrix3f getRotX(double f)
{
  SMatrix3f rs;
  rs.SetIdentity();
  rs(1,1) = cos(f);
  rs(2,2) = cos(f);
  rs(1,2) = -sin(f);
  rs(2,1) = sin(f);
  return rs;
}

SMatrix3f getRotY(double f)
{
  SMatrix3f rs;
  rs.SetIdentity();
  rs(0,0) = cos(f);
  rs(2,2) = cos(f);
  rs(0,2) = -sin(f);
  rs(2,0) = sin(f);
  return rs;
}

SMatrix3f getRotZ(double f)
{
  SMatrix3f rs;
  rs.SetIdentity();
  rs(0,0) = cos(f);
  rs(1,1) = cos(f);
  rs(0,1) = -sin(f);
  rs(1,0) = sin(f);
  return rs;
}

// поворот вокруг произвольной оси
SMatrix3f getArbRot(double f, SVec3f axis)
{
  SMatrix3f rs;
  axis.Normalize();
  SVec3f& a(axis);
  double c = cos(f);
  double s = sin(f);
  double c1 = 1 - c;
  rs(0, 0) = a.x * a.x * c1 + c;
  rs(0, 1) = a.x * a.y * c1 - a.z * s;
  rs(0, 2) = a.x * a.z * c1 + a.y * s;

  rs(1, 0) = a.x * a.y * c1 + a.z * s;
  rs(1, 1) = a.y * a.y * c1 + c;
  rs(1, 2) = a.y * a.z * c1 - a.x * s;

  rs(2, 0) = a.x * a.z * c1 - a.y * s;
  rs(2, 1) = a.x * a.z * c1 + a.x * s;
  rs(2, 2) = a.z * a.z * c1 + c;
  return rs;
}

SMatrix4f getMatrRotTransform(const SVec3f& point, SVec3f axis, double ang)
{
  SMatrix4f rs;
  axis.Normalize();
  SVec3f& a(axis);
  double c = cos(ang);
  double s = sin(ang);
  double c1 = 1 - c;
  rs(0, 0) = a.x * a.x * c1 + c;
  rs(0, 1) = a.x * a.y * c1 - a.z * s;
  rs(0, 2) = a.x * a.z * c1 + a.y * s;

  rs(1, 0) = a.x * a.y * c1 + a.z * s;
  rs(1, 1) = a.y * a.y * c1 + c;
  rs(1, 2) = a.y * a.z * c1 - a.x * s;

  rs(2, 0) = a.x * a.z * c1 - a.y * s;
  rs(2, 1) = a.x * a.z * c1 + a.x * s;
  rs(2, 2) = a.z * a.z * c1 + c;
  rs(3, 3) = 1;
  rs(3,0) = rs(3, 1) = rs(3, 2) = rs(0, 3) = rs(1, 3) = rs(2, 3) = 0; 
  
  rs = getMoveMatr(point) * rs * getMoveMatr(-point);
  
  return rs;        
}

SMatrix4f getMoveMatr(const SVec3f& mv)
{
  SMatrix4f rs;
  rs.SetIdentity();
  rs(0, 3) = mv.x;
  rs(1, 3) = mv.y;
  rs(2, 3) = mv.z;
  return rs;
}

SVec3f intersecton_point(const SVec3f& pA, const SVec3f& pB, const SVec3f& pC,
                           const SVec3f& p_beg, const SVec3f& p_dir)
{
  SVec3f norm((pB - pA) ^ (pC - pA));
  if(norm.length() < eps_precision_compare_double)
    throw SGeomNoPlate();
  return intersecton_point(pA, norm, p_beg, p_dir);
}

double intersecton_point_param(const SVec3f& pA, const SVec3f& pB, const SVec3f& pC,
                               const SVec3f& p_beg, const SVec3f& p_dir)
{
  SVec3f norm((pB - pA) ^ (pC - pA));
  if(norm.length() < eps_precision_compare_double)
    throw SGeomNoPlate();
  return intersecton_point_param(pA, norm, p_beg, p_dir);
}

double angle(const SVec3f& v1, const SVec3f& v2)
{
  double l1(v1.length()), l2(v2.length());
  if(same_val(l1, 0) || same_val(l2, 0 ))
    throw SGeomZeroLength();
  return acos(v1 * v2 / (l1 * l2));
}

bool IsPointInside(const SVec3f& pA, const SVec3f& pB, const SVec3f& pC, const SVec3f& intP)
{
  SVec3f r1(pA - intP);
  SVec3f r2(pB - intP);
  SVec3f r3(pC - intP);
  double val;
  try{
    val = angle(r1, r2) + angle(r2, r3) + angle(r1, r3);
  }
  catch(SGeomZeroLength ex)
  {
    return same_val(r1.length(), 0.0) || same_val(r2.length(), 0.0) || same_val(r3.length(), 0.0);
  }
  return same_val(val, 2 * M_PI);
}

bool does_segment_intersect_box2D(const SVec2f& beg, const SVec2f& end, const SVec2f& bl,  const SVec2f& tr)
{
  SVec2f AB(end - beg);
  if(same_val(AB.x,0))
  {
    double y_t, y_b;
    if(beg.y > end.y)
    {
      y_t = beg.y;
      y_b = end.y;
    }
    else
    {
      y_t = end.y;
      y_b = beg.y;
    }
    return (bl.x <= beg.x) && (tr.x >= beg.x) &&
            (y_b <= tr.y) && (y_t >= bl.y);
  }
  else if(same_val(AB.y, 0))
  {
    double x_t, x_b;
    if(beg.x > end.x)
    {
      x_t = beg.x;
      x_b = end.x;
    }
    else
    {
      x_t = end.x;
      x_b = beg.x;
    }
    return (bl.y <= beg.y) && (tr.y >= beg.y) &&
            (x_b <= tr.x) && (x_t >= bl.x);
  }
  else
  {
    double t, y_, x_;
    // left edge
    t = (bl.x - beg.x) / AB.x;
    if(t >= 0 && t <= 1.0)
    {
      y_ = beg.y + t * AB.y;
      return (y_ <= tr.y) && (y_ >= bl.y);
    }

    // right edge
    t = (tr.x - beg.x) / AB.x;
    if(t >= 0 && t <= 1.0)
    {
      y_ = beg.y + t * AB.y;
      return (y_ <= tr.y) && (y_ >= bl.y);
    }

    // bottom edge
    t = (bl.y - beg.y) / AB.y;
    if(t >= 0 && t <= 1.0)
    {
      x_ = beg.x + t * AB.x;
      return (x_ <= tr.x) && (x_ >= bl.x);
    }

    // top edge
    t = (tr.y - beg.y) / AB.y;
    if(t >= 0 && t <= 1.0)
    {
      x_ = beg.x + t * AB.x;
      return (x_ <= tr.x) && (x_ >= bl.x);
    }
  }

  return false;
}

