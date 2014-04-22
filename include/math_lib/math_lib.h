/******************************************************************************************
Author: Vasiliy Sazonov

******************************************************************************************/
#ifndef H_MATH_LIB_H
#define H_MATH_LIB_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <memory.h>
#include <set>
#include <complex>

#include "exceptions.h"
#define EPS_PRECISION  0.000001
const double eps_precision_compare_double = 1e-09;
const double M_PI = 3.14159265358979323846;

inline double grad2rad(double g)
{
  return g *  M_PI / 180 ;
}

inline double rad2grad(double g)
{
  return 180 * g / M_PI ;
}


inline bool same_val(const double& x, const double& y)
{
  return (fabs(x-y) <= 0.0);// eps_precision_compare_double);
}

class SVec4f;

class SVec3f
{
 public:
   union
   {
	   struct
	   {
		   double x;
		   double y;
		   double z;
	   };
	   double elem[3];
     struct
     {
       double rho;
       double phi;
       double theta;
     };
   };
   SVec3f() : x(0.0), y(0.0), z(0.0) {}
   SVec3f(double ix, double iy, double iz) : x(ix), y(iy), z(iz) {}
   SVec3f(const SVec3f& f) : x(f.x), y(f.y), z(f.z) { }
   SVec3f(const SVec4f& f);
   SVec3f(double * Y) : x(Y[0]), y(Y[1]), z(Y[2]) { }

   const SVec3f& operator=(const SVec3f& f)  
   { 
	 x = f.x; y = f.y; z = f.z; 
	 return *this;
   }

   const SVec3f& operator+=(const SVec3f& f)  
   { 
	   x += f.x; y += f.y; z += f.z; 
	   return *this;
   }

   const SVec3f& operator-=(const SVec3f& f) 
   { 
	 x -= f.x; y -= f.y; z -= f.z; 
	 return *this;
   }

   const SVec3f& operator*=(const double& f) 
   { 
	 x *= f; y *= f; z *= f; 
	 return *this;
   }

   const SVec3f& operator/=(const double& f) 
   { 
	 x /= f; y /= f; z /= f; 
	 return *this;
   }
   
   double sq_length() const
   {
     return x * x + y * y + z * z;
   }

   double length() const 
   {
     return sqrt(sq_length());
   }

   void Normalize()
   {
	 double l(length());
	 if(l > EPS_PRECISION) 
	 {
     double m(1.0 / l);
	   x *= m;
	   y *= m;
	   z *= m;
	 }
   }


   SVec3f get_spherical() const // rho, phi, theta
   {
     double rho = length();
     if(rho == 0) 
       return SVec3f();
     double theta = asin(z / rho);
     double rho2 = sqrt (x * x + y * y);
     if(rho2 == 0)
       return SVec3f (rho, 0, theta);
     else
     {
       double phi = (y >= 0) ? acos(x / rho2) : M_PI - acos(x / rho2); 
       return SVec3f (rho, phi, theta);
     }
   }

   SVec3f operator+(const SVec3f& f) const { return SVec3f(x + f.x, y + f.y, z + f.z); }
   SVec3f operator-(const SVec3f& f) const { return SVec3f(x - f.x, y - f.y, z - f.z); }
   SVec3f operator*(const double& f) const { return SVec3f(f * x, f * y, f * z); }	  
   friend SVec3f operator*(const double& f, const SVec3f& r)  ;
   double& operator[](int i) { return elem[i]; }
   const double& operator[](int i) const { return elem[i]; }
   double operator*(const SVec3f& f) const { return x * f.x + y * f.y + z * f.z; }
   SVec3f operator/(const double& f) const { return SVec3f(x / f, y / f, z / f); }	  
   friend SVec3f operator/(double& f, const SVec3f& r) ;
   bool operator==(const SVec3f& f) const 
   { 
     return same_val(x, f.x) && same_val(y, f.y) && same_val(z, f.z) ; 
   }
   bool operator!=(const SVec3f& f) const 
   { 
     return !((*this) == f); 
   }
   SVec3f operator^(const SVec3f& f) const
   {
	   return SVec3f(y * f.z - f.y * z, f.x * z - x * f.z, x * f.y - f.x * y);	
   }

   SVec3f get_norm_vector() const
   {
     if(x != 0 || y != 0)
     {
       return SVec3f(-y, x, 0);
     }
     return SVec3f(-z, 0, x);
   }

   bool operator<(const SVec3f& r) const
   {
     return (x < r.x) || (x == r.x && y < r.y) ||  (x == r.x && y == r.y && z < r.z); 
   }

   // косинус угла между векторами
   double operator%(const SVec3f& f) const
   {
     return (*this * f) / (length() * f.length());
   }

   
   void WriteToFile(FILE* fp)
   {
     if(fp)
	 {
	   if(fwrite(elem, sizeof(double), 3, fp) != 3)	
		 throw SIOException(IOError);
	 }
	 else
       throw SIOException(WrongFile);
   }

   void ReadFromFile(FILE* fp)
   {
     if(fp)
	 {
	   if(fread(elem, sizeof(double), 3, fp) != 3)	
		 throw SIOException(IOError);
	 }
	 else
       throw SIOException(WrongFile);
   }

   std :: string get_string() const;
   void set_from_string(const std :: string& str);

   friend std :: ostream& operator<<(std :: ostream& st, const SVec3f& vec);
   friend std :: istream& operator>>(std :: istream& st, SVec3f& vec);
   friend double mixed_product(const SVec3f& a, const SVec3f& b, const SVec3f& c) ;

  };
   
  double angle(const SVec3f& v1, const SVec3f& v2);

  inline double intersecton_point_param(const SVec3f& pA, const SVec3f& norm,
                                        const SVec3f& p_beg, const SVec3f& p_dir)
  {
    double divider(norm * p_dir);
    if(fabs(divider) < eps_precision_compare_double)
      throw SGeomNoIntersection();
    return ((pA - p_beg) * norm / divider);
  }

  inline SVec3f reflect_ray(const SVec3f & r, const SVec3f & n) // r - ray, n - norm
  {
    return r - 2 * (r * n) * n ;
  }

  inline SVec3f intersecton_point(const SVec3f& pA, const SVec3f& norm,
                                  const SVec3f& p_beg, const SVec3f& p_dir)
  {
    return p_beg + intersecton_point_param(pA, norm, p_beg, p_dir) * p_dir ;
  }

  SVec3f intersecton_point(const SVec3f& pA, const SVec3f& pB, const SVec3f& pC,
                           const SVec3f& p_beg, const SVec3f& p_dir);
  
  double intersecton_point_param(const SVec3f& pA, const SVec3f& pB, const SVec3f& pC,
                                 const SVec3f& p_beg, const SVec3f& p_dir);
  
  bool IsPointInside(const SVec3f& pA, const SVec3f& pB, const SVec3f& pC, const SVec3f& intP);

  inline SVec3f operator-(const SVec3f& src)
  {
    return SVec3f(-src.x, -src.y, -src.z); 
  }
  
  std :: ostream& operator<<(std :: ostream& st, const SVec3f& vec);   

  double mixed_product(const SVec3f& a, const SVec3f& b, const SVec3f& c) ;

typedef std :: complex<double> complex_double;

class Complex_SVec3f
{
public: 
  SVec3f re;
  SVec3f im;

  const SVec3f& Re() const { return re; }
  const SVec3f& Im() const { return im; }

  Complex_SVec3f() { }
  Complex_SVec3f(const SVec3f& r, const SVec3f& i) : re(r), im(i) { }
  Complex_SVec3f(const SVec3f& r) : re(r), im(SVec3f()) { }
  Complex_SVec3f operator+(const Complex_SVec3f& f) const
  {
    return Complex_SVec3f(re + f.re, im + f.im);
  }

  Complex_SVec3f operator-(const Complex_SVec3f& f) const
  {
    return Complex_SVec3f(re - f.re, im - f.im);
  }

  Complex_SVec3f operator+=(const Complex_SVec3f& f) 
  {
    re += f.re;
    im += f.im;
    return *this;
  }

  Complex_SVec3f operator-=(const Complex_SVec3f& f) 
  {
    re -= f.re;
    im -= f.im;
    return *this;
  }

  complex_double operator*(const Complex_SVec3f& f) const
  {
    return complex_double(re * f.re - im * f.im, re * f.im + im * f.re);
  }
  
  Complex_SVec3f operator^(const Complex_SVec3f& f) const
  {
    return Complex_SVec3f((re ^ f.re) - (im ^ f.im), (re ^ f.im) + (im ^ f.re));
  }

  Complex_SVec3f operator^(const SVec3f& f) const
  {
    return Complex_SVec3f(re ^ f,  im ^ f);
  }

  Complex_SVec3f operator*(const double& f) const
  {
    return Complex_SVec3f(re * f,  im * f);
  }
  
  Complex_SVec3f operator*(const complex_double& f) const
  {
    return Complex_SVec3f(re * f.real() - im * f.imag(),  im * f.real() + re * f.imag());
  }
  
  Complex_SVec3f operator/(const double& f) const
  {
    return Complex_SVec3f(re / f,  im / f);
  }

  Complex_SVec3f conv() const
  {
    return Complex_SVec3f(re, -im);
  }

};

Complex_SVec3f operator*(const double& f, const Complex_SVec3f& vec);
Complex_SVec3f operator^(const SVec3f& f, const Complex_SVec3f& vec);
complex_double operator*(const SVec3f& f, const Complex_SVec3f& vec); 

inline Complex_SVec3f operator*(const complex_double& f, const Complex_SVec3f& vec) 
{
  return Complex_SVec3f(vec.re * f.real() - vec.im * f.imag(),  vec.im * f.real() + vec.re * f.imag());
}

class SVec4f
{
 public:
   union
   {
	   struct
	   {
		   double x;
		   double y;
		   double z;
       double t;
	   };
	   double elem[4];
   };
   SVec4f() : x(0.0), y(0.0), z(0.0), t(0.0) {}
   SVec4f(double ix, double iy, double iz, double it) : x(ix), y(iy), z(iz), t(it) {}
   SVec4f(const SVec4f& f) : x(f.x), y(f.y), z(f.z), t(f.t) { }
   SVec4f(const SVec3f& f) : x(f.x), y(f.y), z(f.z), t(1.0) { }

   const SVec4f& operator=(const SVec4f& f)  
   { 
	   x = f.x; y = f.y; z = f.z; t = f.t;
	   return *this;
   }

   const SVec4f& operator+=(const SVec4f& f)  
   { 
	   x += f.x; y += f.y; z += f.z; t += f.t;
	   return *this;
   }

   const SVec4f& operator-=(const SVec4f& f) 
   { 
	 x -= f.x; y -= f.y; z -= f.z; t -= f.t;
	 return *this;
   }

   const SVec4f& operator*=(const double& f) 
   { 
	 x *= f; y *= f; z *= f; t *= z;
	 return *this;
   }

   const SVec4f& operator/=(const double& f) 
   { 
	 x /= f; y /= f; z /= f; t /= f;
	 return *this;
   }

   SVec3f get_vector() const { return SVec3f(x / t, y / t, z / t); }
   
   double sq_length() const
   {
     return x * x + y * y + z * z + t * t;
   }

   double length() const 
   {
     return sqrt(sq_length());
   }

   void Normalize()
   {
	 double l(length());
	 if(l > EPS_PRECISION) 
	 {
     double m(1.0 / l);
	   x *= m;
	   y *= m;
	   z *= m;
     t *= m;
	 }
   }

   SVec4f operator+(const SVec4f& f) const 
   { 
     return SVec4f(x + f.x, y + f.y, z + f.z, t + f.t); 
   }

   SVec4f operator-(const SVec4f& f) const { return SVec4f(x - f.x, y - f.y, z - f.z, t - f.t); }
   SVec4f operator*(const double& f) const { return SVec4f(f * x, f * y, f * z, f * t); }	  
   friend SVec4f operator*(const double& f, const SVec4f& r)  ;
   double& operator[](int i) { return elem[i]; }
   const double& operator[](int i) const { return elem[i]; }
   double operator*(const SVec4f& f) const { return x * f.x + y * f.y + z * f.z + t * f.t; }
   SVec4f operator/(const double& f) const { return SVec4f(x / f, y / f, z / f, t / f); }	  
   friend SVec4f operator/(double& f, const SVec4f& r) ;
   bool operator==(const SVec4f& f) const 
   { 
     return same_val(x, f.x) && same_val(y, f.y) && same_val(z, f.z) && same_val(t, f.t); 
   }
   bool operator!=(const SVec4f& f) const 
   { 
     return !((*this) == f); 
   }
  
   bool operator<(const SVec4f& r) const
   {
     return (x < r.x) || (x == r.x && y < r.y) ||  (x == r.x && y == r.y && z < r.z) || (x == r.x && y == r.y && z == r.z && t < r.t); 
   }

   // косинус угла между векторами
   double operator%(const SVec4f& f) const
   {
     return (*this * f) / (length() * f.length());
   }

   
   void WriteToFile(FILE* fp)
   {
     if(fp)
	 {
	   if(fwrite(elem, sizeof(double), 4, fp) != 4)	
		 throw SIOException(IOError);
	 }
	 else
       throw SIOException(WrongFile);
   }

   void ReadFromFile(FILE* fp)
   {
     if(fp)
	 {
	   if(fread(elem, sizeof(double), 4, fp) != 4)	
		 throw SIOException(IOError);
	 }
	 else
       throw SIOException(WrongFile);
   }

   friend std :: ostream& operator<<(std :: ostream& st, const SVec4f& vec);

  };

class SMatrix3f 
   {
	   double m_data[9];
 public:
	SMatrix3f()  {  } 	
  
	SMatrix3f(const SMatrix3f& d)  
  { 
    memcpy(m_data, d.m_data, 9 * sizeof(double)); 
  }
  // инициализаци€ столбцами
  SMatrix3f(const SVec3f& col1, const SVec3f& col2, const SVec3f& col3)
  {
    SetByColumns(col1, col2, col3);
  }
  void SetByColumns(const SVec3f& col1, const SVec3f& col2, const SVec3f& col3);
 
	inline   void SetZero() { memset(m_data, 0,  9 * sizeof(double)); }
	inline   void SetIdentity() { SetZero(); m_data[0] = m_data[4] = m_data[8] = 1.0; }
	inline   SMatrix3f operator*(const SMatrix3f& matr) const 
	   {
		 SMatrix3f m;
		 for(size_t i(0); i < 3; i++)
		   for(size_t j(0); j < 3; j++)
		   {
			  double sum(0.0);
			  for(size_t q(0); q < 3; q++)
				  sum += m_data[3 * i + q] * matr.m_data[3 * q + j];	
			  m.m_data[3 * i + j] = sum;   
		   }
	     return m;
	   }

	   double& operator()(size_t a, size_t b)
	   {
		   return m_data[b + a * 3];
	   }

     double& operator[](size_t ind)
     {
       return m_data[ind];
     }

	   const double& operator()(size_t a, size_t b) const 
	   {
		   return m_data[b + a * 3];
	   }

	   SMatrix3f Transpose() const
	   {
		  SMatrix3f va(*this); 
		  double * v(va.m_data);
		  std :: swap(v[1], v[3]);
		  std :: swap(v[2], v[6]);
		  std :: swap(v[5], v[7]);
		  return va;
	   }

     const double * get_data() const { return m_data; }

     double Determinant() const 
     {
       return m_data[0] * m_data[4] * m_data[8] + m_data[3] * m_data[7] * m_data[2] + 
              m_data[6] * m_data[1] * m_data[5] - m_data[6] * m_data[4] * m_data[2] -
              m_data[3] * m_data[1] * m_data[8] - m_data[0] * m_data[7] * m_data[5];
     }

     SMatrix3f Revert() const
     {
       // вычисление обратной матрицы
       SMatrix3f ret;
       double d = Determinant();
       if(d != 0)
       {
          double d1 = 1.0 / d;
          ret.m_data[0] =   d1 * (m_data[4] * m_data[8] - m_data[5] * m_data[7]);
          ret.m_data[1] = - d1 * (m_data[1] * m_data[8] - m_data[2] * m_data[7]);
          ret.m_data[2] =   d1 * (m_data[1] * m_data[5] - m_data[2] * m_data[4]);
          ret.m_data[3] = - d1 * (m_data[3] * m_data[8] - m_data[5] * m_data[6]);
          ret.m_data[4] =   d1 * (m_data[0] * m_data[8] - m_data[2] * m_data[6]);
          ret.m_data[5] = - d1 * (m_data[0] * m_data[5] - m_data[2] * m_data[3]);
          ret.m_data[6] =   d1 * (m_data[3] * m_data[7] - m_data[4] * m_data[6]);
          ret.m_data[7] = - d1 * (m_data[0] * m_data[7] - m_data[1] * m_data[6]);
          ret.m_data[8] =   d1 * (m_data[0] * m_data[4] - m_data[1] * m_data[3]);
       }
       return ret; 
     }
           
     bool Solve(const SVec3f& in, SVec3f& out) const;

	   SVec3f operator*(const SVec3f& vec) const
	   {
		    SVec3f ret_vec;
	   	  for(size_t i = 0; i < 3; i++)
		    {
          double sum(0.0);
			    for(int q(0); q < 3; q++)
			      sum += vec.elem[q] * m_data[i * 3 + q];
			    ret_vec.elem[i] = sum;
		    }
		    return ret_vec;
	   }

 	   friend SVec3f operator*(const SVec3f& vec, const SMatrix3f& matr);
	   friend std :: ostream& operator<<(std :: ostream& st, const SMatrix3f& matr);
   };

   std :: ostream& operator<<(std :: ostream& st, const SMatrix3f& matr);



//=======================Projection operations========================================
   SMatrix3f OrthoProjection(const SVec3f& vec); 
   SMatrix3f OrthoProjection(const SVec3f& vec, const SVec3f& up);
   //=======================Geometric algorithms=========================================
   
///Plate description
class SMatrix4f 
{
  double m_data[16];
public:
	SMatrix4f()  {  } 	
  
	SMatrix4f(const SMatrix4f& d)  
  { 
    memcpy(m_data, d.m_data, 16 * sizeof(double)); 
  }
  
  void SetByColumns(const SVec4f& col1, const SVec4f& col2, const SVec4f& col3, const SVec4f& col4);
  // инициализаци€ столбцами
  SMatrix4f(const SVec4f& col1, const SVec4f& col2, const SVec4f& col3, const SVec4f& col4)
  {
    SetByColumns(col1, col2, col3, col4);
  }
  
  SMatrix4f(const SMatrix3f& d) 
  {
    m_data[0] = d(0,0);
    m_data[1] = d(0,1);
    m_data[2] = d(0,2);
    m_data[3] = 0;
    
    m_data[4] = d(1,0);
    m_data[5] = d(1,1);
    m_data[6] = d(1,2);
    m_data[7] = 0;

    m_data[8] = d(2,0);
    m_data[9] = d(2,1);
    m_data[10] = d(2,2);
    m_data[11] = 0;

    m_data[12] = m_data[13] = m_data[14] = m_data[15] = 0;
  }

	inline   void SetZero() { memset(m_data, 0,  16 * sizeof(double)); }
	inline   void SetIdentity() { SetZero(); m_data[0] = m_data[5] = m_data[10] = m_data[15] = 1.0;  }
	inline   SMatrix4f operator*(const SMatrix4f& matr) const 
	   {
		 SMatrix4f m;
		 for(size_t i(0); i < 4; i++)
		   for(size_t j(0); j < 4; j++)
		   {
			  double sum(0.0);
			  for(size_t q(0); q < 4; q++)
				  sum += m_data[4 * i + q] * matr.m_data[4 * q + j];	
			  m.m_data[4 * i + j] = sum;   
		   }
	     return m;
	   }

	   double& operator()(size_t a, size_t b)
	   {
		   return m_data[b + a * 4];
	   }

	   const double& operator()(size_t a, size_t b) const 
	   {
		   return m_data[b + a * 4];
	   }

	   SMatrix4f Transpose() const
	   {
		  SMatrix4f va(*this); 
		  double * v(va.m_data);
		  std :: swap(v[1], v[4]);
		  std :: swap(v[2], v[8]);
		  std :: swap(v[3], v[12]);
      std :: swap(v[9], v[6]);
      std :: swap(v[13], v[7]);
      std :: swap(v[11], v[14]);
		  return va;
	   }

             
	   SVec4f operator*(const SVec4f& vec) const
	   {
		    SVec4f ret_vec;
	   	  for(size_t i = 0; i < 4; i++)
		    {
          double sum(0.0);
			    for(int q(0); q < 4; q++)
			      sum += vec.elem[q] * m_data[i * 4 + q];
			    ret_vec.elem[i] = sum;
		    }
		    return ret_vec;
	   }

 	   friend SVec4f operator*(const SVec4f& vec, const SMatrix4f& matr);
	   friend std :: ostream& operator<<(std :: ostream& st, const SMatrix4f& matr);
   };

   std :: ostream& operator<<(std :: ostream& st, const SMatrix4f& matr);


class S3DPlate
{
  SVec3f normal;
  double D;
public:
  S3DPlate(const SVec3f& a, const SVec3f& b, const SVec3f& c)
  {
    normal.x = (a.y - b.y) * (b.z - c.z) - (a.z - b.z) * (b.y - c.y);
    normal.y = -(a.x - b.x) * (b.z - c.z) + (a.z - b.z) * (b.x - c.x);
    normal.z = (a.x - b.x) * (b.y - c.y) - (a.y - b.y) * (b.x - c.x);
    D = - normal * a;
  }

  S3DPlate(const SVec3f& norm, const SVec3f& point) : normal(norm)
  {
    D = - normal * point;
  }

  int sign(const SVec3f& probe_point) const
  {
    double d(probe_point * normal + D);
    if(d > 0) 
      return 1;
    else if(d < 0) 
      return -1;
    else  
      return 0;
  }

  double distance_to_plane(const SVec3f& point_from, const SVec3f& way) const
  {
    double z(normal * way);
    if(fabs(z) > 0)
    {
      return (D + point_from * normal) / z;
    }
    return -1.0;
  }
};


class SPlate
{
protected:
  SVec3f v_a, v_b, v_c;
  SVec3f normal;
  double  D;
public:
  SPlate(const SVec3f& a, const SVec3f& b, const SVec3f& c);
  SVec3f GetNormal() const { return normal; }
  double GetD() const { return D; }
  int sign(const SVec3f& v) const 
  { 
    double zn(normal * v + D);
    if(fabs(zn) < EPS_PRECISION)
    {
      return 0;
    }
    else
      return zn > 0.0 ? 1 : -1;
  } 
  double DistanceToPointByVector(const SVec3f& point, const SVec3f& vector)const ;

  bool InfinitePlateIntersection(const SVec3f& beg, const SVec3f& end);
  // 0 - no intersection, 1 - overlaped, 2 - visible projected
  int PlateSegmentIntersections(const SVec3f& beg, const SVec3f& end);
  bool GetPlateIntersection(const SVec3f& beg, const SVec3f& end, SVec3f& res) const;
  double GetRayPlateIntersection(const SVec3f&, const SVec3f&, SVec3f&) const;
  SVec3f GetProjection(const SVec3f& point, const SVec3f& dir) const;
};
  

typedef enum {Left, Right, Beyond, Behind, Between, Begin, End}  t2dclassify;

struct SVec2f
{
  static double org_x;
  static double org_y;
   union
   {
	   struct
	   {
		   double x;
		   double y;
	   };
	   double elem[2];
   };

  SVec2f(double ix = 0, double iy = 0) : x(ix), y(iy) { }
  SVec2f(const SVec2f& src) : x(src.x), y(src.y) { }
  SVec2f(const SVec3f& src) : x(src.x), y(src.y) { }
  static void set_org(double i_x, double i_y) { org_x = i_x; org_y = i_y; }
  static void set_org(const SVec2f& vec)  { org_x = vec.x; org_y = vec.y; }
  double arg() const
  {
    double len_ = sqrt((x - org_x) * (x - org_x) + (y - org_y) * (y - org_y));
    if(len_ != 0.0)
    {
      double cos_fi = (x - org_x) / len_; // * 1.0 + (y - org_y) * 0.0
      double sin_fi = /*- (x - org_x) * 0.0 + */ (y - org_y) ;
      return sin_fi >= 0 ? acos(cos_fi) : 2 * M_PI - acos(cos_fi);
    }
    else
      return 0.0;
  }

  double arg2() const
  {
    double len_ = sqrt(x * x + y * y);
    if(len_ != 0.0)
    {
      double cos_fi = x / len_; // * 1.0 + (y - org_y) * 0.0
      return y >= 0 ? acos(cos_fi) : 2 * M_PI - acos(cos_fi);
    }
    else
      return 0.0;
  }


  double dist_org() const
  {
    return sqrt((x - org_x) * (y - org_y));
  }

  bool operator<(const SVec2f& oth) const
  {
    return (arg() < oth.arg()) || (arg() == oth.arg() && dist_org() < oth.dist_org()) ;
  }

  bool right_turn(const SVec2f& oth) const
  {
    return x * oth.y - oth.x * y > 0.0;
  }
  
  SVec2f& operator=(const SVec3f& src) 
  {
    x = src.x;
    y = src.y;
    return * this;
  }

  SVec2f operator+(const SVec2f& src) const
  {
    return SVec2f(x+src.x, y + src.y);
  }

  SVec2f operator-(const SVec2f& src) const
  {
    return SVec2f(x - src.x, y - src.y);
  }

  double operator*(const SVec2f& src) const
  {
    return double(x * src.x + y * src.y);
  }

  SVec2f operator*(double src) const
  {
    return SVec2f(x * src, y * src);
  }

  SVec2f operator/(double src) const
  {
    double inv(1.0 / src);
    return SVec2f(x * inv, y * inv);
  }

  // Angle between vectors
  double operator^(const SVec2f& f) const
  {
    return acos((*this / length()) * (f / f.length()));
  }

  t2dclassify classify(const SVec2f &p0, const SVec2f &p1) const 
  {
    const SVec2f& p2(*this);
    SVec2f a(p1 - p0);
    SVec2f b(p2 - p0);
    double sa(a. x * b.y - b.x * a.y);
    if (sa > 0.0)
      return Left;
    if (sa < 0.0)
      return Right;
    if ((a.x * b.x < 0.0) || (a.y * b.y < 0.0))
      return Behind;
    if (a.length() < b.length())
      return Beyond;
    if (p0 == p2)
      return Begin;
    if (p1 == p2)
      return End;
    return Between;
  }

  void Normalize()
  {
    if(x == 0.0 && y == 0.0) 
      return;

    if(fabs(x) > fabs(y))
    {
      y = y / fabs(x);
      x = x > 0.0 ? 1.0 : -1.0;
    }
    else
    {
      x = x / fabs(y);
      y = y > 0.0 ? 1.0 : -1.0;
    }

       
    double m = 1.0 / length();
    x *= m;
    y *= m;
    
  }

  const SVec2f& operator*=(double src) 
  {
    x *= src; y *= src;
    return *this;
  }

  double operator%(const SVec2f& oth) const;
  double length() const { return sqrt(x * x + y * y); }
  SVec2f get_polar() const
  {
    return SVec2f(length(), arg2());
  }

  double distance(const SVec2f& oth) const 
  {
    return (oth - *this).length();
  }

  bool operator==(const SVec2f& src) const
  {
    return x == src.x && y == src.y;
  }
  void ReadFromFile(FILE* fp)
  {
    if(fp)
	  {
	    if(fread(elem, sizeof(double), 2, fp) != 2)	
	    throw SIOException(IOError);
	  }
    else
      throw SIOException(WrongFile);
  }

   void WriteToFile(FILE* fp)
   {
     if(fp)
	   {
	     if(fwrite(elem, sizeof(double), 2, fp) != 2)	
		   throw SIOException(IOError);
	   }
	   else
         throw SIOException(WrongFile);
   }
  friend std :: ostream& operator<<(std :: ostream& ost, const SVec2f& vec);
};
  inline std :: ostream& operator<<(std :: ostream& ost, 
                                    const SVec2f& vec)
   {
     ost << "(" << vec.x << "," << vec.y << ")";
     return ost;
   }

 inline SVec2f operator*(double alpha, const SVec2f& vec)
 {
   return SVec2f(alpha * vec.x, alpha * vec.y);
 }

 inline SVec2f operator-(const SVec2f& src)
  {
    return SVec2f(-src.x, -src.y);
  }
 
  bool does_segment_intersect_box2D(const SVec2f& beg, const SVec2f& end,
                                    const SVec2f& bl,  const SVec2f& tr);

bool Solve2x2(const SVec2f& e1, const SVec2f& e2, 
              const SVec2f& right, SVec2f& result);

SMatrix3f getRotX(double f);
SMatrix3f getRotY(double f);
SMatrix3f getRotZ(double f);
SMatrix3f getArbRot(double f, SVec3f axis);
SMatrix4f getMatrRotTransform(const SVec3f& point, SVec3f axis, double ang);

SMatrix4f getMoveMatr(const SVec3f& mv);

// выход дл€ триангул€ции
struct IntTriade
{
  int v1;
  int v2;
  int v3;
  IntTriade(int n1, int n2, int n3): v1(n1), v2(n2), v3(n3) { }
  IntTriade() { v1 = v2 = v3 = -1; }
  IntTriade(const IntTriade& src) : v1(src.v1), v2(src.v2), v3(src.v3) { }
};

inline bool Is_Equal(const double& d1, const double& d2)
{
  return fabs(d1 - d2) < EPS_PRECISION;
}

bool MakeTriangulation(const std :: vector<SVec3f>& points, 
                             std :: vector<IntTriade>& triangles);

bool MakeTriangulation(const std :: vector<SVec3f>& points, 
                             std :: vector<IntTriade>& triangles,
                             const SMatrix3f& project);

bool MakeTriangulation(const std :: vector<SVec2f>& points, 
                             std :: vector<IntTriade>& triangles);

struct S2Dsegment;

struct S2DLine
{
  double a_, b_, c_;
  S2DLine(const double& a, const double& b, const double& c) 
    :   a_(a), b_(b), c_(c) { }
  S2DLine()  :   a_(0.0), b_(0.0), c_(0.0) { }
  S2DLine(const S2DLine& src) 
    :   a_(src.a_), b_(src.b_), c_(src.c_) { }
  S2DLine(const S2Dsegment& seg);
  S2DLine(const SVec2f& beg, const SVec2f& end);
  S2DLine(const SVec3f& beg, const SVec3f& end);
  
  const S2DLine& operator=(const S2DLine& seg)
  {
    a_ = seg.a_; b_ = seg.b_; c_ = seg.c_;
    return *this;
  }

  int sign(const SVec2f& pnt) const 
  { 
    double val(a_ * pnt.x + b_ * pnt.y + c_);
    if(val < -EPS_PRECISION)
      return -1;
    else
    if(val > EPS_PRECISION)
      return 1;
    else
      return 0;
  }
};

class OrthoProjector
{
  SVec3f e1, e2;
public: 
  OrthoProjector(const SVec3f& ld, const SVec3f& up);
  SVec3f GetProjection(const SVec3f& vec) const;
};

template <class T> class STriangle
{
  T A_;
  T B_;
  T C_;

public:
  STriangle() { }
  STriangle(const T& iA, const T& iB, const T& iC) : A_(iA), B_(iB), C_(iC) { }
  T center() const { return 1.0 / 3 * (A_ + B_ + C_); }
  T AB() const { return B_ - A_; } 
  T BA() const { return A_ - B_; } 
  T AC() const { return C_ - A_; } 
  T CA() const { return A_ - C_; } 
  T CB() const { return B_ - C_; } 
  T BC() const { return C_ - B_; } 
  const T& A() const { return A_; }
  const T& B() const { return B_; }
  const T& C() const { return C_; }
  void ChangeOri()
  {
    std :: swap(B_, C_);
  }

};

class  STriangle3d : public STriangle<SVec3f> 
{
public:
  SVec3f Normal() const 
  { 
    SVec3f vec(AB() ^ AC());
    vec.Normalize();
    return vec;
  }
  
  double Area() const
  {
    return 0.5 * (AB() ^ AC()).length(); 
  } 
};
typedef STriangle<SVec2f> STriangle2d;


#endif

