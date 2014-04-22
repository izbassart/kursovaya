#ifndef H_QUAT_H
#define H_QUAT_H

#include "math_lib.h"

struct SQuat
{
  union
  {
    struct
    {
      double t;
      double x;
      double y;
      double z;
    };
    double elem[4];
  };
  SQuat() { }
  SQuat(double it, double ix, double iy, double iz) : t(it), x(ix), y(iy), z(iz)  { }
  SQuat(double it) : t(0), x(0), y(0), z(0)  { }
  SQuat(double it, const SVec3f& vec) : t(it), x(vec.x), y(vec.y), z(vec.z) { }
  SQuat(const SVec3f& vec) : t(0), x(vec.x), y(vec.y), z(vec.z) { }

  SQuat operator+(const SQuat& iq) const
  {
    return SQuat(t + iq.t, x + iq.x, y + iq.y, z + iq.z);
  }

  SQuat operator-(const SQuat& iq) const
  {
    return SQuat(t - iq.t, x - iq.x, y - iq.y, z - iq.z);
  }

  SQuat operator/(double d) const
  {
    return SQuat(t / d, x / d, y / d, z / d);
  }

  double length2() const { return x * x + y * y + z * z + t * t; }
  
  double length() const { return sqrt(length2()); }

  SQuat operator*(const SQuat& iq) const 
  {
    /*return SQuat(x * iq.x - y * iq.y - z * iq.z - t * iq.t,
                 x * iq.y + iq.x * y + z * iq.t - iq.z * t,
                 x * iq.z + iq.x * z - y * iq.t + iq.y * t,
                 x * iq.t + iq.x * t + y * iq.z - iq.y * z );*/

    return SQuat(t * iq.t - x * iq.x - y * iq.y - z * iq.z,
                 t * iq.x + iq.t * x + y * iq.z - iq.y * z,
                 t * iq.y + iq.t * y - x * iq.z + iq.x * z,
                 t * iq.z + iq.t * z + x * iq.y - iq.x * y);
  }

  SQuat conjugate() const { return SQuat(t, - x, - y, - z); }

  SQuat inverse() const { return conjugate() / length2(); }

  SQuat left_quotient(const SQuat& iq) const 
  { 
    return iq.conjugate() * (*this) / iq.length2(); 
  }

  double Real() const { return t; }
  SVec3f vector() const { return SVec3f(x, y, z); }

  SQuat right_quotient(const SQuat& iq) const 
  { 
    return (*this) * iq.conjugate() / iq.length2(); 
  }

  SVec3f apply_to_vec(const SVec3f& vec) const
  {
    SQuat res = (*this) * vec * inverse();
    return res.vector();
  }

  SMatrix3f rot_matrix() const;

};

// кватернион поворота на угол alpha вокруг оси vec
SQuat gen_rot_quat(double alpha, SVec3f vec);

// генерация квартерниона по углам Крылова
SQuat gen_krylov_angles(double Kurs, double Kren, double Tangazh);
inline SQuat gen_krylov_angles_grad(double Kurs, double Kren, double Tangazh)
{
  return gen_krylov_angles(grad2rad(Kurs), grad2rad(Kren), grad2rad(Tangazh));
}

SQuat gen_from_matrix(const SMatrix3f& mat);

std :: ostream& operator<<(std :: ostream& ost, const SQuat& iq);

#endif