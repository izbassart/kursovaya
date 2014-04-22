#include "quat.h"

std :: ostream& operator<<(std :: ostream& ost, const SQuat& iq)
{
  ost << "(" << iq.t << ", " << iq.x << ", " << iq.y << ", " << iq.z << ")";
  return ost;
}

// кватернион поворота на угол alpha вокруг оси vec
SQuat gen_rot_quat(double alpha, SVec3f vec)
{
  double alph2 = alpha / 2;
  vec.Normalize();
  vec *= sin(alph2);
  return SQuat(cos(alph2), vec.x, vec.y, vec.z);
}

// генерация квартерниона по углам Крылова
//                             aY          aX           aZ
SQuat gen_krylov_angles(double Kurs, double Kren, double Tangazh)
{
  double c_psi = cos(Kurs / 2);
  double s_psi = sin(Kurs / 2);
  
  double c_nu = cos(Kren / 2);
  double s_nu = sin(Kren / 2);

  double s_fi = sin(Tangazh / 2);
  double c_fi = cos(Tangazh / 2);

  return SQuat(c_psi * c_nu * c_fi + s_psi * s_nu * s_fi,
               c_psi * s_nu * c_fi + s_psi * c_nu * s_fi,
               s_psi * c_nu * c_fi - c_psi * s_nu * s_fi,
               c_psi * c_nu * s_fi - s_psi * s_nu * c_fi);
}

// возвращает номер максимального элемента
int get_max_num(double a, double b, double c)
{
 /* if(a > b)
  {
    return a > c ? 0 : 2;
  }
  else
  {
    return b > c ? 1 : 2;
  }*/

  return (a > b) ? ((a > c) ? 0 : 2) : ((b > c) ? 1 : 2);
}

// преобразование матрицы вращения в кватернион
// взято с http://www.rossprogrammproduct.com/translations/Matrix%20and%20Quaternion%20FAQ.htm#Q55
SQuat gen_from_matrix(const SMatrix3f& i_matr)
{
  SQuat res_q;
  const double * mat = i_matr.get_data();
  double T = mat[0] + mat[4] + mat[8] + 1;
  if(T > 0)
  {
      double S = 0.5 / sqrt(T);
      res_q.t = 0.25 / S;
      res_q.x = ( mat[7] - mat[5] ) * S;
      res_q.y = ( mat[2] - mat[6] ) * S;
      res_q.z = ( mat[3] - mat[1] ) * S;
  }
  else
  {
    int i = get_max_num(mat[0], mat[4], mat[8]);
    if(i == 0)
    {
      // Столбец 0:
      double S  = sqrt( 1.0 + mat[0] - mat[4] - mat[8] ) * 2;

      res_q.x = 0.5 / S;
      res_q.y = (mat[1] + mat[3] ) / S;
      res_q.z = (mat[2] + mat[6] ) / S;
      res_q.t = (mat[7] + mat[5] ) / S;
    }
    else if(i == 1)
    {
      // Столбец 1:
      double S  = sqrt( 1.0 + mat[4] - mat[0] - mat[8] ) * 2;

      res_q.x = (mat[1] + mat[3] ) / S;
      res_q.y = 0.5 / S;
      res_q.z = (mat[5] + mat[7] ) / S;
      res_q.t = (mat[2] + mat[6] ) / S;
    }
    else //if(i == 2)
    {
      // Столбец 2:
      double S  = sqrt( 1.0 + mat[8] - mat[0] - mat[4] ) * 2;

      res_q.x = (mat[2] + mat[6] ) / S;
      res_q.y = (mat[5] + mat[7] ) / S;
      res_q.z = 0.5 / S;
      res_q.t = (mat[1] + mat[3] ) / S;
    }
  }
  return res_q;
}

SMatrix3f SQuat :: rot_matrix() const
{
  SMatrix3f matr;
  double xx = x * x;
  double xy = x * y;
  double xz = x * z;
  double xt = x * t;

  double yy = y * y;
  double yz = y * z;
  double yt = y * t;

  double zz = z * z;
  double zt = z * t;
  
  matr[0]  = 1 - 2 * ( yy + zz );
  matr[1]  =     2 * ( xy - zt );
  matr[2]  =     2 * ( xz + yt );

  matr[3]  =     2 * ( xy + zt );
  matr[4]  = 1 - 2 * ( xx + zz );
  matr[5]  =     2 * ( yz - xt );

  matr[6]  =     2 * ( xz - yt );
  matr[7]  =     2 * ( yz + xt );
  matr[8] = 1 - 2 * ( xx + yy );
  return matr;
}
