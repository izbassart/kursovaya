#ifndef H_GRAVITY_FIELD_H
#define H_GRAVITY_FIELD_H

#include "math_lib.h"

class GraviationField
{
  int NHarm;                // number of used harmonics
  double GravityConstant;   // Gravity constant
  double Rad;               // Radius of Earth 

  static const double eps; 
  static double EG_coeff[];

public:

  GraviationField(int NH, double GC, double R)
  {
      NHarm = NH;
      GravityConstant = GC;
      Rad = R;

      // Initialization of EarthCoeff
  }

  SVec3f GraviationField :: Compute(double& U, // potential
                                  const SVec3f& ss_coord // position of spaceship
                                  ) const;   
};


#endif