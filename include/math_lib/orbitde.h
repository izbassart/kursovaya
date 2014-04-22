#ifndef H_ORBITDE_H
#define H_ORBITDE_H

#include "math_lib.h"
#include "astronomy.h"
#include <ostream>
#include "quat.h"

struct OrbitPoint
{
  double T; // in thousands of seconds
  SVec3f position;
  SVec3f velocity;
  OrbitPoint() { }
  OrbitPoint(double iT, const SVec3f& pos, const SVec3f& vel) : T(iT), position(pos), velocity(vel) { }
};

inline std::ostream& operator<<(std :: ostream& ost, const OrbitPoint& kd)
{
  ost << "Time: " << kd.T << " pos: " << kd.position << " vel: " << kd.velocity << " height: " 
    << (kd.position.length() - 6.378136) * 1000 << "vel module: " << kd.velocity.length();
  return ost;
}

struct KnotData
{
  SVec3f position;
  SVec3f velocity;
  Astronomy::DateTime time;
  KnotData() { }
  KnotData(const SVec3f& pos, const SVec3f& vel, const Astronomy::DateTime dt) : position(pos), velocity(vel), time(dt) { }
};

inline std::ostream& operator<<(std :: ostream& ost, const KnotData& kd)
{
  ost << "Time: " << kd.time.get_string() << " pos: " << kd.position << " vel: " << kd.velocity ;
  return ost;
}

class OrbitDE
{
  KnotData m_knot;
  std :: vector<OrbitPoint> m_orbit_data;
  double JD_begin;

public:

  OrbitDE(const KnotData& kd);
  void Compute(double Dt, double step);
  const std :: vector<OrbitPoint>& get_data() const { return m_orbit_data; }
};



#endif
