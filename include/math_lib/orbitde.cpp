#include "orbitde.h"
#include "gravity_field.h"
#include "rkdp78.h"

class OrbitDE_output : public RKDP::OutPut 
{
  std :: vector<OrbitPoint>& m_data;
public:
  OrbitDE_output(std :: vector<OrbitPoint>& md) : m_data(md) { }
  virtual void DoOutput(double T, double step, size_t nvar, double* Y)
  {
    // переменных должно быть 6
    if(m_data.empty() || fabs(T - m_data.back().T) >= step)
    {
      // запись
      m_data.push_back(OrbitPoint(T, SVec3f(Y[0], Y[1], Y[2]), SVec3f(Y[3], Y[4], Y[5])));
    }
  }
};

const double GravityConstant = 398.60044;
const double W = 0.072921151467;
const double EarthRadius = 6.378136;

const int number_Harm = 8;


class OrbitDE_rightpart : public RKDP::RightPart
{

  GraviationField m_field;

public:
  OrbitDE_rightpart() : m_field(number_Harm, GravityConstant, EarthRadius) { }
  virtual void Compute(size_t nvar, double T, double *Y, double *G)
  {
    G[0] = Y[3];
    G[1] = Y[4];
    G[2] = Y[5];

    double R2 = Y[0] * Y[0] + Y[1] * Y[1] + Y[2] * Y[2];
    double R = sqrt(R2);
    double H = - GravityConstant / (R * R2);
    
    double A = 0.0;
    SVec3f F =  m_field.Compute(A, SVec3f(Y));
    //F.Normalize();
    A = H + W * W;
    G[3] = A * Y[0] + 2 * W * Y[4] + F[0];
    G[4] = A * Y[1] - 2 * W * Y[3] + F[1];
    G[5] = H * Y[2] + F[2];
  }
};

OrbitDE::OrbitDE(const KnotData& kd) : m_knot(kd) 
{
  JD_begin = Astronomy::DateTimeToJs(kd.time);
}

const double orbit_eps = 1e-10;

void OrbitDE::Compute(double Dt, double step)
{
  m_orbit_data.clear();
  size_t num_steps = size_t(fabs(Dt / step));
  if(fabs(Dt) - num_steps * step > 0)
    num_steps ++;
  if(Dt < 0)
    step = -step;

  OrbitDE_output output(m_orbit_data);
  OrbitDE_rightpart rpart;
  double cur_Beg = 0;//Astronomy::DateTimeToJs(m_knot.time);
  double curY[6];

  curY[0] = m_knot.position.x;
  curY[1] = m_knot.position.y;
  curY[2] = m_knot.position.z;

  curY[3] = m_knot.velocity.x;
  curY[4] = m_knot.velocity.y;
  curY[5] = m_knot.velocity.z;

  for(size_t q(0); q < num_steps; q++)
  {
    double cur_End = cur_Beg + step;
    RKDP::Solver78().integrate(6, curY, cur_Beg, cur_End, step / 2, orbit_eps, false, 
                               &rpart,  &output);
    cur_Beg = cur_End;
  }
}


