#ifndef H_RKDP78_H
#define H_RKDP78_H

#include <vector>


namespace RKDP
{

class RightPart
{
public:
  virtual void Compute(size_t nvar, double T, double *Y, double *G) = 0;
};

class OutPut
{
public:
  virtual void DoOutput(double T, double step, size_t nvar, double* Y) = 0;
};

class Solver78
{
  //double[,] A;
  static double A[14][13]; 
  // 1..13
  static double B[];

     
  // 2..13
  static double C[];
      
    
  double  D[14]; 

  void InitD();
  void InitVar();

public:
  Solver78() 
  {
    InitVar();
  }

  std :: vector<double> check_constants();  
  bool integrate(size_t Nvar, 
                 double * Y, 
                 double t_beg, 
                 double t_end, 
                 double beg_step, 
                 double start_precision,
                 bool control,
                 RightPart * in_rp,
                 OutPut * in_op);
                         
};

};
#endif