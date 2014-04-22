#include "simpson_integral_comp.h"

#include <math.h>

double test_func(double t)
{
  return sin(t);       
}

double simpson_integral(double a, double b, int number, double_funct f)
{
  double res(0);
  if(a != b)
  {
    double step_1 = (b - a) / number;
    double x_odd = a + step_1 / 2;
    double x_even = a + step_1;
    double sum_odd(0); // нечетный
    double sum_even(0);// четный 
    for(int q(0); q < number - 1; q++)
    {
      sum_odd += f(x_odd); 
      sum_even += f(x_even);
      x_odd += step_1;
      x_even += step_1;       
    }
    sum_odd += f(x_odd);
    res = (f(a) + f(b) + 2 * sum_even + 4 * sum_odd) * step_1 / 6;    
  }
  return res;       
}

