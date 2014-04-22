#ifndef H_SIMPSON_H
#define H_SIMPSON_H

typedef double(*double_funct)(double);

double test_func(double t);
double simpson_integral(double a, double b, int number, double_funct f);

#endif
