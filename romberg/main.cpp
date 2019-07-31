#include <iostream>

#include "romberg.hpp"

double fun(const double&);

int main()
{
  double a = 0.0, b = 1.0;
  double I = qromb(fun, a, b);
  std::cout << " I = " << I << std::endl;
  return 0;
}

double fun(const double& x)
{

  return x*x;

}
