#include "integration.hpp"
#include "types.hpp"

#include <iostream>
#include <functional>
#include <cmath>

real_number f (real_number);

int main ()
{

using namespace ev_numeric;

std::function<real_number(real_number)> kernel_function = f;

real_number dist2 = 0;
std::cout << "Insert offset squared:" << std::endl;
std::cin >> dist2;

std::function<real_number(real_number)> psi = [&](const real_number& z)
  -> real_number
{
  return kernel_function( sqrt( dist2 + z*z ) );
};

// NumericalIntegrator<Finite> integrator( psi, 1.0, 2.0 );
NumericalIntegrator<Infinite> integrator( psi, 1.0, 2.0 );

real_number I = integrator.integrate( psi, 1.0, ev_const::pinfty );

std::cout << "I = " << I << std::endl;

return 0;

}

real_number f (real_number x)
{
  return 1.0/(x*x);
}
