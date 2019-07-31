#include "potential.hpp"

SutherlandMie::SutherlandMie
(real_number phi, real_number sigma, real_number gamma):
  pot_well(phi), mol_diam(sigma), exponent(gamma)
  {
    set_potential();
    set_dpotential_dr();
    set_pot_kernel();
  }

void
SutherlandMie::set_potential
(void)
{
  potential = [this](real_number r) -> real_number {
    return -pot_well*std::pow((mol_diam/r), exponent);
  };
}

void
SutherlandMie::set_dpotential_dr
(void)
{
  dpotential_dr = [this](real_number r) -> real_number {
    return -potential(r)*exponent/r;
  };
}

void
SutherlandMie::set_pot_kernel
(void)
{
  pot_kernel = [this](real_number r) -> real_number {
    return dpotential_dr(r)/r;
  };
}

SutherlandMorse::SutherlandMorse
(real_number phi, real_number sigma, real_number a):
  pot_well(phi), mol_diam(sigma), alpha(a)
  {
    set_potential();
    set_dpotential_dr();
    set_pot_kernel();
  }

void
SutherlandMorse::set_potential
(void)
{
  potential = [this](real_number r) -> real_number {
    return -pot_well*exp(-alpha*(r-mol_diam));
  };
}

void
SutherlandMorse::set_dpotential_dr
(void)
{
  dpotential_dr = [this](real_number r) -> real_number {
    return -potential(r)*alpha;
  };
}

void
SutherlandMorse::set_pot_kernel
(void)
{
  pot_kernel = [this](real_number r) -> real_number {
    return dpotential_dr(r)/r;
  };
}
