/*! \file integration.hpp
 *  \brief Header containing classes for numerical integration
 */

#ifndef EV_INTEGRATION_HPP
#define EV_INTEGRATION_HPP

#include "types.hpp"
#include "numerical_integration.hpp"
#include "numerical_rule.hpp"

#define DEFAULT_NINT 1000

#define QUADRATURE_TYPE Simpson

/*! \namespace ev_numeric
 *  \brief Namespace for numerical procedures
 *
 *  In this namespace are defined classes for numerical integration.
 */
namespace ev_numeric
{

using namespace NumericalIntegration;
using namespace Geometry;

enum IntegralType {Finite, Infinite};

/*! \class NumericalIntegrator
 *  \brief A class to perform numerical integration in finite and infinite 1D domains
 */
template<IntegralType dummy_integral_type>
class NumericalIntegrator
{ };

/* Template specialization for finite-domain integrals */
template<>
class NumericalIntegrator<Finite>
{
private:
  std::function<real_number(real_number)> func;
  int nint;
  Domain1D domain;
  Mesh1D mesh;
  Quadrature quad;
public:
  NumericalIntegrator() = default;
  NumericalIntegrator
  (std::function<real_number(real_number)> FUN, real_number A, real_number B, int N = DEFAULT_NINT):
    func(FUN), nint(N), domain(A,B), mesh(domain,N), quad(QUADRATURE_TYPE(),mesh) { }
  real_number integrate ( void )
  {
    return quad.apply(func);
  }
  real_number integrate
  (std::function<real_number(real_number)> FUN, real_number A, real_number B)
  {
    domain.left() = A;
    domain.right() = B;
    Mesh1D dummy_mesh(domain,nint);
    mesh = dummy_mesh;
    return quad.apply(func);
  }
  ~NumericalIntegrator() = default;
};

/* Template specialization for infinite-domain integrals */
template<>
class NumericalIntegrator<Infinite>
{
private:
  std::function<real_number(real_number)> func;                                 // integrand for finite integral
  std::function<real_number(real_number)> funk                                  // integrand for infinite integral
    = [this](real_number x) { return func(1.0/x)/(x*x); };
  int nint;
  Domain1D domain;
  Mesh1D mesh;
  Quadrature quad;
public:
  NumericalIntegrator() = default;
  NumericalIntegrator
  (std::function<real_number(real_number)> FUN, real_number A, real_number B, int N = DEFAULT_NINT):
    func(FUN), nint(N), domain(1.0/B, 1.0/A), mesh(domain,N), quad(QUADRATURE_TYPE(),mesh) { }
  real_number integrate ( void )
  {
    return quad.apply(funk);
  }
  real_number integrate
  (std::function<real_number(real_number)> FUN, real_number A, real_number B)
  {
    domain.left() = 1.0/B;
    domain.right() = 1.0/A;
    Mesh1D dummy_mesh(domain,nint);         mesh = dummy_mesh;
    Quadrature dummy_quad(QUADRATURE_TYPE(),mesh); quad = dummy_quad;
    return quad.apply(funk);
  }
  ~NumericalIntegrator() = default;
};

} /* namespace ev_numeric */

#endif /* EV_INTEGRATION_HPP */
