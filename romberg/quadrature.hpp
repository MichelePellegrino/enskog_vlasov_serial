#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include "recipe_types.hpp"

struct Quadrature
{
  /*
    Abstract base class for elementary quadrature algorithms.
  */
  Int n;                    //Current level of refinement.
  virtual Doub next() = 0;
  /*
    Returns the value of the integral at the nth stage of refinement. The function
    next() must be defined in the derived class.
  */
};

#endif /* QUADRATURE_HPP */
