#ifndef POLY_INTERP_HPP
#define POLY_INTERP_HPP

#include "recipe_types.hpp"
#include "base_interp.hpp"

struct Poly_interp : Base_interp
/*
  Polynomial interpolation object. Construct with x and y vectors, and the number
  M of points to be used locally (polynomial order plus one), then call interp for
  interpolated values.
*/
{
  Doub dy;
  Poly_interp(VecDoub_I &xv, VecDoub_I &yv, Int m)
    : Base_interp(xv,&yv[0],m), dy(0.) {}
  Doub rawinterp(Int jl, Doub x);
};

#endif /* POLY_INTERP_HPP */
