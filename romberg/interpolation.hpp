#ifndef INTERP_HPP
#define INTERP_HPP

#include "recipe_types.hpp"

#include <cmath>
#include <algorithm>

struct Base_interp
/*
  Abstract base class used by all interpolation routines in this chapter.
  Only the routine interp is called directly by the user.
*/
{
  Int n, mm, jsav, cor, dj;
  const Doub *xx, *yy;
  Base_interp(VecDoub_I &x, const Doub *y, Int m)
  /*
    Constructor: Set up for interpolating on a table of x’s and y’s of length m. Normally called by a derived class, not by the user.
  */
  : n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y)
  {
    dj = MIN(1,(int)pow((Doub)n,0.25));
  }
  Doub interp(Doub x)
  {
    /*
      Given a value x, return an interpolated value, using data pointed to by xx and yy.
    */
    Int jlo = cor ? hunt(x) : locate(x);
    return rawinterp(jlo,x);
  }
  Int locate(const Doub x);
  /*
    See definitions below.
  */
  Int hunt(const Doub x);
  Doub virtual rawinterp(Int jlo, Doub x) = 0;
  /*
    Derived classes provide this as the actual interpolation method.
  */
};

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
  Doub rawinterp(Int jl, Doub x) override;
};

#endif /* INTERP_HPP */
