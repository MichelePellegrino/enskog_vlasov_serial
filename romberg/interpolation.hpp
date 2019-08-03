#ifndef INTERP_HPP
#define INTERP_HPP

#include "recipe_types.hpp"

#include <cmath>
#include <algorithm>

struct Poly_interp : Base_interp
{
  Doub dy;
  Poly_interp(VecDoub_I &xv, VecDoub_I &yv, Int m):
    n(xv.size()), mm(m), jsav(0), cor(0), xx(&xv[0]), yy(yv), dy(0.) {}
  
  Doub interp(Doub x)
  {
    Int jlo = cor ? hunt(x) : locate(x);
    return rawinterp(jlo,x);
  }
  Doub rawinterp(Int jl, Doub x);
};

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
  Int locate(const Doub x);
  /*
    See definitions below.
  */
  Int hunt(const Doub x);
};

#endif /* INTERP_HPP */
