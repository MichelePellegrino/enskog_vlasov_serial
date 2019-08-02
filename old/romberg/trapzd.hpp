#ifndef TRAPZD_HPP
#define TRAPZD_HPP

#include "recipe_types.hpp"
#include "quadrature.hpp"

template<class T>
struct Trapzd : Quadrature
{
  /*
    Routine implementing the extended trapezoidal rule.
  */
  /*
    Limits of integration and current value of integral.
  */
  T &func;
  Doub a,b,s;
  Trapzd() {};
  Trapzd(T &funcc, const Doub aa, const Doub bb) :
    func(funcc), a(aa), b(bb) {n=0;}
  /*
    The constructor takes as inputs func, the function or functor to be integrated
    between limits a and b, also input.
  */
  Doub next()
  {
    /*
      Returns the nth stage of refinement of the extended trapezoidal rule. On
      the first call (n=1), the routine returns the crudest estimate of Rab f .x/dx.
      Subsequent calls set n=2,3,... and improve the accuracy by adding 2n-2 additional
      interior points.
    */
    Doub x,tnm,sum,del;
    Int it,j;
    n++;
    if (n == 1)
    {
      return (s=0.5*(b-a)*(func(a)+func(b)));
    }
    else
    {
      for (it=1,j=1;j<n-1;j++) it <<= 1;
      tnm=it;
      del=(b-a)/tnm;
      /*
        This is the spacing of the points to be added.
      */
      x=a+0.5*del;
      for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(x);
      s=0.5*(s+(b-a)*sum/tnm);
      /*
        This replaces s by its refined value.
      */
      return s;
    }
  }
};

#endif /* TRAPZD_HPP */
