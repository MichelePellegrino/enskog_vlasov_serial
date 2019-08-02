#ifndef ROMBERG_HPP
#define ROMBERG_HPP

#include "recipe_types.hpp"
#include "interpolation.hpp"
#include "trapzd.hpp"

#define ITERMAX 25
#define MULT_STEP 0.25

template <class T>
Doub qromb(T &func, Doub a, Doub b, const Doub eps=1.0e-10)
{
  /*
    Returns the integral of the function or functor func from a to b. Integration
    is performed by Romberg’s method of order 2K, where, e.g., K=2 is Simpson’s rule.
  */
  const Int JMAX=ITERMAX, JMAXP=JMAX+1, K=5;
  /*
    Here EPS is the fractional accuracy desired, as determined by the
    extrapolation error estimate; JMAX limits the total number of steps; K is the
    number of points used in the extrapolation.
  */
  VecDoub s(JMAX), h(JMAXP);
  Poly_interp polint(h,s,K);
  h[0]=1.0;
  Trapzd<T> t(func,a,b);
  Doub ss = 0.0;
  for (Int j=1; j<=JMAX; j++)
  {
    s[j-1]=t.next();
    if (j >= K)
    {
      ss=polint.rawinterp(j-K,0.0);
      if (abs(polint.dy) <= eps*abs(ss)) return ss;
    }
    h[j]=MULT_STEP*h[j-1];
    /*
      This is a key step: The factor is 0.25 even though the stepsize is
      decreased by only 0.5. This makes the extrapolation a polynomial in h2 as
      allowed by equation (4.2.1), not just a polynomial in h.
    */
  }
  // std::cout << "WARNING: Too many steps in routine qromb" << std::endl;
  return ss;
}

#endif /* ROMBERG_HPP */
