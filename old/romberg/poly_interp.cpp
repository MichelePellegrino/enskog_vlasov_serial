#include "poly_interp.hpp"

Doub Poly_interp::rawinterp(Int jl, Doub x)
/*
  Given a value x, and using pointers to data xx and yy, this routine returns an
  interpolated value y, and stores an error estimate dy. The returned value is
  obtained by mm-point polynomial interpolation on the subrange xx[jl..jl+mm-1].
*/
{
  Int i,m,ns=0;
  Doub y,den,dif,dift,ho,hp,w;
  const Doub *xa = &xx[jl], *ya = &yy[jl];
  VecDoub c(mm),d(mm);
  dif=abs(x-xa[0]);
  for (i=0;i<mm;i++)
  {
  /*
    Here we find the index ns of the closest table entry,
  */
  if ((dift=abs(x-xa[i])) < dif)
  {
    ns=i;
    dif=dift;
  }
  c[i]=ya[i];
  /*
    and initialize the tableau of c’s and d’s.
    This is the initial approximation to y. For each column of the tableau,
  */
  d[i]=ya[i];
  }
  y=ya[ns--];
  for (m=1;m<mm;m++)
  {
    for (i=0;i<mm-m;i++)
    {
      ho=xa[i]-x;
      /*
        we loop over the current c’s and d’s and update them.
      */
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ((den=ho-hp) == 0.0) throw("Poly_interp error");
      /*
        This error can occur only if two input xa’s are (to within roundoff) identical.
      */
      den=w/den;
      d[i]=hp*den;
      /*
        Here the c’s and d’s are updated.
      */
      c[i]=ho*den;
    }
    y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
    /*
      After each column in the tableau is completed, we decide which correction,
      c or d, we want to add to our accumulating value of y, i.e., which path to
      take through the tableau — forking up or down. We do this in such a way as
      to take the most “straight line” route through the tableau to its apex,
      updating ns accordingly to keep track of where we are. This route keeps the
      partial approximations centered (insofar as possible) on the target x. The
      last dy added is thus the error indication.
    */
  }
  return y;
}
