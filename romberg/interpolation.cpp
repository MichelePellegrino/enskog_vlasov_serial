#include "interpolation.hpp"

Int Base_interp::locate(const Doub x)
/*
  Given a value x, return a value j such that x is (insofar as possible) centered
  in the subrange xx[j..j+mm-1], where xx is the stored pointer. The values in
  xx must be monotonic, either increasing or decreasing. The returned value is
  not less than 0, nor greater than n-1.
*/
{
  Int ju,jm,jl;
  if (n < 2 || mm < 2 || mm > n) throw("locate size error");
  Bool ascnd=(xx[n-1] >= xx[0]);
  jl=0;
  ju=n-1;
  while (ju-jl > 1)
  {
    jm = (ju+jl) >> 1;
    if (x >= xx[jm] == ascnd)
      jl=jm;
    else
      ju=jm;
  }
  cor = abs(jl-jsav) > dj ? 0 : 1;
  jsav = jl;
  return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}

Int Base_interp::hunt(const Doub x)
/*
  Given a value x, return a value j such that x is (insofar as possible) centered
  in the subrange xx[j..j+mm-1], where xx is the stored pointer. The values in
  xx must be monotonic, either increasing or decreasing. The returned value is
  not less than 0, nor greater than n-1.
*/
{
  Int jl=jsav, jm, ju, inc=1;
  if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
  /*
    bisection phase
  */
  Bool ascnd=(xx[n-1] >= xx[0]);
  if (jl < 0 || jl > n-1)
  {
    jl=0;
    ju=n-1;
  }
  else
  {
    if (x >= xx[jl] == ascnd)
    {
        for (;;)
        {
        /*
          True if ascending order of table, false otherwise. Input guess not useful.
          Go immediately to bisec-
        */
        ju = jl + inc;
        if (ju >= n-1) { ju = n-1; break;}
        else if (x < xx[ju] == ascnd) break;
        /*
          Off end of table.
        */
        else
        {
          jl = ju;
          inc += inc;
        }
      }
    }
    else
    {
      ju = jl;
      for (;;)
      {
        /*
          Found bracket. Not done, so double the increment and try again. tion. Hunt up:
          Hunt down:
        */
        jl = jl - inc;
        if (jl <= 0) { jl = 0; break;}
        else if (x >= xx[jl] == ascnd) break;
        else
        {
          /*
            Not done, so double the increment and try again.
          */
          ju = jl;
          inc += inc;
        }
      }
    }
  }
  while (ju-jl > 1)
  {
    jm = (ju+jl) >> 1;
    if (x >= xx[jm] == ascnd)
    jl=jm; else
    ju=jm;
    /*
      Hunt is done, so begin the final bisection phase:
    */
  }
  cor = abs(jl-jsav) > dj ? 0 : 1;
  jsav = jl;
  return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
}

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
