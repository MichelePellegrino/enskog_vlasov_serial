#include "base_interp.hpp"

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
