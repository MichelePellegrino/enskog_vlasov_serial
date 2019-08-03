struct Poly_interp
{

  Int n, mm, jsav, cor, dj;
  const Doub *xx, *yy;
  Doub dy;

  Poly_interp(VecDoub_I &xv, VecDoub_I &yv, Int m):
    n(xv.size()), mm(m), jsav(0), cor(0), xx(&xv[0]), yy(yv), dy(0.)
    {
      dj = MIN(1,(int)pow((Doub)n,0.25));
    }

  Doub interp(Doub x)
  {
    Int jlo = cor ? hunt(x) : locate(x);
    return rawinterp(jlo,x);
  }

  Int locate(const Doub x)
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

  Int hunt(const Doub x)
  {
    Int jl=jsav, jm, ju, inc=1;
    if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
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
          ju = jl + inc;
          if (ju >= n-1) { ju = n-1; break;}
          else if (x < xx[ju] == ascnd) break;
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
          jl = jl - inc;
          if (jl <= 0) { jl = 0; break;}
          else if (x >= xx[jl] == ascnd) break;
          else
          {
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
    }
    cor = abs(jl-jsav) > dj ? 0 : 1;
    jsav = jl;
    return MAX(0,MIN(n-mm,jl-((mm-2)>>1)));
  }

  Doub rawinterp(Int jl, Doub x)
  {
    Int i,m,ns=0;
    Doub y,den,dif,dift,ho,hp,w;
    const Doub *xa = &xx[jl], *ya = &yy[jl];
    VecDoub c(mm),d(mm);
    dif=abs(x-xa[0]);
    for (i=0;i<mm;i++)
    {
      if ((dift=abs(x-xa[i])) < dif)
      {
        ns=i;
        dif=dift;
      }
      c[i]=ya[i];
      d[i]=ya[i];
    }
    y=ya[ns--];
    for (m=1;m<mm;m++)
    {
      for (i=0;i<mm-m;i++)
      {
        ho=xa[i]-x;
        hp=xa[i+m]-x;
        w=c[i+1]-d[i];
        if ((den=ho-hp) == 0.0) throw("Poly_interp error");
        den=w/den;
        d[i]=hp*den;
        c[i]=ho*den;
      }
      y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
    }
    return y;
  }

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
