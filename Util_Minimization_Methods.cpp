// -*- mode: C++; c-basic-offset: 2; -*-
#include <cmath> 
#include <limits>
#include <vector>
#include <iostream>
#include "Util_Misc.hpp"
#include "Util_Abstract_Group.hpp"
#include "Util_Minimization_Methods.hpp"

using namespace std;

namespace DFLib
{
  namespace Util
  {
    double Minimizer::simpleF(double &x,vector<double> &X0, 
                              vector<double> &dir)
    {
      int vecSize=X0.size();
      int i;

      vector<double> X(vecSize);
        
      for (i=0;i<vecSize;++i)
        X[i]=X0[i]+x*dir[i];
      theGroup->setEvaluationPoint(X);
      return(theGroup->getFunctionValue());
    }


    double Minimizer::simpleFandDeriv(double &x,vector<double> &X0, 
                                      vector<double> &dir, double &df)
    {
      int vecSize=X0.size();
      int i;
      double f;

      vector<double> X(vecSize);
      vector<double> grad(vecSize);

      for (i=0;i<vecSize;++i)
        X[i]=X0[i]+x*dir[i];
      theGroup->setEvaluationPoint(X);
      f=theGroup->getFunctionValueAndGradient(grad);
      df=0;
      for (i=0;i<vecSize;++i)
        df += grad[i]*dir[i];
      return f;
    }
    
    void Minimizer::bracketMinimum(double &a, double &b, double &c,
                                   vector<double> &X0, 
                                   vector<double> &direction)
    {

      double fa,fb,fc,fu;
      double temp;
      int i;
      const double GOLD=1.618034;
      const double TINY=1e-20;
      const double GLIMIT=100.0;
      bool done=false;


      fa=simpleF(a,X0,direction);
      fb=simpleF(b,X0,direction);

      // make sure we're going downhill from a to b
      if (fa < fb)        // switch role of a and b 
      {
        temp=a;   a=b;   b=temp;
        temp=fa; fa=fb; fb=temp;
      }
        
      c= b+GOLD*(b-a);
      fc=simpleF(c,X0,direction);

      while (fb > fc)
      {
        double r,q,u,denom,ulim;
        r=(b-a)*(fb-fc);
        q=(b-c)*(fb-fa);
        denom=q-r;
        if (fabs(denom)< TINY)
        {
          denom=TINY*((denom<0)?-1:1);
        }
        // U is abscissa of minimum of parabola that passes through
        // (a,fa),(b,fb),(c,fc)
        u=b-((b-c)*q-(b-a)*r)/(2*denom);
        ulim=b+GLIMIT*(c-b);
          
        if ((b-u)*(u-c) >0)  // u is between b and c
        {
          fu=simpleF(u,X0,direction);
          if (fu<fc)  // minimum is between b and c
          {
            a=b;
            b=u;
            fa=fb;
            fb=fu;
            break;
          }
          else if (fu>fb)  // minimum is between a and u
          {
            c=u;
            fc=fu;
            break;
          }
          // parabolic fit didn't get us anywhere
          u=c+GOLD*(c-b);
          fu=simpleF(u,X0,direction);
        }
        else if ((c-u)*(u-ulim) > 0.0)  // fit is between c and upper limit
        {
          fu=simpleF(u,X0,direction);
          if (fu<fc)
          {
            b=c;
            c=u;
            u=c+GOLD*(c-b);
            fb=fc;
            fc=fu;
            fu=simpleF(u,X0,direction);
          }
        } 
        else if ((u-ulim)*(ulim-c) >=0.0) // limit to maximum value
        {
          u=ulim;
          fu=simpleF(u,X0,direction);
        } 
        else                  // reject parabolic u, default magnification.
        {
          u=c+GOLD*(c-b);
          fu=simpleF(u,X0,direction);
        }

        a=b;
        b=c;
        c=u;
        fa=fb;
        fb=fc;
        fc=fu;
      }
    }

    double Minimizer::brentMinimize(double ax, double bx, double cx, 
                                    vector<double> &X0,
                                    vector<double> &dir,
                                    double &xmin
                                    )
    {

      double tol=sqrt(numeric_limits<double>::epsilon());
      bool ok1,ok2;
      int iter;
      double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
      double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
      const int ITMAX=100;
      const double ZEPS=1e-10;

      a=(ax<cx)?ax:cx;
      b=(ax>cx)?ax:cx;
      x=w=v=bx;
      fx=simpleFandDeriv(x,X0,dir,dx);
      fw=fv=fx;
      dw=dv=dx;
    
      for (iter=1;iter<=ITMAX;++iter)
      {
        xm=0.5*(a+b);
        tol1=tol*fabs(x)+ZEPS;
        tol2=2.0*tol1;
        if (fabs(x-xm) <= (tol2-0.5*(b-a)))
        {
          xmin=x;
          break;
        }
        if (fabs(e) > tol1)
        {
          // initialize to an out-of-bracket value

          d1=2*(b-a);
          d2=d1;
          // Secant method
          if (dw != dx) d1=(w-x)*dx/(dx-dw);
          if (dv != dx) d2=(v-x)*dx/(dx-dv);

          // use d within bracket, and on side pointed to by derivative
          u1=x+d1;
          u2=x+d2;
          ok1=((a-u1)*(u1-b)>0.0 && dx*d1<=0);
          ok2=((a-u2)*(u2-b)>0.0 && dx*d2<=0);
          olde=e;
          e=d;

          if (ok1||ok2) // take only an acceptable d, and smallest if both OK
          {
            if (ok1 && ok2)
              d=(fabs(d1)<fabs(d2)? d1:d2);
            else if (ok1)
              d=d1;
            else
              d=d2;
            if (fabs(d)<fabs(0.5*olde))
            {
              u=x+d;
              if (u-a < tol2 || b-u<tol2)
                d=(xm-x >0.0)?fabs(tol1):-fabs(tol1);
            }
            else
            {
              // do bisection instead, use derivative sign to decide which
              // side to look on.
              d=0.5*(e=((dx>0)?(a-x):(b-x)));
            }
          }
          else
          {
            d=0.5*(e=((dx>0)?(a-x):(b-x)));
          }
        }
        else 
        {
          d=0.5*(e=((dx>0)?(a-x):(b-x)));
        }
        if (fabs(d) >= tol1)
        {
          u=x+d;
          fu=simpleFandDeriv(u,X0,dir,du);
        }
        else
        {
          u=x+(d>0)?fabs(tol1):-fabs(tol1);
          fu=simpleFandDeriv(u,X0,dir,du);
          if (fu>fx) // if the minimum step in downhil direction takes us uphill,
            // we're done.
          {
            xmin=x;
            break;
          }
        }

        if (fu<=fx)
        {
          if (u>=x) 
            a=x;
          else
            b=x;
          v=w; fv=fw; dv=dw;
          w=x; fw=fx; dw=dx;
          x=u; fx=fu; dx=du;
        }
        else
        {
          if (u<x) 
            a=u;
          else
            b=u;
          if (fu<fw||w==x)
          {
            v=w; fv=fw; dv=dw;
            w=u; fw=fu; dw=du;
          }
          else if (fu<fv||v==x||v==w)
          {
            v=u; fv=fu; dv=du;
          }
        }
      }

      if (iter>ITMAX) // we took all our allotted iterations
        throw(Exception("Too many iterations in brentMinimize"));

      return fx;
    }

    double Minimizer::lineSearch(vector<double> &X0,vector<double> &dir)
    {

      double a=0.0,b=1.0,c=2.0;
      double fret;
      double xmin;
      int vecSize=X0.size();
      int j;

      bracketMinimum(a,b,c,X0,dir);
      fret=brentMinimize(a,b,c,X0,dir,xmin);

      for(j=0;j<vecSize;++j)
      {
        dir[j] *= xmin;
        X0[j] += dir[j];
      }
      return fret;
    }

    double Minimizer::conjugateGradientMinimize(vector<double>&X0,
                                                double ftol,
                                                int &iter)
    {
      int j,its;
      double gg,gam,fp,dgg;
      int vecSize=X0.size();
      vector<double> g(vecSize);
      vector<double> h(vecSize);
      vector<double> xi(vecSize);
      const int ITMAX=200;
      const double EPS=1e-10;
      double fret;

      theGroup->setEvaluationPoint(X0);
      fp=theGroup->getFunctionValueAndGradient(xi);
      for (j=0;j<vecSize;++j)
      {
        g[j]=-xi[j];
        xi[j]=h[j]=g[j];
      }
      for (iter=1;iter<=ITMAX;++iter)
      {
        fret=lineSearch(X0,xi);
        if (2*fabs(fret-fp) <= ftol*(fabs(fret)+fabs(fp)+EPS))
        {
          break;
        }
        theGroup->setEvaluationPoint(X0);
        fp=theGroup->getFunctionValueAndGradient(xi);
        dgg=gg=0.0;
        for (j=0;j<vecSize;++j)
        {
          gg += g[j]*g[j];
          dgg += (xi[j]+g[j])*xi[j];
        }

        if (gg == 0.0)   // Unlikely, but if gradient is 0, we're there.
        {
          break;
        }
        gam = dgg/gg;
        for (j=0;j<vecSize;++j)
        {
          g[j]=-xi[j];
          xi[j]=h[j]=g[j]+gam*h[j];
        }
      }
      if (iter>ITMAX)
        throw(Exception("Too many iterations in conjugateGradientMinimize"));
      return fret;
    }
  }    
}            
