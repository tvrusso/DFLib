// -*- mode: C++; c-basic-offset: 2; -*-
//    DFLib: A library of Bearings Only Target Localization algorithms
//    Copyright (C) 2009-2011  Thomas V. Russo
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile$
//
// Purpose        : Provide a collection of minimization algorithms
//
// Special Notes  : To use these methods, one creates a subclass of the 
//                  DFLib::Abstract::Group inteface, providing a function
//                  evaluation method.  This Group object is passed to the
//                  constructor of the Minimization_Methods class, and then
//                  the various methods of the class can be used to find the
//                  minima.
//
// Creator        : 
//
// Creation Date  : 
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision$
//
// Revision Date  : $Date$
//
// Current Owner  : $Author$
//-------------------------------------------------------------------------
#include <cmath> 
#include <limits>
#include <vector>
#include <iostream>

#ifdef DFLIB_GRASS_OUTPUT
#include <fstream>
#endif

#include "Util_Misc.hpp"
#include "Util_Abstract_Group.hpp"
#include "Util_Minimization_Methods.hpp"

namespace DFLib
{
  namespace Util
  {
    double Minimizer::simpleF(double &x,std::vector<double> &X0, 
                              std::vector<double> &dir)
    {
      int vecSize=X0.size();
      int i;

      std::vector<double> X(vecSize);
        
      for (i=0;i<vecSize;++i)
        X[i]=X0[i]+x*dir[i];
      theGroup->setEvaluationPoint(X);
      return(theGroup->getFunctionValue());
    }


    double Minimizer::simpleFandDeriv(double &x,std::vector<double> &X0, 
                                      std::vector<double> &dir, double &df)
    {
      int vecSize=X0.size();
      int i;
      double f;

      std::vector<double> X(vecSize);
      std::vector<double> grad(vecSize);

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
                                   std::vector<double> &X0, 
                                   std::vector<double> &direction)
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

      int iter=0;
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
        iter++;
      }
    }

    double Minimizer::brentMinimize(double ax, double bx, double cx, 
                                    std::vector<double> &X0,
                                    std::vector<double> &dir,
                                    double &xmin
                                    )
    {

      double tol=sqrt(std::numeric_limits<double>::epsilon());
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

    double Minimizer::lineSearch(std::vector<double> &X0,std::vector<double> &dir)
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

    double Minimizer::conjugateGradientMinimize(std::vector<double>&X0,
                                                double ftol,
                                                int &iter)
    {
      int j,its;
      double gg,gam,fp,dgg;
      int vecSize=X0.size();
      std::vector<double> g(vecSize);
      std::vector<double> h(vecSize);
      std::vector<double> xi(vecSize);
      const int ITMAX=200;
      const double EPS=1e-10;
      double fret;

      iter=0;
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

    // returns index of simplex vertex with best function value
    int Minimizer::nelderMeadMinimize(std::vector<std::vector<double> >&Simplex)
    {
      int ndim=Simplex[0].size();
      int npts=Simplex.size();

      if (npts != ndim+1)
        throw(Exception("Number of points in simplex must be one more than number of dimensions of vectors"));

      std::vector<double> fVals(npts);
      // wasteful of space, but who cares.  
      std::vector<double> x0(ndim); // centroid of face through which we reflect
      std::vector<double> xr(ndim); // reflected vertex
      std::vector<double> xe(ndim); // reflected/expanded vertex
      std::vector<double> xc(ndim); // contracted vertex
      
      int indexOfBest, indexOfWorst, indexOfSecondWorst;
      int i;
      bool done=false;

      const double Alpha=1;
      const double Gamma=2;
      const double Rho=0.5;
      const double Sigma=0.5;
      const int maximumIterations=5000;

      double fTestR,fTestE,fTestC;
      int nFunctionEvals=0;
      int niters=0;
      double ftol=sqrt(std::numeric_limits<double>::epsilon());
      double rtol;

      // compute the function values for our simplex corners
      for (i=0; i<npts; i++)
      {
        theGroup->setEvaluationPoint(Simplex[i]);
        fVals[i]=theGroup->getFunctionValue();
        nFunctionEvals++;
      }

#ifdef DFLIB_GRASS_OUTPUT
      ofstream grassVector("testNelderMead.ascii");
      grassVector.precision(16); grassVector.width(20);
      grassVector << "ORGANIZATION: DFLIB Nelder-Mead code"<<std::endl;
      grassVector << "DIGIT DATE:   5/22/2009"<<std::endl;
      grassVector << "DIGIT NAME:   -"<<std::endl;
      grassVector << "MAP NAME:   nelder1"<<std::endl;
      grassVector << "MAP DATE:   2009"<<std::endl;
      grassVector << "MAP SCALE:  24000"<<std::endl;
      grassVector << "OTHER INFO:   Minimization Simplexes"<<std::endl;
      grassVector << "ZONE:   0"<<std::endl;
      grassVector << "MAP_THRESH:   0.500000"<<std::endl;
      grassVector << "VERTI:"<<std::endl;
#endif
      
      while (!done)
      {
#ifdef DFLIB_GRASS_OUTPUT
        // dump the current simplex to a GRASS vector
        grassVector << "B 4" << std::endl;
        for (i=0;i<npts;i++)
        {
          for (int j=0;j<ndim;j++)
            grassVector << " " << Simplex[i][j];
          grassVector<<std::endl;
        }
        for (int j=0;j<ndim;j++)
          grassVector << " " << Simplex[0][j];
        grassVector<<std::endl;
        grassVector<<"C  1 1" << std::endl;
        for (int j=0;j<ndim;j++)
        {
          double temp=0;
          for (int k=0;k<npts;k++)
            temp += Simplex[k][j];
          temp /= npts;
          grassVector << " " << temp;
        }
        grassVector<< std::endl << " 1 " << niters << std::endl;
#endif

        // locate best, worst, and second worst values.  An extraordinarily
        //inefficient way to do it.
        indexOfBest=0;
        if (fVals[0]>fVals[1])
        {
          indexOfWorst=0;
          indexOfSecondWorst=1;
        }
        else
        {
          indexOfWorst=1;
          indexOfSecondWorst=0;
        }
        for (i=0;i<npts;i++)
        {
          if (fVals[i]<fVals[indexOfBest]) 
            indexOfBest=i;
          if (fVals[i]>fVals[indexOfWorst])
          {
            indexOfSecondWorst=indexOfWorst;
            indexOfWorst=i;
          }
          else if (fVals[i]>fVals[indexOfSecondWorst] && i!=indexOfWorst)
          {
            indexOfSecondWorst=i;
          }
        }

#ifdef DFLIB_DEBUG
        std::cout << "nM " << niters++ << " nfuncs=" << nFunctionEvals << " "
             << "fVals["<<indexOfBest<<"]="<<fVals[indexOfBest] << " ";
        for (i=0;i<ndim;i++)
          std::cout << Simplex[indexOfBest][i] << " ";

        std::cout << "fVals["<<indexOfSecondWorst<<"]="<<fVals[indexOfSecondWorst] 
             << " " ;
        for (i=0;i<ndim;i++)
          std::cout << Simplex[indexOfSecondWorst][i] << " ";

        std::cout << "fVals["<<indexOfWorst<<"]="<<fVals[indexOfWorst] << " ";
        for (i=0;i<ndim;i++)
          std::cout << Simplex[indexOfWorst][i] << " ";
        std::cout << std::endl;
#endif
          
        rtol=2*abs(fVals[indexOfWorst]-fVals[indexOfBest])
          /(abs(fVals[indexOfWorst])+abs(fVals[indexOfBest]));
        if (rtol<ftol)
        {
          done=true;
        }
        else
        {
          if (nFunctionEvals>maximumIterations)
            throw(Exception("Maximum function evals exceeded in nelderMead"));
          
          // Now compute the center of mass of the side opposite the worst
          // point:
          x0.assign(ndim,0.0);
          for (i=0;i<npts;i++)
          {
            if (i!=indexOfWorst)
            {
              for(int j=0;j<ndim;j++)
                x0[j]+=Simplex[i][j];
            }
          }
          for(int j=0;j<ndim;j++)
            x0[j]/=ndim;

#ifdef DFLIB_DEBUG
          std::cout << "x0=";
          for (i=0;i<ndim;i++)
            std::cout << x0[i] << " ";
          std::cout << std::endl;
#endif
          
          // Compute the reflection point through the centroid
          for (int j=0;j<ndim;j++)
            xr[j]=x0[j]+Alpha*(x0[j]-Simplex[indexOfWorst][j]);
          
          // evaluate the function at xr
          theGroup->setEvaluationPoint(xr);
          fTestR=theGroup->getFunctionValue();
          nFunctionEvals++;

#ifdef DFLIB_DEBUG
          std::cout << " fTestR = " << fTestR << " Xr=";
          for (i=0;i<ndim;i++)
            std::cout << xr[i] << " ";
          std::cout << std::endl;
#endif

          // If this is the best of all....
          if (fTestR<fVals[indexOfBest])
          {
            // then try to expand it, too
            for (int j=0;j<ndim;j++)
              xe[j]=x0[j]+Gamma*(x0[j]-Simplex[indexOfWorst][j]);
            
            theGroup->setEvaluationPoint(xe);
            fTestE=theGroup->getFunctionValue();
            nFunctionEvals++;
            
            // if this is the best so far, replace the worst with it
            if (fTestE<fVals[indexOfBest])
            {
              fVals[indexOfWorst]=fTestE;
              Simplex[indexOfWorst]=xe;
#ifdef DFLIB_DEBUG
              std::cout << " nM expanded best, " << fTestE << " replacing " 
                   << indexOfWorst << std::endl;
#endif
            }
            else
            {
              // the reflected is the best so far, replace worst with it
              fVals[indexOfWorst]=fTestR;
              Simplex[indexOfWorst]=xr;
#ifdef DFLIB_DEBUG
              std::cout << " nM reflected best, " << fTestR << " replacing " 
                   << indexOfWorst << std::endl;
#endif
            }
          }
          else // reflected is not better than everything
          {
            // is reflected better than second worst?
            if (fTestR<fVals[indexOfSecondWorst])
            {
              // yes, toss the worst and use this one
              fVals[indexOfWorst]=fTestR;
              Simplex[indexOfWorst]=xr;
#ifdef DFLIB_DEBUG
              std::cout << " nM reflected better than second worst, " 
                   << fTestR << " replacing " 
                   << indexOfWorst << std::endl;
#endif
            }
            else
            {
              // no, worst than second worst
              // try contraction:
              // this is a point part way along the line connecting the
              // worst and the centroid.
              for (int j=0;j<ndim;j++)
                xc[j]=Simplex[indexOfWorst][j]
                  +Rho*(x0[j]-Simplex[indexOfWorst][j]);
              
              theGroup->setEvaluationPoint(xc);
              fTestC=theGroup->getFunctionValue();
              nFunctionEvals++;
              
#ifdef DFLIB_DEBUG
              std::cout << " fTestC = " << fTestC << " Xc=";
              for (i=0;i<ndim;i++)
                std::cout << xc[i] << " ";
              std::cout << std::endl;
#endif
              // is this better than the worst point?
              if (fTestC<=fVals[indexOfWorst])
              {
                // then toss the worst and replace with contracted
                fVals[indexOfWorst]=fTestC;
                Simplex[indexOfWorst]=xc;
#ifdef DFLIB_DEBUG
                std::cout << " nM contracted better than worst, " 
                   << fTestC << " replacing " 
                   << indexOfWorst << std::endl;
#endif
              }
              else
              {
                // we really can't win, can we?  Reduce the whole thing
                // toward the best point
#ifdef DFLIB_DEBUG
                std::cout << " nM reducing the whole deal " << std::endl;
                std::cout << "   best one is " << indexOfBest 
                     << " with function value " << fVals[indexOfBest] << std::endl;
#endif

                for (int vertex=0;vertex<npts;vertex++)
                  if (vertex != indexOfBest)
                  {
                    for (int component=0;component<ndim;component++)
                    {
                      Simplex[vertex][component]=
                        Simplex[indexOfBest][component]+
                        Sigma*(Simplex[vertex][component]-
                               Simplex[indexOfBest][component]);
                    }
                    theGroup->setEvaluationPoint(Simplex[vertex]);
                    fVals[vertex]=theGroup->getFunctionValue();
#ifdef DFLIB_DEBUG
                    std::cout << " Just changed vertex " << vertex << " value to " 
                         << fVals[vertex]<<std::endl;
                    std::cout << "  value of best is still fVals["
                         <<indexOfBest<<"]=" << fVals[indexOfBest] 
                         << std::endl;
#endif
                    nFunctionEvals++;
                  }
#ifdef DFLIB_DEBUG
                std::cout << " Finished reducing... " << std::endl;
                std::cout << "nM  nfuncs=" << nFunctionEvals << " "
                     << "fVals["<<indexOfBest<<"]="<<fVals[indexOfBest] << " " 
                     << "fVals["<<indexOfSecondWorst<<"]="
                     <<fVals[indexOfSecondWorst] << " " 
                     << "fVals["<<indexOfWorst<<"]="<<fVals[indexOfWorst] 
                     << std::endl;
#endif
              }
            }
          }
        }
      }
      return(indexOfBest);
    }

  }                                    
}
