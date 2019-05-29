
/*
// ROUTINES FROM NUMERICAL RECIPES IN C
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "numrec.h"

#define TOL 2.0e-4

#define ITMAX 200
#define CGOLD 0.3819660
#define ZEPS 1.0e-10

/* Here ITMAX is the maximum allowed number of iterations; CGOLD is 
   the golden ratio; ZEPS is a small number that protects against 
   trying to achieve fractional accuracy for a minimum that happens 
   to be exactly zero.*/

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-25

#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

/* Here GOLD is the default ratio by which successive intervals are 
   magnified; GLIMIT is the maximum magnification allowed for a 
   parabolic-fit step. */

int ncom;

/* global variables to communicate with n1dim */
double *pcom,*xicom,(*nrfunc)(double []);

/* defines for lnsrch */
#define ALF 1.0e-4 /* Ensures sucient decrease in function value. */
#define TOLX 1.0e-7 /* Convergence criterion on x.*/

/* defines and vars for newt */
#define MAXITS 200
#define TOLF 1.0e-4
#define TOLMIN 1.0e-6
#define STPMX 100.0
/*Here MAXITS is the maximum number of iterations; TOLF sets the convergence criterion on
function values; TOLMIN sets the criterion for deciding whether spurious convergence to a
minimum of fmin has occurred; TOLX is the convergence criterion on x; STPMX is the scaled
maximum step length allowed in line searches. */
int nn; /* Global variables to communicate with fmin.*/
double *fvec;
void (*nrfuncv)(int n, double v[], double f[]);
#define FREERETURN {free_dvector(fvec,1,n);free_dvector(xold,1,n);\
free_dvector(p,1,n);free_dvector(g,1,n);free_dmatrix(fjac,1,n,1,n);\
free_ivector(indx,1,n);return;}


#define EPS 1.0e-8 /*Approximate square root of the machine double precision.*/

void spline(double x[], double y[], int n, 
             double yp1, double ypn, double y2[])
/*
  Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., 
  y i = f(xi), with x1 < x2 < : : : < xN , and given values yp1 and ypn 
  for the first derivative of the interpolating function at points 1 and 
  n, respectively, this routine returns an array y2[1..n] that contains
  the second derivatives of the interpolating function at the tabulated 
  points xi. If yp1 and/or ypn are equal to 1 x 10^30 or larger, the 
  routine is signaled to set the corresponding boundary condition for a 
  natural spline, with zero second derivative on that boundary.
*/
{
  int i,k;
  double p,qn,sig,un,*u;

  u=dvector(1,n-1);
  if (yp1 > 0.99e30) 
  /* The lower boundary condition is set either to be natural
     or else to have a specified first derivative. */ 

    y2[1]=u[1]=0.0; 
  else { 
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }

  for (i=2;i<=n-1;i++) { 
  /* This is the decomposition loop of the tridiagonal algorithm. 
     y2 and u are used for temporary storage of the decomposed
     factors. */
    
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30) 
  /* The upper boundary condition is set either to be natural
     or else to have a specified first derivative. */

    qn=un=0.0;
  else { 
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--) 
  /* This is the backsubstitution loop of the tridiagonal
     algorithm. */
  y2[k]=y2[k]*y2[k+1]+u[k];
  free_dvector(u,1,n-1);
}


void splint(double xa[], double ya[], double y2a[], 
             int n, double x, double *y)
/* Given the arrays xa[1..n] and ya[1..n], which tabulate a function 
   (with the xai's in order), and given the array y2a[1..n], which is 
   the output from spline above, and given a value of x, this routine 
   returns a cubic-spline interpolated value y. */
{
  void nrerror(char error_text[]);
  int klo,khi,k;
  double h,b,a;
  klo=1; 
  /* We will find the right place in the table by means of
     bisection. This is optimal if sequential calls to this
     routine are at random values of x. If sequential calls
     are in order, and closely spaced, one would do better
     to store previous values of klo and khi and test if
     they remain appropriate on the next call.
  */

  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  } 
  /* klo and khi now bracket the input value of x. */

  h=xa[khi]-xa[klo];
  if (h == 0.0) nrerror("Bad xa input to routine splint"); 
  /* The xa's must be distinct. */

  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h; 
  /* Cubic spline polynomial is now evaluated. */

  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void splie2(double x1a[], double x2a[], double **ya, 
             int m, int n, double **y2a)
/*
  Given an m by n tabulated function ya[1..m][1..n], and tabulated 
  independent variables x2a[1..n], this routine constructs 
  one-dimensional natural cubic splines of the rows of ya and returns 
  the second-derivatives in the array y2a[1..m][1..n]. (The array 
  x1a[1..m] is included in the argument list merely for consistency 
  with routine splin2.)  Typically execute this once.
*/
{
  int j;
  for (j=1;j<=m;j++)
    spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]); 

  /* Values 1 x 10 ^ 30 signal a natural spline. */

}



void splin2(double x1a[], double x2a[], 
            double **ya, double **y2a, 
            int m, int n,
            double x1, double x2, 
            double *y)
/*
  Given x1a, x2a, ya, m, n as described in splie2 and y2a as produced 
  by that routine; and given a desired interpolating point x1,x2; this 
  routine returns an interpolated function value y by bicubic spline 
  interpolation.
*/
{
  int j;
  double *ytmp,*yytmp;

  ytmp=dvector(1,m);
  yytmp=dvector(1,m); 
  /* Perform m evaluations of the row splines constructed by splie2, 
     using the one-dimensional spline evaluator splint.*/

  for (j=1;j<=m;j++)
    splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);

  spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp); 
  /* Construct the one-dimensional column spline and evaluate it. */
  splint(x1a,yytmp,ytmp,m,x1,y);

  free_dvector(yytmp,1,m);
  free_dvector(ytmp,1,m);
}


void powell(double p[], double **xi, int n, double ftol, int *iter, 
            double *fret,
            double (*func)(double []))

/* Minimization of a function func of n variables. Input consists of 
   an initial starting point p[1..n]; an initial matrix xi[1..n][1..n],
   whose columns contain the initial set of directions (usually the n 
   unit vectors); and ftol, the fractional tolerance in the function 
   value such that failure to decrease by more than this amount on one
   iteration signals doneness. On output, p is set to the best point 
   found, xi is the then-current direction set, fret is the returned
   function value at p, and iter is the number of iterations taken. 
   The routine linmin is used. */
{
  int i,ibig,j;
  double del,fp,fptt,t,*pt,*ptt,*xit;
  pt=dvector(1,n);
  ptt=dvector(1,n);
  xit=dvector(1,n);
  *fret=(*func)(p);
  for (j=1;j<=n;j++) pt[j]=p[j]; /* Save the initial point. */
  for (*iter=1;;++(*iter)) {
    fp=(*fret);
    ibig=0;
    del=0.0; /* Will be the biggest function decrease. */
    for (i=1;i<=n;i++) { 
      /*In each iteration, loop over all directions in the set.*/
      for (j=1;j<=n;j++) xit[j]=xi[j][i];
      fptt=(*fret);
      linmin(p,xit,n,fret,func);  
      /*Copy the direction, minimize along it, and record it if it 
        is the largest decrease so far. */
      if (fptt-(*fret) > del) { 
        del=fptt-(*fret);
        ibig=i;
      }
    }
    if (2.0*(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))+TINY) {
      free_dvector(xit,1,n); /* Termination criterion.*/
      free_dvector(ptt,1,n);
      free_dvector(pt,1,n);
      return;
    }
    if (*iter == ITMAX) nrerror("powell exceeding maximum iterations.");
    for (j=1;j<=n;j++) {
    /* Construct the extrapolated point and the average direction 
       moved. Save the old starting point.*/
      ptt[j]=2.0*p[j]-pt[j];
      xit[j]=p[j]-pt[j];
      pt[j]=p[j];
    }
    fptt=(*func)(ptt); /* Function value at extrapolated point. */
    if (fptt < fp) {
      t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
      if (t < 0.0) {
        linmin(p,xit,n,fret,func);
       /* Move to the minimum of the new direction, and save the new 
          direction.*/
        for (j=1;j<=n;j++) {
          xi[j][ibig]=xi[j][n];
          xi[j][n]=xit[j];
        }
      }
    }
  }
}



void linmin(double p[], double xi[], int n, double *fret, 
            double (*func)(double []))
/* Given an n-dimensional point p[1..n] and an n-dimensional 
   direction xi[1..n], moves and resets p to where the function 
   func(p) takes on a minimum along the direction xi from p,
   and replaces xi by the actual vector displacement that p was 
   moved. Also returns as fret the value of func at the returned 
   location p. This is actually all accomplished by calling the
   routines mnbrak and brent. */
{
  int j;
  double xx,xmin,fx,fb,fa,bx,ax;
  /* define the global variables */
  ncom=n; 
  pcom=dvector(1,n);
  xicom=dvector(1,n);
  nrfunc=func;
  for (j=1;j<=n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0; /* Initial guess for brackets. */
  xx=1.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
  for (j=1;j<=n;j++) { 
    /* Construct the vector results to return. */
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  free_dvector(xicom,1,n);
  free_dvector(pcom,1,n);
}

double f1dim(double x)
/* Must accompany linmin. */
{
  int j;
  double f,*xt;
  xt=dvector(1,ncom);
  for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt);
  free_dvector(xt,1,ncom);
  return f;
}


double brent(double ax, double bx, double cx, 
            double (*f)(double), double tol, double *xmin)

/* Given a function f, and given a bracketing triplet of abscissas ax, 
   bx, cx (such that bx is between ax and cx, andf(bx) is less than 
   both f(ax) and f(cx)), this routine isolates the minimum to a 
   fractional precision of about tol using Brent's method. The 
   abscissa of the minimum is returned as xmin, and the minimum 
   function value is returned as brent, the returned function value. */
{
  int iter;
  double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0; 
  /* This will be the distance moved on the step before last. */
  a=(ax < cx ? ax : cx); 

  /* a and b must be in ascending order, but input abscissas need 
     not be. */
  b=(ax > cx ? ax : cx);

  x=w=v=bx; /* Initializations... */
  fw=fv=fx=(*f)(x);
  for (iter=1;iter<=ITMAX;iter++) {  /* Main program loop. */
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) { /* Test for done here. */
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) { /*Construct a trial parabolic fit.*/
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
        d=CGOLD*(e=(x >= xm ? a-x : b-x));
      /* The above conditions determine the acceptability of the 
         parabolic fit. Here we take the golden section step into 
         the larger of the two segments. */
      else {
        d=p/q; /* Take the parabolic step. */
        u=x+d;
        if (u-a < tol2 || b-u < tol2)
          d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u);
    /* This is the one function evaluation per iteration. */
    if (fu <= fx) { 
      /* Now decide what to do with our function evaluation.*/
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u) /* Housekeeping follows: */
      SHFT(fv,fw,fx,fu)
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      } else if (fu <= fv || v == x || v == w) {
        v=u;
        fv=fu;
      }
    } /* Done with housekeeping. Back for
         another iteration. */
  }
  nrerror("Too many iterations in brent");
  *xmin=x; /* Never get here. */
  return fx;
}


void mnbrak(double *ax, double *bx, double *cx, 
            double *fa, double *fb, double *fc,
            double (*func)(double))
/* Given a function func, and given distinct initial points ax and 
   bx, this routine searches in the downhill direction (defined by 
   the function as evaluated at the initial points) and returns
   new points ax, bx, cx that bracket a minimum of the function. 
   Also returned are the function values at the three points, fa, 
   fb, and fc. */
{
  double ulim,u,r,q,fu,dum;
  *fa=(*func)(*ax);
  *fb=(*func)(*bx);
  if (*fb > *fa) { 
    /* Switch roles of a and b so that we can go downhill in the 
       direction from a to b. */
    SHFT(dum,*ax,*bx,dum)
    SHFT(dum,*fb,*fa,dum)
  }
  *cx=(*bx)+GOLD*(*bx-*ax); /* First guess for c. */
  *fc=(*func)(*cx);
  while (*fb > *fc) { /* Keep returning here until we bracket. */
    r=(*bx-*ax)*(*fb-*fc); 
    /* Compute u by parabolic extrapolation from a; b; c. 
       TINY is used to prevent any possible division by zero. */
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(MAXD(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    /* We won't go farther than this. Test various possibilities: */
    if ((*bx-u)*(u-*cx) > 0.0) { 
      /* Parabolic u is between b and c: try it. */
      fu=(*func)(u);
      if (fu < *fc) { 
        /* Got a minimum between b and c. */
        *ax=(*bx);
        *bx=u;
        *fa=(*fb);
        *fb=fu;
        return;
      } else if (fu > *fb) { 
        /* Got a minimum between between a and u. */
        *cx=u;
        *fc=fu;
        return;
      }
      u=(*cx)+GOLD*(*cx-*bx); 
      /* Parabolic fit was no use. Use default magnification. */ 
      fu=(*func)(u);
    } else if ((*cx-u)*(u-ulim) > 0.0) { 
      /* Parabolic fit is between c and its allowed limit. */
      fu=(*func)(u);
      if (fu < *fc) {
        SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
        SHFT(*fb,*fc,fu,(*func)(u))
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) { 
      /* Limit parabolic u to maximum allowed value. */
      u=ulim;
      fu=(*func)(u);
    } else { 
      /* Reject parabolic u, use default magnification. */ 
      u=(*cx)+GOLD*(*cx-*bx); 
      fu=(*func)(u);
    }
    SHFT(*ax,*bx,*cx,u) /* Eliminate oldest point and continue. */
    SHFT(*fa,*fb,*fc,fu)
  }
}


double zbrent(double (*func)(double), double x1, double x2, double tol)
/* Using Brent's method, find the root of a function func known to 
   lie between x1 and x2. The root, returned as zbrent, will be 
   refined until its accuracy is tol. */
{
  int iter;
  double a=x1,b=x2,c=x2,d,e,min1,min2;
  double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    nrerror("Root must be bracketed in zbrent");
  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a; /* Rename a, b, c and adjust bounding interval d. */
      fc=fa;
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*ZEPS*fabs(b)+0.5*tol;  /* Convergence check. */
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa; /* Attempt inverse quadratic interpolation. */
      if (a == c) {
        p=2.0*xm*s;
        q=1.0-s;
      } else {
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q; /* Check whether in bounds. */
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=d; /* Accept interpolation. */
        d=p/q;
      } else {
        d=xm; /* Interpolation failed, use bisection. */
        e=d;
      }
    } else { /* Bounds decreasing too slowly, use bisection. */
      d=xm;
      e=d;
    }
    a=b; /* Move last best guess to a. */
    fa=fb;
    if (fabs(d) > tol1) /* Evaluate new trial root. */
      b +=d;
    else
      b += SIGN(tol1,xm);
    fb=(*func)(b);
  }
  nrerror("Maximum number of iterations exceeded in zbrent");
  return 0.0; /* Never get here. */
}

void bcucof(double y[], double y1[], double y2[], double y12[], 
            double d1, double d2, double **c)

/* Given arrays y[1..4], y1[1..4], y2[1..4], and y12[1..4], containing 
   the function, gradients, and cross derivative at the four grid 
   points of a rectangular grid cell (numbered counterclockwise from 
   the lower left), and given d1 and d2, the length of the grid cell 
   in the 1- and 2-directions, this routine returns the table 
   c[1..4][1..4] that is used by routine bcuint for bicubic 
   interpolation. 
*/
{
  static int wt[16][16]=
  { 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
    -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,
    2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,
    0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
    0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,
    0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,
    -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
    9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,
    -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,
    2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
    -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,
    4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1};

  int l,k,j,i;
  double xx,d1d2,cl[16],x[16];
  d1d2=d1*d2;
  for (i=1;i<=4;i++) { /* Pack a temporary vector x. */
    x[i-1]=y[i];
    x[i+3]=y1[i]*d1;
    x[i+7]=y2[i]*d2;
    x[i+11]=y12[i]*d1d2;
  }
  for (i=0;i<=15;i++) { /* Matrix multiply by the stored table. */
    xx=0.0;
    for (k=0;k<=15;k++) xx += wt[i][k]*x[k];
    cl[i]=xx;
  }
  l=0;
  for (i=1;i<=4;i++) /* Unpack the result into the output table. */
  for (j=1;j<=4;j++) c[i][j]=cl[l++];
}

void bcuint(double y[], double y1[], double y2[], double y12[], 
            double x1l, double x1u, double x2l, double x2u, 
            double x1, double x2, 
            double *ansy, double *ansy1, double *ansy2)
/* Bicubic interpolation within a grid square. Input quantities 
   are y,y1,y2,y12 (as described in bcucof); x1l and x1u, the lower 
   and upper coordinates of the grid square in the 1-direction;
   x2l and x2u likewise for the 2-direction; and x1,x2, the 
   coordinates of the desired point for the interpolation. The 
   interpolated function value is returned as ansy, and the 
   interpolated gradient values as ansy1 and ansy2. This routine 
   calls bcucof.
*/
{
  int i;
  double t,u,d1,d2,**c;
  c=dmatrix(1,4,1,4);
  d1=x1u-x1l;
  d2=x2u-x2l;
  bcucof(y,y1,y2,y12,d1,d2,c); /* Get the c's.*/
  if (x1u == x1l || x2u == x2l) nrerror("Bad input in routine bcuint");
  t=(x1-x1l)/d1; /* Equation (3.6.4). */
  u=(x2-x2l)/d2;
  *ansy=(*ansy2)=(*ansy1)=0.0;
  for (i=4;i>=1;i--) { /* Equation (3.6.6). */
    *ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
    *ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
    *ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];
  }
  *ansy1 /= d1;
  *ansy2 /= d2;

  free_dmatrix(c,1,4,1,4);
}



void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
            double *f, double stpmax, int *check, double (*func)(double []))
/*
Given an n-dimensional point xold[1..n], the value of the function and gradient there, fold
and g[1..n], and a direction p[1..n], finds a new point x[1..n] along the direction p from
xold where the function func has decreased "sufficiently." The new function value is returned
in f. stpmax is an input quantity that limits the length of the steps so that you do not try to
evaluate the function in regions where it is undefined or subject to overflow. p is usually the
Newton direction. The output quantity check is false (0) on a normal exit. It is true (1) when
x is too close to xold. In a minimization algorithm, this usually signals convergence and can
be ignored. However, in a zero-finding algorithm the calling program should check whether the
convergence is spurious. Some "difficult" problems may require double precision in this routine.
*/
{
  int i;
  double a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,
         test,tmplam;

  *check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if (sum > stpmax) for (i=1;i<=n;i++) p[i] *= stpmax/sum; /*Scale if attempted step is too big.*/
  for (slope=0.0,i=1;i<=n;i++) slope += g[i]*p[i];
  if (slope >= 0.0) nrerror("Roundoff problem in lnsrch.");
  test=0.0; /* Compute lambda min. */
  for (i=1;i<=n;i++) {
    temp=fabs(p[i])/MAXD(fabs(xold[i]),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;
  alam=1.0; /* Always try full Newton step first.*/
  for (;;) { /* Start of iteration loop. */
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
    *f=(*func)(x);
    if (alam < alamin) { /* Convergence on delta x. For zero finding,
                            the calling program should verify the convergence. */
      for (i=1;i<=n;i++) x[i]=xold[i];
      *check=1;
      return;
    } else if (*f <= fold+ALF*alam*slope) return; /* Sufficient function decrease.*/
    else { /* Backtrack.*/
      if (alam == 1.0) tmplam = -slope/(2.0*(*f-fold-slope)); /* First time.*/
      else { /* Subsequent backtracks. */
        rhs1 = *f-fold-alam*slope;
        rhs2=f2-fold-alam2*slope;
        a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
        b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
        if (a == 0.0) tmplam = -slope/(2.0*b);
        else {
          disc=b*b-3.0*a*slope;
          if (disc < 0.0) tmplam=0.5*alam;
          else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
          else tmplam=-slope/(b+sqrt(disc));
        }
        if (tmplam > 0.5*alam) tmplam=0.5*alam;   /* lamda <= 0.5 lamda_1.*/
      }
    }
    alam2=alam;
    f2 = *f;
    alam=MAXD(tmplam,0.1*alam);  /* lamda >= 0.1 lamda_1.*/
  } /* Try again. */
}




void newt(double x[], int n, int *check, void (*vecfunc)(int, double [], double []))
/*
Given an initial guess x[1..n] for a root in n dimensions, find the root by a globally convergent
Newton's method. The vector of functions to be zeroed, called fvec[1..n] in the routine
below, is returned by the user-supplied routine vecfunc(n,x,fvec). The output quantity
check is false (0) on a normal return and true (1) if the routine has converged to a local
minimum of the function fmin defined below. In this case try restarting from a different initial
guess.
*/
{
  void fdjac(int n, double x[], double fvec[], double **df,
  void (*vecfunc)(int, double [], double []));
  double fmin(double x[]);
  void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
              double *f, double stpmax, int *check, double (*func)(double []));
  void lubksb(double **a, int n, int *indx, double b[]);
  void ludcmp(double **a, int n, int *indx, double *d);

  int i,its,j,*indx;
  double d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;
  indx=ivector(1,n);
  fjac=dmatrix(1,n,1,n);
  g=dvector(1,n);
  p=dvector(1,n);
  xold=dvector(1,n);
  fvec=dvector(1,n); /* Define global variables.*/
  nn=n;
  nrfuncv=vecfunc;
  f=fmin(x); /*fvec is also computed by this call.*/
  test=0.0; /*Test for initial guess being a root. Use
              more stringent test than simply TOLF.*/
  for (i=1;i<=n;i++)
    if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
  if (test < 0.01*TOLF) {
    *check=0;
    FREERETURN
  }
  for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]); /* Calculate stpmax for line searches. */
  stpmax=STPMX*MAXD(sqrt(sum),(float)n);
  for (its=1;its<=MAXITS;its++) { /* Start of iteration loop.*/
    fdjac(n,x,fvec,fjac,vecfunc);
    /* If analytic Jacobian is available, you can replace the routine fdjac below with your own routine. */
    for (i=1;i<=n;i++) { /*Compute grad(f) for the line search.*/
      for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
      g[i]=sum;
    }
    for (i=1;i<=n;i++) xold[i]=x[i]; /*Store x,*/
    fold=f; /*and f.*/
    for (i=1;i<=n;i++) p[i] = -fvec[i]; /*Right-hand side for linear equations.*/
    ludcmp(fjac,n,indx,&d); /*Solve linear equations by LU decomposition.*/
    lubksb(fjac,n,indx,p);
    lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin);
    /*lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.*/
    test=0.0; /*Test for convergence on function values.*/
    for (i=1;i<=n;i++) if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < TOLF) {
      *check=0;
       FREERETURN
    }
    if (*check) { /* Check for gradient of f zero, i.e., spurious convergence.*/
      test=0.0;
      den=MAXD(f,0.5*n);
      for (i=1;i<=n;i++) {
        temp=fabs(g[i])*MAXD(fabs(x[i]),1.0)/den;
        if (temp > test) test=temp;
      }
      *check=(test < TOLMIN ? 1 : 0);
      FREERETURN
    }
    test=0.0; /*Test for convergence on dx.*/
    for (i=1;i<=n;i++) {
      temp=(fabs(x[i]-xold[i]))/MAXD(fabs(x[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX) FREERETURN
  }
  nrerror("MAXITS exceeded in newt");
}



void fdjac(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double []))
/*
Computes forward-difference approximation to Jacobian. On input, x[1..n] is the point at
which the Jacobian is to be evaluated, fvec[1..n] is the vector of function values at the
point, and vecfunc(n,x,f) is a user-supplied routine that returns the vector of functions at
x. On output, df[1..n][1..n] is the Jacobian array.
*/
{
  int i,j;
  double h,temp,*f;
  f=dvector(1,n);
  for (j=1;j<=n;j++) {
    temp=x[j];
    h=EPS*fabs(temp);
    if (h == 0.0) h=EPS;
    x[j]=temp+h; /* Trick to reduce finite precision error. */
    h=x[j]-temp;
    (*vecfunc)(n,x,f);
    x[j]=temp;
    for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h; /* Forward difference formula.*/
  }
  free_dvector(f,1,n);
}



double fmin(double x[])
/*
Returns f = 1
2 F  F at x. The global pointer *nrfuncv points to a routine that returns the
vector of functions at x. It is set to point to a user-supplied routine in the calling program.
Global variables also communicate the function values back to the calling program.
*/
{
  int i;
  double sum;
  (*nrfuncv)(nn,x,fvec);
  for (sum=0.0,i=1;i<=nn;i++) sum += DSQR(fvec[i]);
  return 0.5*sum;
}



void ludcmp(double **a, int n, int *indx, double *d)
/*
Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
indx[1..n] is an output vector that records the row permutation effected by the partial
pivoting; d is output as  1 depending on whether the number of row interchanges was even
or odd, respectively. This routine is used in combination with lubksb to solve linear equations
or invert a matrix.
*/
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv; /* vv stores the implicit scaling of each row. */
  vv=dvector(1,n);
  *d=1.0; /* No row interchanges yet. */
  for (i=1;i<=n;i++) { /* Loop over rows to get the implicit scaling information. */
    big=0.0;
    for (j=1;j<=n;j++) if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
    /* No nonzero largest element. */
    vv[i]=1.0/big; /* Save the scaling. */
  }
  for (j=1;j<=n;j++) { /* This is the loop over columns of Crout's method. */
    for (i=1;i<j;i++) { /* This is equation (2.3.12) except for i = j. */
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0; /* Initialize for the search for largest pivot element. */
    for (i=j;i<=n;i++) { /* This is i = j of equation (2.3.12) and i = j +1: ::N
                            of equation (2.3.13).*/
      sum=a[i][j];
      for (k=1;k<j;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
        /* Is the gure of merit for the pivot better than the best so far? */
        big=dum;
        imax=i;
      }
    }
    if (j != imax) { /* Do we need to interchange rows?  */
      for (k=1;k<=n;k++) { /* Yes, do so... */
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d); /* ...and change the parity of d.  */
      vv[imax]=vv[j]; /* Also interchange the scale factor. */
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    /* If the pivot element is zero the matrix is singular (at least to the precision of the
       algorithm). For some applications on singular matrices, it is desirable to substitute
       TINY for zero. */
    if (j != n) { /* Now, finally, divide by the pivot element. */
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  } /* Go back for the next column in the reduction. */
  free_dvector(vv,1,n);
}


void lubksb(double **a, int n, int *indx, double b[])
/*
Solves the set of n linear equations A  X = B. Here a[1..n][1..n] is input, not as the matrix
A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
B, and returns with the solution vector X. a, n, and indx are not modified by this routine
and can be left in place for successive calls with different right-hand sides b. This routine takes
into account the possibility that b will begin with many zero elements, so it is efficient for use
in matrix inversion.
*/
{
  int i,ii=0,ip,j;
  double sum;
  for (i=1;i<=n;i++) { /* When ii is set to a positive value, it will become the
                          index of the first nonvanishing element of b. We now
                          do the forward substitution, equation (2.3.6). The
                          only new wrinkle is to unscramble the permutation
                          as we go. */
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii) for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i; /* A nonzero element was encountered, so from now on we
                           will have to do the sums in the loop above.*/
    b[i]=sum;
  }
  for (i=n;i>=1;i--) { /* Now we do the backsubstitution, equation (2.3.7). */
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i]; /* Store a component of the solution vector X. */
  } /* All done! */
}

double maxvec(double *vec,int svec) {
  int i;
  double res;
  res=vec[1];
  for (i=2;i<=svec;i++) {
    res=MAXD(vec[i],res);
  }
  return res;
}

double minvec(double *vec,int svec) {
  int i;
  double res;
  res=vec[1];
  for (i=2;i<=svec;i++) {
    res=MIND(vec[i],res);
  }
  return res;
}

void nrerror(char error_text[])
{
  fprintf(stderr," *** Numerical Recipes RUN-TIME ERROR: %s\n",error_text);
  fprintf(stderr," *** NOW EXITING TO SYSTEM ... FUUUUUUCKK!!\n");
  scanf("?");
  exit(1);
}

/* subscript ranges go v[nl...nh] etc...*/

float *vector(long nl, long nh)
{
  float *v;

  v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
  if (!v) nrerror("MEMORY ALLOCATION FAILURE IN vector()");
  return v-nl+NR_END;
}

int *ivector(long nl, long nh)
{
  int *v;

  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) nrerror("MEMORY ALLOCATION FAILURE IN ivector()");
  return v-nl+NR_END;
}

unsigned int *uivector(long nl, long nh)
{
  unsigned int *v;

  v=(unsigned int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned int)));
  if (!v) nrerror("MEMORY ALLOCATION FAILURE IN uivector()");
  return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
{
  unsigned char *v;

  v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
  if (!v) nrerror("MEMORY ALLOCATION FAILURE IN cvector()");
  return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
{
  unsigned long *v;

  v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
  if (!v) nrerror("MEMORY ALLOCATION FAILURE IN lvector()");
  return v-nl+NR_END;
}

double *dvector(long nl, long nh)
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) nrerror("MEMORY ALLOCATION FAILURE IN dvector()");
  return v-nl+NR_END;
}

real *rvector(long nl, long nh)
{
  real *v;
  v=(real *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(real)));
  if (!v) nrerror("MEMORY ALLOCATION FAILURE IN rvector()");
  return v-nl+NR_END;
}

fcomplex *complex_vector(long nl, long nh)
{
  fcomplex *v;

  v=(fcomplex *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(fcomplex)));
  if (!v) nrerror("MEMORY ALLOCATION FAILURE IN complex_vector()");
  return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
{
  long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("MEMORY ALLOCATION FAILURE 1 IN matrix()");

  m+=NR_END;
  m-=nrl;

  m[nrl]=(float*)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl]) nrerror("MEMORY ALLOCATION FAILURE 2 IN matrix()");

  m[nrl]+=NR_END;
  m[nrl]-=ncl;

  for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  return m;
}

real **rmatrix(long nrl, long nrh, long ncl, long nch)
{
  long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  real **m;

  m=(real **)malloc((size_t)((nrow+NR_END)*sizeof(real*)));
  if (!m) nrerror("MEMORY ALLOCATION FAILURE 1 IN rmatrix()");

  m+=NR_END;
  m-=nrl;

  m[nrl]=(real*)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("MEMORY ALLOCATION FAILURE 2 IN rmatrix()");
 
  m[nrl]+=NR_END;
  m[nrl]-=ncl;

  for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
  long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  m=(double **)malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) nrerror("MEMORY ALLOCATION FAILURE 1 IN dmatrix()");

  m+=NR_END;
  m-=nrl;

  m[nrl]=(double*)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) nrerror("MEMORY ALLOCATION FAILURE 2 IN dmatrix()");

  m[nrl]+=NR_END;
  m[nrl]-=ncl;

  for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
{
  long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  m=(int **)malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) nrerror("MEMORY ALLOCATION FAILURE 1 IN matrix()");

  m+=NR_END;
  m-=nrl;

  m[nrl]=(int*)malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) nrerror("MEMORY ALLOCATION FAILURE 2 IN matrix()");

  m[nrl]+=NR_END;
  m[nrl]-=ncl;

  for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch]*/
{
  long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
  float **m;

  m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("MEMORY ALLOCATION FAILURE IN submatrix()");

  m+=NR_END;
  m-=newrl;

  for (i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

  return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard c manner as a[nrow][ncol].  The routine should be 
called with the address &a[0][0] as the first argument. */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m) nrerror("MEMORY ALLOCATION FAILURE IN convert_matrix()");

  m+=NR_END;
  m-=nrl;

  m[nrl]=a-ncl;
  for (i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;

  return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;

  t=(float ***)malloc((size_t)((nrow+NR_END)*sizeof(float**)));
  if (!t) nrerror("MEMORY ALLOCATION FAILURE 1 IN f3tensor()");
  t+=NR_END;
  t-=nrl;

  t[nrl]=(float**)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
  if (!t[nrl]) nrerror("MEMORY ALLOCATION FAILURE 2 IN f3tensor()");
  t[nrl]+=NR_END;
  t[nrl]-=ncl;

  t[nrl][ncl]=(float*)malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
  if (!t[nrl][ncl]) nrerror("MEMORY ALLOCATION FAILURE 3 IN f3tensor()");
  t[nrl][ncl]+=NR_END;
  t[nrl][ncl]-=ndl;

  for (j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for (i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
	t[i][ncl]=t[i-1][ncl]+ncol*ndep;
	for (j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  return t;
}


void free_vector(float *v, long nl, long nh)
{
  free((FREE_ARG)(v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
{
  free((FREE_ARG)(v+nl-NR_END));
}

void free_uivector(unsigned int *v, long nl, long nh)
{
  free((FREE_ARG)(v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
{
  free((FREE_ARG)(v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
{
  free((FREE_ARG)(v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
{
  free((FREE_ARG)(v+nl-NR_END));
}

void free_rvector(real *v, long nl, long nh)
{
  free((FREE_ARG)(v+nl-NR_END));
}

void free_complex_vector(fcomplex *v, long nl, long nh)
{
  free((FREE_ARG)(v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG)(m[nrl]+ncl-NR_END));
  free((FREE_ARG)(m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG)(m[nrl]+ncl-NR_END));
  free((FREE_ARG)(m+nrl-NR_END));
}

void free_rmatrix(real **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG)(m[nrl]+ncl-NR_END));
  free((FREE_ARG)(m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG)(m[nrl]+ncl-NR_END));
  free((FREE_ARG)(m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG)(b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG)(b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  free((FREE_ARG)(t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG)(t[nrl]+ncl-NR_END));
  free((FREE_ARG)(t+nrl-NR_END));
}

/*changed to double precision*/

fcomplex Cadd(fcomplex a, fcomplex b) 
{
  fcomplex c;
  c.r=a.r+b.r;
  c.i=a.i+b.i;
  return c;
}

fcomplex Csub(fcomplex a, fcomplex b)
{
  fcomplex c;
  c.r=a.r-b.r;
  c.i=a.i-b.i;
  return c;
}

fcomplex Cmul(fcomplex a, fcomplex b)
{
  fcomplex c;
  c.r=a.r*b.r-a.i*b.i;
  c.i=a.i*b.r+a.r*b.i;
  return c;
}

fcomplex Complex(double re, double im)
{
  fcomplex c;
  c.r=re;
  c.i=im;
  return c;
}

fcomplex Conjg(fcomplex z)
{
  fcomplex c;
  c.r=z.r;
  c.i=-z.i;
  return c;
}

fcomplex Cdiv(fcomplex a, fcomplex b) 
{
  fcomplex c;
  double r, den;
  if (fabs(b.r) >= fabs(b.i)) {
	r=b.i/b.r;
	den=b.r+r*b.i;
	c.r=(a.r+r*a.i)/den;
	c.i=(a.i-r*a.r)/den;
  } else {
    r=b.r/b.i;
	den=b.i+r*b.r;
	c.r=(a.r*r+a.i)/den;
	c.i=(a.i*r-a.r)/den;
  }
  return c;
}

double Cabs(fcomplex z)
{
  double x,y,ans,temp;
  x=fabs(z.r);
  y=fabs(z.i);
  if (x==0.0) ans=y;
  else if (y==0.0) ans=x;
  else if (x>y) {
	temp=y/x;
	ans=x*sqrt(1.0+temp*temp);
  } else {
	temp=x/y;
	ans=y*sqrt(1.0+temp*temp);
  }
  return ans;
}

fcomplex Csqrt(fcomplex z)
{
  fcomplex c;
  double x,y,w,r;
  if ((z.r==0.0) && (z.i==0.0)) {
	c.r=0.0;
	c.i=0.0;
	return c;
  } else {
    x=fabs(z.r);
	y=fabs(z.i);
	if (x>=y) {
      r=y/x;
	  w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
    } else {
	  r=x/y;
	  w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
	}
    if (z.r >= 0.0) {
      c.r=w;
	  c.i=z.i/(2.0*w);
	} else {
	  c.i=(z.i>=0) ? w : -w;
	  c.r=z.i/(2.0*c.i);
    }
	return c;
  }
}

fcomplex RCmul(double x, fcomplex a) {
  fcomplex c;
  c.r=x*a.r;
  c.i=x*a.i;
  return c;
}

