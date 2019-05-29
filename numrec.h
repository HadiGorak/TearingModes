#include <math.h>
#include <stdlib.h>

#ifndef _COMPLEX_H_
#define _COMPLEX_H_

typedef struct FCOMPLEX {double r,i;} fcomplex;

fcomplex Cadd(fcomplex a, fcomplex b);
fcomplex Csub(fcomplex a, fcomplex b);
fcomplex Cmul(fcomplex a, fcomplex b);
fcomplex Complex(double re, double im);
fcomplex Conjg(fcomplex z);
fcomplex Cdiv(fcomplex a, fcomplex b) ;
double Cabs(fcomplex z);
fcomplex Csqrt(fcomplex z);
fcomplex RCmul(double x, fcomplex a);

#endif

void spline(double x[], double y[], int n, 
             double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], 
             int n, double x, double *y);
void splie2(double x1a[], double x2a[], double **ya, 
             int m, int n, double **y2a);
void splin2(double x1a[], double x2a[], 
            double **ya, double **y2a, 
            int m, int n,
            double x1, double x2, 
            double *y);
void powell(double p[], double **xi, int n, double ftol, int *iter,
            double *fret, double (*func)(double []));
void linmin(double p[], double xi[], int n, double *fret,
            double (*func)(double []));
double brent(double ax, double bx, double cx, 
            double (*f)(double), double tol, double *xmin);
void mnbrak(double *ax, double *bx, double *cx, 
            double *fa, double *fb, double *fc,
            double (*func)(double));
double f1dim(double x);
double zbrent(double (*func)(double), double x1, double x2, double tol);
void bcucof(double y[], double y1[], double y2[], double y12[],
            double d1, double d2, double **c);
void bcuint(double y[], double y1[], double y2[], double y12[],
            double x1l, double x1u, double x2l, double x2u,
            double x1, double x2,
            double *ansy, double *ansy1, double *ansy2);
void lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
            double *f, double stpmax, int *check, double (*func)(double []));
void newt(double x[], int n, int *check,void (*vecfunc)(int, double [], double []));
void fdjac(int n, double x[], double fvec[], double **df, void (*vecfunc)(int, double [], double []));
double fmin(double x[]);
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);
double maxvec(double *vec,int svec);
double minvec(double *vec,int svec);

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

typedef double real;

#define NR_END 1
#define FREE_ARG char*

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0: dsqrarg*dsqrarg)

static real rmaxarg1,rmaxarg2;
#define MAXR(a,b) (rmaxarg1=(a),rmaxarg2=(b),(rmaxarg1) > (rmaxarg2) ? (rmaxarg1) : (rmaxarg2))

static real rminarg1,rminarg2;
#define MINR(a,b) (rminarg1=(a),rminarg2=(b),(rminarg1) < (rminarg2) ? (rminarg1) : (rminarg2))

static double dmaxarg1,dmaxarg2;
#define MAXD(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define MIND(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define MAXF(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define MINF(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define MAXL(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ? (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define MINL(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ? (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define MAXI(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define MINI(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


void   nrerror(char error_text[]);

float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned int *uivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
real *rvector(long nl, long nh);
fcomplex *complex_vector(long nl, long nh);

float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
real **rmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch, 
                  long newrl, long newcl);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_uivector(unsigned int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_rvector(real *v, long nl, long nh);
void free_complex_vector(fcomplex *v, long nl, long nh);

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_rmatrix(real **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

#endif 


