
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "numrec.h"


struct PtEnv {
 //unsigned int loc : 1;         /* Location - INSIDE or OUTSIDE the f.c.  */
 //unsigned int pos : 1;         /* Position - INTERNAL or near BOUNDARY.  */
 //unsigned int vi  : 1;         /* Very Internal flag */
 double dwdt[2][2], dpol[2][2], dnc[2][2];
 double w[2];
};
void makepole(double nrm, double offset, double maxy,
              double betamin, double betamax,
              int pdn, double *dp, double *bn);
void read_dp(double *bmin, double *bmax, unsigned int *pdn, 
             double *dp, double *bn);
double dwdt(double t_in, double w);
double dwdt_func(double w);
double D_pol_term(double w);
double D_tot_term(double w);
double D_prm_term(double t_in, double w);
double getbetan(double betanpre, double dt, double dbetadt, double w,
                double betamin, double betamax);
void store_dwdt_vs_w();
