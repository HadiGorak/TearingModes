
#include "dwdt.h"

unsigned int iprev,dpt,bnt,dpol_mod,dtot_mod,i,j,dp_w_pow,dp_w_tan,
             dp_w_w,numsteps,numsav,dpnum,bnnum,wnum,tmax;
double *betan,*dp,*dpp2,*t,*betan_t,*betan_t_f,*dp_t,*dp_t_f,*w_t,*dwdt_t,
       *wthrs_t,*bn_file,*t_file,*dpb,*dpb2;
double D_nc,D_I,D_R,D_pol,D_tot,w_d,k_0,eta_star,d_star,w_0,w_b,
       alpha_l,alpha_s,dt,tau,dbetadt,betamin,betamax,nrm,offset,maxy,dcy,
       bnf_min,bnf_max,dbwfac,w_sat,H,betan_sat,alphadp,dp_bnsat,m;
FILE *fpw;

int main (int argc, char **argv) {

  FILE *fp;
  char oneline[1024];
  double k1,k2,k3,k4,dy1,dyn,f1,f2,alpha,beta,deltadp;

  /* set defaults */
  D_nc=0.0;
  D_I=-0.5;
  D_R=-0.25;
  D_tot=0.0;
  D_pol=0.0;
  eta_star=1.0;
  k_0=1.0;
  w_d=0.01;
  w_b=0.02;
  w_sat=0.25;
  dpt=1;
  bnt=0;
  dpol_mod=1;
  dtot_mod=1;
  dpnum=100;
  wnum=100;
  w_0=0.001;
  dbetadt=1.0;
  betamin=1.0;
  betamax=2.0;
  betan_sat=betamin+0.75*(betamax-betamin);
  dp_w_pow=0;
  dp_w_tan=0;
  dp_w_w=0;
  nrm=1.0;
  maxy=100.;
  dcy=0.0;
  offset=0.0;
  dt=0.001; /* take 1 ms timesteps */
  tau=1.0; /* integrate for 1 sec */
  numsteps=tau/dt;
  if (numsteps<tau/dt) numsteps+=1;
  numsav=numsteps/1000;
  H=0.;
  m=2;

  /* eta_star is dummy, actually in conversations with Tom G. we decided it's 
     probably good to insert the t_r/r^2 factor in place of the k0/eta_star.
     Since r^2/t_r for the 2/1 in 098549.02030 is about 1.6 eta_star is
     set to that in that case, and k_0 is set to 1 */

  /* parse and store command line input */
  for (i=1;i<=argc-1;i++) {
    if (!strcmp(argv[i],"-dpt")) {
      i++;
      sscanf(argv[i],"%i",&dpt);
    } else if (!strcmp(argv[i],"-bnt")) {
      i++;
      sscanf(argv[i],"%i",&bnt);
    } else if (!strcmp(argv[i],"-dpolmod")) {
      i++;
      sscanf(argv[i],"%i",&dpol_mod);
    } else if (!strcmp(argv[i],"-dtotmod")) {
      i++;
      sscanf(argv[i],"%i",&dtot_mod);
    } else if (!strcmp(argv[i],"-dpn")) {
      i++;
      sscanf(argv[i],"%i",&dpnum);
    } else if (!strcmp(argv[i],"-wn")) {
      i++;
      sscanf(argv[i],"%i",&wnum);
    } else if (!strcmp(argv[i],"-w0")) {
      i++;
      sscanf(argv[i],"%lf",&w_0);
    } else if (!strcmp(argv[i],"-wb")) {
      i++;
      sscanf(argv[i],"%lf",&w_b);
    } else if (!strcmp(argv[i],"-wsat")) {
      i++;
      sscanf(argv[i],"%lf",&w_sat);
    } else if (!strcmp(argv[i],"-betansat")) {
      i++;
      sscanf(argv[i],"%lf",&betan_sat);
    } else if (!strcmp(argv[i],"-dp_w_pow")) {
      i++;
      sscanf(argv[i],"%i",&dp_w_pow);
    } else if (!strcmp(argv[i],"-dp_w_tan")) {
      i++;
      sscanf(argv[i],"%i",&dp_w_tan);
    } else if (!strcmp(argv[i],"-dp_w_w")) {
      i++;
      sscanf(argv[i],"%i",&dp_w_w);
    } else if (!strcmp(argv[i],"-m")) {
      i++;
      sscanf(argv[i],"%lf",&m);
    } else if (!strcmp(argv[i],"-nrm")) {
      i++;
      sscanf(argv[i],"%lf",&nrm);
    } else if (!strcmp(argv[i],"-ofst")) {
      i++;
      sscanf(argv[i],"%lf",&offset);
    } else if (!strcmp(argv[i],"-maxy")) {
      i++;
      sscanf(argv[i],"%lf",&maxy);
    } else if (!strcmp(argv[i],"-dcy")) {
      i++;
      sscanf(argv[i],"%lf",&dcy);
    } else if (!strcmp(argv[i],"-dtt")) {
      i++;
      sscanf(argv[i],"%lf",&dt);
      i++;
      sscanf(argv[i],"%lf",&tau);
      numsteps=tau/dt;
      if (numsteps<tau/dt) numsteps+=1;
      numsav=numsteps/1000;
    } else if (!strcmp(argv[i],"-db")) {
      i++;
      sscanf(argv[i],"%lf",&dbetadt);
    } else if (!strcmp(argv[i],"-bmn")) {
      i++;
      sscanf(argv[i],"%lf",&betamin);
    } else if (!strcmp(argv[i],"-bmx")) {
      i++;
      sscanf(argv[i],"%lf",&betamax);
    } else if (!strcmp(argv[i],"-dnc")) {
      i++;
      sscanf(argv[i],"%lf",&D_nc);
    } else if (!strcmp(argv[i],"-di")) {
      i++;
      sscanf(argv[i],"%lf",&D_I);
    } else if (!strcmp(argv[i],"-dr")) {
      i++;
      sscanf(argv[i],"%lf",&D_R);
    } else if (!strcmp(argv[i],"-dtot")) {
      i++;
      sscanf(argv[i],"%lf",&D_tot);
    } else if (!strcmp(argv[i],"-dpol")) {
      i++;
      sscanf(argv[i],"%lf",&D_pol);
    } else if (!strcmp(argv[i],"-hfac")) {
      i++;
      sscanf(argv[i],"%lf",&H);
    } else if (!strcmp(argv[i],"-k0")) {
      i++;
      sscanf(argv[i],"%lf",&k_0);
    } else if (!strcmp(argv[i],"-wd")) {
      i++;
      sscanf(argv[i],"%lf",&w_d);
    } else if (!strcmp(argv[i],"-es")) {
      i++;
      sscanf(argv[i],"%lf",&eta_star);
    } else if (!strcmp(argv[i],"-h")) {
      printf("\nIsland Evolution Integrator with time dependent Delta Prime\n");
      printf("\n Input: the parameter values to the modified Rutherford");
      printf("\n        model and a delta_prime as a function of time");
      printf("\nOutput: dw/dt and w as functions of time");
      printf("\n\nSwitch     Choices and Meaning");
      printf("\n----------------------------------------");
      printf("\n-dpt       delta_prime(t) model 0=file pole.dat");
      printf("\n           with the number of points as first line");
      printf("\n           and two columns: betan delta_prime;");
      printf("\n           1=analytic with parameters set by input");
      printf("\n           arguments with delta_prime linear decay");
      printf("\n           after pole imposable; 2=analytic with");
      printf("\n           params set by input with betan linear decay");
      printf("\n           and delta_prime determined from betan (not");
      printf("\n           yet implemented, useful?); 3=analytic with");
      printf("\n           params set by input with dp(w) after betan");
      printf("\n           reaches max");
      printf("\n-bnt       betan(t) model 0=linear ramp up at db from");
      printf("\n           bmn to bmx; 1=from file with two columns");
      printf("\n           and first line with how many rows");
      printf("\n           first column is time in ms and second");
      printf("\n           column is betan");
      printf("\n-dpn       number of points in the analytical pole");
      printf("\n-wn        number of points in w space");
      printf("\n-w0        initial island width");
      printf("\n-nrm       delta_prime peak value");
      printf("\n-ofst      the offset of the pole curve");
      printf("\n-maxy      the maxy of the pole curve (See code)");
      printf("\n-dcy       the decay rate delta-prime/s after pole (0 direct)");
      printf("\n-dpolmod   model used for Dpol term (see code)");
      printf("\n-dtotmod   model used for Dtot term (see code)");
      printf("\n-dp_w_pow  set to 1 for delta-star form of dp(w)");
      printf("\n-dp_w_tan  set to 1 for tanh form dp(w)");
      printf("\n-dp_w_w    set to 1 for dp-a*w form of dp(w)");
      printf("\n-m         m poloidal mode number");
      printf("\n-dtt       timestep and total time in seconds");
      printf("\n-db        dbetan/dt initial ramp rate");
      printf("\n-bmn       beta min");
      printf("\n-bmx       beta max");
      printf("\n-dnc       D_nc");
      printf("\n-di        D_I");
      printf("\n-dr        D_R");
      printf("\n-dtot      D_tot");
      printf("\n-dpol      D_pol");
      printf("\n-hfac      H");
      printf("\n-k0        k_0");
      printf("\n-wd        w_d thermal diffusivity scale length");
      printf("\n-wb        w_b banana width");
      printf("\n-wsat      saturated island width apriori estimate");
      printf("\n-betansat  betan at saturated island width apriori estimate");
      printf("\n-es        eta_star");
      printf("\n\n");
      exit(0);
    }
  }

  /***************************************************/
  /* Set Up Input Dependent Variables                */
  /***************************************************/


  /*printf("%i %i %lf %lf\n",numsteps,numsav,tau,dt);*/

  /* allocate time dependent arrays */
  t=dvector(1,numsteps);
  dp_t=dvector(1,numsteps);
  dp_t_f=dvector(1,numsteps);
  betan_t=dvector(1,numsteps);
  betan_t_f=dvector(1,numsteps);
  w_t=dvector(1,numsteps);
  dwdt_t=dvector(1,numsteps);
  wthrs_t=dvector(1,numsteps);

  /* create the delta-prime(beta) function */
  if (dpt==0) {

    /* read_dp(&betamin,&betamax,&dpnum,dp,betan); */

    /* open and read the pole and equilibrium data */

    fp=fopen("pole.dat","r");
    fgets(oneline,1024,fp);
    sscanf(oneline,"%i\n",&dpnum);

    dp=dvector(1,dpnum);
    betan=dvector(1,dpnum);

    i=1;
    while (fgets(oneline,1024,fp)!=NULL) {
      sscanf(oneline,"%lf %lf\n",&betan[i],&dp[i]);
      dp[i]=nrm*dp[i];
      i++;
    }
    fclose(fp);

    betamin=minvec(betan,dpnum);
    betamax=maxvec(betan,dpnum);

  } else {
    dp=dvector(1,dpnum);
    betan=dvector(1,dpnum);
    makepole(nrm,offset,maxy,betamin,betamax,dpnum,dp,betan);
  }

  /* setup the Mercier parameters */
  alpha_l=0.5-sqrt(-D_I);
  alpha_s=0.5+sqrt(-D_I);

  /* if the betan time dependence is taken from an input file read it in */
  if (bnt==1) {
    /* read betan(t) from a file */
    fp=fopen("betan.dat","r");
    fgets(oneline,1024,fp);
    sscanf(oneline,"%i\n",&bnnum);
    t_file=dvector(1,bnnum);
    bn_file=dvector(1,bnnum);

    i=1;
    while (fgets(oneline,1024,fp)!=NULL) {
      sscanf(oneline,"%lf %lf\n",&t_file[i],&bn_file[i]);
      t_file[i]/=1000.0;
      i++;
    }
    fclose(fp);

    bnf_min=minvec(bn_file,bnnum);
    bnf_max=maxvec(bn_file,bnnum);

    if (bnf_max<betamax) {

      if (dt*(numsteps-1)+t_file[1]<=t_file[bnnum]) {

        /* setup and spline the values of betan to the integration times */
        dpp2=dvector(1,bnnum);
        dy1=(2*bn_file[4]-9*bn_file[3]+18*bn_file[2]-11*bn_file[1])/
            (6*(t_file[2]-t_file[1]));
        dyn=(-2*bn_file[bnnum-3]+9*bn_file[bnnum-2]
             -18*bn_file[bnnum-1]+11*bn_file[bnnum])/
            (6*(t_file[2]-t_file[1]));
        spline(t_file,bn_file,bnnum,dy1,dyn,dpp2); 
    
        for (i=1;i<=numsteps;i++) {
          t[i]=t_file[1]+dt*(i-1);
          splint(t_file,bn_file,dpp2,bnnum,t[i],&betan_t[i]);
        }

        free_dvector(dpp2,1,bnnum);

      } else {
        printf("\n\nERROR: file betan.dat contains time values out of range");
        printf("\n\n");
        exit(0);
      }

    } else {
      printf("\n\nERROR: file betan.dat contains beta values out of range\n");
      exit(0);
    }
  }


  /* set up a spline of dp(beta) using finite diff derivs at endpoints */
  dpb2=dvector(1,dpnum);
  dy1=(2*dp[4]-9*dp[3]+18*dp[2]-11*dp[1])/(6*(betan[2]-betan[1]));
  dyn=(-2*dp[dpnum-3]+9*dp[dpnum-2]
       -18*dp[dpnum-1]+11*dp[dpnum])/
      (6*(betan[2]-betan[1]));
  spline(betan,dp,dpnum,dy1,dyn,dpb2); 

  /* set up a spline of beta(dp) for evaluation of betan_sat estimate */
  dpp2=dvector(1,dpnum);
  if (dp[2]==dp[1]) dy1=1e30;
  else dy1=(betan[2]-betan[1])/(dp[2]-dp[1]);
  if (dp[dpnum]==dp[dpnum-1]) dyn=1e30;
  else dyn=(betan[dpnum]-betan[dpnum-1])/(dp[dpnum]-dp[dpnum-1]);
  spline(dp,betan,dpnum,dy1,dyn,dpp2); 

  /* note that choosing 
  
  delta'(betan_sat)-alphadp*w_sat=-m 

  and 

  dwdt(w_sat)=0 

  gives a third order polynomial for w_s: 

  w_s^3 - D_tot/m w_s^2 - D_pol/m = 0

  then use the root to get betan_sat.  For now estimate with w_sat set apriori.

  Compare to standard w_sat method

  */

  /* splint(dp,betan,dpp2,dpnum,alphadp*w_sat-m,&betan_sat); */

  /* get the raw dp(beta) at betan_sat */
  splint(betan,dp,dpb2,dpnum,betan_sat,&dp_bnsat);

  /* get alphadp for delta-prime(w) equation dp(betan_sat)-alphadp*w_sat=-m */
  alphadp=(dp_bnsat+m)/w_sat;

  /* alphadp gets used in D_prm_term */

  /* open a file for writing frames of dwdt vs. w for each time */
  fpw=fopen("dwdt_w.out","w");
  fprintf(fpw,"%i\n",wnum); /* the first line is the frame length */

  /*******************************************/ 
  /* integrate the island evolution equation */ 
  /*******************************************/

  t[1]=0.;
  w_t[1]=w_0;
  dwdt_t[1]=0.0;
  fp=fopen("w.out","w");
  iprev=0;
  for (i=2;i<=numsteps;i++) {

    t[i]=t[i-1]+dt;

    /* use 4th order runge kutta method to integrate dwdt(w,t) */

    k1=dt*dwdt(t[i-1],w_t[i-1]);
    k2=dt*dwdt(t[i-1]+0.5*dt,w_t[i-1]+0.5*k1);
    k3=dt*dwdt(t[i-1]+0.5*dt,w_t[i-1]+0.5*k2);
    k4=dt*dwdt(t[i-1]+dt,w_t[i-1]+k3);
    dwdt_t[i]=(k1+2*k2+2*k3+k4)/6;      /* this is actually dw here */
    w_t[i]=w_t[i-1]+dwdt_t[i];
    dwdt_t[i]/=dt;                      /* change this to dwdt from dw */

    if (w_t[i]<w_0) w_t[i]=w_0; /* we don't want to go below w_0 */

    if (i-iprev==numsav) {
      iprev=i;

      fprintf(fp,"%lf %lf %lf %lf %lf %lf %lf\n",
              t[i],w_t[i],wthrs_t[i],dwdt_t[i],dp_t_f[i],betan_t_f[i],
              (betan_t_f[i] - betan_t_f[i-1])/dt);

      store_dwdt_vs_w();

    /* printf("%lf %lf %lf %lf %lf %lf\n",
                t[i],w_t[i],wthrs_t[i],dwdt_t[i],dp_t_f[i],betan_t[i]);
      exit(0);*/
    }

  }

  fclose(fp);
  fclose(fpw);

  /* clean up and quit */
  free_dvector(dp,1,dpnum);
  free_dvector(betan,1,dpnum);
  free_dvector(dpp2,1,dpnum);
  free_dvector(t,1,numsteps);
  free_dvector(betan_t,1,numsteps);
  free_dvector(dp_t,1,numsteps);
  free_dvector(dp_t_f,1,numsteps);
  free_dvector(w_t,1,numsteps);
  free_dvector(dwdt_t,1,numsteps);
  free_dvector(wthrs_t,1,numsteps);

  return 0;
}

double dwdt(double t_in, double w) {
 
  return (eta_star/k_0)*(D_prm_term(t_in,w) + D_tot_term(w) + D_pol_term(w));

}

double D_prm_term(double t_in, double w) {
  double dpval,dpw,f1,f2,alpha,beta,deltadp,betan_t_tmp,w_trns,d2betadwdt; 

  /* get the betan at time t_in */
  betan_t_tmp=getbetan(betan_t_f[i-1],t_in-t[i-1],dbetadt,w,betamin,betamax);

  /* get the delta-prime(betan) for that betan */
  splint(betan,dp,dpb2,dpnum,betan_t_tmp,&dpval);

  /* prevent errors for w going below w_0 */
  if (w<w_0) w=w_0;

  /* Now modify the raw delta-prime(betan) with a choice of w functions */

  /* review the paper on this*/
  if (dp_w_pow==1) dpval=dpval*pow(w/2.0,-2.0*alpha_l)*sqrt(-4.0*D_I);

  /* the tanh function is used to transition from the dp(beta) to
     dp(w) */ 
  /* if (dp_w_tan==1) {
    f1=0.0;
    f2=1.0;
    f2=(f1-f2)/2.0;
    alpha=10;
    beta=1;
    w_trns=0.1;
    deltadp=f1+f2*(tanh(alpha*(w_trns-w))-beta);

    dpval=dpval+deltadp*(dpw-dpval); 
  } */

  /* standard linear decrease with w */
  if (dp_w_w==1) dpval=dpval - alphadp*w;

  /* store the final value of dp and betan*/
  dp_t_f[i]=dpval;
  betan_t_f[i]=betan_t_tmp;

  return dpval;

}

double D_tot_term(double w) {

  switch (dtot_mod) {
    case 1: return D_tot*w/(w*w + w_d*w_d);
    case 2: return D_tot*w/(w*w + w_d*w_d) +
                   D_R/(alpha_s-H)/sqrt(w*w + w_d*w_d);
  }
}

double D_pol_term(double w) {
  switch (dpol_mod) {
    case 1: return D_pol/pow(w,3.0);
    case 2: return w*D_pol/(pow(w_b,4.0)+pow(w,4.0));
  }
}

double dwdt_func(double w) {
  return dwdt(t[i],w);
}

double getbetan(double betanpre, double dt, double dbetadt, double w,
                double betamin, double betamax) {
  double betanval,d2betadwdt;

  d2betadwdt=2.0/w_sat;  /* ~2/0.25=dbetadt(fastish)/w_sat(typical) */
  /* d2betadwdt=dbetadt/w_sat; */

  betanval=betanpre+(dbetadt-d2betadwdt*w)*dt;

  if (betanval>betamax) betanval=betamax;
  if (betanval<betamin) betanval=betamin;

  return betanval;
}

void makepole(double nrm, double offset, double maxy, 
              double betamin, double betamax, 
              int pdn, double *dp, double *bn) {
  int i;
  double *d2,yp1,ypn,xp;
  FILE *dpf;

  d2=dvector(1,pdn);

  /* set up the pole to run between 0 and 1, offset from them by one step
     because of problems caused at 0 and 1 */
  for (i=1;i<=pdn;i++) {
    bn[i]=1.*i/(pdn+1); /* the 1. is required to convert to double precision 
                           before the division operation */
    dp[i]=-M_PI*bn[i]/tan(M_PI*bn[i]);
  }

  /* determine the point in x where the maximum pole value desired occurs. 
     this is done with a spline evaluation, which must be set up with the
     first derivatives at the end points. */

  /* for simplicity here I take the dy/dx and invert it since the points
     in x are equally spaced. (we're interpolating x on y and need dx/dy) */
  yp1= 6/((pdn+1)*(2*dp[4]-9*dp[3]+18*dp[2]-11*dp[1]));
  ypn=-6/((pdn+1)*(2*dp[pdn-3]-9*dp[pdn-2]+
                           18*dp[pdn-1]-11*dp[pdn]));

  spline(dp,bn,pdn,yp1,ypn,d2);
  splint(dp,bn,d2,pdn,maxy,&xp);

  dpf=fopen("pole_m.dat","w");

  /* set up the pole to run between an adjusted offset from zero and xp */
  for (i=1;i<=pdn;i++) {
    bn[i]=xp*i/pdn;
    dp[i]=-M_PI*bn[i]/tan(M_PI*bn[i]);

    /* now set the x dimension to be the beta scale */
    bn[i]=betamin+(betamax-betamin)*(i-1)/(pdn-1);

    /* the pole naturally starts at -1, re-set it to start at 0 */
    dp[i]=dp[i]+1.0;

    /* the maxy value was chosen for the shape, not the value, 
       normalize it out */
    dp[i]=dp[i]/(maxy+1.);

    /* now we set the pole height to nrm, with offset value inserted for
       the next step */
    dp[i]=(nrm-offset)*dp[i];

    /* finally, shift the pole by the offset, and we arrive at a pole shaped
       by maxy, starting at offset, and with a maximum of nrm */
    dp[i]=dp[i]+offset;

    fprintf(dpf,"%lf %lf\n",bn[i],dp[i]);
  }

  fclose(dpf);

  free_dvector(d2,1,pdn);

}

void read_dp(double *bmin, double *bmax, unsigned int *pdn,
             double *dp_in, double *bn_in) {
  FILE *fp;
  char oneline[1024];

  /* open and read the pole and equilibrium data */

  fp=fopen("pole.dat","r");
  fgets(oneline,1024,fp);
  sscanf(oneline,"%i\n",pdn);

  dp_in=dvector(1,*pdn);
  bn_in=dvector(1,*pdn);

  i=1;
  while (fgets(oneline,1024,fp)!=NULL) {
    sscanf(oneline,"%lf %lf\n",&bn_in[i],&dp_in[i]);
    dp_in[i]=nrm*dp_in[i]; 
    i++;
  }
  fclose(fp);

  *bmin=minvec(bn_in,*pdn);
  *bmax=maxvec(bn_in,*pdn);

}

void store_dwdt_vs_w() {

  /* find the neoclassical threshold island by bracketing dwdt=0 and finding
     the zero crossing.  use dpnum as an index of the accuracy with which to 
     bracket the root in dwdt and record dwdt vs. w */
  double br;

  /* record the current time's dwdt vs. w */
  for (j=0;j<=wnum-1;j++) 
    fprintf(fpw,"%lf %lf\n",(1.*j)/(wnum-1.0),dwdt(t[i],(1.*j)/(wnum-1.0)));

  /* reset the bracket point and walk up in w to find the bracket point */
  br=1.0/(wnum-1.0);
  for (j=1;j<=wnum-1;j++) {
    if(dwdt(t[i],(1.*j)/(wnum-1.0))>0.0) {
      br=j/(wnum-1.0);
      j=wnum;
    }
  }
  /*if (br>1.0/(wnum-1.0)) {
    printf("pre\n");
    wthrs_t[i]=zbrent(dwdt_func,1.0/wnum,br,1e-12);
    printf("post\n");
  } else wthrs_t[i]=0.0;*/
}
