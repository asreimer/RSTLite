
/* bffit.c
   ==========
*/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
#include "rtypes.h"
#include "rmath.h" 
#include "nrfit.h" 
#include "rprm.h"
#include "rawdata.h" 
#include "fitdata.h"
#include "fitblk.h"
#include "bffit.h"
#include "badsmp.h"
#include "badlags.h"
#include "mpfit.h"
#include "selfclutter.h"
#include "error_estimates.h"

/* dbl_cmp
   =========
   Author: R.J.Barnes & K. Baker
*/

struct one_param_fit {
  float param;
  float min_chi;
  float left_p;
  float right_p;
};

struct two_param_fit {
  float param1;
  float param2;
  float min_chi;
  float left_p1;
  float right_p1;
  float left_p2;
  float right_p2;
};

int lm_dbl_cmp(const void *x,const void *y) 
{
  double *a,*b;
  a=(double *) x;
  b=(double *) y;
  if (*a > *b) return 1;
  else if (*a == *b) return 0;
  else return -1;
}


void tx_flag_lags(int range,int *badlag,struct FitACFBadSample *bptr,
	       struct FitPrm *ptr) {
/* Does the same thing as FitACFCkRng in badlags.c EXCEPT
   that is doesn't use the Ponomarenko and Waters 2006 CRI
   badlag detection algorithm */

  int sam1, sam2, i, j;
  for (i=0; i<ptr->mplgs; i++) {
	badlag[i] = 0;
	sam1 = ptr->lag[0][i]*(ptr->mpinc/ptr->smsep)
			+ range - 1;
	sam2 = ptr->lag[1][i]*(ptr->mpinc/ptr->smsep)
			+ range - 1;

	for (j=0; j<bptr->nbad; j++) {
      if ((sam1 == bptr->badsmp[j]) || (sam2 == bptr->badsmp[j])) 
        badlag[i] = 1;
	  if (sam2 < bptr->badsmp[j]) break;
    }
  }

  /* This section of code is only of use to fitacf for reprocessing old
	data that used the 16 lag, 7 pulse code.  */
	
  if ((ptr->mplgs == 17) && (ptr->old !=0) )
	badlag[13] = 1;

  return;
}

/* Function to determine at what range lag0 measurements are "bad" due to 
more than one pulse "in the air". */
int bad_lag0(struct RadarParm *prm,int mplgs,int16 *lagtable[2]) {
  /* taken from badlag.c originally ACFBadLagZero */
  int i;
  int badrng=-1;
  int sampleunit;
  int numrng;
  int numfrng;
  int mpincdiff;
  int pulserisetime;
  int halftxpl;

  /* search which pulse in the pulse-pattern that is used to compute lag 0
      power  */
 
  for(i=0; i < prm->mppul;i++)
    if(lagtable[0][0] == prm->pulse[i]) break;
    if(i >= (prm->mppul - 1)) badrng = -1;  
    else {   
      sampleunit = prm->mpinc / prm->smsep;
      mpincdiff = prm->pulse[i+1] - 
                   prm->pulse[i];
      numrng = mpincdiff * sampleunit;
   
      /* compute the pulse rise time effect */

      halftxpl = prm->txpl / 2;
      pulserisetime = halftxpl /  prm->smsep;
      if ((halftxpl % prm->smsep) != 0)
         pulserisetime++;   /* add one because of the int div */
      numrng = numrng - pulserisetime;  /* subtract pulse rise time */
      numfrng = prm->lagfr / prm->smsep;

      /* now compute the start of the bad range */

      badrng = numrng - numfrng;
      if (badrng < 0) badrng = 0;
   }      
   return badrng;
}


void setup_fblk(struct RadarParm *prm, struct RawData *raw,struct FitBlock *input)
{
  int i,j,n;
  void *tmp=NULL;

  if (prm->time.yr < 1993) input->prm.old=1;

  input->prm.xcf=prm->xcf;
  input->prm.tfreq=prm->tfreq;
  input->prm.noise=prm->noise.search;
  input->prm.nrang=prm->nrang;
  input->prm.smsep=prm->smsep;
  input->prm.nave=prm->nave;
  input->prm.mplgs=prm->mplgs;
  input->prm.mpinc=prm->mpinc;
  input->prm.txpl=prm->txpl;
  input->prm.lagfr=prm->lagfr;
  input->prm.mppul=prm->mppul;
  input->prm.bmnum=prm->bmnum;
  input->prm.cp=prm->cp;
  input->prm.channel=prm->channel;
  input->prm.offset=prm->offset; /* stereo offset */


  /* need to incorporate Sessai's code for setting the offset
     for legacy data here.
  */


  if (input->prm.pulse==NULL) tmp=malloc(sizeof(int)*input->prm.mppul);
  else tmp=realloc(input->prm.pulse,sizeof(int)*input->prm.mppul);
  if (tmp==NULL) return;
  input->prm.pulse=tmp;
  for (i=0;i<input->prm.mppul;i++) input->prm.pulse[i]=prm->pulse[i];

  for (n=0;n<2;n++) {
    if (input->prm.lag[n]==NULL) tmp=malloc(sizeof(int)*(input->prm.mplgs+1));
    else tmp=realloc(input->prm.lag[n],sizeof(int)*(input->prm.mplgs+1));
    if (tmp==NULL) return;
    input->prm.lag[n]=tmp;
    for (i=0;i<=input->prm.mplgs;i++) input->prm.lag[n][i]=prm->lag[n][i];
  }



  if (input->prm.pwr0==NULL) tmp=malloc(sizeof(int)*input->prm.nrang);
  else tmp=realloc(input->prm.pwr0,sizeof(int)*input->prm.nrang);
  if (tmp==NULL) return;
  input->prm.pwr0=tmp;

  if (input->acfd==NULL) tmp=malloc(sizeof(struct complex)*input->prm.nrang*
                                    input->prm.mplgs);
  else tmp=realloc(input->acfd,sizeof(struct complex)*input->prm.nrang*
                                   input->prm.mplgs);
  if (tmp==NULL) return;
  input->acfd=tmp;

  if (input->xcfd==NULL) tmp=malloc(sizeof(struct complex)*input->prm.nrang*
                                    input->prm.mplgs);
  else tmp=realloc(input->xcfd,sizeof(struct complex)*input->prm.nrang*
                                   input->prm.mplgs);
  if (tmp==NULL) return;
  input->xcfd=tmp;

  memset(input->acfd,0,sizeof(struct complex)*input->prm.nrang*
                                   input->prm.mplgs);
  memset(input->xcfd,0,sizeof(struct complex)*input->prm.nrang*
                                   input->prm.mplgs);



  for (i=0;i<input->prm.nrang;i++) {
    input->prm.pwr0[i]=raw->pwr0[i];

    if (raw->acfd[0] !=NULL) {
      for (j=0;j<input->prm.mplgs;j++) {
        input->acfd[i*input->prm.mplgs+j].x=raw->acfd[0][i*input->prm.mplgs+j];
        input->acfd[i*input->prm.mplgs+j].y=raw->acfd[1][i*input->prm.mplgs+j];
      }
    }
    if (raw->xcfd[0] !=NULL) {
      for (j=0;j<input->prm.mplgs;j++) {
        input->xcfd[i*input->prm.mplgs+j].x=raw->xcfd[0][i*input->prm.mplgs+j];
        input->xcfd[i*input->prm.mplgs+j].y=raw->xcfd[1][i*input->prm.mplgs+j];
      }
    }
  }

  return;
}


double calc_phi0(float *x,float *y, float m, int n)
{
	/*this is old*/
  double sum_x = 0, sum_y = 0, b;
  int i;

  for(i=0;i<n;i++)
  {
    sum_x += x[i];
    sum_y += y[i];
  }/*

  b = (sum_y/((float)n))-((m*sum_x)/((float)n));*/


  b = ((m*sum_x)/((float)n));

	/*this is new*/
	double * phases = malloc(n*sizeof(double));
  for(i=0;i<n;i++)
    phases[i] = y[i];

	qsort(phases, n, sizeof(double), lm_dbl_cmp);

	double median_phase = phases[(int)(n/2.)];
	b = median_phase - b;
	free(phases);
  return b;
}


/*function to calculate noise level*/
void lm_noise_stat(struct RadarParm *prm, struct FitBlock * fblk,
                double * skynoise)
{
  int j=0, R, i=0;
  double * pwrd = malloc(prm->nrang*sizeof(double));
  if(prm->cp != -3310 && prm->cp != 3310 && prm->cp != 503 && prm->cp != -503)
  {
    j=0;
    for(R=0;R<prm->nrang;R++)
    {
      pwrd[j] = fblk->prm.pwr0[R];
      if(pwrd[j] > 0.)
        j++;
    }
    qsort(pwrd, j, sizeof(double), lm_dbl_cmp);
    if(j >= 10)
    {
      for(i=0;i<10;i++)
        *skynoise += pwrd[i];
      *skynoise /= 10.;
    }
    else
    {
      for(i=0;i<j;i++)
        *skynoise += pwrd[i];
      *skynoise /= (double)j;
    }
    if(*skynoise <= 1.) *skynoise = prm->noise.search;
  }
  else
    *skynoise = prm->noise.search;

  return;
}


/* Function that calculates the "best fit" velocity for the decaying complex
sinusoidal model of the SuperDARN ACF. This function requires previous 
knowledge of the spectral width as it only fits for the velocity. */
struct one_param_fit get_v_brute(struct RadarParm *prm, double *good_lags, int goodcnt, 
                                 int *lag_inds, double *repwr, double *impwr, float *error, double width)
{
  
  int i,j,k,lag,L;
  int min_ind;
  double F1,F2;
  const double delta_chi = 9.0; /* delta chi to calculate error bars 3sigma gives delta_chi = 9.0 */
  const int num_f = 1001;
  const double nyquist_f = 1.0/(2.0 * prm->mpinc * 1.e-6);
  const double f_step = (nyquist_f - (-nyquist_f)) / ((double)(num_f) - 1);  
  const double tau = prm->mpinc*1.e-6;
  const double pi = 3.1415926;
  const double lamda = 299792458.0/(prm->tfreq*1000.0);

  /* temporary variables for output */
  struct one_param_fit vfit;

  /* frequency and chi-squared arrays */
  double *fs = malloc(num_f*sizeof(double));
  double *chi = malloc(num_f*sizeof(double));
  double *times = malloc(goodcnt*sizeof(double));
  double *D1 = malloc(goodcnt*sizeof(double));
  double *D2 = malloc(goodcnt*sizeof(double));
  double *exp_factors = malloc(goodcnt*sizeof(double));
  double *errors = malloc(goodcnt*sizeof(double));

  for(j=0;j<goodcnt;j++)
  {
    lag = (int)(good_lags[j]);
    times[j] = tau*lag;
    exp_factors[j] = repwr[0]*exp(-times[j]*2.*pi*width/lamda);
    D1[j] = repwr[lag];
    D2[j] = impwr[lag];
    errors[j] = (double)(error[lag_inds[j]]);
  }

  /* initialize chi array and generate an array of spectral widths */
  for (i=0;i<num_f;i++)
  {
    fs[i] = -nyquist_f + i*f_step;
    chi[i] = 0.0;
  }

  /* calculate chi2 at each velocity */
  min_ind = 0;
  vfit.min_chi = 10000000000000000.;
  vfit.param = fs[0];

  for (i=0;i<num_f;i++)
  {
      for(j=0;j<goodcnt;j++)
      {
        F1 = exp_factors[j]*cos(2.*pi*fs[i]*times[j]);
        F2 = exp_factors[j]*sin(2.*pi*fs[i]*times[j]);
        chi[i] = chi[i] + ((D1[j]-F1)*(D1[j]-F1) + (D2[j]-F2)*(D2[j]-F2))/(errors[j]*errors[j]);
      }

      if (chi[i] < vfit.min_chi)
      {
        vfit.min_chi = chi[i];
        vfit.param = fs[i];
        min_ind = i;
      }
  }

  /* using 3-sigma, get error bars for spectral width (assymetric) */
  for (i=min_ind;i<num_f;i++)
  {
    if (chi[i] <= vfit.min_chi + delta_chi)
    {
      vfit.right_p = fs[i];
    } /* else {
      break;
    } */
  }

  for (i=min_ind;i>-1;i--)
  {
    if (chi[i] <= vfit.min_chi + delta_chi)
    {
      vfit.left_p = fs[i];
    } /*else {
      break;
    }*/
  }

  vfit.param *= lamda/2.;
  vfit.left_p *= lamda/2.;
  vfit.right_p *= lamda/2.;

  free(times);
  free(D1);
  free(D2);
  free(exp_factors);
  free(errors);
  free(fs);
  free(chi);

  return vfit;
}

/* Function that calculates the "best fit" spectral width for the 
decaying complex sinusoidal model of the SuperDARN ACF. Function 
assumes a decaying exponential ACF envelope. */
struct one_param_fit get_w_brute(struct RadarParm *prm, double *good_lags, int goodcnt, 
                                 int *lag_inds, double *lagpwr, float *error)
{

  int i,j,k,lag,L;
  double F;
  int min_ind;
  const double delta_chi = 9.0; /* delta chi to calculate error bars 3sigma gives delta_chi = 9.0 */
  const int num_w = 1001;
  const double max_w = 1000.0;
  const double min_w = -100.0;
  const double ws_step = (max_w - min_w) / ((double)(num_w) - 1);  
  const double tau = prm->mpinc*1.e-6;
  const double pi = 3.1415926;
  const double lamda = 299792458.0/(prm->tfreq*1000.0);

  /* temporary variables for output */
  struct one_param_fit wfit;

  /* spectral width and chi-squared arrays */
  double *ws = malloc(num_w*sizeof(double));
  double *chi = malloc(num_w*sizeof(double));
  double *times = malloc(goodcnt*sizeof(double));
  double *D = malloc(goodcnt*sizeof(double));
  double *errors = malloc(goodcnt*sizeof(double));

  for(j=0;j<goodcnt;j++)
  {
    lag = (int)(good_lags[j]);
    times[j] = tau*lag;
    D[j] = lagpwr[lag];
    errors[j] = (double)(error[lag_inds[j]]);
  }

  /* initialize chi array and generate an array of spectral widths */
  for (i=0;i<num_w;i++)
  {
    ws[i] = min_w + i*ws_step;
    chi[i] = 0.0;
  }
  /* calculate chi2 at each spectral width */
  min_ind = 0;
  wfit.min_chi = 10000000000000000.;
  wfit.param = ws[0];

  for (i=0;i<num_w;i++)
  {
      for(j=0;j<goodcnt;j++)
      {
        F = lagpwr[0]*exp(-times[j]*2.*pi*ws[i]/lamda);
        chi[i] += (D[j]-F)*(D[j]-F)/(errors[j]*errors[j]);
      }

      if (chi[i] < wfit.min_chi)
      {
        wfit.min_chi = chi[i];
        wfit.param = ws[i];
        min_ind = i;
      }
  }
  
  /* using 3-sigma, get error bars for spectral width (assymetric) */
  for (i=min_ind;i<num_w;i++)
  {
    if (chi[i] <= wfit.min_chi + delta_chi)
    {
      wfit.right_p = ws[i];
    } else {
      break;
    } 
  }

  for (i=min_ind;i>-1;i--)
  {
    if (chi[i] <= wfit.min_chi + delta_chi)
    {
      wfit.left_p = ws[i];
    } else {
      break;
    } 
  }

  free(times);
  free(D);
  free(errors);
  free(ws);
  free(chi);

  return wfit;
}

/* Use JP Villain power envelope */
/* Function that calculates the "best fit" spectral width and diffusion coefficient 
for the decaying complex sinusoidal model of the SuperDARN ACF. Function assumes a 
decaying exponential ACF envelope. */
struct two_param_fit get_w_and_d_brute(struct RadarParm *prm, double *good_lags, int goodcnt, 
                                       int *lag_inds, double *lagpwr, float *error)
{

  int i,j,k,lag,L;
  double F;
  int min_ind1,min_ind2;
  const double delta_chi = 9.0; /* delta chi to calculate error bars 3sigma gives delta_chi = 9.0 */
  const int num_w = 1000;
  const double max_w = 1000.0;
  const double min_w = -100.0;
  const double ws_step = (max_w - min_w) / ((double)(num_w) - 1);  

  const int num_tl = 1000;
  const double max_tl = 0.1;
  const double min_tl = 1e-6;
  const double tl_step = (max_tl - min_tl) / ((double)(num_tl) - 1);  

  const double tau = prm->mpinc*1.e-6;
  const double pi = 3.1415926;
  const double lamda = 299792458.0/(prm->tfreq*1000.0);
  const double kmag = 2.*pi/lamda;

  /* temporary variables for output */
  struct two_param_fit wtlfit;

  /* spectral width and chi-squared arrays */
  double *ws = malloc(num_w*sizeof(double));
  double *tl = malloc(num_tl*sizeof(double));
  double **chi = malloc(num_w*sizeof(double *));
  for (i=0;i<num_w;i++)
  {
    chi[i] = malloc(num_tl*sizeof(double));
  }
  double *times = malloc(goodcnt*sizeof(double));
  double *D = malloc(goodcnt*sizeof(double));
  double *errors = malloc(goodcnt*sizeof(double));

  for(k=0;k<goodcnt;k++)
  {
    lag = (int)(good_lags[k]);
    times[k] = tau*lag;
    D[k] = lagpwr[lag];
    errors[k] = (double)(error[lag_inds[k]]);
  }

  /* initialize chi array and generate an array of spectral widths */
  for (i=0;i<num_w;i++)
  {
    ws[i] = min_w + i*ws_step;
  }
  for (j=0;j<num_tl;j++)
  {
    tl[j] = min_tl + j*tl_step;
  }
  for (i=0;i<num_w;i++)
  {
    for (j=0;j<num_tl;j++)
    {
      chi[i][j] = 0.0;
    }
  }

  /* calculate chi2 at each spectral width */
  min_ind1 = 0;
  min_ind2 = 0;
  wtlfit.min_chi = 1000000000.;
  wtlfit.param1 = ws[0];
  wtlfit.param2 = tl[0];

  for (i=0;i<num_w;i++)
  {
    for(j=0;j<num_tl;j++)
    {
      for(k=0;k<goodcnt;k++)
      {
        F = lagpwr[0]*exp(-times[k]*kmag*ws[i] - kmag*ws[i]*tl[j]*(exp(-times[k]/tl[j]) - 1));
        chi[i][j] += (D[k]-F)*(D[k]-F)/(errors[k]*error[k]);
      }
      
      if (chi[i][j] < wtlfit.min_chi)
      {
        wtlfit.min_chi = chi[i][j];
        wtlfit.param1 = ws[i];
        wtlfit.param2 = tl[j];
        min_ind1 = i;
        min_ind2 = j;
      }
    }
  }
  
  /* using 3-sigma, get error bars for spectral width (assymetric) */
  for (i=min_ind1;i<num_w;i++)
  {
    if (chi[i][min_ind2] < wtlfit.min_chi + delta_chi)
    {
      wtlfit.right_p1 = ws[i];
    } else {
      break;
    }
  }

  for (i=min_ind1;i>-1;i--)
  {
    if (chi[i][min_ind2] < wtlfit.min_chi + delta_chi)
    {
      wtlfit.left_p1 = ws[i];
    } else {
      break;
    }
  }

  for (j=min_ind2;j<num_tl;j++)
  {
    if (chi[min_ind1][j] < wtlfit.min_chi + delta_chi)
    {
      wtlfit.right_p2 = tl[j];
    } else {
      break;
    }
  }

  for (j=min_ind2;j>-1;j--)
  {
    if (chi[min_ind1][j] < wtlfit.min_chi + delta_chi)
    {
      wtlfit.left_p2 = tl[j];
    } else {
      break;
    }
  }

  free(times);
  free(D);
  free(errors);
  free(ws);
  for (i=0;i<num_w;i++)
  {
    free(chi[i]);
  }
  free(chi);

  return wtlfit;
}

void bffit(struct RadarParm *prm,struct RawData *raw,
              struct FitData *fit, struct FitBlock *fblk, int print)
{

  /*fprintf(stderr,"Starting bffit function...\n");
  fprintf(stderr,"%d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %lf\n",prm->stid,prm->time.yr,prm->time.mo,
            prm->time.dy,prm->time.hr,prm->time.mt,(int)prm->time.sc,prm->bmnum,
            prm->nave,prm->cp,prm->lagfr,prm->smsep,fblk->prm.vdir);*/
  int k;
  double minpwr  = 3.0;
  double skynoise = 0.;
  int   minlag  = 6;
  int pwr_flg,fit_flg;
  int *lag_inds=NULL;
  double *lagpwr=NULL,*good_lags=NULL;
  double *repwr=NULL, *impwr=NULL;
  double lag0pwr,re,im;
  double dog;
  int   *lag_avail=NULL,availcnt=0,goodcnt=0;
  int   lastlag,lag,i,j,R,L,mplgs,tauflg = 0;
  int status;

  struct one_param_fit vfit;
  struct one_param_fit wfit;

  int *badlag = malloc(prm->mplgs * sizeof(int));
  struct FitACFBadSample badsmp;

  /* definitions for calculating Cmpse self-clutter */
  int badrng;
  float *selfclutter = malloc(sizeof(float)*prm->mplgs);
  float *pwrd = malloc(sizeof(float)*prm->nrang);

  /* definitions for error estimates */
  float *error = malloc(sizeof(float)*prm->mplgs);
  float *lag0error = malloc(sizeof(float)*prm->nrang);

  /*check for tauscan*/
  if(prm->cp == -3310 || prm->cp == 3310 || prm->cp == 503 || prm->cp == -503)
    tauflg = 1;

  /* Find the highest lag, and allocate memory */
  if(!(tauflg && prm->mplgs == 18))
  {
    lastlag = 0;
    for (j=0;j<prm->mplgs;j++)
    {
      if (abs(prm->lag[0][j]-prm->lag[1][j])>lastlag)
        lastlag = abs(prm->lag[0][j]-prm->lag[1][j]);
    }
  } else {
    lastlag=prm->mplgs-1;
  }

  /*define some stuctures using # of lags*/
  lagpwr       = malloc(sizeof(double)*(lastlag+1));
  repwr       = malloc(sizeof(double)*(lastlag+1));
  impwr       = malloc(sizeof(double)*(lastlag+1));
  lag_avail    = malloc(sizeof(int)*(lastlag+1));
  good_lags    = malloc(sizeof(double)*(lastlag+1));
  lag_inds = malloc((lastlag+1)*sizeof(int));

  /*setup fitblock parameter*/
  setup_fblk(prm, raw, fblk);

  FitSetRng(fit,fblk->prm.nrang);
  if(fblk->prm.xcf)
  {
    FitSetXrng(fit,fblk->prm.nrang);
    FitSetElv(fit,fblk->prm.nrang);
  }

  /*calculate noise levels*/
  lm_noise_stat(prm,fblk,&skynoise);
  if(!tauflg)
  {
    /*check for stereo operation*/
    if(fblk->prm.channel==0) 
    {
      FitACFBadlags(&fblk->prm,&badsmp);
    } else {
      FitACFBadlagsStereo(&fblk->prm,&badsmp);
    }
  }

  if(prm->cp == 153) 
  {
    mplgs = prm->mplgs - 1;
  } else {
    mplgs = prm->mplgs;
  }

  prm->noise.mean = skynoise;

  if(print)
  {
    fprintf(stderr,"%d  %d  %lf  %d  %lf\n",prm->nrang,mplgs,skynoise,prm->tfreq,prm->mpinc*1.e-6);
    fprintf(stderr,"%d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %lf\n",prm->stid,prm->time.yr,prm->time.mo,
            prm->time.dy,prm->time.hr,prm->time.mt,(int)prm->time.sc,prm->bmnum,
            prm->nave,prm->cp,prm->lagfr,prm->smsep,fblk->prm.vdir);
  }

  /* determine the "bad range" range where lag0 pwr estimates become bad */
  badrng = bad_lag0(prm,prm->mplgs,prm->lag);

  /* need an array of lag0 power without noise for selfclutter calculation */
  for (R=0;R<prm->nrang;R++)
  {
    pwrd[R] = fblk->prm.pwr0[R] - skynoise;
    if (pwrd[R] < skynoise*0.00001)
      pwrd[R] = skynoise*0.00001;
  }

  status = lag0_error(prm->nrang, pwrd, skynoise, prm->nave, prm->nave, lag0error);

  /* Loop every range gate and calculate parameters */
  for (R=0;R<prm->nrang;R++)
  {
    /*initialize parameters*/
    fit->rng[R].v        = 0.;
    fit->rng[R].v_err    = HUGE_VAL;
    fit->rng[R].p_0      = 0.0;
    fit->rng[R].w_l      = 0.0;
    fit->rng[R].w_l_err  = 0.0;
    fit->rng[R].p_l      = 0.0;
    fit->rng[R].p_l_err  = 0.0;
    fit->rng[R].w_s      = 0.0;
    fit->rng[R].w_s_err  = 0.0;
    fit->rng[R].p_s      = 0.0;
    fit->rng[R].p_s_err  = 0.0;
    fit->rng[R].sdev_l   = 0.0;
    fit->rng[R].sdev_s   = 0.0;
    fit->rng[R].sdev_phi = 0.0;
    fit->rng[R].qflg     = 0;
    fit->rng[R].nump     = 0;
    fit->rng[R].gsct     = 0;
    availcnt = 0;

    /* calculate the self_clutter in each lag */
    status = Cmpse(&fblk->prm, prm->lag, R, selfclutter, badrng, pwrd);

    /*calculate SNR of lag0power*/
    lag0pwr  = 10.0*log10(pwrd[R]/skynoise);

    /*not tauscan, check for badlags*/
    if(!tauflg) tx_flag_lags(R+1,badlag,&badsmp,&fblk->prm);

    /*Preliminaries, badlag checking, power level checking*/

    for(L=0;L<mplgs;L++)
    {
      if(tauflg && prm->mplgs == 18) /*tauscan , new ROS*/
      { 
        lag = L;
      } else {                       /*old ROS*/
        lag = abs(prm->lag[0][L] - prm->lag[1][L]);
      }

      re  = fblk->acfd[R*fblk->prm.mplgs+L].x;
      im  = fblk->acfd[R*fblk->prm.mplgs+L].y;
      lagpwr[lag] = sqrt(re*re + im*im);
      repwr[lag] = re;
      impwr[lag] = im;

      if((tauflg || badlag[L] == 0))/* && lagpwr[lag]>pwrd[R]/sqrt(1.0*prm->nave))*/
      {
        lag_avail[availcnt] = lag;
        availcnt++;
      } else {
        lagpwr[lag] = 0.0;
        repwr[lag] = 0.0;
        impwr[lag] = 0.0;
      }
    }

    /* set pwr_flg and minlag based on tauscan or not */
    /* basically, pwr_flag == 1 if SNR >= 1 at lag0 */
    if(tauflg) /* for tauscan */
    {
      pwr_flg = (lag0pwr >= minpwr);
    } else {   /* for everything else */
      pwr_flg = (pwrd[R] >= skynoise);
    }

    pwr_flg = 1;

    /* if SNR is high enough and we have ge minlag "good" lags */
    if((pwr_flg) && (availcnt>=minlag))
    {
      /* build array of "good" lags */
      goodcnt = 0;
      for(i=0;i<availcnt;i++)
      {
        lag = lag_avail[i];
        good_lags[goodcnt] = lag;
        goodcnt++;
      }

      if((tauflg) && prm->mplgs == 18)
      {
        /*tauscan AND new ROS*/
        for(j=0;j<goodcnt;j++)
        {
          lag_inds[j] = (int)(good_lags[j]);
        }
      } else {   
        /*non-tauscan OR old ROS*/
        for(j=0;j<goodcnt;j++)
        {
          lag = (int)(good_lags[j]);       
          for(k=0;k<prm->mplgs;k++)
          {
            if(abs(prm->lag[0][k]-prm->lag[1][k])==lag)
            {
              lag_inds[j] = k;
              break;
            }
          }
        }
      }

      acf_error(prm->mplgs, pwrd[R], skynoise, selfclutter, prm->nave, error);
      error[0] = lag0error[R];

      /**********************/
      /*single component fit*/
      /**********************/

      /* First fit for the spectral width, assuming exponential ACF envelope */
      wfit = get_w_brute(prm,good_lags,goodcnt,lag_inds,lagpwr,error);

      /* Next fit for the velocity, using fitted spectral width, assuming exponential ACF envelope */
      vfit = get_v_brute(prm,good_lags,goodcnt,lag_inds,repwr,impwr,error,wfit.param);

      /* Were the fits good? Use reduced chi**2 to determine yes/no */
      /*degrees of freedom */
      /*dog = number of data points - # of fitted parameters - 1 */
      dog = ((float)(goodcnt)) - 1 - 1;
      fit_flg = 1; /*((wfit.min_chi/dog < 1) && (vfit.min_chi/dog < 1)); */

      /* Now save parameters to datafile */
      fit->rng[R].p_0   = lag0pwr;
      fit->rng[R].p_l   = lag0pwr;
      fit->rng[R].p_l_err = 10.0*log10((lag0error[R]+pwrd[R])/skynoise)-lag0pwr;
      fit->rng[R].v     = vfit.param;
      fit->rng[R].w_l   = wfit.param;
      /* What are the error bars on the fitted parameters? */
      /* width - generally very assymetric. For now store in w_l_err and w_s_err */
      fit->rng[R].w_l_err = fabs(wfit.right_p - wfit.param);
      fit->rng[R].w_s_err = fabs(wfit.param - wfit.left_p);
      /* velocity - generally very close to symetric. Return largest of two. */
      if (fabs(vfit.param - vfit.left_p) >= fabs(vfit.right_p - vfit.param))
      {
          fit->rng[R].v_err = fabs(vfit.param - vfit.left_p);
      } else {
          fit->rng[R].v_err = fabs(vfit.right_p - vfit.param);
      }
      fit->rng[R].nump  = goodcnt;
      fit->noise.skynoise = skynoise;
      fit->rng[R].w_s   = vfit.left_p;
      fit->rng[R].p_s   = vfit.right_p;
      fit->rng[R].qflg  = fit_flg;
      fit->rng[R].gsct = (fabs(vfit.param)-(30-1./3.*fabs(wfit.param)) < 0);
    }
  }

  free(lag_inds);
  free(lagpwr);
  free(repwr);
  free(impwr);
  free(lag_avail);
  free(good_lags);
  free(badlag);
  free(selfclutter);
  free(pwrd);
  free(error);
  free(lag0error);

  return;
}
