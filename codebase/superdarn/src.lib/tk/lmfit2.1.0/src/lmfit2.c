
/* lmfit.c
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
#include "lmfit2.h"
#include "badsmp.h"
#include "badlags.h"
#include "mpfit.h"
#include "selfclutter.h"
#include "error_estimates.h"

/* dbl_cmp
   =========
   Author: R.J.Barnes & K. Baker
*/
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
void lm_noise_stat(struct RadarParm *prm, struct RawData * raw,
                double * skynoise)
{
  int j=0, R, i=0;
  double * pwrd = malloc(prm->nrang*sizeof(double));
  if(prm->cp != -3310 && prm->cp != 3310 && prm->cp != 503 && prm->cp != -503)
  {
    j=0;
    for(R=0;R<prm->nrang;R++)
    {
      pwrd[j] = raw->pwr0[R];
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
double get_v_brute(struct RadarParm *prm, float *good_lags, int goodcnt, 
                   float *repwr, float *impwr, float *error, float width, 
                   float *min_chi, float *min_chi_v, float *left_chi_v, 
                   float *right_chi_v)
{
  
  int i,j,k,lag,L;
  float time;
  float delta_chi = 9.0; /* delta chi to calculate error bars 3sigma gives delta_chi = 9.0 */
  float num_f = 1000.0;
  float nyquist_f = 1.0/(2.0 * prm->mpinc * 1.e-6);
  float f_step = (nyquist_f - (-nyquist_f)) / (num_f - 1);
  float F1,F2,D1,D2;  
  float tau = prm->mpinc*1.e-6;
  float pi = 3.1415926;
  float lamda = 299792458.0/(prm->tfreq*1000.0);

  /* temporary variables for output */
  float chi_min = -1;
  float chi_min_f = -1;
  float left_f = -1;
  float right_f = -1;

  /* frequency and chi-squared arrays */
  float *fs = malloc(num_f*sizeof(float));
  float *chi = malloc(num_f*sizeof(float));

  /* initialize chi array and generate an array of spectral widths */
  fs[0] = -nyquist_f;
  chi[0] = 0.0;
  for (i=1;i<num_f;i++)
  {
      fs[i] = fs[i-1] + f_step;
      chi[i] = 0.0;
  }

  /* calculate chi2 at each spectral width */
  for (i=0;i<num_f;i++)
  {
      for(j=0;j<goodcnt;j++)
      {
        lag = (int)(good_lags[j]);
        /*tauscan AND new ROS*/
        if((prm->cp == -3310 || prm->cp == 3310 || prm->cp == 503 || prm->cp == -503) && prm->mplgs == 18)
        {
          L = lag;
        } else {  /*non-tauscan OR old ROS*/
          for(k=0;k<prm->mplgs;k++)
          {
            if(abs(prm->lag[0][k]-prm->lag[1][k])==lag)
            {
              L = k;
            }
          }
        }
        time = tau*lag;
        D1 = repwr[lag];
        D2 = impwr[lag];
        F1 = repwr[0]*exp(-time*2.*pi*width/lamda)*cos(2.*pi*fs[i]*time);
        F2 = repwr[0]*exp(-time*2.*pi*width/lamda)*sin(2.*pi*fs[i]*time);

        chi[i] += ((D1-F1)*(D1-F1) + (D2-F2)*(D2-F2))/(error[L]*error[L]);
      }
  }

  /* determine the minimum chi2 and the spectral width at it */
  chi_min = chi[0];
  chi_min_f = fs[0];
  for (i=0;i<num_f;i++)
  {
    if (chi[i] < chi_min)
    {
      chi_min = chi[i];
      chi_min_f = fs[i];
    }
  }

  /* using 3-sigma, get error bars for spectral width (assymetric) */
  for (i=0;i<num_f;i++)
  {
    if ((chi[i] > chi_min + delta_chi) && (chi_min_f > fs[i]))
    {
      left_f = fs[i];
    }
    if ((chi[i] < chi_min + delta_chi) && (chi_min_f < fs[i]))
    {
      right_f = fs[i];
    }
  }

  *min_chi = chi_min;
  *min_chi_v = chi_min_f * lamda/2.;
  *left_chi_v = left_f * lamda/2.;
  *right_chi_v = right_f * lamda/2.;

  free(fs);
  free(chi);

  return 0;
}


/* Function that calculates the "best fit" spectral width for the 
decaying complex sinusoidal model of the SuperDARN ACF. Function 
assumes a decaying exponential ACF envelope. */
double get_w_brute(struct RadarParm *prm,
                      float *good_lags, int goodcnt, float *lagpwr,
                      float *error, float *min_chi, float *min_chi_w, float *left_chi_w, float *right_chi_w)
{

  int i,j,k,lag,L;
  float time;
  float delta_chi = 9.0; /* delta chi to calculate error bars 3sigma gives delta_chi = 9.0 */
  float num_w = 1000.0;
  float max_w = 1000.0;
  float min_w = 0.0;
  float ws_step = (max_w - min_w) / (num_w - 1);
  float F,D;  
  float tau = prm->mpinc*1.e-6;
  float pi = 3.1415926;
  float lamda = 299792458.0/(prm->tfreq*1000.0);

  /* temporary variables for output */
  float chi_min = -1;
  float chi_min_w = -1;
  float left_w = -1;
  float right_w = -1;

  /* spectral width and chi-squared arrays */
  float *ws = malloc(num_w*sizeof(float));
  float *chi = malloc(num_w*sizeof(float));


  /* initialize chi array and generate an array of spectral widths */
  ws[0] = min_w;
  chi[0] = 0.0;
  for (i=1;i<num_w;i++)
  {
      ws[i] = ws[i-1] + ws_step;
      chi[i] = 0.0;
  }

  /* calculate chi2 at each spectral width */
  for (i=0;i<num_w;i++)
  {
      for(j=0;j<goodcnt;j++)
      {
        lag = (int)(good_lags[j]);
        /*tauscan AND new ROS*/
        if((prm->cp == -3310 || prm->cp == 3310 || prm->cp == 503 || prm->cp == -503) && prm->mplgs == 18)
        {
          L = lag;
        } else {  /*non-tauscan OR old ROS*/
          for(k=0;k<prm->mplgs;k++)
          {
            if(abs(prm->lag[0][k]-prm->lag[1][k])==lag)
            {
              L = k;
            }
          }
        }
        time = tau*lag;
        D = lagpwr[lag];
        F = lagpwr[0]*exp(-time*2.*pi*ws[i]/lamda);
        chi[i] += (D-F)*(D-F)/(error[L]*error[L]);
      }
  }
  
  /* determine the minimum chi2 and the spectral width at it */
  chi_min = chi[0];
  chi_min_w = ws[0];
  for (i=0;i<num_w;i++)
  {
    if (chi[i] < chi_min)
    {
      chi_min = chi[i];
      chi_min_w = ws[i];
    }
  }

  /* using 3-sigma, get error bars for spectral width (assymetric) */
  for (i=0;i<num_w;i++)
  {
    if ((chi[i] > chi_min + delta_chi) && (chi_min_w > ws[i]))
    {
      left_w = ws[i];
    }
    if ((chi[i] < chi_min + delta_chi) && (chi_min_w < ws[i]))
    {
      right_w = ws[i];
    }
  }

  *min_chi = chi_min;
  *min_chi_w = chi_min_w;
  *left_chi_w = left_w;
  *right_chi_w = right_w;

  free(ws);
  free(chi);

  return 0;
}

/* Use JP Villain power envelope and extract spectral width and diffusion coefficient 

NOT READY YET. TO BE IMPLEMENTED.

double get_w_and_d_brute(struct RadarParm *prm,struct RawData *raw,
                      struct FitData *fit, struct FitBlock *fblk,
		      int rang, double skynoise)
{
  int i,j,lastlag;
  float num_w = 1000.0;

  return 0;
}
*/

void lmfit2(struct RadarParm *prm,struct RawData *raw,
              struct FitData *fit, struct FitBlock *fblk, int print)
{
  float minpwr  = 3.0;
  double skynoise = 0.;
  int   minlag  = 6;
  int pwr_flg,fit_flg;
  float *lagpwr=NULL,*good_lags=NULL;
  float *repwr=NULL, *impwr=NULL;
  float lag0pwr,re,im;
  float fitted_width=0.0,fitted_velocity=0.0;
  float dog;
  int   *lag_avail=NULL,availcnt=0,goodcnt=0;
  int   lastlag,lag,i,j,R,L,mplgs,tauflg = 0;
  int status;

  float min_chi_w = -1;
  float left_chi_w = -1;
  float right_chi_w = -1;
  float min_chi_v = -1;
  float left_chi_v = -1;
  float right_chi_v = -1;

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
  lagpwr       = malloc(sizeof(float)*(lastlag+1));
  repwr       = malloc(sizeof(float)*(lastlag+1));
  impwr       = malloc(sizeof(float)*(lastlag+1));
  lag_avail    = malloc(sizeof(int)*(lastlag+1));
  good_lags    = malloc(sizeof(float)*(lastlag+1));

  /*setup fitblock parameter*/
  setup_fblk(prm, raw, fblk);

  FitSetRng(fit,fblk->prm.nrang);
  if(fblk->prm.xcf)
  {
    FitSetXrng(fit,fblk->prm.nrang);
    FitSetElv(fit,fblk->prm.nrang);
  }

  /*calculate noise levels*/
  lm_noise_stat(prm,raw,&skynoise);
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
    pwrd[R] = raw->acfd[0][R*prm->mplgs] - skynoise;
    if (pwrd[R] < 0)
      pwrd[R] = 0;
  }

  status = lag0_error(prm->nrang, pwrd, skynoise, prm->nave, prm->nave, lag0error);

  /* Loop every range gate and calculate parameters */
  for (R=0;R<prm->nrang;R++)
  {
    /* subtract noise level from lag 0 if tauscan */
    if(!(tauflg && prm->mplgs == 18))
    {
      raw->acfd[0][R*prm->mplgs] -= skynoise;
      if (raw->acfd[0][R*prm->mplgs] < 0)
        raw->acfd[0][R*prm->mplgs] = 0;
    }

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
    lag0pwr  = 10.0*log10((raw->acfd[0][R*prm->mplgs])/skynoise);

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

      re  = raw->acfd[0][R*prm->mplgs+L];
      im  = raw->acfd[1][R*prm->mplgs+L];
      lagpwr[lag] = sqrt(re*re + im*im);
      repwr[lag] = re;
      impwr[lag] = im;

      if((tauflg || badlag[L] == 0) && lagpwr[lag]>raw->acfd[0][R*prm->mplgs]/sqrt(1.0*prm->nave))
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
      pwr_flg = (lag0pwr>=minpwr);
    } else {   /* for everything else */
      pwr_flg = raw->acfd[0][R*prm->mplgs] >= skynoise;
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

      acf_error(prm->mplgs, pwrd[R], skynoise, selfclutter, prm->nave, error);
      error[0] = lag0error[R];

      /**********************/
      /*single component fit*/
      /**********************/

      /* First fit for the spectral width, assuming exponential ACF envelope */
      status = get_w_brute(prm,good_lags,goodcnt,lagpwr,error,&min_chi_w,&fitted_width,&left_chi_w,&right_chi_w);

      /* Next fit for the velocity, using fitted spectral width, assuming exponential ACF envelope */
      status = get_v_brute(prm,good_lags,goodcnt,repwr,impwr,error,fitted_width,&min_chi_v,&fitted_velocity,&left_chi_v,&right_chi_v);

      /* Were the fits good? Use reduced chi**2 to determine yes/no */
      /*degrees of freedom */
      /*dog = number of data points - # of fitted parameters - 1 */
      dog = ((float)(goodcnt)) - 1 - 1;
      fit_flg = ((min_chi_w/dog < 1) && (min_chi_v/dog < 1));      

      /* Now save parameters to datafile */
      fit->rng[R].p_l   = lag0pwr;
      fit->rng[R].p_l_err = lag0pwr/log(10.0)*lag0error[R]/skynoise; /* propagated from error in amplitude */
      fit->rng[R].v     = fitted_velocity;
      fit->rng[R].w_l   = fitted_width;
      /* What are the error bars on the fitted parameters? */
      /* width - generally very assymetric. For now store in w_l_err and w_s_err */
      fit->rng[R].w_l_err = fabs(right_chi_w - fitted_width);
      fit->rng[R].w_s_err = fabs(fitted_width - left_chi_w);
      /* velocity - generally very close to symetric. Return largest of two. */
      if (fabs(fitted_velocity - left_chi_v) > fabs(right_chi_v - fitted_velocity))
      {
          fit->rng[R].v_err = fabs(fitted_velocity - left_chi_v);
      } else {
          fit->rng[R].v_err = fabs(right_chi_v - fitted_velocity);
      }
      fit->rng[R].nump  = goodcnt;
      fit->noise.skynoise = skynoise;
      fit->rng[R].w_s   = 0;
      fit->rng[R].qflg  = fit_flg;
      fit->rng[R].gsct = (fabs(min_chi_v)-(30-1./3.*fabs(fitted_width)) < 0);
    }
  }


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
