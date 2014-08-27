 /*COPYRIGHT:
Copyright (C) 2011 by Virginia Tech

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


 MODIFICATION HISTORY:
 Written by AJ Ribeiro 06/16/2011
 Based on code orginally written by Pasha Ponomarenko
*/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <sys/time.h>
#include <zlib.h>
#include "rtypes.h"
#include "dmap.h"
#include "sim_data.h"
#include "rtypes.h"
#include "rprm.h"
#include "rawdata.h"
#include "rawwrite.h"
#include "radar.h"
#include "iq.h"
#include "iqwrite.h"
#include "fitdata.h"
#include "fitread.h"

struct gates
{
  double vel;
	double wid;
	double pow;
};

/* Build the radar parameter structure */
void makeRadarParm2(struct RadarParm * prm, char * argv[], int argc, int cpid, int nave,
                    int lagfr, double smsep, double noise_lev, double amp0, int n_samples,
                    double dt, int n_pul, int n_lags, int nrang, double rngsep, double freq,
                    int * pulse_t,int yr, int mo, int dy, int hr, int mt, int sc, int bmnum)
{
  int i;
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = gmtime(&rawtime);

  if (prm->origin.time !=NULL) free(prm->origin.time);
  if (prm->origin.command !=NULL) free(prm->origin.command);
  if (prm->pulse !=NULL) free(prm->pulse);
  for (i=0;i<2;i++) if (prm->lag[i] !=NULL) free(prm->lag[i]);

  memset(prm,0,sizeof(struct RadarParm));
  prm->origin.time=NULL;
  prm->origin.command=NULL;
  prm->pulse=NULL;
  prm->lag[0]=NULL;
  prm->lag[1]=NULL;
  prm->combf=NULL;


  prm->revision.major = 1;
  prm->revision.minor = 0;
  prm->origin.code = 0;

  RadarParmSetOriginTime(prm,asctime(timeinfo));
  RadarParmSetOriginCommand(prm,"sim_real");


  prm->cp = (int16)cpid;
  prm->stid = 0;

  prm->time.yr = (int16)yr;
  prm->time.mo = (int16)mo;
  prm->time.dy = (int16)dy;
  prm->time.hr = (int16)hr;
  prm->time.mt = (int16)mt;
  prm->time.sc = (int16)sc;
	prm->bmnum = bmnum;

	
  prm->time.us = 0;

  prm->txpow = 9000;                              /* Transmitted power in kW */
  prm->nave = (int16)nave;                        /* Number of samples per average */
  prm->atten = 0;
  prm->lagfr = (int16)(lagfr*smsep*1e6);          /* Lag in time to first range gate (us) */
  prm->smsep = (int16)(smsep*1.e6);               /* Sample separation in time (us) */
  prm->ercod = 0;
  prm->stat.agc = 8192;
  prm->stat.lopwr = 8192;
  prm->noise.search = noise_lev;
  prm->noise.mean = noise_lev;
  prm->channel = 0;
  prm->bmazm = 0;
  prm->scan = 1;
  prm->rxrise = 100;
  prm->intt.sc = (int16)(smsep*n_samples*nave);
  prm->intt.us = (int)(((smsep*n_samples*nave)-(int)(smsep*n_samples*nave))*1e6);
  prm->txpl = 300;
  prm->mpinc = (int16)(dt*1e6);
  prm->mppul = (int16)n_pul;
  prm->mplgs = (int16)n_lags;
  prm->mplgexs = (int16)n_lags;
  prm->nrang = (int16)nrang;
  prm->frang = (int16)(rngsep*lagfr*1e-3);
  prm->rsep = (int16)(rngsep*1e-3);
  prm->xcf = 0;
  prm->tfreq = (int16)(freq*1e-3);
  prm->offset = 0;
  prm->mxpwr = 1070000000;
  prm->lvmax = 20000;

  int16 temp_pul[n_pul];

  for(i=0;i<n_pul;i++)
    temp_pul[i] = (int16)pulse_t[i];

  RadarParmSetPulse(prm,n_pul,temp_pul);


  if(cpid == 1)
  { /* Schiffler */
    int16 temp_lag[100] = {0,0,26,27,20,22,9,12,22,26,22,27,20,26,20,27,12,20,0,9,
                                12,22,9,20,0,12,9,22,12,26,12,27,9,26,9,27};
    RadarParmSetLag(prm,n_lags,temp_lag);
  }
  else if(cpid == 503)
  { 
    int16 temp_lag[100] = {0,0,15,16,27,29,29,32,23,27,27,32,23,29,16,23,15,23,
                            23,32,16,27,15,27,16,29,15,29,32,47,16,32,15,32};
    RadarParmSetLag(prm,n_lags,temp_lag);
  }
  else
  { /* Katscan */
    int16 temp_lag[100] = {0,0,42,43,22,24,24,27,27,31,22,27,24,31,14,22,22,
                                31,14,24,31,42,31,43,14,27,0,14,27,42,27,43,14,31,
                                24,42,24,43,22,42,22,43,0,22,0,24};
    RadarParmSetLag(prm,n_lags,temp_lag);
  }

  RadarParmSetCombf(prm,"sim_real");

}
/* this is a driver program from the data simulator */

int main(int argc,char *argv[])
{

  /********************************************************
  ** definitions of variables needed for data generation **
  ********************************************************/
  double t_d = .04;                         /*Irregualrity decay time s*/
  double w = -9999.;                        /*spectral width*/
  double t_g = 1e-6;                        /*irregularity growth time*/
  double t_c = 1000.;                       /*precipitation time constant*/
  double v_dop =450.;                       /*Background velocity (m/s)*/
  double velo = 0;
  double c = 3.e8;                          /*Speed of light (m/s)*/
  double lambda;
  double freq;                              /*transmit frequency*/
  double amp0 = 1.;                         /*amplitude scaling factor*/
  int noise_flg = 1;                        /*flag to indicate whether white noise is included*/
  double noise_lev = 0.;                    /*white noise level (ratio)*/
  int lagfr = 4;                            /*lag to first range*/
  int life_dist = 0;                        /*lifetime distribution*/
  int cpid = 150;                             /*control program ID number*/

  double dt;                                /*basic lag time*/
  int cri_flg = 1;                          /*cross-range interference flag*/
  int smp_flg = 0;                          /*output raw samples flag*/
  int decayflg = 0;

  int rt = 0; 				    /* variable to catch return values of fscanf, 
                                               prevents compile warnings which resulted in 
                                               bizzare terminal behaviour in Ubuntu 14.04 LTS */

  /*other variables*/
  int s,lag;
  long i,j,in_tau;
  double taus;
  struct RadarParm * prm;
  struct FitData *fitacf;
  int nave;
  int nrang;

  int n_samples;                            /*Number of datapoints in a single pulse sequence*/
  int n_pul,n_lags,*pulse_t,*tau;
  double smsep;                     /*sample spearation*/
  double rngsep;                    /*range gate spearation*/

 
char helpstr[] =
  "\nmake_sim:  generates simulated single-component lorentzian ACFs\n\n"
  "Calling Sequence:  ./sim_fitacf [-options] > output.txt\n"
  "Options:\n"
  "--help: show this information\n"
  "-v_spread v: set gaussian Doppler velocity spread (standard devation) to v\n"
  "         default is 0\n"
  "-t_g t: set growth time to t (in microseconds)\n"
  "         default is 1 us (negligible)\n"
  "-t_c t: set precipitation time constant (lifetime) to t (in microseconds)\n"
  "         default is 1e6 ms (negligible)\n"
  "-noise n: add in white noise level to produce SNR n (in dB)\n"
  "         default is no noise\n"
  "-nave n: set number of averages in the integration period to n\n"
  "         default is 70/50/20 for oldscan/katscan/tauscan\n"
  "-nocri: remove cross range interference from the ACFs\n"
  "         default is CRI on\n"
  "         WARNING: removing cross-range interference will make\n"
  "         the raw samples unuseable, since each range gate will\n"
  "         have to be integrated seperately\n"
  "-samples: output raw samples (to iqdat file) instead of ACFs\n"
  "         default is output ACFs (to rawacf file)\n"
  "\nNOTE: all option inputs must be integers\n";




  /*process command line options*/
  for (i = 1; i < argc; i++)
  {
    /*command line velocity spread*/
    if (strcmp(argv[i], "-v_spread") == 0)
      velo = (double)atoi(argv[i+1]);
    /*command line decorrelation time*/
    else if (strcmp(argv[i], "-t_d") == 0)
      t_d = 1e-3*atoi(argv[i+1]);
    /*command line growth time*/
    else if (strcmp(argv[i], "-t_g") == 0)
      t_g = 1e-6*atoi(argv[i+1]);
    /*command line precipitation time constant*/
    else if (strcmp(argv[i], "-t_c") == 0)
      t_c = 1e-3*atoi(argv[i+1]);
    /*command line noise*/
    else if (strcmp(argv[i], "-noise") == 0)
    {
      noise_flg = 1;
      noise_lev = 1./pow(10.,((double)atoi(argv[i+1]))/10.);
    }
    /*command line CRI flag*/
    else if (strcmp(argv[i], "-nocri") == 0)
      cri_flg = 0;
    /*command line output samples*/
    else if (strcmp(argv[i], "-samples") == 0)
      smp_flg = 1;
    /*display help*/
    else if (strcmp(argv[i], "--help") == 0)
    {
      fprintf(stderr,"%s\n",helpstr);
      exit(0);
    }
  }

  /* #############
     FILE IO STUFF 
     ############# */
  /*fit file to recreate*/
  char * filename = argv[argc-1];

  /* Open the fit file */
  FILE * fitfp=fopen(filename,"r");
  fprintf(stderr,"%s\n",filename);
  if(fitfp==NULL)
  {
    fprintf(stderr,"File %s not found.\n",filename);
    exit(-1);
  }
 


  /* ###################################
     INITIALIZE STRUCTURES AND VARIABLES
     ################################### */

  /* Initialize the parameter structure*/
  prm = RadarParmMake();

  /* Initialize the fit structure*/
  fitacf=FitMake();

  fprintf(stderr,"Reading first line of fit file.\n");
  /* Read the first line of the fit file */
  s=FitFread(fitfp,prm,fitacf);
  fprintf(stderr,"Status %d.\n",s);
  if (s==-1) {
    fprintf(stderr,"Error reading file.\n");
    exit(-1);
  }


  do
  {
        fprintf(stderr,"Loop.\n");
        /* set some variables for ease of programming/readability */
        nrang = prm->nrang;
        nave = prm->nave;
        n_lags = prm->mplgs;
        freq = prm->tfreq*1.e3;
        lambda = c/freq;
        cpid = prm->cp;
        n_pul = prm->mppul;
	rngsep = prm->rsep*1.e3;
	smsep = prm->smsep*1.e-6;
	dt = prm->mpinc*1.e-6;
        lagfr = prm->lagfr/prm->smsep;


        fprintf(stderr,"Pulse table.\n");
        /* GET THE PULSE TABLE */
        pulse_t = malloc(n_pul*sizeof(int));
	for (i=0;i<prm->mppul;i++) pulse_t[i]=prm->pulse[i];

        fprintf(stderr,"Lags.\n");
        /* FIGURE OUT WHICH LAGS ARE IN THE MULTIPLE PULSE SEQUENCE */
        tau = malloc(n_lags*sizeof(int));
        for (i=0;i<=prm->mplgs;i++) {
            lag = prm->lag[1][i] - prm->lag[0][i];
            in_tau = 0;

            for (j=0;j<=prm->mplgs;j++) {
              if (tau[j] == lag) {
                in_tau=1;
              }
            }

            if (!in_tau) {
              tau[i]=lag;
              fprintf(stderr,"%d \n",lag);
            }
        }

        /*control program dependent variables*/
	taus = dt/smsep;                                      /*lag time in samples*/
	n_samples = (pulse_t[n_pul-1]*taus+nrang+lagfr);      /*number of samples in 1 pulse sequence*/

        fprintf(stderr,"Other arrays.\n");

	
	double * t_d_arr = malloc(nrang*sizeof(double));
	double * t_g_arr = malloc(nrang*sizeof(double));
	double * t_c_arr = malloc(nrang*sizeof(double));
	double * v_dop_arr = malloc(nrang*sizeof(double));
	double * velo_arr = malloc(nrang*sizeof(double));
	double * amp0_arr = malloc(nrang*sizeof(double));
	int * qflg = malloc(nrang*sizeof(int));
	complex double ** acfs = malloc(nrang*sizeof(complex double *));

	for(i=0;i<nrang;i++)
        {

                /*array with the irregularity decay time for each range gate*/
                if (fitacf->rng[i].qflg) {
			t_d_arr[i] = lambda/(fitacf->rng[i].w_l*2.*PI);
		} else {
			t_d_arr[i] = 1000.0;
		}
	        /*array with the irregularity growth time for each range gate*/
		t_g_arr[i] = t_g;
	        /*array with the irregularity lifetime for each range gate*/
		t_c_arr[i] = t_c;
        	/*array with the irregularity doppler velocity for each range gate*/
		v_dop_arr[i] = fitacf->rng[i].v;
		/*array with the irregularity doppler velocity for each range gate*/	
		velo_arr[i] = velo;
	        /* white noise level */
	        noise_lev = fitacf->noise.skynoise;
		/*array with the ACF amplitude for each range gate*/
                if (fitacf->rng[i].qflg) {
		    amp0_arr[i] = pow(10.0,fitacf->rng[i].p_0/10.0)*noise_lev; /* THIS IS IN LOG */
                } else {
                    amp0_arr[i] = 0;
                }
		/*flags to tell which range gates contain scatter*/
		qflg[i] = fitacf->rng[i].qflg;
		/*Creating the output array for ACFs*/
                acfs[i] = malloc(n_lags*sizeof(complex double));
		for(j=0;j<n_lags;j++) {
			acfs[i][j] = 0.+I*0.;
                }

                fprintf(stderr,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %d  %lf  %lf\n", t_d_arr[i],
                        t_g_arr[i], t_c_arr[i], v_dop_arr[i], velo_arr[i], noise_lev, amp0_arr[i], 
                        qflg[i], creal(acfs[i][0]), cimag(acfs[i][0]));
        }


		fprintf(stderr,"%hd  %hd  %hd  %hd  %hd  %hd\n",prm->time.yr,prm->time.mo,
                        prm->time.dy,prm->time.hr,prm->time.mt,prm->time.sc);
                fprintf(stderr,"%d  %lf %lf  %hd  %lf  %d  %d  %lf  %lf  %lf  %d\n",cpid,freq, lambda,
                        prm->bmnum,noise_lev,nave,lagfr,dt,smsep,rngsep,nrang);


		/*create a structure to store the raw samples from each pulse sequence*/
		complex double * raw_samples = malloc(n_samples*nave*sizeof(complex double));
                /*fprintf(stderr, "makeRadarParm2 \n");
		makeRadarParm2(prm, argv, argc, cpid, nave, lagfr, smsep, noise_lev, amp0, n_samples,
                    dt, n_pul, n_lags, nrang, rngsep, freq, pulse_t,prm->time.yr,prm->time.mo,
		    prm->time.dy,prm->time.hr,prm->time.mt,prm->time.sc,prm->bmnum);*/

		/*call the simulation function*/
                fprintf(stderr, "sim_data \n");
                cpid = 150;
		sim_data(t_d_arr, t_g_arr, t_c_arr, v_dop_arr, qflg, velo_arr, amp0_arr, freq, noise_lev,
							noise_flg, nave, nrang, lagfr, smsep, cpid, life_dist,
							n_pul, cri_flg, n_lags, pulse_t, tau, dt, raw_samples, acfs, decayflg);

                for(i=0;i<nrang;i++)
		{
                	fprintf(stderr,"%lf %lf \n", creal(acfs[i][0]),cimag(acfs[i][0]));
                }

		if(!smp_flg)
		{
			/*fill the rawdata structure*/
			struct RawData * raw;
			raw = RawMake();

			raw->revision.major = 1;
			raw->revision.minor = 1;
			raw->thr=0.0;
			int * slist = malloc(nrang*sizeof(int));
			float * pwr0 = malloc(nrang*sizeof(float));
			float * acfd = malloc(nrang*n_lags*2*sizeof(float));
			float * xcfd = malloc(nrang*n_lags*2*sizeof(float));
			for(i=0;i<nrang;i++)
			{
				slist[i] = i;
				pwr0[i] = creal(acfs[i][0]);
				for(j=0;j<n_lags;j++)
				{
					acfd[i*n_lags*2+j*2] = creal(acfs[i][j]);
					acfd[i*n_lags*2+j*2+1] = cimag(acfs[i][j]);
					xcfd[i*n_lags*2+j*2] = 0.;
					xcfd[i*n_lags*2+j*2+1] = 0.;
				}
			}

			RawSetPwr(raw,nrang,pwr0,nrang,slist);
			RawSetACF(raw,nrang,n_lags,acfd,nrang,slist);
			RawSetXCF(raw,nrang,n_lags,xcfd,nrang,slist);
			i=RawFwrite(stdout,prm,raw);
			free(slist);
			free(pwr0);
			free(acfd);
			free(xcfd);
		}
		else
		{
			/*fill the iqdata structure*/
			struct IQ *iq;
			iq=IQMake();

			int16 * samples = malloc(n_samples*nave*2*2*sizeof(int16));
			for(i=0;i<nave;i++)
			{
				/*main array samples*/
				for(j=0;j<n_samples;j++)
				{
					samples[i*n_samples*2*2+j*2] = (int16)(creal(raw_samples[i*n_samples+j]));
					samples[i*n_samples*2*2+j*2+1] = (int16)(cimag(raw_samples[i*n_samples+j]));
				}
				/*interferometer array samples*/
				for(j=0;j<n_samples;j++)
				{
					samples[i*n_samples*2*2+j*2+n_samples] = 0;
					samples[i*n_samples*2*2+j*2+1+n_samples] = 0;
				}
			}

			int * badtr = malloc(nave*n_pul*2*sizeof(int));

			IQFwrite(stdout,prm,iq,badtr,samples);
			free(samples);
			free(badtr);
		}

		free(pulse_t);
		free(tau);
		free(raw_samples);

        /*free dynamically allocated memory*/
        for(i=0;i<nrang;i++)
          free(acfs[i]);
        free(acfs);
        free(qflg);
        free(t_d_arr);
        free(t_g_arr);
        free(t_c_arr);
        free(v_dop_arr);
        free(velo_arr);
        free(amp0_arr);
 

  } while ((s=FitFread(fitfp,prm,fitacf)) !=-1);

  FitFree(fitacf);
  fclose(fitfp);


  return 0;
}
