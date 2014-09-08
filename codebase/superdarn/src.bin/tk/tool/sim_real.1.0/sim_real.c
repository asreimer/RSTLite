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


/* this is a driver program from the data simulator */
/* It primes the simulator with parameters read from fitacf, fitex2, or lmfit data files. */

int main(int argc,char *argv[])
{

  /********************************************************
  ** definitions of variables needed for data generation **
  ********************************************************/

  double t_g = 1e-6;                        /*irregularity growth time*/
  double t_c = 1000.;                       /*precipitation time constant*/
  double velo = 0;                          /*velocity distribution width*/
  double c = 3.e8;                          /*Speed of light (m/s)*/
  double lambda;                            /*transmit wavelength*/
  double freq;                              /*transmit frequency*/
  int noise_flg = 1;                        /*flag to indicate whether white noise is included*/
  double noise_lev = 0.;                    /*white noise level (ratio)*/
  int lagfr = 4;                            /*lag to first range*/
  int life_dist = 0;                        /*lifetime distribution*/
  int cpid = 150;                           /*control program ID number*/
  int n = 20000;
  double dt;                                /*basic lag time*/
  int cri_flg = 1;                          /*cross-range interference flag*/
  int smp_flg = 0;                          /*output raw samples flag*/
  int decayflg = 0;
  int scflg = 0;

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
  char * filename1;
  char * filename2;

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
  "-nonoise: do not add noise to the simulation\n"
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
    /*command line growth time*/
    else if (strcmp(argv[i], "-t_g") == 0)
      t_g = 1e-6*atoi(argv[i+1]);
    else if (strcmp(argv[i], "-n") == 0)
      n = atoi(argv[i+1]);
    /*command line precipitation time constant*/
    else if (strcmp(argv[i], "-t_c") == 0)
      t_c = 1e-3*atoi(argv[i+1]);
    /*command line noise*/
    else if (strcmp(argv[i], "-noise") == 0)
    {
      noise_flg = 0;
    }
    /*command line CRI flag*/
    else if (strcmp(argv[i], "-nocri") == 0)
      cri_flg = 0;
    /*command line output samples*/
    else if (strcmp(argv[i], "-samples") == 0)
      smp_flg = 1;
    else if (strcmp(argv[i], "-scf") == 0)
    {
      scflg = 1;
      filename1 = argv[i+1];
      filename2 = argv[i+2];
    }
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

  FILE * fitfp1=NULL;
  FILE * fitfp2=NULL;
  /* Open the fit file */
  FILE * fitfp=fopen(filename,"r");
  if(fitfp==NULL)
  {
    fprintf(stderr,"File %s not found.\n",filename);
    exit(-1);
  }

  fprintf(stderr,"Reading file: %s\n",filename);
  if (scflg)
  {
    fitfp1=fopen(filename1,"w");
    fitfp2=fopen(filename2,"w");
    fprintf(stderr,"Saving to files: %s and %s\n",filename1,filename2);
  
    if(fitfp1==NULL)
    {
      fprintf(stderr,"Cannot open file %s.\n",filename1);
      exit(-1);
    }
    if(fitfp2==NULL)
    {
      fprintf(stderr,"Cannot open file %s.\n",filename2);
      exit(-1);
    }
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
        RadarParmSetOriginCommand(prm,"sim_real");
        RadarParmSetCombf(prm,"sim_real");

        if (scflg)
        {
          prm->scf=1;
        } else {
          prm->scf=0;
        }

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

        /* GET THE PULSE TABLE */
        pulse_t = malloc(n_pul*sizeof(int));
	for (i=0;i<prm->mppul;i++) pulse_t[i]=prm->pulse[i];

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
            }
        }

        /*control program dependent variables*/
	taus = dt/smsep;                                      /*lag time in samples*/
	n_samples = (pulse_t[n_pul-1]*taus+nrang+lagfr);      /*number of samples in 1 pulse sequence*/
	
	double * t_d_arr = malloc(nrang*sizeof(double));
	double * t_g_arr = malloc(nrang*sizeof(double));
	double * t_c_arr = malloc(nrang*sizeof(double));
	double * v_dop_arr = malloc(nrang*sizeof(double));
	double * velo_arr = malloc(nrang*sizeof(double));
	double * amp0_arr = malloc(nrang*sizeof(double));
	int * qflg = malloc(nrang*sizeof(int));
	complex double ** acfs = malloc(nrang*sizeof(complex double *));
	complex double ** scfs = malloc(nrang*sizeof(complex double *));
        /*create a structure to store the raw samples from each pulse sequence*/
	complex double * raw_samples = malloc(n_samples*nave*sizeof(complex double));
        for(i=0; i < n_samples*nave; i++)
            raw_samples[i] = 0. + I*0.;

        /* Since iqdat and rawacf files do not calculate the mean noise, let's clear it here */
        prm->noise.mean = 0.0;

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
	        noise_lev = sqrt(fitacf->noise.skynoise); /* Need the noise from the clear frequency search! */
		/*array with the ACF amplitude for each range gate*/
                if (fitacf->rng[i].qflg) {
                    /* Since we are extracting the amplitude from the lag0 from a fit file,
                       which was derived using the fitacf->noise.skynoise (noise.mean currently 
                       isn't properly set by fitacf... but noise.skynoise and noise.mean are equal
                       for fitex files and noise.skynoise for fitex and fitacf are equal so use
                       that for now */
                       
		    amp0_arr[i] = sqrt(pow(10.0,fitacf->rng[i].p_0/10.0)*fitacf->noise.skynoise + fitacf->noise.skynoise);
                } else {
                    amp0_arr[i] = 0;
                }
		/*flags to tell which range gates contain scatter*/
		qflg[i] = fitacf->rng[i].qflg;
		/*Creating the output array for ACFs*/
                acfs[i] = malloc(n_lags*sizeof(complex double));
                scfs[i] = malloc(n_lags*sizeof(complex double));
		for(j=0;j<n_lags;j++) {
			acfs[i][j] = 0.+I*0.;
                        scfs[i][j] = 0.+I*0.;
                }

                /* fprintf(stderr,"%lf  %lf  %lf  %lf  %lf  %lf  %lf  %d  %lf  %lf\n", t_d_arr[i],
                        t_g_arr[i], t_c_arr[i], v_dop_arr[i], velo_arr[i], noise_lev, amp0_arr[i], 
                        qflg[i], creal(acfs[i][0]), cimag(acfs[i][0])); */
        }


		fprintf(stderr,"%hd  %hd  %hd  %hd  %hd  %hd\n",prm->time.yr,prm->time.mo,
                        prm->time.dy,prm->time.hr,prm->time.mt,prm->time.sc);
                fprintf(stderr,"%d  %lf %lf  %hd  %lf  %d  %d  %lf  %lf  %lf  %d\n",cpid,freq, lambda,
                        prm->bmnum,noise_lev,nave,lagfr,dt,smsep,rngsep,nrang);



		/*call the simulation function*/
                cpid = 150;
                fprintf(stderr,"SIMDATA");
		sim_data(t_d_arr, t_g_arr, t_c_arr, v_dop_arr, qflg, velo_arr, amp0_arr, freq, noise_lev,
							noise_flg, nave, nrang, lagfr, smsep, cpid, life_dist,
							n_pul, cri_flg, n_lags, pulse_t, tau, dt, raw_samples, 
                                                        acfs, scfs, scflg, n, decayflg);


		if(!smp_flg || scflg)
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
                        float * scfd = malloc(nrang*n_lags*2*sizeof(float));
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
					scfd[i*n_lags*2+j*2] = creal(scfs[i][j]);
					scfd[i*n_lags*2+j*2+1] = cimag(scfs[i][j]);
				}
			}

			RawSetPwr(raw,nrang,pwr0,0,NULL);
			RawSetACF(raw,nrang,n_lags,acfd,0,NULL);
			RawSetXCF(raw,nrang,n_lags,xcfd,0,NULL);
			RawSetSCF(raw,nrang,n_lags,scfd,0,NULL);

                        if (scflg) {
			    i=RawFwrite(fitfp1,prm,raw);
                        } else {
			    i=RawFwrite(stdout,prm,raw);
                        }


			free(slist);
			free(pwr0);
			free(acfd);
			free(xcfd);
                        free(scfd);
                        RawFree(raw);
		}
		if(smp_flg || scflg)
		{
			/*fill the iqdata structure*/
			struct IQ *iq;
			iq=IQMake();

			iq->revision.major = 1;
			iq->revision.minor = 1;
                        iq->seqnum = nave;
                        iq->chnnum = 1;
                        iq->smpnum = n_samples;
                        iq->skpnum = lagfr;

 
 			int16 * samples = malloc(n_samples*nave*2*2*sizeof(int16));
                        int * offset = malloc(nave*sizeof(int));
                        float * noise = malloc(nave*sizeof(float));
                        int * atten = malloc(nave*sizeof(int));
                        struct timespec * time = malloc(nave*sizeof(struct timespec));
                        int * size = malloc(nave*sizeof(int));

			for(i=0;i<nave;i++)
			{
                                offset[i] = i*n_samples*2*2;
                                noise[i] = noise_lev;
                                atten[i] = 0;

                                time[i].tv_sec = 0;
                                time[i].tv_nsec = 0;
                                size[i] = n_samples*2*2;
				/*main array samples*/
				for(j=0;j<n_samples;j++)
				{
					samples[i*n_samples*2*2+j*2] = (int16)(cimag(raw_samples[i*n_samples+j]));
					samples[i*n_samples*2*2+j*2+1] = (int16)(creal(raw_samples[i*n_samples+j]));
				}
				/*interferometer array samples*/
				for(j=0;j<n_samples;j++)
				{
					samples[i*n_samples*2*2+j*2+n_samples*2] = 0;
					samples[i*n_samples*2*2+j*2+1+n_samples*2] = 0;
				}
			}

			unsigned int * badtr = malloc(nave*n_pul*2*sizeof(unsigned int));

                        IQSetAtten(iq,nave,atten);
                        IQSetNoise(iq,nave,noise);
                        IQSetOffset(iq,nave,offset); /* offset into the sample buffer */
                        IQSetSize(iq,nave,size);
                        IQSetTime(iq,nave,time);
                        fprintf(stderr,"WRITING IQ\n");
                        if (scflg)
                        {
			    IQFwrite(fitfp2,prm,iq,badtr,samples);
                        } else {
			    IQFwrite(stdout,prm,iq,badtr,samples);
                        }

                        free(time);
                        free(size);
                        free(offset);
                        free(noise);
                        free(atten);
			free(samples);
			free(badtr);
                        IQFree(iq);
		}

		free(pulse_t);
		free(tau);
		free(raw_samples);

        /*free dynamically allocated memory*/
        for(i=0;i<nrang;i++)
          free(acfs[i]);
        free(acfs);
        for(i=0;i<nrang;i++)
          free(scfs[i]);
        free(scfs);
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

  if (scflg)
  {
      fclose(fitfp1);
      fclose(fitfp2);
  }


  return 0;
}
