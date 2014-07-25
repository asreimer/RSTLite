/* make_raw.c
    ==========
   Author: R.J.Barnes
*/

/*
 LICENSE AND DISCLAIMER
 
 Copyright (c) 2012 The Johns Hopkins University/Applied Physics Laboratory
 
 This file is part of the Radar Software Toolkit (RST).
 
 RST is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 any later version.
 
 RST is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with RST.  If not, see <http://www.gnu.org/licenses/>.
 
 
**************************************
VERSION 1.9 ASREIMER
EDITS to make_raw so IQ data from digital receivers
can be properly processed into rawacf files by this 
code. Issue was that xcfoff was hard-coded in such a
way as to set it to 2, even for digital receivers. See below.
**************************************

ADDED FROM ACFCALCULATE DOCUMENTATION:
  * For an analogue receiver, samples are multiplexed,
  * eg. sample 1 & 2 are the I & Q for the main array,
  * end sample 3 & 4 are the I & Q for the interferometer array.
  * For a digital receiver samples are not multiplexed,
  * so for the APL and Saskatoon implmentation the first X
  * samples are the I & Q samples from the main array,
  * followed by the I & Q samples from the secondary array.
  * The Alaskan implementation is slightly different again
  * using two halves of the total DMA buffer. 
  *
  * The differences in implementation are handled by a combination
  * of the rngoff and xcfoff arguments.
  *
  * Analogue receiver:    rngoff        xcfoff    (IQ samples for main and interferometer arrays are interleved)
  *                                               
  * No XCFs                 2             0       (this code expects iq->chnnum = 1)
  * With XCFs               4             2       (this code expects iq->chnnum = 2)
  *
  * Digital Receiver      rngoff        xcfoff    (All IQ for main then all IQ for interferometer)
  * (APL & Saskatoon)     
  * No XCFs                 2             0       (this code expects iq->chnnum = 1)
  * With XCFs               2             nsamp   (this code expects iq->chnnum = 1)
  *                                               (this code expects nsamp = 2*iq->smpnum)
  *
  * Digital Receiver      rngoff        xcfoff
  * (Alaska)
  * No XCFs                 2             0
  * With XCFs               2             8192 (half DMA buffer size)
  *
  *
  *
  * IMPORTANT THINGS TO NOTE:
  * -- iq->chnnum should NEVER be 2 unless IQ is from an analogue receiever.
  * -- Each I and Q values are stored as 16-bit numbers whereas iq->smpnum counts IQ pairs.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <zlib.h>
#include "rtypes.h"
#include "dmap.h"
#include "option.h"
#include "rtime.h"
#include "dmap.h"
#include "rprm.h"
#include "iq.h"
#include "rawdata.h"

#include "iqread.h"
#include "rawwrite.h"

#include "tsg.h"
#include "acf.h"
#include "radar.h"

#include "errstr.h"
#include "hlpstr.h"



#define ACF_PART 0
#define XCF_PART 1
#define REAL_BUF_OFFSET 0
#define IMAG_BUF_OFFSET 1

struct RadarNetwork *network;  
struct Radar *radar;
struct RadarSite *site;

struct OptionData opt;

float *pwr0=NULL;
float *acfd=NULL;
float *xcfd=NULL;

int loadlag(FILE *fp,int *lag[2]) {
  /* load the lag table */
  char txt[256];
  int lagnum=0,c,real,imag;
  while(fgets(txt,255,fp) !=NULL) {
    for (c=0;(txt[c] !=0) && (c<256);c++) 
      if (txt[c]=='#') break;
    if (txt[c]=='#') continue;
    if (sscanf(txt,"%d %d\n",&real,&imag) !=2) continue;
    lag[0][lagnum]=real;
    lag[1][lagnum]=imag; 
    lagnum++;
  }
  return lagnum-1;
}


int main (int argc,char *argv[]) {

  int arg;
  char *envstr;
  unsigned char help=0;
  unsigned char option=0;
  unsigned char vb=0;

  FILE *fp;
  struct RadarParm *prm;
  struct IQ *iq;
  struct RawData *raw;
  unsigned int *badtr=NULL;
  int16 *samples=NULL;  
  int16 *ptr;

  void *tmp=NULL;

  int i,j,n;
  int badrng=0;
  int aflg,abflg;
  int thr=0,lmt=0;
  int atstp=0;

  struct TSGprm tprm;
  int *lag[2]={NULL,NULL};
  int *pulse=NULL;


  int mplgs=0;

  int roff=REAL_BUF_OFFSET;
  int ioff=IMAG_BUF_OFFSET;
  
  float thrsh=0;
  int qiflg=0;
  int skpval=0;

  int s=0;

  /* ASREIMER 23 July 2014
     Need to initialize new variables */
  int digital;
  int rngoff;
  int xcfoff;
  /* End ASREIMER */

  prm=RadarParmMake();
  iq=IQMake();
  raw=RawMake();

  envstr=getenv("SD_RADAR");
  if (envstr==NULL) {
    fprintf(stderr,"Environment variable 'SD_RADAR' must be defined.\n");
    exit(-1);
  }

  fp=fopen(envstr,"r");

  if (fp==NULL) {
    fprintf(stderr,"Could not locate radar information file.\n");
    exit(-1);
  }

  network=RadarLoad(fp);
  fclose(fp); 
  if (network==NULL) {
    fprintf(stderr,"Failed to read radar information.\n");
    exit(-1);
  }

  envstr=getenv("SD_HDWPATH");
  if (envstr==NULL) {
    fprintf(stderr,"Environment variable 'SD_HDWPATH' must be defined.\n");
    exit(-1);
  }

  RadarLoadHardware(envstr,network);

  OptionAdd(&opt,"-help",'x',&help);
  OptionAdd(&opt,"-option",'x',&option);
  OptionAdd(&opt,"vb",'x',&vb);
  OptionAdd(&opt,"t",'f',&thrsh);
  OptionAdd(&opt,"qi",'x',&qiflg);
  OptionAdd(&opt,"skip",'i',&skpval);
  OptionAdd(&opt,"d",'x',&digital);  /* Added d option so we can process IQ data from digital receivers ASREIMER */

  arg=OptionProcess(1,argc,argv,&opt,NULL);

  if (help==1) {
    OptionPrintInfo(stdout,hlpstr);
    exit(0);
  }

  if (option==1) {
    OptionDump(stdout,&opt);
    exit(0);
  }

  if (qiflg) {
    roff=IMAG_BUF_OFFSET;
    ioff=REAL_BUF_OFFSET;
  }

  if (arg==argc) fp=stdin;
  else fp=fopen(argv[arg],"r");

  if (fp==NULL) {
    fprintf(stderr,"File not found.\n");
    exit(-1);
  }

  while (IQFread(fp,prm,iq,&badtr,&samples) !=-1) {
    if (vb) 
      fprintf(stderr,"%d-%d-%d %d:%d:%d beam=%d\n",prm->time.yr,prm->time.mo,
	     prm->time.dy,prm->time.hr,prm->time.mt,prm->time.sc,prm->bmnum);

    
    /* get the hardware info */

     radar=RadarGetRadar(network,prm->stid);
     if (radar==NULL) {
      fprintf(stderr,"Failed to get radar information.\n");
      break;
    }

    site=RadarYMDHMSGetSite(radar,prm->time.yr,prm->time.mo,
	  	          prm->time.dy,prm->time.hr,prm->time.mt,
                          prm->time.sc);

    if (site==NULL) {
      fprintf(stderr,"Failed to get site information.\n");
      break;
    }

    atstp=site->atten;

    /* zero our arrays */

    if (pwr0==NULL) tmp=malloc(sizeof(float)*prm->nrang);
    else tmp=realloc(pwr0,sizeof(float)*prm->nrang);
    if (tmp==NULL) {
      s=-1;
      break;
    }
    memset(tmp,0,sizeof(float)*prm->nrang);
    pwr0=tmp;

    if (acfd==NULL) tmp=malloc(sizeof(float)*2*prm->nrang*prm->mplgs);
    else tmp=realloc(acfd,sizeof(float)*2*prm->nrang*prm->mplgs);
    if (tmp==NULL) {
      s=-1;
      break;
    }
    memset(tmp,0,sizeof(float)*2*prm->nrang*prm->mplgs);
    acfd=tmp;

    if (xcfd==NULL) tmp=malloc(sizeof(float)*2*prm->nrang*prm->mplgs);
    else tmp=realloc(xcfd,sizeof(float)*2*prm->nrang*prm->mplgs);
    if (tmp==NULL) {
      s=-1;
      break;
    }
    memset(tmp,0,sizeof(float)*2*prm->nrang*prm->mplgs);
    xcfd=tmp;

  
    for (i=0;i<2;i++) {
      if (lag[i]==NULL) tmp=malloc(sizeof(int)*(prm->mplgs+1));
      else tmp=realloc(lag[i],sizeof(int)*(prm->mplgs+1));
      if (tmp==NULL) break;
      lag[i]=tmp;
      for (j=0;j<=prm->mplgs;j++) lag[i][j]=prm->lag[i][j];
   }
   if (i !=2) {
     s=-1;
     break;
   }
   mplgs=prm->mplgs;
 

   if (skpval==0) skpval=iq->skpnum;   

    if (pulse==NULL) tmp=malloc(sizeof(int)*prm->mppul);
    else tmp=realloc(pulse,sizeof(int)*prm->mppul);
    if (tmp==NULL) {
       s=-1;
       break;
    }
    pulse=tmp;

    for (i=0;i<prm->mppul;i++) pulse[i]=prm->pulse[i];
 
    tprm.nrang=prm->nrang;
    tprm.lagfr=prm->lagfr;
    tprm.smsep=prm->smsep;
    tprm.mpinc=prm->mpinc;
    tprm.txpl=prm->txpl;
    tprm.mppul=prm->mppul;
    tprm.smdelay=skpval;
    tprm.pat=pulse;
    
    badrng=ACFBadLagZero(&tprm,prm->mplgs,lag);

    /* ASREIMER 23 July 2014
       If processing IQ from a digital receiver we
       have a different range offset than if analogue */

    rngoff = 2*iq->chnnum;

    /* Properly process xcf for IQ from digital receivers.*/
    if (digital) {
      if (prm->xcf==1) {
        xcfoff = 2*iq->smpnum;
      } else {
        xcfoff = 0;
      }
    } else {
      if (prm->xcf==1) {
        xcfoff = 2;
      } else {
        xcfoff = 0;
      }
    }

    /* End ASREIMER */

    for (n=0;n<iq->seqnum;n++) {

      ptr=samples+iq->offset[n];


      aflg=ACFSumPower(&tprm,mplgs,lag,pwr0,
		       ptr,rngoff,skpval !=0, /* rngoff used to be 2*iq->chnum ASREIMER */
                       roff,ioff,badrng,
                       iq->noise[n],prm->mxpwr,prm->atten*atstp,
                       thr,lmt,&abflg);
       
      
      ACFCalculate(&tprm,ptr,rngoff,skpval !=0, /* rngoff used to be 2*iq->chnum ASREIMER */
		   roff,ioff,mplgs,
	  	   lag,acfd,ACF_PART,xcfoff,badrng,iq->atten[n]*atstp,NULL);

      if (prm->xcf==1) ACFCalculate(&tprm,ptr,
				     rngoff,skpval !=0, /* rngoff used to be 2*iq->chnum ASREIMER */
                                     roff,ioff,prm->mplgs,
	  	                     lag,xcfd,XCF_PART,xcfoff,badrng,
				     iq->atten[n]*atstp,NULL);

      if ((n>0) && (iq->atten[n] != iq->atten[n-1]))
              ACFNormalize(pwr0,acfd,xcfd,prm->nrang,prm->mplgs,atstp); 
          

    }   
   
     
    ACFAverage(pwr0,acfd,xcfd,prm->nave,prm->nrang,prm->mplgs);
    
    RawSetPwr(raw,prm->nrang,pwr0,0,NULL);
    RawSetACF(raw,prm->nrang,prm->mplgs,acfd,0,NULL);
    RawSetXCF(raw,prm->nrang,prm->mplgs,xcfd,0,NULL);
     
    raw->thr=thrsh;
    RawFwrite(stdout,prm,raw);
   

  }
  if (fp !=stdin) fclose(fp);
  if (s !=0) {
    fprintf(stderr,"Error allocating memory.\n");
    return -1;
  }
  return 0;
} 






















