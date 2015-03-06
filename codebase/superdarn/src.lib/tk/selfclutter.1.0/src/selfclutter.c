/* selfclutter.c
   ==============
   Author: A.S.Reimer
*/

/*
 LICENSE AND DISCLAIMER
 
 Copyright (c) 2014 The Institute for Space and Atmospheric Study at 
 the University of Saskatchewan
 
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
 
 
 
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "rtypes.h"
#include "tsg.h"
#include "acf.h"
#include "fitblk.h"


/***********************************************************************************************/
/***********************************************************************************************/

/* IQ Voltage Based Self-Clutter Estimator */
int Cvse(struct TSGprm *prm,
 		 int16 *inbuf,int rngoff,int dflg,
		 int roffset,int ioffset,
		 int mplgs,int *lagtable[2],
  	         float *scbuf,
	         int xcf,int xcfoff,
                 int badrange,float atten,float *dco, int cleaned_acf) {
   /* 

      *prm          - pointer to the radar parameter structure
      *inbuf        - pointer to the complex voltage data
      rngoff        - offset needed to get to the next complex voltage sample from main antenna array (see below)
      dflg          - whether or not data was recorded using a digital receiver
      roffset       - offset for the real component of the complex numbers pointed to by inbuf
      ioffset       - offset for the imaginary component of the complex numbers stored to by inbuf
      mplgs         - number of lags
      *lagtable[2]  - the lag table
      *scbuf        - pointer to the output "self-clutter" array (same as acfbuf or xcfbuf in make_raw)
      xcf           - flag for processing voltage data using interfermeter voltages or not
      xcfoff        - offset needed to get to the next complex voltage sample from the interferometer antenna array (see below)
      badrange      - at which range gates, lag0 power measurement is "bad" (ie. more than one pulse has been TXed)
      atten         - an attenuation correction (always 0 in my experience)
      *dco          - ??? always NULL in my experience
      cleaned_acf   - a flag to specify if the acf minus the self-clutter should be returned in scbuf, or just the self-clutter


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
  * Analogue receiver:    rngoff        xcfoff
  *
  * No XCFs                 2             0
  * With XCFs               4             2
  *
  * Digital Receiver      rngoff        xcfoff
  * (APL & Saskatoon)     
  * No XCFs                 2             0
  * With XCFs               2             nsamp
  *
  * Digital Receiver      rngoff        xcfoff
  * (Alaska)
  * No XCFs                 2             0
  * With XCFs               2             8192 (half DMA buffer size)
  *
  *
  */

   int sdelay=0;
   int range;
   int sampleunit;
   int offset1;
   float real;
   float imag;
   int lag;
   int sample1;
   int sample2;
   int nrang;
   float temp1;
   float temp2;
   int offset2;

   int smpfr;
   int tp_in_tau;
   int S1,S2;
   int pul,pul1,pul2;
   int mppul;
   int re_term1,re_term2,re_term3,re_term4;
   int im_term1,im_term2,im_term3,im_term4;
   int r1[prm->mppul],r2[prm->mppul];           /*initialize integer arrays*/
   int sample3,sample4;
   int offset3,offset4;
   int temp;

   float dcor1=0;
   float dcor2=0;
   float dcoi1=0;
   float dcoi2=0;

   if ((dco !=NULL)) {
     if (xcf==ACF_PART) {
       dcor1=dco[0];
       dcor2=dco[0];
       dcoi1=dco[1];
       dcoi2=dco[1];       
     } else {
       dcor1=dco[0];
       dcor2=dco[2];
       dcoi1=dco[1];
       dcoi2=dco[3];
     }
   }

   nrang = prm->nrang;
   smpfr = (prm->lagfr / prm->smsep); /*number of voltage samples to the first range*/
   mppul = prm->mppul;                /*number of pulses in the pulse sequence*/


   if (dflg) sdelay=prm->smdelay; /* digital receiver delay term equal to skpnum */
   sampleunit = (prm->mpinc / prm->smsep) * rngoff; /* equal to number of tx pulses */
                                                    /* that can fit in one lag time times rngoff (tau/tp)*rngoff */
   tp_in_tau = (prm->mpinc / prm->smsep);

   for(range=0;range < nrang ; range++) {

         /*Calculate offsets for the voltage samples for the current range gate*/
         offset1 = (range+sdelay) * rngoff;
	 if (xcf == ACF_PART) offset2 = offset1;
	 else offset2 = ((range+sdelay) * rngoff) + xcfoff;

	 for(lag=0;lag < mplgs; lag++) {

             /* Determine which ranges are intefering in the current lag */
             /* samples are we using for the current lag*/
             S1=tp_in_tau*lagtable[0][lag] + range + smpfr;
             S2=tp_in_tau*lagtable[1][lag] + range + smpfr;

             for (pul=0; pul < mppul; pul++) {
                 /*Find the pulses that were transmitted before the samples were recorded
                   and then save which range gates each pulse is coming from. */
                 if (prm->pat[pul]*tp_in_tau <= S1){
                     temp = (S1 - prm->pat[pul]*tp_in_tau - smpfr);
                     /*Also we need to check and make sure we only save interfering range 
                       gates where we have valid lag0 power.*/
                     if ((temp != range) && (temp >= 0) && (temp < nrang) && (temp < badrange)) {
                         r1[pul]= temp;
                     } else {
                         r1[pul]=-1000;
                     }
                 } else {
                     r1[pul]=-1000;
                 }
                 if (prm->pat[pul]*tp_in_tau <= S2){
                     temp = (S2 - prm->pat[pul]*tp_in_tau - smpfr);
                     if ((temp != range) && (temp >= 0) && (temp < nrang) && (temp < badrange)) {
                         r2[pul]= temp;
                     } else {
                         r2[pul]=-1000;
                     }
                 } else {
                     r2[pul]=-1000;
                 }
             }

             re_term1 = 0;
             re_term2 = 0;
             im_term1 = 0;
             im_term2 = 0;
             re_term3 = 0;
             re_term4 = 0;
             im_term3 = 0;
             im_term4 = 0;


             if ((range >= badrange) && (lag == 0)) {
                 sample1 =lagtable[0][mplgs]*sampleunit + offset1;        
                 sample2 =lagtable[1][mplgs]*sampleunit + offset2;
             } else { 
                 sample1 =lagtable[0][lag]*sampleunit + offset1;        
                 sample2 =lagtable[1][lag]*sampleunit + offset2;
             }

             re_term1 = (float) (inbuf[sample1+ roffset]);
             im_term1 = (float) (inbuf[sample1+ ioffset]);
             re_term2 = (float) (inbuf[sample2+ roffset]);
             im_term2 = (float) (inbuf[sample2+ ioffset]);

             re_term3 = (float) (inbuf[sample1+ roffset]);
             im_term3 = (float) (inbuf[sample1+ ioffset]);
             re_term4 = (float) (inbuf[sample2+ roffset]);
             im_term4 = (float) (inbuf[sample2+ ioffset]);

             /*Third term in the summation for the self-clutter estimate (V_i*V_j^*)*/
             for (pul1=0; pul1 < mppul; pul1++) {

                 if (r1[pul1] != -1000){
                     offset3 = (r1[pul1]+sdelay) * rngoff;  /*for the interfering voltages from sample #2*/
                     sample3 = lagtable[0][0]*sampleunit + offset3;        

                     re_term1 -= (float) (inbuf[sample3+ roffset]);
                     im_term1 -= (float) (inbuf[sample3+ ioffset]);

                 }
             }

             for (pul2=0; pul2 < mppul; pul2++) {
                 if (r2[pul2] != -1000) {
                     offset4 = (r2[pul2]+sdelay) * rngoff;  /*for the interfering voltages from sample #1*/
                     sample4 = lagtable[0][0]*sampleunit + offset4;

                     re_term2 -= (float) (inbuf[sample4+ roffset]);
                     im_term2 -= (float) (inbuf[sample4+ ioffset]);

                 }
             }
             
             /* sum the real and imaginary acfs */

             if (cleaned_acf == 1) {
                 real = (re_term1*re_term2 + im_term1*im_term2); 
                 imag = (re_term1*im_term2 - im_term1*re_term2); 
             } else {
                 real = re_term3*re_term4 + im_term3*im_term4 - (re_term1*re_term2 + im_term1*im_term2); 
                 imag = re_term3*im_term4 - im_term3*re_term4 - (re_term1*im_term2 - im_term1*re_term2); 
             }

             if (atten !=0) {
               real=real/atten;
               imag=imag/atten;
             }
     
             /* Save the clutter estimate to the scbuf (self-clutter buffer) */
             scbuf[range*(2*mplgs)+2*lag]   = real + scbuf[range*(2*mplgs)+2*lag];
             scbuf[range*(2*mplgs)+2*lag+1] = imag + scbuf[range*(2*mplgs)+2*lag+1];
       } 
   }

   return 0;
}  

/***********************************************************************************************/
/***********************************************************************************************/

/* Maximal IQ Voltage Based Self-Clutter Estimator */
int Cmvse(struct TSGprm *prm,
 		 int16 *inbuf,int rngoff,int dflg,
		 int roffset,int ioffset,
		 int mplgs,int *lagtable[2],
  	         float *scbuf,
	         int xcf,int xcfoff,
                 int badrange,float atten,float *dco, int cleaned_acf) {
   /* 

      *prm          - pointer to the radar parameter structure
      *inbuf        - pointer to the complex voltage data
      rngoff        - offset needed to get to the next complex voltage sample from main antenna array (see below)
      dflg          - whether or not data was recorded using a digital receiver
      roffset       - offset for the real component of the complex numbers pointed to by inbuf
      ioffset       - offset for the imaginary component of the complex numbers stored to by inbuf
      mplgs         - number of lags
      *lagtable[2]  - the lag table
      *scbuf        - pointer to the output "self-clutter" array (same as acfbuf or xcfbuf in make_raw)
      xcf           - flag for processing voltage data using interfermeter voltages or not
      xcfoff        - offset needed to get to the next complex voltage sample from the interferometer antenna array (see below)
      badrange      - at which range gates, lag0 power measurement is "bad" (ie. more than one pulse has been TXed)
      atten         - an attenuation correction (always 0 in my experience)
      *dco          - ??? always NULL in my experience
      cleaned_acf   - a flag to specify if the acf minus the self-clutter should be returned in scbuf, or just the self-clutter


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
  * Analogue receiver:    rngoff        xcfoff
  *
  * No XCFs                 2             0
  * With XCFs               4             2
  *
  * Digital Receiver      rngoff        xcfoff
  * (APL & Saskatoon)     
  * No XCFs                 2             0
  * With XCFs               2             nsamp
  *
  * Digital Receiver      rngoff        xcfoff
  * (Alaska)
  * No XCFs                 2             0
  * With XCFs               2             8192 (half DMA buffer size)
  *
  *
  */


   int sdelay=0;
   int range;
   int sampleunit;
   int offset1;
   float real;
   float imag;
   int lag;
   int sample1;
   int sample2;
   int nrang;
   float temp1;
   float temp2;
   int offset2;

   int smpfr;
   int tp_in_tau;
   int S1,S2;
   int pul,pul1,pul2;
   int mppul;
   float re_term1,re_term2,re_term3;
   float im_term1,im_term2,im_term3;
   float term1, term2, term3;
   int r1[prm->mppul],r2[prm->mppul];           /*initialize integer arrays*/
   int sample3,sample4;
   int offset3,offset4;
   int temp;

   float dcor1=0;
   float dcor2=0;
   float dcoi1=0;
   float dcoi2=0;

   if ((dco !=NULL)) {
     if (xcf==ACF_PART) {
       dcor1=dco[0];
       dcor2=dco[0];
       dcoi1=dco[1];
       dcoi2=dco[1];       
     } else {
       dcor1=dco[0];
       dcor2=dco[2];
       dcoi1=dco[1];
       dcoi2=dco[3];
     }
   }

   nrang = prm->nrang;
   smpfr = (prm->lagfr / prm->smsep); /*number of voltage samples to the first range*/
   mppul = prm->mppul;                /*number of pulses in the pulse sequence*/


   if (dflg) sdelay=prm->smdelay; /* digital receiver delay term equal to skpnum */
   sampleunit = (prm->mpinc / prm->smsep) * rngoff; /* equal to number of tx pulses */
                                                    /* that can fit in one lag time times rngoff (tau/tp)*rngoff */
   tp_in_tau = (prm->mpinc / prm->smsep);

   for(range=0;range < nrang ; range++) {

         /*Calculate offsets for the voltage samples for the current range gate*/
         offset1 = (range+sdelay) * rngoff;
	 if (xcf == ACF_PART) offset2 = offset1;
	 else offset2 = ((range+sdelay) * rngoff) + xcfoff;

	 for(lag=0;lag < mplgs; lag++) {

             /* Determine which ranges are intefering in the current lag */
             /* samples are we using for the current lag*/
             S1=tp_in_tau*lagtable[0][lag] + range + smpfr;
             S2=tp_in_tau*lagtable[1][lag] + range + smpfr;

             for (pul=0; pul < mppul; pul++) {
                 /*Find the pulses that were transmitted before the samples were recorded
                   and then save which range gates each pulse is coming from. */
                 if (prm->pat[pul]*tp_in_tau <= S1){
                     temp = (S1 - prm->pat[pul]*tp_in_tau - smpfr);
                     /*Also we need to check and make sure we only save interfering range 
                       gates where we have valid lag0 power.*/
                     if ((temp != range) && (temp >= 0) && (temp < nrang) && (temp < badrange)) {
                         r1[pul]= temp;
                     } else {
                         r1[pul]=-1000;
                     }
                 } else {
                     r1[pul]=-1000;
                 }
                 if (prm->pat[pul]*tp_in_tau <= S2){
                     temp = (S2 - prm->pat[pul]*tp_in_tau - smpfr);
                     if ((temp != range) && (temp >= 0) && (temp < nrang) && (temp < badrange)) {
                         r2[pul]= temp;
                     } else {
                         r2[pul]=-1000;
                     }
                 } else {
                     r2[pul]=-1000;
                 }
             }

             re_term1 = 0;
             re_term2 = 0;
             re_term3 = 0;
             im_term1 = 0;
             im_term2 = 0;
             im_term3 = 0;
             term1 = 0;
             term2 = 0;
             term3 = 0;

             if ((range >= badrange) && (lag == 0)) {
                 sample1 =lagtable[0][mplgs]*sampleunit + offset1;        
                 sample2 =lagtable[1][mplgs]*sampleunit + offset2;
             } else { 
                 sample1 =lagtable[0][0]*sampleunit + offset1;        
                 sample2 =lagtable[1][0]*sampleunit + offset2;
             }

             /*First term in the summation for the self-clutter estimate (V_r*V_j^*)*/
             for (pul=0; pul < mppul; pul++) {
                 if (r2[pul] != -1000) {

                     /*offset3 = (r1[pul]+sdelay) * rngoff; */ /*for the interfering voltages from sample #1*/
                     offset4 = (r2[pul]+sdelay) * rngoff;  /*for the interfering voltages from sample #2*/

                    /* sample1 = lagtable[0][0]*sampleunit + offset1;        
                     sample2 = lagtable[1][0]*sampleunit + offset2; */

                     /*sample3 = lagtable[0][0]*sampleunit + offset3;        */
                     sample4 = lagtable[0][0]*sampleunit + offset4;

                     /*power 1*/
                     temp1 = (float) (inbuf[sample1+ roffset]) * 
                             (float) (inbuf[sample1+ roffset]);

                     temp2 = (float) (inbuf[sample1 + ioffset]) * 
                             (float) (inbuf[sample1 + ioffset]);

                     re_term1 = temp1 + temp2;

                     /*power2*/
                     temp1 = (float) (inbuf[sample4 + roffset]) *
                             (float) (inbuf[sample4 + roffset]);
                     temp2 = (float) (inbuf[sample4 + ioffset]) * 
	              	     (float) (inbuf[sample4 + ioffset]); 

                     im_term1 = temp1 + temp2;

                     term1 += sqrt(re_term1*im_term1);

                 }
             }

             /*Second term in the summation for the self-clutter estimate (V_i*V_r^*)*/
             for (pul=0; pul < mppul; pul++) {
                 if (r1[pul] != -1000) {

                     offset3 = (r1[pul]+sdelay) * rngoff;  /*for the interfering voltages from sample #1*/
                     /*offset4 = (r2[pul]+sdelay) * rngoff; */  /*for the interfering voltages from sample #2*/

                     /*sample1 = lagtable[0][0]*sampleunit + offset1;        
                     sample2 = lagtable[1][0]*sampleunit + offset2;*/
	
                     sample3 = lagtable[0][0]*sampleunit + offset3;        
                     /*sample4 = lagtable[0][0]*sampleunit + offset4;    */

                     /*Real*/
                     temp1 = (float) (inbuf[sample3+ roffset]) * 
                             (float) (inbuf[sample3+ roffset]);

                     temp2 = (float) (inbuf[sample3 + ioffset]) * 
                             (float) (inbuf[sample3 + ioffset]);

                     re_term2 = temp1 + temp2;

                     /*Imaginary*/
                     temp1 = (float) (inbuf[sample2 + roffset]) *
                             (float) (inbuf[sample2 + roffset]);
                     temp2 = (float) (inbuf[sample2 + ioffset]) * 
	              	     (float) (inbuf[sample2 + ioffset]); 

                     im_term2 = temp1 + temp2;
 
                     term2 += sqrt(re_term2*im_term2);

                 }
             }

             /*Third term in the summation for the self-clutter estimate (V_i*V_j^*)*/
             for (pul1=0; pul1 < mppul; pul1++) {
                 for (pul2=0; pul2 < mppul; pul2++) {
                     if ((r1[pul1] != -1000) && (r2[pul2] != -1000)) {

                         offset3 = (r1[pul1]+sdelay) * rngoff;  /*for the interfering voltages from sample #2*/
                         offset4 = (r2[pul2]+sdelay) * rngoff;  /*for the interfering voltages from sample #1*/

                         /*sample1 = lagtable[0][lag]*sampleunit + offset1;        
                           sample2 = lagtable[1][lag]*sampleunit + offset2;*/

                         sample3 = lagtable[0][0]*sampleunit + offset3;        
                         sample4 = lagtable[0][0]*sampleunit + offset4;

                         /*Real*/
                         temp1 = (float) (inbuf[sample3+ roffset]) * 
                                 (float) (inbuf[sample3+ roffset]);
    
                         temp2 = (float) (inbuf[sample3 + ioffset]) * 
                                 (float) (inbuf[sample3 + ioffset]);

                         re_term3 = temp1 + temp2;

                         /*Imaginary*/
                         temp1 = (float) (inbuf[sample4 + roffset]) *
                                 (float) (inbuf[sample4 + roffset]);
                         temp2 = (float) (inbuf[sample4 + ioffset]) * 
	                         (float) (inbuf[sample4 + ioffset]); 

                         im_term3 = temp1 + temp2;

                         term3 += sqrt(re_term3*im_term3);

                     }
                 }
             }
             

             /* sum the real and imaginary acfs */

             real = term1 + term2 + term3;
             imag = 0;


             if (atten !=0) {
               real=real/atten;
               imag=imag/atten;
             }
     
             /* Save the clutter estimate to the scbuf (self-clutter buffer) */
             scbuf[range*(2*mplgs)+2*lag]   = real + scbuf[range*(2*mplgs)+2*lag];
             scbuf[range*(2*mplgs)+2*lag+1] = imag + scbuf[range*(2*mplgs)+2*lag+1];
       } 
   }

   return 0;
}  

/***********************************************************************************************/
/***********************************************************************************************/

/* Maximal lag0 Power Based Self-Clutter Estimator */
int Cmpse(struct FitPrm *prm, int16 *lagtable[2], int gate,
  	  float *self_clutter, int badrange, float *pwr0) {
   /* 

      *prm          - pointer to the radar parameter structure
      *lagtable[2]  - the lag table
      *self_clutter - pointer to array to store output maximal self-clutter of each lag
      gate          - the range gate to calculate the self-clutter for
      badrange      - at which range gates, lag0 power measurement is "bad" (ie. more than one pulse has been TXed)
      *pwr0         - pointer to array of lag0 power (in mag units, NOT DB) at each range

      Example usage in FITACF (in the do_fit.c code after ptr[i].p_0 is set to SNR in units of DB around line 140):
      ***********************************************
      int status;
      int r;
      self_clutter = malloc(sizeof(float)*iptr->prm.mplgs);
      
      badrng=ACFBadLagZero(&iptr->prm,iptr->prm.mplgs,&iptr->prm.lag);

      for (r=0;r<iptr->prm.nrang;r++) {
          status = Cmpse(&iptr->prm, &iptr->prm.lag, r, &self_clutter, badrng, &pwrd);
          SOME_FUNCTION_THAT_USES_SELF-CLUTTER (ie. to make error bars)
      }
      free(self_clutter);
   */


   /* constants from prm structure */
   int nrang = prm->nrang;
   int smpfr = (prm->lagfr / prm->smsep); /*number of voltage samples to the first range*/
   int mppul = prm->mppul;                /*number of pulses in the pulse sequence*/
   int mplgs = prm->mplgs;
   int tp_in_tau = (prm->mpinc / prm->smsep);

   /* for loop indicies */
   int lag;
   int pul,pul1,pul2;

   /* temporary storage variables */
   int S1,S2;
   float term1, term2, term3;
   int r1[prm->mppul],r2[prm->mppul];           /*initialize integer arrays*/
   int temp;

   if (badrange < 0)
     badrange = nrang;

   for(lag=0;lag < mplgs; lag++) {
       fprintf(stderr,"Gate, lag:%d %d %d %d %d \n", gate, lag, lagtable[0][lag], lagtable[1][lag], smpfr);
       /* First, initialize self_clutter power to 0 */
       self_clutter[lag] = 0;

       /* Determine which ranges are intefering in the current lag */
       /* samples are we using for the current lag*/
       S1=tp_in_tau*lagtable[0][lag] + gate + smpfr;
       S2=tp_in_tau*lagtable[1][lag] + gate + smpfr;

       for (pul=0; pul < mppul; pul++) {
           /*Find the pulses that were transmitted before the samples were recorded
             and then save which range gates each pulse is coming from. */
           if (prm->pulse[pul]*tp_in_tau <= S1){
               temp = (S1 - prm->pulse[pul]*tp_in_tau - smpfr);
               /*Also we need to check and make sure we only save interfering range 
                 gates where we have valid lag0 power.*/
               if ((temp != gate) && (temp >= 0) && (temp < nrang) && (temp < badrange)) {
                   r1[pul]= temp;
               } else {
                   r1[pul]=-1000;
               }
           } else {
               r1[pul]=-1000;
           }
           if (prm->pulse[pul]*tp_in_tau <= S2){
               temp = (S2 - prm->pulse[pul]*tp_in_tau - smpfr);
               if ((temp != gate) && (temp >= 0) && (temp < nrang) && (temp < badrange)) {
                   r2[pul]= temp;
               } else {
                   r2[pul]=-1000;
               }
           } else {
               r2[pul]=-1000;
           }
       }

       /* initialize temporary variables */
       term1 = 0;
       term2 = 0;
       term3 = 0;

       /*First term in the summation for the self-clutter estimate (P_r*P_j^*)*/
       for (pul=0; pul < mppul; pul++) {
           if (r2[pul] != -1000) {
               term1 += sqrt(pwr0[gate]*pwr0[r2[pul]]);
           }
       }

       /*Second term in the summation for the self-clutter estimate (P_i*P_r^*)*/
       for (pul=0; pul < mppul; pul++) {
           if (r1[pul] != -1000) {
               term2 += sqrt(pwr0[gate]*pwr0[r1[pul]]);
           }
       }

       /*Third term in the summation for the self-clutter estimate (P_i*P_j^*)*/
       for (pul1=0; pul1 < mppul; pul1++) {
           for (pul2=0; pul2 < mppul; pul2++) {
               if ((r1[pul1] != -1000) && (r2[pul2] != -1000)) {
                   term3 += sqrt(pwr0[r1[pul1]]*pwr0[r2[pul2]]);
               }
           }
       }


       self_clutter[lag] = term1 + term2 + term3;


   } 

   return 0;
} 


