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
#include "rtypes.h"
#include "tsg.h"
#include "acf.h"




 /* For an analogue receiver, samples are multiplexed,
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

 

int EstimateSelfClutter(struct TSGprm *prm,
 		 int16 *inbuf,int rngoff,int dflg,
		 int roffset,int ioffset,
		 int mplgs,int *lagtable[2],
  	         float *scbuf,
	         int xcf,int xcfoff,
                 int badrange,float atten,float *dco) {

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

             real = re_term3*re_term4 + im_term3*im_term4 - (re_term1*re_term2 + im_term1*im_term2);
             imag = re_term3*im_term4 - im_term3*re_term4 - (re_term1*im_term2 - im_term1*re_term2);

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



