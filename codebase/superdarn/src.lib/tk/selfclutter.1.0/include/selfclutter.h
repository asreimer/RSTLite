/* selfclutter.h
   =====
   Author: A.S. Reimer
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




int EstimateSelfClutter(struct TSGprm *prm,
 		 int16 *inbuf,int rngoff,int dflg,
		 int roffset,int ioffset,
		 int mplgs,int *lagtable[2],
  	         float *scbuf,
	         int xcf,int xcfoff,
                 int badrange,float atten,float *dco, int cleaned_acf, int slow_way);




