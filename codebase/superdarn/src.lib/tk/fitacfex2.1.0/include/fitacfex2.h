/* fitacfex2.h
   ==========
*/

#ifndef _fitacfex2_H
#define _fitacfex2_H
#include "fitblock.h"

void fitacfex2(struct RadarParm *prm,struct RawData *raw,
                struct FitData *fit,struct FitBlock *fblk, int print, int more_badlags);

#endif
