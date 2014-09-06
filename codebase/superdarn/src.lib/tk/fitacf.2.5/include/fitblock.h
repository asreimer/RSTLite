/* fitblock.h
   ========
   Author: A.S.Reimer
*/

#ifndef _FITBLOCK_H
#define _FITBLOCK_H

struct FitBlock {
  struct FitPrm prm;
  struct complex *acfd;
  struct complex *xcfd;
};

#endif
