# Makefile for sim_test
# =====================
# Author: R.J.Barnes
# by R.J.Barnes
#
#
include $(MAKECFG).$(SYSTEM)

INCLUDE=-I$(IPATH)/base -I$(IPATH)/general -I$(IPATH)/superdarn
OBJS = sim_test.o
SRC=sim_test.c
DSTPATH = $(BINPATH)
OUTPUT = sim_test
LIBS=-lsim_data.1 -lradar.1 -lfitacf.1 -lraw.1 -lfit.1 -ldmap.1 -lrtime.1 -lrcnv.1 -lrmath.1 -liqdata.1
SLIB=-lm -lz
 

include $(MAKEBIN).$(SYSTEM)
