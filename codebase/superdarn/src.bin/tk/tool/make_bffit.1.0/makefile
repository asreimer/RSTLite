# Makefile for make_bffit
# =====================
#
#
include $(MAKECFG).$(SYSTEM)

INCLUDE=-I$(IPATH)/base -I$(IPATH)/general -I$(IPATH)/superdarn
OBJS = make_bffit.o

SRC=make_bffit.c

DSTPATH = $(BINPATH)
OUTPUT = make_bffit
LIBS=-loldraw.1 -loldfit.1 -lcfit.1 -lrscan.1 -lradar.1 -lfitacf.1 -lraw.1 -lfit.1 -ldmap.1 -lopt.1 -lrtime.1 -lrcnv.1 -lbffit.1 -lmpfit  -lrmath.1 -lacf.1 -lselfclutter.1 -lerror_estimates.1
SLIB=-lm -lz


include $(MAKEBIN).$(SYSTEM)
