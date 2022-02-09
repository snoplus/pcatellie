CFLAGS = $(RATROOT)/build/linuxx8664gcc/calib/pmtcal/PCABits.o -fPIC -Wall -Wno-write-strings -DSWAP_BYTES -g -fdiagnostics-show-option
LFLAGS = -L$(ROOTSYS)/lib -L$(RATROOT)/lib

INCLUDES= `root-config --cflags --glibs` \
-I$(RATROOT)/include \
-I$(RATROOT)/include/external \
-lRATEvent_Linux \
-I$(RATROOT)/src/libpq/include
#-lRATEvent_Darwin \

CC= g++

.PHONY: all

all: Compare

Compare: compare.cc
	$(CC) $(CFLAGS) $(LFLAGS) -L/usr/lib -o $@ $^ $(INCLUDES)  -lSpectrum

clean:
	rm Compare
