FC=/usr/lib64/openmpi/bin/mpif90
FFLAGS=-lblas -llapack
LIB=/shared/lib/libscalapack.a
SRCS=$(wildcard *.f90)
PROGS=$(patsubst %.f90,%,$(SRCS))

all: $(PROGS)

%: %.f90
	$(FC) $(FFLAGS) -o $@ $< $(LIB)
