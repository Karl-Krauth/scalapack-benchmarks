FC=mpiifort
FFLAGS=-i8 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
LINK=-L${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -lpthread -lm -ldl
SRCS=$(wildcard *.f90)
PROGS=$(patsubst %.f90,%,$(SRCS))

all: $(PROGS)

%: %.f90
	$(FC) $(FFLAGS) -o $@ $< $(LINK)
