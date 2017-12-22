#include $(CS2_HOME)/make.inc

FC=gfortran -g -Wall -fbacktrace -ffree-line-length-0 -fcheck=all

SOURCES_F90 =\
module_mesh.f90\
module_discretization.f90\
module_solver.f90\
dagmg.f90\
dagmg_mumps.f90

OBJECTS=$(SOURCES_F90:.f90=.o) $(SOURCES_F77:.f=.o)



all: 


.PHONY: clean
clean:
	rm -f *.o *.mod lib*.a *.x

.PHONY: realclean
realclean: clean
	rm -f *~ \#*\#

ns2d: $(OBJECTS)
	$(FC) -o $@.x $@.f90 $(OBJECTS) /usr/lib/lapack/liblapack.so /usr/lib/libblas/libblas.so

%.o: %.f90
	$(FC) $(INCS) $(FCFLAGS) -c $< -o $@

%.o: %.f
	$(FC) $(INCS) $(FCFLAGS) -c $< -o $@
