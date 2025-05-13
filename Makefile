#
FC       = gfortran
FFLAGS   = -O -w
LIBS     = -lblas -llapack
#

cfiles: hoecs.f90
	$(FC) $(FFLAGS) hoecs.f90 $(LIBS) -o hoecs.x
	-/bin/rm *.mod

clean:
	-/bin/rm *.mod *.x
