FFLAGS = -O2
FORTRAN = gfortran -ffixed-line-length-none -frecord-marker=4

SOURCE = \
CONVLV.f  REALFT.f TWOFFT.f FOUR1.f\

OBJECTS = $(SOURCE:.f=.o)

all:	$(OBJECTS)
	$(FORTRAN) $(FFLAGS) -c gen_iso.f 
	$(FORTRAN) $(FFLAGS) $(OBJECTS) gen_iso.o -o gen_iso
	$(FORTRAN) $(FFLAGS) -c isothermal.f 
	$(FORTRAN) $(FFLAGS) $(OBJECTS) isothermal.o -o isothermal
	rm -rf *.o

