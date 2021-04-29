
FC = gfortran
FFLAGS = -O0 -fopenmp  -g  -fbacktrace   -finit-real=snan   -fcheck=all -finit-integer=-88888888 
# -ffpe-trap=invalid,zero,overflow


OBJS = plume_rise_driver.o plume_rise_IS4FIRES.o plumeRise_PRMv1.o _plume_rise_test.o rconstants.o

PROG = plume_rise_test

all: $(PROG)

plume_rise_driver.o: plume_rise_IS4FIRES.o plumeRise_PRMv1.o rconstants.o
plumeRise_PRMv1.o: rconstants.o

%.o: %.f90
	$(FC) -c $(FFLAGS)  -o $@ $<


$(PROG): $(OBJS)
	$(FC) $(OBJS) -fopenmp -leccodes -leccodes_f90 -o $(PROG)
	mv *.o ../obj/
	mv *.mod ../obj/
	mv $(PROG) ../bin/

clean:
	rm -rf $(OBJS) $(PROG) 


ctags:
	ctags -R
