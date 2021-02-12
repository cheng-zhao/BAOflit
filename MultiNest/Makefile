AR = ar r  
LINKLIB = ld -shared

OBJS = utils.o utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^

%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $^

default: libnest3.a

libnest3.so: $(OBJS) 
	$(LINKLIB) -o $(FLIBS) $@ $^ 
 
libnest3.a: $(OBJS) 
	$(AR) $@ $^ 
 
clean: 
	rm libnest3.*  *.o *.mod
