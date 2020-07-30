#Makefile using new version of p3dfft
FC      = mpif90
FFLAGS  = -fast -Mfree -Mextend -Mbounds -Minfo -acc ta=tesla,managed
#FFLAGS  = -fast -Mfree -Mextend -Mbounds -Minfo 
P3DFFT  = -L ./p3dfft_summit/lib -I ./p3dfft_summit/include -lp3dfft -L $(OLCF_FFTW_ROOT)/lib -l fftw3 -mp
LINPLIB = /ccs/home/radha/linpack/linpack_summit_dev/liblinpack.a

OBJECTS = ModGparams.mod ModFunc.mod ModInitials.mod ModElastic.mod

all: $(OBJECTS) ran_2.o Multi_MEUMAPPS.o Multi_MEUMAPPS.x

# 
ModGparams.mod: ModGparams.f
	$(FC) -c -o $@ $< $(FFLAGS) $(P3DFFT)

ModFunc.mod: ModFunc.f
	$(FC) -c -o $@ $< $(FFLAGS) $(P3DFFT)

ran_2.o: ran_2.f
	$(FC) -c -o $@ $< $(FFLAGS) $(P3DFFT)

ModInitials.mod: ModInitials.f
	$(FC) -c -o $@ $< $(FFLAGS) $(P3DFFT)

ModElastic.mod: ModElastic.f
	$(FC) -c -o $@ $< $(FFLAGS) $(P3DFFT)

Multi_MEUMAPPS.o: Multi_MEUMAPPS.f
	$(FC) -c -o $@ $< $(FFLAGS) $(P3DFFT)

Multi_MEUMAPPS.x : Multi_MEUMAPPS.o $(OBJECTS) ran_2.o $(LINPLIB)
	$(FC) -o $@ $(FFLAGS) Multi_MEUMAPPS.o $(OBJECTS) ran_2.o $(LINPLIB) $(P3DFFT)
