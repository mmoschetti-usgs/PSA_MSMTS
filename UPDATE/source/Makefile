rsobjects = responseSpectral.o \
	 rscalc_interp_acc.o
cobjects = lowcornfreq.o get_amp_phase.o \
	 smooth_interpolate.o 
mcobjects = modiflowcornfreq.o get_amp_phase.o \
	 smooth_interpolate.o 
rotobjects = gmrotd50.o rotate.o \
	 rscalc_interp_acc.o trapint.o
rotSobjects = rotd50.o rotate.o \
	 rscalc_interp_acc.o trapint.o

FC = gfortran
CC = gcc
CCFLAGS = -I/home/mmoschetti/OK_KS_database/include -L/home/mmoschetti/OK_KS_database/lib -lgsl -lgslcblas -lm 

All : RScalc LowCornFreq Within \
	RotD50 Mlcf \
	 GetZeros CompGT GmRotD50 PolySub

RScalc : $(rsobjects)
	$(FC) -o RScalc $(rsobjects)
LowCornFreq: $(cobjects)
	$(FC) -o LowCornFreq $(cobjects)
Mlcf: $(mcobjects)
	$(FC) -o Mlcf $(mcobjects)
Within: within.c
	$(CC) -o Within within.c
GetZeros: getZeros.c
	$(CC) -o GetZeros getZeros.c
CompGT : comp.c
	$(CC) -o CompGT comp.c
GmRotD50: $(rotobjects)
	$(FC) -o GmRotD50 $(rotobjects)
RotD50: $(rotSobjects)
	$(FC) -o RotD50 $(rotSobjects)
PolySub: polysub.c
	$(CC) $(CCFLAGS) -o PolySub polysub.c

fork.o : fork.f
get_amp_phase.o : get_amp_phase.f
smooth_interpolate.o : smooth_interpolate.f
rscalc_interp_acc.o : rscalc_interp_acc.f
rotate.o : rotate.f
rscalc_ts.o : rscalc_ts.f
trapint.o : trapint.f

clean :
	rm $(cobjects)
