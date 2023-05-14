#============================================
FC      = gfortran
DBG     = 
#DBG    = -pg
#FFLAGS  = -O -C -traceback $(DBG)
FFLAGS  = -O3  $(DBG)
#FFLAGS  = $(DBG)
PROGRAM = xtLoop

LDR     = f90
LDFLAGS =
 
SRCS    = constants.f90 rng.f90 read_input.f90 metropolis.f90 mag.f90  \
	tLoop.f90 energy.f90 
OBJS    = $(SRCS:.f90=.o)
 
.SUFFIXES: .o .f90

.f90.o:
		$(FC) -c $(FFLAGS) $<
 

$(PROGRAM): $(OBJS) 
		$(FC) $(OBJS) $(FFLAGS) -o $@
 
clean:
		rm *.o *.mod
#================================================
