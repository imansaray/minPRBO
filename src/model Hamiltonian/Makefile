.SUFFIXES: .o .f90

#--- compilers, CF got serial, PCF for parallel ----
CF = gfortran
# CF = ifort
#
#--- MPI ----
#PCF = mpif90
PCF = mpiifort
#PCF = mpiicc

OPT = 
#OPT = -fdefault-real-8 
POPT = -g -traceback -fp-stack-check -check bounds

#---- Libraries, LIBS for serial, PBLIBS for parallel
LIBS = -lm -lgfortran -lblas -llapack  
PLIBS = -lm -lgfortran -lscalapack -lblas -llapack  

RM = rm -f

# Serial codes are listed first, then parallel codes
all: lanczos planczos_full planczos_minprbo	


clean:
	$(RM) *.o

# Serial Code
lanczos: lanczos.f90
	$(CF) -o lanczos $(OPT) lanczos.f90 $(LIBS)


# Parallel code

planczos_full: planczos_full.f90
	$(PCF) -o planczos_full $(POPT) planczos_full.f90 $(PLIBS)

planczos_prbo: planczos_minprbo.f90
	$(PCF) -o planczos_prbo $(POPT) planczos_prbo.f90 $(PLIBS)
