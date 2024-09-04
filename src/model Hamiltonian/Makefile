.SUFFIXES: .o .f90
OPT = -g -fbacktrace -fcheck=all 
#-Wall
#-fdump-tree-original  

#--- compilers ----
CF = gfortran
#CF = ifort
#
#--- MPI ----
#CF = mpif90
#CF = mpiicc

#LIBS =  /share/sw/lapack/3.9.0/gcc-9/lib/liblapack.a
CFT = $(CF) -c $(OPT)

RM = rm -f

all: lanczos 

clean:
	$(RM) *.o

lanczos: lanczos.f90
	$(CF) -o lanczos $(OPT) lanczos.f90 -lm -llapack -lblas  
