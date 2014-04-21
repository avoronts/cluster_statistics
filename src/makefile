SHELL = /bin/sh

# Path to LAMMPS extraction directory
LAMMPS_ROOT = $(HOME)/moldin
LAMMPS_SRC = ./lmp #$(LAMMPS_ROOT)/mylammps/

# Uncomment the line below if using the MPI stubs library
MPI_STUBS = #-I$(LAMMPS_SRC)/STUBS

FC = mpiifort    # replace with your Fortran compiler
CXX = mpiicc   # replace with your C++ compiler

# Flags for Fortran compiler, C++ compiler, and C preprocessor, respectively
FFLAGS = -O2 -fPIC
CXXFLAGS = -O2 -fPIC 
CPPFLAGS =
LDFLAGS = -L. -llammps -L$(HOME)/lib

all : program.x move

program.x : LAMMPS.o LAMMPS-wrapper.o program.o
	$(FC) $(LDFLAGS) $+ -o $@

LAMMPS.o lammps.mod : LAMMPS.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $<

LAMMPS-wrapper.o : LAMMPS-wrapper.cpp LAMMPS-wrapper.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -I$(LAMMPS_SRC) -I.

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean :
	$(RM) *.o *.mod 

move :
	cp program.x ../test

dist :
	tar -czf Fortran-interface.tar.gz LAMMPS-wrapper.h LAMMPS-wrapper.cpp LAMMPS.F90 makefile README
