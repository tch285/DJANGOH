pwd    = ${PWD}
SHELL  = /bin/sh
AR     = ar
ARFLAGS= -cr
CC     = gcc
CXX    = gfortran
GFF    = gfortran
CXX    = gfortran

PROGRAM=djangoh
LHAPDF=${LHAPDF5}
#PATH=${LHAPDF}/bin:/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:./
PATH=${LHAPDF}:/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:./

#LHAPDF_LIB=$(LHAPDF)/lib
LHAPDF_LIB=$(LHAPDF)

LIBS=-L$(LHAPDF_LIB) -lLHAPDF 	

CFFLAGS = -O -Wall -pedantic -fno-automatic  -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero
FCFLAGS = -O -Wall -pedantic -fno-automatic  -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero
FGFLAGS = -O -Wall -pedantic -fno-automatic  -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero

INCDIR = -I$(LHAPDF)/include 
