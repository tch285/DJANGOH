pwd    = ${PWD}
SHELL  = /bin/sh
AR     = ar
ARFLAGS= -cr
CC     = gcc
CXX    = gfortran
GFF    = gfortran
CXX    = gfortran

PROGRAM=djangoh

LIBS = $(shell lhapdf-config --ldflags)

CFFLAGS = -O -Wall -std=legacy -fno-automatic -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero
FCFLAGS = -O -Wall -std=legacy -fno-automatic -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero
FGFLAGS = -O -Wall -std=legacy -fno-automatic -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero

LHAPDFINC = $(shell lhapdf-config --incdir)
INCDIR = -I$(LHAPDFINC)
