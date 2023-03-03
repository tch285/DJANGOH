pwd    = $(PWD)
SHELL  = /bin/sh
AR     = ar
ARFLAGS= -cr
GFF    = gfortran

PROGRAM=djangoh
# Place the environment variable containing the LHAPDF libraries within the braces.
LHAPDF=$(LHAPDF6)

# could help if we used lhapdf-config but we don't
# PATH=$(LHAPDF)/bin:$PATH
LIBS=-L$(LHAPDF)/lib -lLHAPDF
INCDIR = -I$(LHAPDF)/include

GFFLAGS = -O -Wall -std=legacy -fno-automatic -fcheck=all -g -fbacktrace -ffpe-trap=invalid,zero

