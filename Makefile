###################################################
#                                                 #
#   Makefile for package djangoh 4.6.17           #
#   interfaced via PYTHIA to LHAPDF               #
#                                                 #
#   HS: 24.02.2011                                #
#   last mod: 23.07.2021                          #
#                                                 #
###################################################


PATH=/Users/spiesber/fortran/LHAPDF6/bin:/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin
CXX=gfortran

PROGRAM=djangoh
DJANGOH=Users/spiesber/fortran/django/djangoh-4.6.17
LHAPDF=/Users/spiesber/fortran/LHAPDF6
LHAPDF_LIB=$(LHAPDF)/lib
LIBS=-L$(LHAPDF_LIB) -lLHAPDF -lmathlib -lkernlib

CFFLAGS = -Wall -pedantic -fno-automatic -g -fbacktrace -ffpe-trap=invalid,zero

OBJS = $(addsuffix .o, $(basename $(SRCS)))

SRCS = djangoh_h.f djangoh_l.f djangoh_u.f djangoh_t.f \
	djangoh_p.f sophia-dj.f pythia-6.4.24-dj.f jetset7409.f

%.o: %.f
	$(CXX) $(CFFLAGS) -o $@ -c $<

gmc_random.o: gmc_random.f
	$(CXX) $(CFFLAGS) -ffree-form -c gmc_random.f -o gmc_random.o

djangoh: $(OBJS) gmc_random.o
	$(CXX) $(CFFLAGS) -o $@ $(OBJS) gmc_random.o $(LIBS)



clean:
	rm -f $(OBJS) $(PROGRAM)
	rm -f gmc_random.o

