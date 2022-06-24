# DJANGOH - version 4.6.20-eic (24 June 2022)
Monte Carlo simulation for deep inelastic lepton nucleon scattering

## Overview (as given by the original developer)

DJANGOH performs event simulation of neutral and charged current 
lepton nucleon scattering. The program was developed originally 
for deep inelastic electron proton scattering at HERA, but has 
been extended and includes options for muon scattering, heavy 
nuclear targets, elastic scattering and polarized protons. 
The emphasis is put on the inclusion of radiative corrections, 
comprising single soft and hard photon emission and the complete 
set of electroweak 1-loop corrections in the Standard Model. 
Large-mass hadronic final states are generated by an interface 
to LEPTO which simulates QCD effects. Low-mass hadronic final 
states are included by an interface to SOPHIA. 

Please send complaints, observations, suggestions to 
spiesber@uni-mainz.de

## Additional notes for EIC/BNL use:

This repository is a fork of [the original DJANGOH]{https://github.com/spiesber/DJANGOH}
code.  The main simulation code is untouched, but modifications 
have been made to the user routine `HSUSER` such that the event 
record contains the necessary information to make `EICTree`s.  In
addition, the `config.mk` file has been changed to direct the 
system to the correct location of the LHAPDF libraries.

To install, first check where your LHAPDF libraries are.  This may 
be in a couple different environment variables:

```bash
echo $LHAPDF
echo $LHAPDF5
echo $LHAPDF6
```

Whichever one gives the LHAPDF libraries, go to `config.mk` and make 
the necessary modifications to line 11.  For example, if LHAPDF5 
contains the libraries, then change line 11 to

```bash
LHAPDF=${LHAPDF5}
```

Next, `cd` into this local repository, rename `makefile-sample` to 
`makefile` and run `make`.  In my experience, there were a ton of 
warnings, mostly yelling about obsolescent features in the Fortran 
code, but no errors.  The executable will be built in this directory.
Finally, add the following line of code to your startup script(s) 
such that it can find the right DJANGOH executable (if there are 
multiple):

```bash
export PATH="/path/to/DJANGOH:${PATH}"
```

Verify by trying

```bash
which djangoh
```

which should yield `/path/to/DJANGOH`.