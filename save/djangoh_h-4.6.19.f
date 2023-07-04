C
      PROGRAM HSMAIN
C
C
C***********************************************************************
C
C
C                     H E R A C L E S
C
C
C     EVENT GENERATOR FOR DEEP-INELASTIC E-P SCATTERING
C     INCLUDING RADIATIVE CORRECTIONS
C
C***********************************************************************
C
C     AUTHORS:
C              H. SPIESBERGER, 
C                               Johannes Gutenberg Univseritaet 
C                               Mainz
C                               Email: spiesber@uni-mainz.de
C
C              A. KWIATKOWSKI
C                               II. INSTITUT FUER THEORETISCHE PHYSIK,
C                               UNIVERSITAET HAMBURG
C                               FRG
C
C              H.-J.MOEHRING
C                               SEKTION PHYSIK
C                               UNIVERSITAET LEIPZIG
C                               FRG
C
C     with extensions by T. Martini (for hadron polarization, 2012+13)
C
C     VERSION 4.6.19,  Sep 01 2021
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      EXTERNAL HSNCG1,HSNCG2,HSTSK1,HSTSK2,HSK1TS,HSK1K3
      EXTERNAL HSCCG1,HSCCG2,HSCCKL,HSCCQI,HSCCQF
      EXTERNAL HSELG1,HSELG2,HSELK1,HSELK2,HSELCO
      CHARACTER FNAME*19,FNAMET*19
      PARAMETER (NCHN2=1,NCHC2=2,NCHE2=3)
      PARAMETER (NCHN31=6,NCHN32=7,NCHN33=8,NCHN34=9)
      PARAMETER (NCHC31=12,NCHC32=13,NCHC33=14)
      PARAMETER (NCHE31=15,NCHE32=16,NCHE33=17)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSISGM/ TCUTQ,TCUTQS
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD,IDPVR
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSELEP/ IDIPOL
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSWGTC/ IWEIGS
      COMMON /HSONLY/ IHSONL
      COMMON /HSLTYP/ LEPIN1
      COMMON /HSASY0/ PONULL,CHNULL,ICPASY
C---
      PARAMETER(NDIM2=2,NBIN2=50)
      PARAMETER(NREG2N=2500)
C---                     NREG2N=NBIN2**NDIM2
      LOGICAL LGLO2,LLOC2
      double precision ntot2
      COMMON /HSSNC2/ SIG2,SIG2E,T2GGMA,T2GMAX(NREG2N),
     +                XX2(50,2),
     +                FFGO2,DNCG2,FFLO2,DNCL2,GOLD2,
     +                NTOT2,NM2(NREG2N),NDO2,
     +                NCAL2,NCA12,NCA22,IBIM2,JCOR2,
     +                LGLO2,LLOC2
      PARAMETER(NREG2C=2500)
C---                     NREG2C=NBIN2**NDIM2
      LOGICAL LGLO2C,LLOC2C
      double precision ntot2c
      COMMON /HSSCC2/ SIG2C,SIG2EC,T2GMAC,T2MAXC(NREG2C),
     +                XX2C(50,2),
     +                FFGO2C,DNCG2C,FFLO2C,DNCL2C,GOLD2C,
     +                NTOT2C,NM2C(NREG2C),NDO2C,
     +                NCAL2C,NCA12C,NCA22C,IBIM2C,JCOR2C,
     +                LGLO2C,LLOC2C
      PARAMETER(NREG2E=50)
      PARAMETER(NDIM1=1)
C---                     NREG2E=NBIN2**NDIM1
      LOGICAL LGLO2E,LLOC2E
      double precision ntot2e
      COMMON /HSSEL2/ SIG2L,SIG2EE,T2GMAE,T2MAXE(NREG2E),
     +                XX2E(50,1),
     +                FFGO2E,DNCG2E,FFLO2E,DNCL2E,GOLD2E,
     +                NTOT2E,NM2E(NREG2E),NDO2E,
     +                NCAL2E,NCA12E,NCA22E,IBIM2E,JCOR2E,
     +                LGLO2E,LLOC2E
      PARAMETER(NDIM31=5,NBIN31=6,NREG31=7776)
C---                     NREG31=NBIN31**NDIM31
      LOGICAL LGLO31,LLOC31
      double precision ntot31
      COMMON /HSSN31/ SIG31,SIG31E,T31GMA,T31MAX(NREG31),
     +                XX31(50,NDIM31),
     +                FFGO31,DNCG31,FFLO31,DNCL31,GOLD31,
     +                SI31,SI2N31,SWGT31,SCHI31,IT31,
     +                NM31(NREG31),NDO31,
     +                NTOT31,NCAL31,NCA131,NCA231,IBIM31,JCOR31,
     +                LGLO31,LLOC31
      PARAMETER(NDIM32=5,NBIN32=6,NREG32=7776)
C---                     NREG32=NBIN32**NDIM32
      LOGICAL LGLO32,LLOC32
      double precision ntot32
      COMMON /HSSN32/ SIG32,SIG32E,T32GMA,T32MAX(NREG32),
     +                XX32(50,NDIM32),
     +                FFGO32,DNCG32,FFLO32,DNCL32,GOLD32,
     +                SI32,SI2N32,SWGT32,SCHI32,IT32,
     +                NM32(NREG32),NDO32,
     +                NTOT32,NCAL32,NCA132,NCA232,IBIM32,JCOR32,
     +                LGLO32,LLOC32
      PARAMETER(NDIM33=5,NBIN33=6,NREG33=7776)
C---                     NREG33=NBIN33**NDIM33
      LOGICAL LGLO33,LLOC33
      double precision ntot33
      COMMON /HSSN33/ SIG33,SIG33E,T33GMA,T33MAX(NREG33),
     +                XX33(50,NDIM33),
     +                FFGO33,DNCG33,FFLO33,DNCL33,GOLD33,
     +                SI33,SI2N33,SWGT33,SCHI33,IT33,
     +                NM33(NREG33),NDO33,
     +                NTOT33,NCAL33,NCA133,NCA233,IBIM33,JCOR33,
     +                LGLO33,LLOC33
      PARAMETER(NDIM34=5,NBIN34=6,NREG34=7776)
C---                     NREG34=NBIN34**NDIM34
      LOGICAL LGLO34,LLOC34
      double precision ntot34
      COMMON /HSSN34/ SIG34,SIG34E,T34GMA,T34MAX(NREG34),
     +                XX34(50,NDIM34),
     +                FFGO34,DNCG34,FFLO34,DNCL34,GOLD34,
     +                SI34,SI2N34,SWGT34,SCHI34,IT34,
     +                NM34(NREG34),NDO34,
     +                NTOT34,NCAL34,NCA134,NCA234,IBIM34,JCOR34,
     +                LGLO34,LLOC34
      PARAMETER(NDM3CC=5)
      PARAMETER(NBN31C=6,NRG31C=7776)
C---                     NRG31C=NBN31C**NDM3CC
      LOGICAL LGL31C,LLC31C
      double precision ntt31c
      COMMON /HSSC31/ SIG31C,SG31EC,T31GMC,T31MXC(NRG31C),
     +                XX31C(50,5),
     +                FFG31C,DNG31C,FFL31C,DNL31C,GLD31C,
     +                SI31C,S2N31C,SWT31C,SCH31C,IT31C,
     +                NM31C(NRG31C),NDO31C,
     +                NTT31C,NCL31C,NC131C,NC231C,IBM31C,JCR31C,
     +                LGL31C,LLC31C
      PARAMETER(NBN32C=6,NRG32C=7776)
C---                     NRG32C=NBN32C**NDM3CC
      LOGICAL LGL32C,LLC32C
      double precision ntt32c
      COMMON /HSSC32/ SIG32C,SG32EC,T32GMC,T32MXC(NRG32C),
     +                XX32C(50,5),
     +                FFG32C,DNG32C,FFL32C,DNL32C,GLD32C,
     +                SI32C,S2N32C,SWT32C,SCH32C,IT32C,
     +                NM32C(NRG32C),NDO32C,
     +                NTT32C,NCL32C,NC132C,NC232C,IBM32C,JCR32C,
     +                LGL32C,LLC32C
      PARAMETER(NBN33C=6,NRG33C=7776)
C---                     NRG33C=NBN33C**NDM3CC
      LOGICAL LGL33C,LLC33C
      double precision ntt33c
      COMMON /HSSC33/ SIG33C,SG33EC,T33GMC,T33MXC(NRG33C),
     +                XX33C(50,5),
     +                FFG33C,DNG33C,FFL33C,DNL33C,GLD33C,
     +                SI33C,S2N33C,SWT33C,SCH33C,IT33C,
     +                NM33C(NRG33C),NDO33C,
     +                NTT33C,NCL33C,NC133C,NC233C,IBM33C,JCR33C,
     +                LGL33C,LLC33C
      PARAMETER(NDM3EL=4)
      PARAMETER(NBN31E=10,NRG31E=10000)
C---                     NRG31E=NBN31E**NDM3EL
      LOGICAL LGL31E,LLC31E
      double precision ntt31e
      COMMON /HSSE31/ SIG31L,SG31EE,T31GME,T31MXE(NRG31E),
     +                XX31E(50,5),
     +                FFG31E,DNG31E,FFL31E,DNL31E,GLD31E,
     +                SI31E,S2N31E,SWT31E,SCH31E,IT31E,
     +                NM31E(NRG31E),NDO31E,
     +                NTT31E,NCL31E,NC131E,NC231E,IBM31E,JCR31E,
     +                LGL31E,LLC31E
      PARAMETER(NBN32E=10,NRG32E=10000)
C---                     NRG32E=NBN32E**NDM3EL
      LOGICAL LGL32E,LLC32E
      double precision ntt32e
      COMMON /HSSE32/ SIG32L,SG32EE,T32GME,T32MXE(NRG32E),
     +                XX32E(50,5),
     +                FFG32E,DNG32E,FFL32E,DNL32E,GLD32E,
     +                SI32E,S2N32E,SWT32E,SCH32E,IT32E,
     +                NM32E(NRG32E),NDO32E,
     +                NTT32E,NCL32E,NC132E,NC232E,IBM32E,JCR32E,
     +                LGL32E,LLC32E
      PARAMETER(NBN33E=10,NRG33E=10000)
C---                     NRG33E=NBN33E**NDM3EL
      LOGICAL LGL33E,LLC33E
      double precision ntt33e
      COMMON /HSSE33/ SIG33L,SG33EE,T33GME,T33MXE(NRG33E),
     +                XX33E(50,5),
     +                FFG33E,DNG33E,FFL33E,DNL33E,GLD33E,
     +                SI33E,S2N33E,SWT33E,SCH33E,IT33E,
     +                NM33E(NRG33E),NDO33E,
     +                NTT33E,NCL33E,NC133E,NC233E,IBM33E,JCR33E,
     +                LGL33E,LLC33E
C-------------------------
      CHARACTER CHNAME*45
      COMMON /HSNAMC/ CHNAME(20)
      COMMON /HSNUME/ SIGTOT,SIGTRR,SIGG(20),SIGGRR(20),NEVENT,NEVE(20)
      COMMON /SPPASS/ NSOPH,NSPOUT,NFAILP,NSPACC
      COMMON /DJPASS/ NTOTDJ,NPASDJ,NQELAS,NFAILL,NFAILQ
      COMMON /HSPASS/ NREJCW
      COMMON /HSVGLP/ NPOVEG,NUMINT,NPHYP
      COMMON /HSRDIO/ ISDINP,ISDOUT
      COMMON /VGRES/  S1,S2,S3,S4
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      REAL            PYSTOP,PYSLAM
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMAX,NPYMIN
      COMMON /HSALFS/ PAR111,PAR112,PARL11,PARL19,MST111,MST115
      REAL            PAR111,PAR112,PARL11,PARL19
      INTEGER         MST111,MST115
      CHARACTER LHAPATHI*232
      CHARACTER LHAPATH*232
      COMMON /LHAPDFC/ LHAPATH
      SAVE /LHAPDFC/
      CHARACTER OUTFILENAM*80
      COMMON /HSOUTF/ OUTFILENAM,ICH
      COMMON /LEPTOU/ CUTDJ(14),LST(40),PARL(30),XDJ,YDJ,W2DJ,Q2DJ,UDJ
      REAL            CUTDJ,            PARL,    XDJ,YDJ,W2DJ,Q2DJ,UDJ
      CHARACTER TITLE*80
      CHARACTER CODE*10,CODEWD*10
      DIMENSION CODE(40)
      DIMENSION UIO(97)
      DIMENSION INT2C(5),ISAM2C(5),INT3C(15),ISAM3C(15)
C---for date and time, used as seeds for random number generator 
      integer today(3), now(3)
C
C-------------------------------
C---INITIALIZE DEFAULTS NOT YET INITIALIZED IN BLOCK DATA
C---AND SET STARTING VALUES FOR VARIOUS COUNTERS
      DATA FNAME /' HERACLES DATA FILE'/
      DATA INT2C,ISAM2C,INT3C,ISAM3C /40*0/
      DATA CODE/
     1'TITLE     ','EL-BEAM   ','PR-BEAM   ','KINEM-CUTS','EGAM-MIN  ',
     2'INT-OPT-NC','INT-OPT-CC','INT-POINTS','HYP-CUBES ','GSW-PARAM ',
     3'STRUCTFUNC','NFLAVORS  ','SAM-OPT-NC','SAM-OPT-CC','RNDM-SEEDS',
     4'GSW-MASS  ','THMIN-QRAD','FLONG     ','ALFAS     ','          ',
     5'EP-DIPOLE ','NUCLEUS   ','NUCL-MOD  ','LHAPATH   ','OUTFILENAM',
     6'THETA-CUT ','PT-CUT    ','POLPDF    ','          ','          ',
     7'WEIGHTS   ','ASYM-ONLY ','          ','          ','          ',
     8'INT-ONLY  ','TEST-OPT  ','IOUNITS   ','START     ','STOP      '/
      DATA TITLE/' '/
C
C---ACCURACY FOR INTEGRATING THE NON-RADIATIVE CONTRIBUTION
      DATA EPSO /1D-4/
C
C---DEFINE VEGAS PARAMETERS
      DATA ACCVEG,IPRVEG,IGRVEG   /1D-5,    3,     1/
      DATA NPOIN /20/
      DATA IWEIGR/0D0/
C
C---INITIALIZE AND TEST THE RANDOM NUMBER GENERATOR
      CALL HSRNST(12,34,56,78)
C     CALL HSRNTE(1)
C
C---For the moment, fix output units to standard output
      LUNOUT=6
      LUNTES=6
C---Move output file definitions to where 
C   file name is known, only set default name here, 
C   see code word OUTFILENAM
      OUTFILENAM='djangoh-default-output'
      ICH=INDEX(OUTFILENAM,' ')-1
      IODEF=0

C---PRINT THE TITLE
      WRITE(LUNOUT,9)
 9    FORMAT(
     1'**************************************************',
     2'*****************************',
     3//,10X,'                         HERACLES '
     4//,10X,'     Event generator for deep-inelastic e-P collisions '
     5 /,10X,'              including radiative corrections  '
     6//,10X,'                 VERSION 4.6.19, 01.09.2021 '//
     8//,10X,'                      H. Spiesberger '//
     9' **************************************************',
     1'****************************',//)
C
C***********************************************************************
C               READ INPUT DATA
C
C     STRUCTURE OF INPUT:
C                         1)  CODEWD  (A10)
C                         2)  CORRESPONDING DATA (FORMAT FREE)
C***********************************************************************
C
 1    CONTINUE
      READ(LUNIN,90,END=4) CODEWD
      WRITE(LUNOUT,91) CODEWD
      DO 2 ISW=1,40
      IF(CODEWD.EQ.CODE(ISW)) GOTO 3
 2    CONTINUE
      WRITE(LUNOUT,92)
      GOTO 1
 3    CONTINUE
C   TITLE
      IF (ISW.EQ.1) THEN
         GOTO 100
C   EL-BEAM
      ELSEIF (ISW.EQ.2) THEN
         GOTO 200
C   PR-BEAM
      ELSEIF (ISW.EQ.3) THEN
         GOTO 300
C   KINEM-CUTS
      ELSEIF (ISW.EQ.4) THEN
         GOTO 400
C   EGAM-MIN
      ELSEIF (ISW.EQ.5) THEN
         GOTO 500
C   INT-OPT-NC
      ELSEIF (ISW.EQ.6) THEN
         GOTO 600
C   INT-OPT-CC
      ELSEIF (ISW.EQ.7) THEN
         GOTO 700
C   INT-POINTS
      ELSEIF (ISW.EQ.8) THEN
         GOTO 800
C   HYP-CUBES
      ELSEIF (ISW.EQ.9) THEN
         GOTO 900
C   GSW-PARAM
      ELSEIF (ISW.EQ.10) THEN
         GOTO 1000
C   STRUCTFUNC
      ELSEIF (ISW.EQ.11) THEN
         GOTO 1100
C   NFLAVORS
      ELSEIF (ISW.EQ.12) THEN
         GOTO 1200
C   SAM-OPT-NC
      ELSEIF (ISW.EQ.13) THEN
         GOTO 1300
C   SAM-OPT-CC
      ELSEIF (ISW.EQ.14) THEN
         GOTO 1400
C   RNDM-SEEDS
      ELSEIF (ISW.EQ.15) THEN
         GOTO 1500
C   GSW-MASS
      ELSEIF (ISW.EQ.16) THEN
         GOTO 1600
C   THMIN-QRAD
      ELSEIF (ISW.EQ.17) THEN
         GOTO 1700
C   FLONG
      ELSEIF (ISW.EQ.18) THEN
         GOTO 1800
C   ALFAS
      ELSEIF (ISW.EQ.19) THEN
         GOTO 1900
C   not used
      ELSEIF (ISW.EQ.20) THEN
         GOTO 2000
C   EP-DIPOLE
      ELSEIF (ISW.EQ.21) THEN
         GOTO 2100
C   NUCLEUS
      ELSEIF (ISW.EQ.22) THEN
         GOTO 2200
C   NUCL-MOD
      ELSEIF (ISW.EQ.23) THEN
         GOTO 2300
C   LHAPATH
      ELSEIF (ISW.EQ.24) THEN
         GOTO 2400
C   OUTFILENAM
      ELSEIF (ISW.EQ.25) THEN
         GOTO 2500
C   THETA-CUT
      ELSEIF (ISW.EQ.26) THEN
         GOTO 2600
C   PT-CUT
      ELSEIF (ISW.EQ.27) THEN
         GOTO 2700
C   POLPDF
      ELSEIF (ISW.EQ.28) THEN
         GOTO 2800
C   not used
      ELSEIF (ISW.EQ.29) THEN
         GOTO 2900
C   not used
      ELSEIF (ISW.EQ.30) THEN
         GOTO 3000
C   WEIGHTS
      ELSEIF (ISW.EQ.31) THEN
         GOTO 3100
C   not used
      ELSEIF (ISW.EQ.32) THEN
         GOTO 3200
C   not used
      ELSEIF (ISW.EQ.33) THEN
         GOTO 3300
C   not used
      ELSEIF (ISW.EQ.34) THEN
         GOTO 3400
C   not used
      ELSEIF (ISW.EQ.35) THEN
         GOTO 3500
C   INT-ONLY
      ELSEIF (ISW.EQ.36) THEN
         GOTO 3600
C   TEST-OPT
      ELSEIF (ISW.EQ.37) THEN
         GOTO 3700
C   IOUNITS
      ELSEIF (ISW.EQ.38) THEN
         GOTO 3800
C   START
      ELSEIF (ISW.EQ.39) THEN
         GOTO 3900
C   STOP
      ELSEIF (ISW.EQ.40) THEN
         GOTO 4000
      ENDIF
      GOTO 1
 4    CONTINUE
      WRITE(LUNOUT,93)
      GOTO 4000
C
 90   FORMAT(A10)
 91   FORMAT(//' *****NEXT CONTROL CARD ***** ',A10/)
 92   FORMAT(/,' UNKNOWN CODEWORD - CONTROL CARD IGNORED')
 93   FORMAT(/,' UNEXPECTED END OF INPUT - STOP ASSUMED.')
C
C***********************************************************************
C               CONTROL CARD: CODEWD = TITLE
C               DEFINES THE TITLE OF THE JOB
C***********************************************************************
 100  CONTINUE
      READ(LUNIN,190) TITLE
      WRITE(LUNOUT,191) TITLE
      GO TO 1
 190  FORMAT(A80)
 191  FORMAT(/,6X,A80,/)
C
C***********************************************************************
C               CONTROL CARD: CODEWD = EL-BEAM
C               DEFINES THE PROPERTIES OF THE ELECTRON BEAM
C
C     EELE    =  ENERGY OF THE ELECTRON BEAM
C     POLARI  =  DEGREE OF ELECTRON BEAM POLARIZATION
C     LLEPT   =  -1  ELECTRON BEAM
C             =  +1  POSITRON BEAM
C     LLEPT   =  -3  muon- BEAM
C             =  +3  muon+ BEAM
C***********************************************************************
 200  CONTINUE
      READ(LUNIN,*) EELE, POLARI, LLEPT
      WRITE(LUNOUT,'(5X,2(A,1PE12.3,5X),A,I3)')
     *        ' EELE=',EELE,'POLARI=',POLARI,'LLEPT=',LLEPT
C...27.05.2016: new option for muon scattering via this input
      IF (LLEPT.EQ.-1) THEN
        LEPIN1=11
      ELSEIF (LLEPT.EQ.1) THEN
        LEPIN1=-11
      ELSEIF (LLEPT.EQ.-3) THEN
        LEPIN1=13
      ELSEIF (LLEPT.EQ.3) THEN
        LEPIN1=-13
      ENDIF
      LEPIN=LEPIN1
      IF (LLEPT.EQ.-3) LLEPT=-1
      IF (LLEPT.EQ.3) LLEPT=1
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = PR-BEAM
C               DEFINES THE PROPERTIES OF THE PROTON BEAM
C
C     EPRO    =  ENERGY OF THE PROTON BEAM
C     HPOLAR  =  DEGREE OF PROTON BEAM POLARIZATION
C***********************************************************************
 300  CONTINUE
      READ(LUNIN,*) EPRO,HPOLAR
      WRITE(LUNOUT,'(5X,2(A,1PE12.3,5X))') 
     *        ' EPRO=',EPRO,'HPOLAR=',HPOLAR
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = KINEM-CUTS
C
C     ICUT      :  DEFINITION OF THE KINEMATICAL CUTS FOR INTEGRATION
C        ICUT=1 :      LOWER CUT IN Q**2, Y-LIMITS IGNORED
C            =2 :      LOWER CUT IN W,    Y-LIMITS IGNORED
C            =3 :      CUT IN (YMIN,YMAX)
C
C     NOTE :  FINAL SETTING OF CUTS
C             ACCORDING TO THE MOST RESTRICTIVE CONDITIONS
C             IN SUBROUTINE HSPRLG
C***********************************************************************
 400  CONTINUE
      READ(LUNIN,*) ICUT, XMIN,XMAX, YMIN,YMAX, Q2MIN,Q2MAX, WMIN
      WRITE(LUNOUT,'(5X,A/4X,I3,2X,4(1PE13.4))')
     &       ' ICUT, XMIN,        XMAX,        YMIN,        YMAX ',
     &         ICUT, XMIN,XMAX, YMIN,YMAX
      WRITE(LUNOUT,'(11X,A/9X,3(1PE13.4))')
     &             ' Q2MIN,       Q2MAX,       WMIN ',
     &               Q2MIN,Q2MAX, WMIN
      IF (XMIN.GT.XMAX) THEN
        WRITE(LUNOUT,'(5X,A)') ' ****** XMIN > XMAX: EXECUTION STOPPED'
        STOP
      ENDIF
      IF (YMIN.GT.YMAX) THEN
        WRITE(LUNOUT,'(5X,A)') ' ****** YMIN > YMAX: EXECUTION STOPPED'
        STOP
      ENDIF
      IF (Q2MIN.GT.Q2MAX) THEN
        WRITE(LUNOUT,'(5X,A)') 
     &  ' ****** Q2MIN > Q2MAX: EXECUTION STOPPED'
        STOP
      ENDIF
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = EGAM-MIN
C***********************************************************************
  500 CONTINUE
      READ(LUNIN,*) EGMIN
      WRITE(LUNOUT,'(5X,A,3X,1PE13.4)') ' EGMIN = ',EGMIN
      IF(EGMIN.GT.0D0) IOPEGM=1
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = INT-OPT-NC
C
C               DEFINES THE CONTRIBUTIONS TO THE NEUTRAL CURRENT
C               CROSS SECTION FOR WHICH THE INTEGRATED CROSS SECTION
C               HAS TO BE CALCULATED AND THE SAMPLING PROCEDURE
C               PREPARED
C
C               INC2  :  NEUTRAL CURRENT NON-RADIATIVE CONTRIBUTION
C               INC31 :  NEUTRAL CURRENT CHANNEL 1 (LEPTONIC INITIAL
C                        STATE RADIATION)
C               INC32 :  NEUTRAL CURRENT CHANNEL 2 (LEPTONIC FINAL
C                        STATE RADIATION)
C               INC33 :  NEUTRAL CURRENT CHANNEL 3 (COMPTON PART)
C               INC34 :  NEUTRAL CURRENT CHANNEL 4 (QUARKONIC RADIATION)
C               IEL2  :  ELASTIC EP NON-RADIATIVE CONTRIBUTION
C               IEL31 :  ELASTIC TAIL, CHANNEL 15 (LEPTONIC INITIAL
C                        STATE RADIATION)
C               IEL32 :  ELASTIC TAIL, CHANNEL 16 (LEPTONIC FINAL
C                        STATE RADIATION)
C               IEL33 :  ELASTIC TAIL, CHANNEL 17 (COMPTON PART)
C***********************************************************************
 600  CONTINUE
      READ(LUNIN,*) INC2,INC31,INC32,INC33,INC34,
     +              IEL2,IEL31,IEL32,IEL33
      WRITE(LUNOUT,'(5X,A,3X,I5)')
     +               ' INC2  = ' ,INC2
      WRITE(LUNOUT,'(2(5X,A,3X,I5))')
     +               ' INC31 = ',INC31,' INC32 = ',INC32
      WRITE(LUNOUT,'(2(5X,A,3X,I5))')
     +               ' INC33 = ',INC33,' INC34 = ' ,INC34
      WRITE(LUNOUT,'(5X,A,3X,I5)')
     +               ' IEL2  = ' ,IEL2
      WRITE(LUNOUT,'(2(5X,A,3X,I5))')
     +               ' IEL31 = ',IEL31,' IEL32 = ',IEL32
      WRITE(LUNOUT,'(5X,A,3X,I5)')
     +               ' IEL33 = ',IEL33
      INT2(1)=INC2
      INT2(3)=IEL2
      INT3(1)=INC31
      INT3(2)=INC32
      INT3(3)=INC33
      INT3(4)=INC34
      INT3(10)=IEL31
      INT3(11)=IEL32
      INT3(12)=IEL33
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = INT-OPT-CC
C
C               DEFINES THE CONTRIBUTIONS TO THE CHARGED CURRENT
C               CROSS SECTION FOR WHICH THE INTEGRATED CROSS SECTION
C               HAS TO BE CALCULATED AND THE SAMPLING PROCEDURE
C               PREPARED
C
C               ICC2  :  CHARGED CURRENT NON-RADIATIVE CONTRIBUTION
C               ICC31 :  CHARGED CURRENT CHANNEL 1 (KP)
C               ICC32 :  CHARGED CURRENT CHANNEL 2 (KQ)
C               ICC33 :  CHARGED CURRENT CHANNEL 3 (KQS)
C***********************************************************************
 700  CONTINUE
      READ(LUNIN,*) ICC2,ICC31,ICC32,ICC33
      WRITE(LUNOUT,'(2(5X,A,3X,I5))')
     +               ' ICC2  = ' ,ICC2,  ' ICC31 = ',ICC31
     +             , ' ICC32 = ' ,ICC32, ' ICC33 = ',ICC33
      INT2(2)=ICC2
      INT3(7)=ICC31
      INT3(8)=ICC32
      INT3(9)=ICC33
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = INT-POINTS
C
C   NUMBER OF INTEGRATION POINTS FOR VEGAS
C   DEFAULT: 1000
C***********************************************************************
 800  CONTINUE
      READ(LUNIN,*) NPOVEG
      WRITE(LUNOUT,'(5X,A,5X,4I6)') ' NPOVEG', NPOVEG
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = HYP-CUBES
C
C       NPHYP  =  NUMBER OF POINTS TO BE SAMPLED PER HYPERCUBE
C                 FOR ESTIMATION OF THE LOCAL MAXIMA
C***********************************************************************
 900  CONTINUE
      READ(LUNIN,*) NPHYP
      IF(NPHYP.LT.3) NPHYP=3
      WRITE(LUNOUT,'(5X,A/5X,4I6)') ' NPHYP', NPHYP
      NPOIN=NPHYP
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = GSW-PARAM
C
C      INPUT FOR DEFINITION OF GSW-PARAMETERS AND
C                DEFINITION OF THOSE CONTRIBUTIONS
C                      TO BE INCLUDED INTO THE ACTUAL CALCULATION
C                      OF SOFT AND VIRTUAL CORRECTIONS
C
C               COMPARE SUBROUTINE HSSETP( 0/1 = NO/YES)
C
C      LPARIN( 1) = 1 :  ELWEAK PARAMETERS SET WITH FIXED W-, Z-MASS
C                   2 :  ELWEAK PARAMETERS SET WITH FIXED G/MU, Z-MASS
C      LPARIN( 2) = 0 :  ONLY BORN TERM WITHOUT ANY CORRECTIONS
C                   1 :  BORN TERM INCLUDING CORRECTIONS
C                        ACCORDING TO THE FOLLOWING PARAMETERS
C                           ( 0/1 = NO/YES)
C      LPARIN( 3) :  SOFT PHOTON EXPONENTIATION AND HIGHER ORDER TERMS
C                    IN DELTA_R
C      LPARIN( 4) :  LEPTONIC QED CORRECTIONS (VIRT. + REAL BREMSSTR.)
C      LPARIN( 5) :  QUARKONIC QED CORRECTIONS (VIRT. + REAL BREMSSTR.)
C      LPARIN( 6) :  LEPTON-QUARK INTERFERENCE (VIRT. + REAL BREMSSTR.)
C      LPARIN( 7) :  FERMIONIC CONTRIBUTIONS TO PHOTON SELF ENERGY
C      LPARIN( 8) :  FERMIONIC CONTRIBUTION TO GAM-Z MIXING TERM
C                    OF THE SELF ENERGY
C      LPARIN( 9) :  FERMIONIC CONTRIBUTIONS TO THE Z-SELF-ENERGY
C      LPARIN(10) :  FERMIONIC CONTRIBUTIONS TO THE W-SELF-ENERGY
C      LPARIN(11) :  PURELY WEAK CONTRIBUTIONS TO THE SELF ENERGIES,
C                    VERTEX CORRECTIONS AND BOXES
C      LPARIN(12) :  Z-EXCHANGE INCLUDED
C***********************************************************************
 1000 CONTINUE
      READ(LUNIN,*) (LPARIN(I),I=1,11)
      WRITE(LUNOUT,'(5X,20I2)') (LPARIN(I),I=1,11)
      LPARIN(12)=1
C---REDEFINITION FOR INTERNAL USE
      LPAR(1)=1
      LPAR(2)=LPARIN(2)
      LPAR(3)=LPARIN(3)
      LPAR(4)=LPARIN(1)
      LPAR(7)=LPARIN(7)
      LPAR(8)=LPARIN(8)
      LPAR(9)=LPARIN(9)
      LPAR(10)=LPARIN(10)
      LPAR(11)=0
      IF(LPARIN(4).EQ.1 .OR. LPARIN(5).EQ.1 .OR. LPARIN(6).EQ.1)
     &   LPAR(11)=1
      LPAR(12)=LPARIN(4)
      LPAR(13)=LPARIN(5)
      LPAR(14)=LPARIN(6)
      LPAR(15)=LPARIN(11)
      LPAR(16)=LPARIN(11)
      LPAR(17)=LPARIN(12)
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = STRUCTFUNC
C
C               DEFINES THE PARAMETRIZATION OF PARTON DENSITIES
C               OR STRUCTURE FUNCTIONS
C               APPLIED IN THE ACTUAL CALCULATION
C
C***********************************************************************
 1100 CONTINUE
      READ(LUNIN,*) ILQMOD,ILIB,ICODE
      IPDFOP=0
      IF (ILQMOD.LE.1) IPDFOP=1
      LPAR(6)=ICODE
      WRITE(LUNOUT,'(5X,A,I7)') ' ILQMOD = ',ILQMOD
      WRITE(LUNOUT,'(5X,A,I7)') ' ILIB   = ',ILIB
      WRITE(LUNOUT,'(5X,A,I7)') ' ICODE  = ',ICODE
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = NFLAVORS
C***********************************************************************
 1200 CONTINUE
      READ(LUNIN,*) NPYMIN,NPYMAX
      IF(NPYMIN.GT.6) NPYMIN=6
      IF(NPYMAX.LT.NPYMIN) NPYMAX=NPYMIN
      IF(NPYMAX.LE.0) NPYMAX=6
      WRITE(LUNOUT,'(5X,A,I6)') ' NPYMIN = ',NPYMIN
      WRITE(LUNOUT,'(5X,A,I6)') ' NPYMAX = ',NPYMAX
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = SAM-OPT-NC
C
C        DEFINES THE  NEUTRAL CURRENT CONTRIBUTIONS REQUESTED IN
C        EVENT SAMPLING
C               ISNC2 :  NEUTRAL CURRENT NON-RADIATIVE CONTRIBUTION
C               ISNC31:  NEUTRAL CURRENT CHANNEL 1 (LEPTONIC INITIAL
C                        STATE RADIATION)
C               ISNC32:  NEUTRAL CURRENT CHANNEL 2 (LEPTONIC FINAL
C                        STATE RADIATION)
C               ISNC33:  NEUTRAL CURRENT CHANNEL 3 (COMPTON PART)
C               ISNC34:  NEUTRAL CURRENT CHANNEL 4 (QUARKONIC RADIATION)
C               ISEL2 :  ELASTIC EP NON-RADIATIVE CONTRIBUTION
C               ISEL31:  ELASTIC TAIL, CHANNEL 1 (LEPTONIC INITIAL
C                        STATE RADIATION)
C               ISEL32:  ELASTIC TAIL, CHANNEL 2 (LEPTONIC FINAL
C                        STATE RADIATION)
C               ISEL33:  ELASTIC TAIL, CHANNEL 3 (COMPTON PART)
C***********************************************************************
1300  CONTINUE
      READ(LUNIN,*) ISNC2,ISNC31,ISNC32,ISNC33,ISNC34
     +             ,ISEL2,ISEL31,ISEL32,ISEL33
      WRITE(LUNOUT,'(5X,A,3X,I5)')
     +               ' ISNC2  = ' ,ISNC2
      WRITE(LUNOUT,'(2(5X,A,3X,I5))')
     +               ' ISNC31 = ',ISNC31,' ISNC32 = ' ,ISNC32
      WRITE(LUNOUT,'(2(5X,A,3X,I5))')
     +               ' ISNC33 = ',ISNC33,' ISNC34 = ' ,ISNC34
      WRITE(LUNOUT,'(5X,A,3X,I5)')
     +               ' ISEL2  = ' ,ISEL2
      WRITE(LUNOUT,'(2(5X,A,3X,I5))')
     +               ' ISEL31 = ',ISEL31,' ISEL32 = ' ,ISEL32
      WRITE(LUNOUT,'(5X,A,3X,I5)')
     +               ' ISEL33 = ',ISEL33
      ISAM2(1)=ISNC2
      ISAM2(3)=ISEL2
      ISAM3(1)=ISNC31
      ISAM3(2)=ISNC32
      ISAM3(3)=ISNC33
      ISAM3(4)=ISNC34
      ISAM3(10)=ISEL31
      ISAM3(11)=ISEL32
      ISAM3(12)=ISEL33
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = SAM-OPT-CC
C
C        DEFINES THE  CHARGED CURRENT CONTRIBUTIONS REQUESTED IN
C        EVENT SAMPLING
C
C              ISCC2  :  CHARGED CURRENT NON-RADIATIVE CONTRIBUTION
C              ISCC31 :  CHARGED CURRENT CHANNEL 1 (KP)
C              ISCC32 :  CHARGED CURRENT CHANNEL 2 (KQ)
C              ISCC33 :  CHARGED CURRENT CHANNEL 3 (KQS)
C***********************************************************************
1400  CONTINUE
      READ(LUNIN,*) ISCC2,ISCC31,ISCC32,ISCC33
      WRITE(LUNOUT,'(2(5X,A,3X,I5))')
     +               ' ISCC2  = ' ,ISCC2,  ' ISCC31 = ',ISCC31
     +             , ' ISCC32 = ' ,ISCC32, ' ISCC33 = ',ISCC33
      ISAM2(2)=ISCC2
      ISAM3(7)=ISCC31
      ISAM3(8)=ISCC32
      ISAM3(9)=ISCC33
C---CHANNELS NOT YET DEFINED
      ISAM3(8)=0
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = RNDM-SEEDS
C
C   INPUT / OUTPUT OF ACTUAL RANDOM NUMBER SEEDS
C***********************************************************************
 1500 CONTINUE
      READ(LUNIN,*) ISDINP,ISDOUT
      WRITE(LUNOUT,'(5X,A,I6)') ' ISDINP = ',ISDINP
      WRITE(LUNOUT,'(5X,A,I6)') ' ISDOUT = ',ISDOUT
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = GSW-MASS
C
C   ELECTROWEAK MASS PARAMETERS
C***********************************************************************
 1600 CONTINUE
      READ(LUNIN,*) MW,MZ,MH,MT
      LPAR(5)=1
      WRITE(LUNOUT,'(5X,A,4F12.4)') ' MW, MZ, MH, MT = ',MW,MZ,MH,MT
      GOTO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = THMIN-QRAD
C
C        DEFINES THE ANGULAR CUTS FOR HARD QUARKONIC BREMSSTRAHLUNG
C
C              TCUTQ  :  INITIAL STATE RADIATION
C              TCUTQS :  FINAL STATE RADIATION
C
C***********************************************************************
1700  CONTINUE
      READ(LUNIN,*) TCUTQ,TCUTQS
      WRITE(LUNOUT,'(2(5X,A,F10.4,A))')
     &               ' TCUTQ =  ', TCUTQ, ' RAD',
     &               ' TCUTQS =  ', TCUTQS, ' RAD'
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = FLONG
C
C    INCLUDE THE LONGITUDINAL STRUCTURE FUNCTION (FOR IPART < 1000)
C
C***********************************************************************
 1800 CONTINUE
      READ(LUNIN,*) IFLOPT,PARL11,PARL19
      WRITE(LUNOUT,'(5X,A,I5,A,F10.4,A,F10.4)')
     &               ' IFLOPT =  ', IFLOPT,
     &               '   PARL11 =', PARL11,
     &               '   PARL19 =', PARL19
      CALL DIFLOP
      GO TO 1
 1900 CONTINUE
C
C***********************************************************************
C               CONTROL CARD: CODEWD = ALFAS
C
C    DEFINITION OF ALPHA_S IN THE CALCULATION OF THE LONGITUDINAL
C    STRUCTURE FUNCTION
C
C***********************************************************************
      READ(LUNIN,*) MST111,MST115,PAR111,PAR112
      WRITE(LUNOUT,'(5X,A,I5)') ' MST111 =  ', MST111
      WRITE(LUNOUT,'(5X,A,I5)') ' MST115 =  ', MST115
      WRITE(LUNOUT,'(5X,A,F10.4)') ' PAR111 = ', PAR111
      WRITE(LUNOUT,'(5X,A,F10.4)') ' PAR112 = ', PAR112
      CALL DIALFS
      GO TO 1
 2000 CONTINUE
      GO TO 1
 2100 CONTINUE
C
C***********************************************************************
C               CONTROL CARD: CODEWD = EP-DIPOLE
C
C    INCLUDE PARAMETRIZATION OF STEIN ET AL. FOR THE DEVIATION
C    OF THE DIPOLE FORM FACTOR FOR ELASTIC EP SCATTERING
C
C***********************************************************************
      READ(LUNIN,*) IDIPOL
      WRITE(LUNOUT,'(5X,A,I5)') ' IDIPOL =  ', IDIPOL
      GO TO 1
 2200 CONTINUE
C
C***********************************************************************
C               CONTROL CARD: CODEWD = NUCLEUS
C
C    SCATTERING OFF HEAVY NUCLEI
C
C***********************************************************************
      READ(LUNIN,*) EPRO,HNA,HNZ
      WRITE(LUNOUT,'(5X,A,1PE13.3)') ' E PER NUCLEON=',EPRO
      WRITE(LUNOUT,'(5X,A,F5.0)')    ' A-NUCLEUS=',HNA
      WRITE(LUNOUT,'(5X,A,F5.0)')    ' Z-NUCLEUS=',HNZ
      GO TO 1
 2300 CONTINUE
C
C***********************************************************************
C               CONTROL CARD: CODEWD = NUCL-MOD
C
C    MODEL FOR NUCLEAR PARTON DISTRIBUTIONS
C
C***********************************************************************
      READ(LUNIN,*) INUMOD
      WRITE(LUNOUT,'(5X,A,I5)')    ' INUMOD=',INUMOD
      GO TO 1
 2400 CONTINUE
C
C***********************************************************************
C               CONTROL CARD: CODEWD = LHAPATH
C
C    PATH NAME FOR GRID FILES OF LHAPDF LIBRARY
C
C***********************************************************************
      READ(LUNIN,2590) LHAPATHI
      LHAPATH=LHAPATHI
      WRITE(LUNOUT,'(5X,A,A)')    ' LHAPATH=',LHAPATHI
      GOTO 1
 2590 FORMAT(A232)
 2500 CONTINUE
C
C***********************************************************************
C               CONTROL CARD: CODEWD = OUTFILENAM
C
C    FILE NAME TO REDIRECT STANDARD OUTPUT VIA LUNOUT (D=6)
C
C***********************************************************************
      READ(LUNIN,2590) OUTFILENAM
      ICH=INDEX(OUTFILENAM,' ')-1
      IODEF=1
      WRITE(LUNOUT,'(5X,A,A)')    ' OUTFILENAM=',OUTFILENAM
C...Open files when reading input is finished
      GOTO 1
 2600 CONTINUE
C
C***********************************************************************
C               CONTROL CARD: CODEWD = THETA-CUT
C
C    CUT ON MINIMUM ELECTRON SCATTERING ANGLE (MASSES NEGLECTED)
C
C***********************************************************************
      READ(LUNIN,*) THEMIN,THEMAX
      WRITE(LUNOUT,'(5X,A,F12.4)') ' THETA-MIN (deg) =  ', THEMIN
      WRITE(LUNOUT,'(5X,A,F12.4)') ' THETA-MAX (deg) =  ', THEMAX
      GO TO 1
 2700 CONTINUE
C
C***********************************************************************
C               CONTROL CARD: CODEWD = PT-CUT
C
C    CUT ON MINIMUM ELECTRON SCATTERING ANGLE (MASSES NEGLECTED)
C
C***********************************************************************
      READ(LUNIN,*) PTMIN
      WRITE(LUNOUT,'(5X,A,F12.4)') ' PT-MIN =  ', PTMIN
      GO TO 1
 2800 CONTINUE
C
C***********************************************************************
C               CONTROL CARD: CODEWD = POLPDF
C
C    CODE FOR PARAMETRIZATION OF POLARIZED PARTON DISTRIBUTIONS
C
C***********************************************************************
      READ(LUNIN,*) IDPVR
      WRITE(LUNOUT,'(5X,A,I8)') ' IDPVR =  ', IDPVR
      GO TO 1
 2900 CONTINUE
      GO TO 1
 3000 CONTINUE
      GO TO 1
 3100 CONTINUE
C
C***********************************************************************
C               CONTROL CARD: CODEWD = WEIGHTS
C
C    USE EXTERNALLY DEFINED WEIGHT IN EVENT GENERATION
C
C***********************************************************************
      READ(LUNIN,*) IWEIGR
      WRITE(LUNOUT,'(5X,A,I5)') ' IWEIGS =  ', IWEIGR
C...for initialization keep
      IWEIGS=0
      GO TO 1
 3200 CONTINUE
C
C***********************************************************************
C               CONTROL CARD: CODEWD = ASYM-ONLY
C
C    CACLULATE ONLY POLAROZATION OR CHARGE ODD PARTS
C
C***********************************************************************
      READ(LUNIN,*) ICPASY
      WRITE(LUNOUT,'(5X,A,I5)') ' ICPASY =  ', ICPASY
      GO TO 1
 3300 CONTINUE
      GO TO 1
 3400 CONTINUE
      GO TO 1
 3500 CONTINUE
C
C***********************************************************************
C               CONTROL CARD: CODEWD = INT-ONLY
C
C               PERFORM ONLY INTEGRATION, NO ETIMATION OF MAXIMA
C               IOPLOT < 0: NO CALL TO HSESTM
C               IOPLOT >= 0: CALL TO HSESTM
C***********************************************************************
 3600 CONTINUE
      READ(LUNIN,*) IOPLOT
      WRITE(LUNOUT,'(5X,A,5X,I3)') ' INTOPT', IOPLOT
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = TEST-OPT
C
C               DEFINE AMOUNT OUTPUT FOR TEST RUNS
C
C     IPRINT :  DIFFERENT QUANTITY OF TEST OUTPUT FOR IPRINT .GT. 0
C
C***********************************************************************
 3700 CONTINUE
      READ(LUNIN,*) IPRINT
      WRITE(LUNOUT,'(5X,A,5X,I3)') ' IPRINT',
     *                               IPRINT
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = IOUNITS
C
C               I / O UNITS
C
C     LUNOUT: FOR STANDARD INPUT	
C     LUNRND: IN-/OUTPUT FOR RANDOM NUMBER SEED
C     LUNDAT: IN-OUTPUT OF SAMPLING INFORMATION
C
C     disabled as of version 4.6.16 onwards
C
C***********************************************************************
 3800 CONTINUE
c      LUOOLD=LUNOUT
      READ(LUNIN,*) LDUMY1,LDUMY2,LDUMY3
c      IF (LUNOUT.NE.LUOOLD) THEN
c        WRITE(LUOOLD,'(5X,A,I3,A)')
c     *  ' ******* WARNING: LOGICAL UNIT FOR STANDARD OUTPUT CHANGED TO '
c     * ,LUNOUT
c     * ,'         THE FOLLOWING OUTPUT IS WRITTEN TO THIS UNIT'
c      ENDIF
c      WRITE(LUNOUT,'(5X,A,5X,5I3)')
c     *            ' LUNOUT,LUNRND,LUNDAT'
c     *             ,LUNOUT,LUNRND,LUNDAT
c      OPEN(LUNDAT,FILE='djh.dat',STATUS='UNKNOWN',FORM='UNFORMATTED')
c      OPEN(LUNRND,FILE='djhrnd.dat',STATUS='UNKNOWN',FORM='FORMATTED')
        WRITE(LUNOUT,'(5X,A)')
     *  ' ******* WARNING: CODEWORD IOUNITS IS DISABLED '
      GO TO 1
C
C***********************************************************************
C               CONTROL CARD: CODEWD = START
C               STARTS THE SAMPLING OF EVENTS AND STANDARD HISTOGRAM
C               OUTPUT
C
C     NEVENT :  NUMBER OF EVENTS TO BE SAMPLED
C***********************************************************************
 3900 CONTINUE
      READ(LUNIN,*) NEVENT
      WRITE(LUNOUT,'(5X,A,I12,/)') ' NEVENT =',NEVENT
      INFOCA=0
      DO 3921 I=1,5
        IF(NEVENT.LE.0) ISAM2(I)=0
        IF(INT2(I).GT.0) INFOCA=1
 3921   CONTINUE
      DO 3922 I=1,15
        IF(NEVENT.LE.0) ISAM3(I)=0
        IF(INT3(I).GT.0) INFOCA=1
 3922 CONTINUE
C
C---Now we can open output files
      IF (IODEF.EQ.1) THEN
        WRITE(6,'(5X,A)')
     *   '***** All output is redirected to file '
        WRITE(6,'(5X,A,A)')    '      OUTFILENAM = '
     *        ,OUTFILENAM(1:ICH)//'_out.dat'
        WRITE(6,'(5X,A)')
     *   '***** from now on '
        LUNOUT=NextUn()
        OPEN(LUNOUT,FILE=OUTFILENAM(1:ICH)//'_out.dat',STATUS='NEW')
      ENDIF
      LUNDAT=NextUn()
      OPEN(LUNDAT,FILE=OUTFILENAM(1:ICH)//'_smp.dat',STATUS='UNKNOWN',
     *     FORM='UNFORMATTED')
      LUNRND=NextUn()
      OPEN(LUNRND,FILE=OUTFILENAM(1:ICH)//'_rnd.dat'
     *     ,STATUS='UNKNOWN', FORM='FORMATTED')
C---Could also open a separate file for test output, but we keep
C   one single file for all output for the moment
      LUNTES=LUNOUT
C---  PRINT THE TITLE ALSO TO FILE
      WRITE(LUNOUT,9)

C---Initialize random number generator
      IF(ISDINP.GT.0) THEN
        READ(LUNRND,*) UIO
        READ(LUNRND,*) CIO,CDIO,CMIO,IIO,JIO
        CALL HSRNIN(UIO,CIO,CDIO,CMIO,IIO,JIO)
      elseif (isdinp.lt.0) then
        write(lunout,'(5x,A)') ' Use ranlux with date and time as seed'
        call idate(today)   ! today(1)=day, (2)=month, (3)=year
        call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
        initseed = today(1)+10*today(2)+today(3)+now(1)+5*now(3)
        write(lunout,'(5X,A,I12)') ' SEED = ', initseed
        call rndmq (idum1,idum2,initseed,' ')
      elseif (isdinp.eq.0) then
        write(lunout,'(5x,A)') ' Use ranlux with fixed seed'
        initseed = 1213141567
        write(lunout,'(5X,A,I12)') ' SEED = ', initseed
        call rndmq (idum1,idum2,initseed,' ')
      ENDIF

C---Check asymmetry only
      IF (ICPASY.EQ.1) THEN
        IOPLOT=-1
        PONULL=0D0
        POLARI=+1.0D0
ctest        POLARI=-1.0D0
        CHNULL=1.0D0
        WRITE(LUNOUT,'(5X,A)')
     *  ' *** Calculate only polarization-dependent part '
        WRITE(LUNOUT,'(5X,A)')
     *  '     of the cross section: No event generation '
      ELSEIF (ICPASY.EQ.2) THEN
        IOPLOT=-1
        CHNULL=0.0D0
ctest        LLEPT=+1
        LLEPT=-1
        WRITE(LUNOUT,'(5X,A)')
     *  ' *** Calculate only charge-odd part '
        WRITE(LUNOUT,'(5X,A)')
     *  '     of the cross section: No event generation '
      ELSE
        PONULL=1.0D0
        CHNULL=1.0D0
      ENDIF

C---DO THE NECESSARY INITIALIZATION / PRINT PARAMETERS
      CALL HSPRLG

C---CHECK CONSISTENCY OF INPUT DATA FILE
      IF(INFOCA.EQ.1) THEN
        READ(LUNDAT,ERR=3999,END=3999,IOSTAT=IOS) FNAMET
        IF(FNAME.EQ.FNAMET) THEN
          READ(LUNDAT,ERR=3999,END=3999,IOSTAT=IOS)
     &         INT2C,INT3C,ISAM2C,ISAM3C
          CALL HSDTIN
          GOTO 3998
        ENDIF
 3999   CONTINUE
        WRITE(LUNOUT,'(/A,I3/A)')
     &  ' *** NO STANDARD HERACLES INPUT FILE ASSIGNED TO UNIT',LUNDAT,
     &  ' *** COLD START OF HERACLES ASSUMED'
        DO 3901 I=1,5
          INT2C(I)=0
          ISAM2C(I)=0
 3901   CONTINUE
        DO 3902 I=1,15
          ISAM3C(I)=0
          INT3C(I)=0
 3902   CONTINUE
C
 3998   CONTINUE
C
        INTEST=0
        DO 3903 I=1,15
          IF(INT3(I).GT.100.AND.INT3C(I).EQ.0) INTEST=1
 3903   CONTINUE
        IF(INTEST.EQ.1) THEN
          WRITE(LUNOUT,'(A/A,15I4/A,15I4/A)')
     &        ' *** INCONSISTENT INPUT DATA FOR INTEGRATION (INT3) ***',
     &        ' *** INT3(I) : ',INT3,
     &        ' *** INT3C(I): ',INT3C,
     &        ' *** EXECUTION STOPPED ***'
          STOP
        ENDIF
      ENDIF
C
C***********************************************************************
C
C               TEST RUN 
C
C***********************************************************************
c
c      call tstrrr
c      stop
c
C
C***********************************************************************
C               INTEGRATION / INITIALIZATION FOR EVENT SAMPLING
C               DETERMINATION OF GLOBAL/LOCAL MAXIMA
C***********************************************************************
      INFOSA=0
C
C---NEUTRAL CURRENT
C---BORN TERM + SOFT & VIRTUAL CORRECTIONS------------------------------
C
      IF(INT2(1).GE.1) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',CHNAME(1)
C---INTEGRATION
        INCCC=1
        CALL HSINIT(INCCC,EPSO,NBIN2,NDO2,SIG2,SIG2E,XX2)
        WRITE(LUNOUT,'(//A/5X,A/5X,1PE12.4,A,1PE12.4,A)')
     *       ' CROSS SECTION (WITH ERROR ESTIMATE) FOR THE CHANNEL:',
     *       CHNAME(1), SIG2, ' +/- ', SIG2E, '  NB'
        WRITE(LUNTES,'(5X,A,1PD10.1)') ' RELATIVE ACCURACY REQUIRED:',
     &      EPSO
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSNCG2,NCHN2,NDIM2,NPOIN,NDO2,NBIN2,
     *              T2GGMA,T2GMAX,XX2,IBIM2,NREG2N)
C---SET OPTION TO SAVE INFORMATION FROM INTEGRATION ONTO UNIT LUNDAT
        INFOSA=1
        INT2C(1)=1
      ENDIF
C
C---NEUTRAL CURRENT
C---LEPTONIC BREMSSTRAHLUNG: INITIAL STATE------------------------------
      IF(INT3(1).GT.0) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',CHNAME(6)
C---INTEGRATION
        IF(INT3(1).LT.100) THEN
          ITMX31=INT3(1)
          CALL VEGAS(HSTSK1,ACCVEG,NDIM31,NPOVEG,ITMX31,IPRVEG,IGRVEG,
     &               NDO31,IT31,SI31,SI2N31,SWGT31,SCHI31,XX31)
        ELSEIF(INT3(1).LT.200) THEN
          ITMX31=INT3(1) - 100
          CALL VEGAS1(HSTSK1,ACCVEG,NDIM31,NPOVEG,ITMX31,IPRVEG,IGRVEG,
     &                NDO31,IT31,SI31,SI2N31,SWGT31,SCHI31,XX31)
        ELSE
          ITMX31=INT3(1) - 200
          CALL VEGAS2(HSTSK1,ACCVEG,NDIM31,NPOVEG,ITMX31,IPRVEG,IGRVEG,
     &                NDO31,IT31,SI31,SI2N31,SWGT31,SCHI31,XX31)
        ENDIF
        SIG31=S1
        SIG31E=S2
C
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSTSK1,NCHN31,NDIM31,NPOIN,NDO31,NBIN31,
     *              T31GMA,T31MAX,XX31,IBIM31,NREG31)
        INFOSA=1
        INT3C(1)=INT3C(1) + ITMX31
      ENDIF
C
C---NEUTRAL CURRENT
C---LEPTONIC BREMSSTRAHLUNG: FINAL STATE -------------------------------
      IF(INT3(2).GT.0) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',CHNAME(7)
C---INTEGRATION
        IF(INT3(2).LT.100) THEN
          ITMX32=INT3(2)
          CALL VEGAS(HSTSK2,ACCVEG,NDIM32,NPOVEG,ITMX32,IPRVEG,IGRVEG,
     &               NDO32,IT32,SI32,SI2N32,SWGT32,SCHI32,XX32)
        ELSEIF(INT3(2).LT.200) THEN
          ITMX32=INT3(2) - 100
          CALL VEGAS1(HSTSK2,ACCVEG,NDIM32,NPOVEG,ITMX32,IPRVEG,IGRVEG,
     &                NDO32,IT32,SI32,SI2N32,SWGT32,SCHI32,XX32)
        ELSE
          ITMX32=INT3(2) - 200
          CALL VEGAS2(HSTSK2,ACCVEG,NDIM32,NPOVEG,ITMX32,IPRVEG,IGRVEG,
     &                NDO32,IT32,SI32,SI2N32,SWGT32,SCHI32,XX32)
        ENDIF
        SIG32=S1
        SIG32E=S2
C
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSTSK2,NCHN32,NDIM32,NPOIN,NDO32,NBIN32,
     *              T32GMA,T32MAX,XX32,IBIM32,NREG32)
        INFOSA=1
        INT3C(2)=INT3C(2) + ITMX32
      ENDIF
C
C---NEUTRAL CURRENT
C---LEPTONIC BREMSSTRAHLUNG: COMPTON CONTRIBUTION ----------------------
      IF(INT3(3).GT.0) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',CHNAME(8)
C---INTEGRATION
        IF(INT3(3).LT.100) THEN
          ITMX33=INT3(3)
          CALL VEGAS(HSK1TS,ACCVEG,NDIM33,NPOVEG,ITMX33,IPRVEG,IGRVEG,
     &               NDO33,IT33,SI33,SI2N33,SWGT33,SCHI33,XX33)
        ELSEIF(INT3(3).LT.200) THEN
          ITMX33=INT3(3) - 100
          CALL VEGAS1(HSK1TS,ACCVEG,NDIM33,NPOVEG,ITMX33,IPRVEG,IGRVEG,
     &                NDO33,IT33,SI33,SI2N33,SWGT33,SCHI33,XX33)
        ELSE
          ITMX33=INT3(3) - 200
          CALL VEGAS2(HSK1TS,ACCVEG,NDIM33,NPOVEG,ITMX33,IPRVEG,IGRVEG,
     &                NDO33,IT33,SI33,SI2N33,SWGT33,SCHI33,XX33)
        ENDIF
        SIG33=S1
        SIG33E=S2
C
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSK1TS,NCHN33,NDIM33,NPOIN,NDO33,NBIN33,
     *              T33GMA,T33MAX,XX33,IBIM33,NREG33)
        INFOSA=1
        INT3C(3)=INT3C(3) + ITMX33
      ENDIF
C
C---NEUTRAL CURRENT
C---LEPTONIC BREMSSTRAHLUNG: QUARKONIC RADIATION -----------------------
      IF(INT3(4).GT.0) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',CHNAME(9)
C---INTEGRATION
        IF(INT3(4).LT.100) THEN
          ITMX34=INT3(4)
          CALL VEGAS(HSK1K3,ACCVEG,NDIM34,NPOVEG,ITMX34,IPRVEG,IGRVEG,
     &               NDO34,IT34,SI34,SI2N34,SWGT34,SCHI34,XX34)
        ELSEIF(INT3(4).LT.200) THEN
          ITMX34=INT3(4) - 100
          CALL VEGAS1(HSK1K3,ACCVEG,NDIM31,NPOVEG,ITMX34,IPRVEG,IGRVEG,
     &                NDO34,IT34,SI34,SI2N34,SWGT34,SCHI34,XX34)
        ELSE
          ITMX34=INT3(4) - 200
          CALL VEGAS2(HSK1K3,ACCVEG,NDIM34,NPOVEG,ITMX34,IPRVEG,IGRVEG,
     &                NDO34,IT34,SI34,SI2N34,SWGT34,SCHI34,XX34)
        ENDIF
        SIG34=S1
        SIG34E=S2
C
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSK1K3,NCHN34,NDIM34,NPOIN,NDO34,NBIN34,
     *              T34GMA,T34MAX,XX34,IBIM34,NREG34)
        INFOSA=1
        INT3C(4)=INT3C(4) + ITMX34
      ENDIF
C
C-----------------------------------------------------------------------
C
C---CHARGED CURRENT
C---BORN TERM + SOFT & VIRTUAL CORRECTIONS------------------------------
C
      IF(INT2(2).GE.1) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',CHNAME(2)
C---INTEGRATION
        INCCC=2
        CALL HSINIT(INCCC,EPSO,NBIN2,NDO2C,SIG2C,SIG2EC,XX2C)
        WRITE(LUNOUT,'(//A/5X,A/5X,1PE12.4,A,1PE12.4,A)')
     *       ' CROSS SECTION (WITH ERROR ESTIMATE) FOR THE CHANNEL:',
     *       CHNAME(2), SIG2C, ' +/- ', SIG2EC, '  NB'
        WRITE(LUNTES,'(5X,A,1PD10.1)') ' RELATIVE ACCURACY REQUIRED:',
     &      EPSO
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSCCG2,NCHC2,NDIM2,NPOIN,NDO2C,NBIN2,
     *              T2GMAC,T2MAXC,XX2C,IBIM2C,NREG2C)
C---SET OPTION TO SAVE INFORMATION FROM INTEGRATION ONTO UNIT LUNDAT
        INFOSA=1
        INT2C(2)=1
      ENDIF
C
C---CHARGED CURRENT
C---PHOTON BREMSSTRAHLUNG: KP - CHANNEL-------------------------------
      IF(INT3(7).GT.0) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',
     &CHNAME(12)
C---INTEGRATION
        IF(INT3(7).LT.100) THEN
          ITM31C=INT3(7)
          CALL VEGAS(HSCCKL,ACCVEG,NDM3CC,NPOVEG,ITM31C,IPRVEG,IGRVEG,
     &               NDO31C,IT31C,SI31C,S2N31C,SWT31C,SCH31C,XX31C)
        ELSEIF(INT3(7).LT.200) THEN
          ITM31C=INT3(7) - 100
          CALL VEGAS1(HSCCKL,ACCVEG,NDM3CC,NPOVEG,ITM31C,IPRVEG,IGRVEG,
     &                NDO31C,IT31C,SI31C,S2N31C,SWT31C,SCH31C,XX31C)
        ELSE
          ITM31C=INT3(7) - 200
          CALL VEGAS2(HSCCKL,ACCVEG,NDM3CC,NPOVEG,ITM31C,IPRVEG,IGRVEG,
     &                NDO31C,IT31C,SI31C,S2N31C,SWT31C,SCH31C,XX31C)
        ENDIF
        SIG31C=S1
        SG31EC=S2
C
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSCCKL,NCHC31,NDM3CC,NPOIN,NDO31C,NBN31C,
     *              T31GMC,T31MXC,XX31C,IBM31C,NRG31C)
        INFOSA=1
        INT3C(7)=INT3C(7) + ITM31C
      ENDIF
C
C---CHARGED CURRENT
C---PHOTON BREMSSTRAHLUNG: KQ - CHANNEL-------------------------------
      IF(INT3(8).GT.0) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',
     &CHNAME(13)
C---INTEGRATION
        IF(INT3(8).LT.100) THEN
          ITM32C=INT3(8)
          CALL VEGAS(HSCCQI,ACCVEG,NDM3CC,NPOVEG,ITM32C,IPRVEG,IGRVEG,
     &               NDO32C,IT32C,SI32C,S2N32C,SWT32C,SCH32C,XX32C)
        ELSEIF(INT3(8).LT.200) THEN
          ITM32C=INT3(8) - 100
          CALL VEGAS1(HSCCQI,ACCVEG,NDM3CC,NPOVEG,ITM32C,IPRVEG,IGRVEG,
     &                NDO32C,IT32C,SI32C,S2N32C,SWT32C,SCH32C,XX32C)
        ELSE
          ITM32C=INT3(8) - 200
          CALL VEGAS2(HSCCQI,ACCVEG,NDM3CC,NPOVEG,ITM32C,IPRVEG,IGRVEG,
     &                NDO32C,IT32C,SI32C,S2N32C,SWT32C,SCH32C,XX32C)
        ENDIF
        SIG32C=S1
        SG32EC=S2
C
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSCCQI,NCHC32,NDM3CC,NPOIN,NDO32C,NBN32C,
     *              T32GMC,T32MXC,XX32C,IBM32C,NRG32C)
        INFOSA=1
        INT3C(8)=INT3C(8) + ITM32C
      ENDIF
C
C---CHARGED CURRENT
C---PHOTON BREMSSTRAHLUNG: KQS - CHANNEL------------------------------
      IF(INT3(9).GT.0) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',
     &CHNAME(14)
C---INTEGRATION
        IF(INT3(9).LT.100) THEN
          ITM33C=INT3(9)
          CALL VEGAS(HSCCQF,ACCVEG,NDM3CC,NPOVEG,ITM33C,IPRVEG,IGRVEG,
     &               NDO33C,IT33C,SI33C,S2N33C,SWT33C,SCH33C,XX33C)
        ELSEIF(INT3(9).LT.200) THEN
          ITM33C=INT3(9) - 100
          CALL VEGAS1(HSCCQF,ACCVEG,NDM3CC,NPOVEG,ITM33C,IPRVEG,IGRVEG,
     &                NDO33C,IT33C,SI33C,S2N33C,SWT33C,SCH33C,XX33C)
        ELSE
          ITM33C=INT3(9) - 200
          CALL VEGAS2(HSCCQF,ACCVEG,NDM3CC,NPOVEG,ITM33C,IPRVEG,IGRVEG,
     &                NDO33C,IT33C,SI33C,S2N33C,SWT33C,SCH33C,XX33C)
        ENDIF
        SIG33C=S1
        SG33EC=S2
C
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSCCQF,NCHC33,NDM3CC,NPOIN,NDO33C,NBN33C,
     *              T33GMC,T33MXC,XX33C,IBM33C,NRG33C)
        INFOSA=1
        INT3C(9)=INT3C(9) + ITM33C
      ENDIF
C
C-----------------------------------------------------------------------
C
C---ELASTIC EP SCATTERING
C---BORN TERM + SOFT & VIRTUAL CORRECTIONS------------------------------
C
      IF(INT2(3).GE.1) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',CHNAME(3)
C---INTEGRATION
        CALL HSINIL(HSELG1,EPSO,NBIN2,NDO2E,SIG2L,SIG2EE,XX2E)
        WRITE(LUNOUT,'(//A/5X,A/5X,1PE12.4,A,1PE12.4,A)')
     *       ' CROSS SECTION (WITH ERROR ESTIMATE) FOR THE CHANNEL:',
     *       CHNAME(3), SIG2L, ' +/- ', SIG2EE, '  NB'
        WRITE(LUNTES,'(5X,A,1PD10.1)') ' RELATIVE ACCURACY REQUIRED:',
     &      EPSO
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSELG2,NCHE2,NDIM1,NPOIN,NDO2E,NBIN2,
     *              T2GMAE,T2MAXE,XX2E,IBIM2E,NREG2E)
C---SET OPTION TO SAVE INFORMATION FROM INTEGRATION ONTO UNIT LUNDAT
        INFOSA=1
        INT2C(3)=1
      ENDIF
C
C---ELASTIC TAIL
C---PHOTON BREMSSTRAHLUNG: INITIAL STATE RADIATION -------------------
      IF(INT3(10).GT.0) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',
     &CHNAME(15)
C---INTEGRATION
        IF(INT3(10).LT.100) THEN
          ITM31E=INT3(10)
          CALL VEGAS(HSELK1,ACCVEG,NDM3EL,NPOVEG,ITM31E,IPRVEG,IGRVEG,
     &               NDO31E,IT31E,SI31E,S2N31E,SWT31E,SCH31E,XX31E)
        ELSEIF(INT3(10).LT.200) THEN
          ITM31E=INT3(10) - 100
          CALL VEGAS1(HSELK1,ACCVEG,NDM3EL,NPOVEG,ITM31E,IPRVEG,IGRVEG,
     &                NDO31E,IT31E,SI31E,S2N31E,SWT31E,SCH31E,XX31E)
        ELSE
          ITM31E=INT3(10) - 200
          CALL VEGAS2(HSELK1,ACCVEG,NDM3EL,NPOVEG,ITM31E,IPRVEG,IGRVEG,
     &                NDO31E,IT31E,SI31E,S2N31E,SWT31E,SCH31E,XX31E)
        ENDIF
        SIG31L=S1
        SG31EE=S2
C
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSELK1,NCHE31,NDM3EL,NPOIN,NDO31E,NBN31E,
     *              T31GME,T31MXE,XX31E,IBM31E,NRG31E)
        INFOSA=1
        INT3C(10)=INT3C(10) + ITM31E
      ENDIF
C
C---ELASTIC TAIL
C---PHOTON BREMSSTRAHLUNG: FINAL STATE RADIATOIN ---------------------
      IF(INT3(11).GT.0) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',
     &CHNAME(16)
C---INTEGRATION
        IF(INT3(11).LT.100) THEN
          ITM32E=INT3(11)
          CALL VEGAS(HSELK2,ACCVEG,NDM3EL,NPOVEG,ITM32E,IPRVEG,IGRVEG,
     &               NDO32E,IT32E,SI32E,S2N32E,SWT32E,SCH32E,XX32E)
        ELSEIF(INT3(11).LT.200) THEN
          ITM32E=INT3(11) - 100
          CALL VEGAS1(HSELK2,ACCVEG,NDM3EL,NPOVEG,ITM32E,IPRVEG,IGRVEG,
     &                NDO32E,IT32E,SI32E,S2N32E,SWT32E,SCH32E,XX32E)
        ELSE
          ITM32E=INT3(11) - 200
          CALL VEGAS2(HSELK2,ACCVEG,NDM3EL,NPOVEG,ITM32E,IPRVEG,IGRVEG,
     &                NDO32E,IT32E,SI32E,S2N32E,SWT32E,SCH32E,XX32E)
        ENDIF
        SIG32L=S1
        SG32EE=S2
C
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSELK2,NCHE32,NDM3EL,NPOIN,NDO32E,NBN32E,
     *              T32GME,T32MXE,XX32E,IBM32E,NRG32E)
        INFOSA=1
        INT3C(11)=INT3C(11)+ITM32E
      ENDIF
C
C---ELASTIC TAIL
C---PHOTON BREMSSTRAHLUNG: COMPTON PART ------------------------------
      IF(INT3(12).GT.0) THEN
        WRITE(LUNOUT,'(///,5X,2A)') ' START INTEGRATION FOR ',
     &CHNAME(17)
C---INTEGRATION
        IF(INT3(12).LT.100) THEN
          ITM33E=INT3(12)
          CALL VEGAS(HSELCO,ACCVEG,NDM3EL,NPOVEG,ITM33E,IPRVEG,IGRVEG,
     &               NDO33E,IT33E,SI33E,S2N33E,SWT33E,SCH33E,XX33E)
        ELSEIF(INT3(12).LT.200) THEN
          ITM33E=INT3(12) - 100
          CALL VEGAS1(HSELCO,ACCVEG,NDM3EL,NPOVEG,ITM33E,IPRVEG,IGRVEG,
     &                NDO33E,IT33E,SI33E,S2N33E,SWT33E,SCH33E,XX33E)
        ELSE
          ITM33E=INT3(12) - 200
          CALL VEGAS2(HSELCO,ACCVEG,NDM3EL,NPOVEG,ITM33E,IPRVEG,IGRVEG,
     &                NDO33E,IT33E,SI33E,S2N33E,SWT33E,SCH33E,XX33E)
        ENDIF
        SIG33L=S1
        SG33EE=S2
C
C---LOCAL/GLOBAL MAXIMA FOR THE MODIFIED SAMPLING FUNCTION
        CALL HSESTM(HSELCO,NCHE33,NDM3EL,NPOIN,NDO33E,NBN33E,
     *              T33GME,T33MXE,XX33E,IBM33E,NRG33E)
        INFOSA=1
        INT3C(12)=INT3C(12) + ITM33E
      ENDIF
C--------------------------------------------------------------------
      IF(INFOSA.EQ.1) THEN
        WRITE(LUNTES,'(///2A/)')
     *       ' ACTUAL CROSS SECTION VALUES (WITH ERROR ESTIMATES)',
     *       ' AFTER INTEGRATION'
        WRITE(LUNTES,'(/A/)')
     *       ' ** NEUTRAL CURRENT: '
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' VIRTUAL & SOFT CONTRIBUTIONS  SIG2   ',
     *       SIG2, ' +/- ', SIG2E, '  NB'
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' INITIAL STATE LEPTONIC RADIATION     ',
     *       SIG31, ' +/- ', SIG31E, '  NB'
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' FINAL STATE LEPTONIC RADIATION       ',
     *       SIG32, ' +/- ', SIG32E, '  NB'
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' COMPTON CONTRIBUTION                 ',
     *       SIG33, ' +/- ', SIG33E, '  NB'
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' QUARKONIC RADIATION                  ',
     *       SIG34, ' +/- ', SIG34E, '  NB'
        WRITE(LUNTES,'(/A/)')
     *       ' ** CHARGED CURRENT: '
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' VIRTUAL & SOFT CONTRIBUTIONS  SIG2C  ',
     *       SIG2C, ' +/- ', SIG2EC, '  NB'
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' LEPTONIC RADIATION (1/k.p)   SIG31C  ',
     *       SIG31C, ' +/- ', SG31EC, '  NB'
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' QUARKONIC RADIATION (1/k.q)  SIG32C  ',
     *       SIG32C, ' +/- ', SG32EC, '  NB'
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' LEPTONIC RADIATION (1/k.qs)  SIG33C  ',
     *       SIG33C, ' +/- ', SG33EC, '  NB'
        WRITE(LUNTES,'(/A/)')
     *       ' ** ELASTIC TAIL: '
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' ELASTIC EP WITH CORRECTIONS          ',
     *       SIG2L, ' +/- ', SIG2EE, '  NB'
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' INITIAL STATE LEPTONIC RADIATION     ',
     *       SIG31L, ' +/- ', SG31EE, '  NB'
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' FINAL STATE LEPTONIC RADIATION       ',
     *       SIG32L, ' +/- ', SG32EE, '  NB'
        WRITE(LUNTES,'(3X,A,1PE12.4,A,1PE12.4,A)')
     *     ' COMPTON CONTRIBUTION                 ',
     *       SIG33L, ' +/- ', SG33EE, '  NB'
        SIGTW=SIG2+SIG31+SIG32+SIG33+SIG34
     *       +SIG2C+SIG31C+SIG32C+SIG33C
     *       +SIG2L+SIG31L+SIG32L+SIG33L
        SIGEW=DSQRT(SIG2E**2+SIG31E**2+SIG32E**2+SIG33E**2+SIG34E**2
     *       +SIG2EC**2+SG31EC**2+SG32EC**2+SG33EC**2
     *       +SIG2EE**2+SG31EE**2+SG32EE**2+SG33EE**2)
        WRITE(LUNTES,'(/A,1PE12.4,A,1PE12.4,A)')
     *     ' ** TOTAL CROSS SECTION:                 ',
     *       SIGTW, ' +/- ', SIGEW, '  NB'
        IMOD=1
        CALL HSCPDF(IMOD)
      ENDIF
C
C***********************************************************************
C   SAVE INFORMATION FROM INTEGRATION
C***********************************************************************
C
      IF(INFOSA.EQ.1) THEN
        REWIND (LUNDAT)
        WRITE(LUNDAT,ERR=3996,IOSTAT=IOS) FNAME
        WRITE(LUNDAT,ERR=3996,IOSTAT=IOS)
     &       INT2C,INT3C,ISAM2C,ISAM3C
        CALL HSDOUT
      ENDIF
      INFOSA=0
C
C***********************************************************************
C   INITIALIZATION FOR SAMPLING
C***********************************************************************

C---  Check that initialization is complete
      if (ioplot.lt.0) then
         write(lunout,'(/A,I3/A/A)')
     &   ' *** Event sampling not possible: INT-ONLY = ',ioplot,
     &   ' *** Only integration was requested ',
     &   ' *** Execution stopped ***'
         stop
      endif
      
      DO 3904 I=1,5
        IF(ISAM2(I).GT.0) INFOSA=1
        IF(ISAM2C(I).LT.ISAM2(I)) ISAM2C(I)=ISAM2(I)
 3904 CONTINUE
      DO 3905 I=1,15
        IF(ISAM3(I).GT.0) INFOSA=1
        IF(ISAM3C(I).LT.ISAM3(I)) ISAM3C(I)=ISAM3(I)
 3905 CONTINUE
      IF(INFOSA.EQ.1) THEN
C
C---CHECK CONSISTENCY OF INPUT DATA FILE
C---READ SAMPLING INFORMATION FOR ALL CONTRIBUTIONS
        REWIND (LUNDAT)
        READ(LUNDAT,ERR=3994,END=3994,IOSTAT=IOS) FNAMET
        IF(FNAME.EQ.FNAMET) THEN
          READ(LUNDAT,ERR=3994,END=3994,IOSTAT=IOS)
     &         INT2C,INT3C,ISAM2C,ISAM3C
          CALL HSDTIN
          GOTO 3993
        ENDIF
 3994   CONTINUE
        WRITE(LUNOUT,'(/A,I3/A)')
     &  ' *** NO STANDARD HERACLES INPUT FILE ASSIGNED TO UNIT',LUNDAT,
     &  ' *** EXECUTION STOPPED ***'
        STOP
C
 3993   CONTINUE
        INFOSA=0
        DO 3906 I=1,5
          IF(ISAM2(I).GT.0.AND.INT2C(I).EQ.0) INFOSA=1
 3906   CONTINUE
        DO 3907 I=1,15
          IF(ISAM3(I).GT.0.AND.INT3C(I).EQ.0) INFOSA=1
 3907   CONTINUE
        IF(INFOSA.EQ.1) THEN
          WRITE(LUNOUT,'(A/A/2(A,5I4/),2(A,15I4/),A)')
     &      ' *** INCONSISTENT INPUT DATA FOR SAMPLING (INT/ISAM) ***',
     &      ' *** ANY CHANNEL(S) REQUESTED NOT YET INTEGRATED     ***',
     &      ' *** INT2C(I): ',INT2C,
     &      ' *** ISAM2(I): ',ISAM2,
     &      ' *** INT3C(I): ',INT3C,
     &      ' *** ISAM3(I): ',ISAM3,
     &      ' *** EXECUTION STOPPED ***'
          STOP
        ENDIF
C
C---DETERMINE TOTAL CROSS SECTIONS AND ERRORS
        IF(ISAM2(1).GT.0) THEN
          SIGG(1)=SIG2
          SIGGRR(1)=SIG2E
        ENDIF
        IF(ISAM2(2).GT.0) THEN
          SIGG(2)=SIG2C
          SIGGRR(2)=SIG2EC
        ENDIF
        IF(ISAM2(3).GT.0) THEN
          SIGG(3)=SIG2L
          SIGGRR(3)=SIG2EE
        ENDIF
        IF(ISAM3(1).GT.0) THEN
          SIGG(6)=SIG31
          SIGGRR(6)=SIG31E
        ENDIF
        IF(ISAM3(2).GT.0) THEN
          SIGG(7)=SIG32
          SIGGRR(7)=SIG32E
        ENDIF
        IF(ISAM3(3).GT.0) THEN
          SIGG(8)=SIG33
          SIGGRR(8)=SIG33E
        ENDIF
        IF(ISAM3(4).GT.0) THEN
          SIGG(9)=SIG34
          SIGGRR(9)=SIG34E
        ENDIF
        IF(ISAM3(7).GT.0) THEN
          SIGG(12)=SIG31C
          SIGGRR(12)=SG31EC
        ENDIF
        IF(ISAM3(8).GT.0) THEN
          SIGG(13)=SIG32C
          SIGGRR(13)=SG32EC
        ENDIF
        IF(ISAM3(9).GT.0) THEN
          SIGG(14)=SIG33C
          SIGGRR(14)=SG33EC
        ENDIF
        IF(ISAM3(10).GT.0) THEN
          SIGG(15)=SIG31L
          SIGGRR(15)=SG31EE
        ENDIF
        IF(ISAM3(11).GT.0) THEN
          SIGG(16)=SIG32L
          SIGGRR(16)=SG32EE
        ENDIF
        IF(ISAM3(12).GT.0) THEN
          SIGG(17)=SIG33L
          SIGGRR(17)=SG33EE
        ENDIF
        DO 3908 I=1,20
          SIGTOT=SIGTOT + SIGG(I)
          SIGTRR=SIGTRR + SIGGRR(I)**2
 3908   CONTINUE
        SIGTRR=SQRT(SIGTRR)
C
C---READ INPUT FOR DJANGO6
        CALL DJGCHC(IHSONL)
C---INITIALIZATION OF USER ROUTINES---------
        ICALL=1
        CALL HSUSER(ICALL,0D0,0D0,0D0)
        IF (IHSONL.EQ.0) THEN
          IF (ICC32.NE.0.OR.ISCC32.NE.0) THEN
            WRITE(LUNOUT,*)' '
            WRITE(LUNOUT,*)' *** WARNING: channel CC32 not active, '
            WRITE(LUNOUT,*)'              set to ICC32=ISCC32=0'
            ICC32=0
            ISCC32=0
          ENDIF
          IF (ICC33.NE.0.OR.ISCC33.NE.0) THEN
            WRITE(LUNOUT,*)' '
            WRITE(LUNOUT,*)' *** WARNING: channel CC33 not active, '
            WRITE(LUNOUT,*)'              set to ICC33=ISCC33=0'
            ICC33=0
            ISCC33=0
          ENDIF
          IF (ILQMOD.GT.1) THEN
            WRITE(LUNOUT,'(A/A/A)')
     &      ' *** INCONSISTENT INPUT DATA: ILQMOD and FRAG ***',
     &      ' *** ILQMOD > 1 not allowed for FRAG .ne. -1  ***',
     &      ' *** EXECUTION STOPPED                        ***'
            STOP
          ENDIF
          INTER=4
          IF (ISAM2(2).GT.0.OR.ISAM3(7).GT.0
     &        .OR. ISAM3(8).GT.0 .OR. ISAM3(9).GT.0) THEN
            INTER=2
            IF (ISAM2(1).GT.0.OR.ISAM3(1).GT.0
     &          .OR. ISAM3(2).GT.0 .OR. ISAM3(3).GT.0
     &          .OR. ISAM3(4).GT.0) THEN
              INTER=4
              WRITE(LUNOUT,'(//10X,A/10X,2A/10X,2A//)')
     &      ' ******* WARNING: DJANGO INITIALISATION IS FOR NC EVENTS',
     &      ' *******          NC & CC NOT POSSIBLE AT THE SAME TIME',
     &      ' IN DJANGO/LEPTO',
     &      ' ******* HADRONIZATION WILL BE TREATED INCORRECTLY',
     &      ' FOR CC EVENTS'
            ENDIF
          ENDIF
          CALL DJGINIT(LEPIN,PELE,-PPRO,INTER)
        ENDIF
C
C***********************************************************************
C   EVENT GENERATION
C***********************************************************************
C
C---CONTROL OF EVENT GENERATION IN HSEVTG
        IWEIGS=IWEIGR
        NREJCW=0
        CALL HSEVTG
C
C---FINAL CALL OF USER TO GENERATE USER MONITORED OUTPUT
Chs-Correct total cross sections for failed hadronization 
        NEVTRN=NPASDJ+NSOPH+NQELAS+NREJCW
        PARL(24)=REAL(SIGTOT)
        IF (IHSONL.EQ.0) 
     &    PARL(24)=REAL(SIGTOT*(1D0-DFLOAT(NREJCW)/NEVTRN))
        ICALL=3
        CALL HSUSER(ICALL,0D0,0D0,0D0)
        IMOD=2
        CALL HSCPDF(IMOD)
C
C---SAVE INFORMATION FOR FURTHER SAMPLING
        REWIND (LUNDAT)
        WRITE(LUNDAT,ERR=3996,IOSTAT=IOS) FNAME
        WRITE(LUNDAT,ERR=3996,IOSTAT=IOS)
     &       INT2C,INT3C,ISAM2C,ISAM3C
        CALL HSDOUT
      ENDIF
C
C---SAVE CURRENT RANDOM NUMBER STATUS IF REQUESTED
      IF(ISDOUT.GT.0) THEN
        REWIND LUNRND
        CALL HSRNOU(UIO,CIO,CDIO,CMIO,IIO,JIO)
        WRITE(LUNRND,*) UIO
        WRITE(LUNRND,*) CIO,CDIO,CMIO,IIO,JIO
        WRITE(LUNOUT,'(/A,I2)')
     *  ' *** ACTUAL RANDOM NUMBER SEEDS WRITTEN TO UNIT LUNRND=',LUNRND
      ENDIF
      STOP
C
C---ERROR WHILE WRITING HERACLES DATA FILE
 3996 CONTINUE
      WRITE(LUNOUT,'(A,I4/A)')
     &   ' *** ERROR WRITING HERACLES DATA FILE ON UNIT',LUNDAT,
     &   ' *** EXECUTION STOPPED ***'
      STOP
C
C***********************************************************************
C               CONTROL CARD: CODEWD = STOP
C
C  STOPS THE EXECUTION OF THE PROGRAM
C***********************************************************************
 4000 CONTINUE
      STOP
C
C***********************************************************************
C               ERROR END
C***********************************************************************
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      BLOCK DATA HSBLKD
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /VGASIO/ NINP,NOUTP
      COMMON /HSRDIO/ ISDINP,ISDOUT
      COMMON /HSVGLP/ NPOVEG,NUMINT,NPHYP
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSIRCX/ XIRDEL
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSISGM/ TCUTQ,TCUTQS
      COMMON /HSWGTC/ IWEIGS
      COMMON /HSASY0/ PONULL,CHNULL,ICPASY
      COMMON /HSONLY/ IHSONL
C---------------------------------------------------------------------
      PARAMETER(NDIM2=2,NBIN2=50)
      PARAMETER(NREG2N=2500)
C---                     NREG2N=NBIN2**NDIM2
      LOGICAL LGLO2,LLOC2
      double precision ntot2
      COMMON /HSSNC2/ SIG2,SIG2E,T2GGMA,T2GMAX(NREG2N),
     +                XX2(50,2),
     +                FFGO2,DNCG2,FFLO2,DNCL2,GOLD2,
     +                NTOT2,NM2(NREG2N),NDO2,
     +                NCAL2,NCA12,NCA22,IBIM2,JCOR2,
     +                LGLO2,LLOC2
      PARAMETER(NREG2C=2500)
C---                     NREG2C=NBIN2**NDIM2
      LOGICAL LGLO2C,LLOC2C
      double precision ntot2c
      COMMON /HSSCC2/ SIG2C,SIG2EC,T2GMAC,T2MAXC(NREG2C),
     +                XX2C(50,2),
     +                FFGO2C,DNCG2C,FFLO2C,DNCL2C,GOLD2C,
     +                NTOT2C,NM2C(NREG2C),NDO2C,
     +                NCAL2C,NCA12C,NCA22C,IBIM2C,JCOR2C,
     +                LGLO2C,LLOC2C
      PARAMETER(NREG2E=50)
      PARAMETER(NDIM1=1)
C---                     NREG2E=NBIN2**NDIM1
      LOGICAL LGLO2E,LLOC2E
      double precision ntot2e
      COMMON /HSSEL2/ SIG2L,SIG2EE,T2GMAE,T2MAXE(NREG2E),
     +                XX2E(50,1),
     +                FFGO2E,DNCG2E,FFLO2E,DNCL2E,GOLD2E,
     +                NTOT2E,NM2E(NREG2E),NDO2E,
     +                NCAL2E,NCA12E,NCA22E,IBIM2E,JCOR2E,
     +                LGLO2E,LLOC2E
C-----------------------------
      PARAMETER(NDIM31=5,NBIN31=6,NREG31=7776)
C---                     NREG31=NBIN31**NDIM31
      LOGICAL LGLO31,LLOC31
      double precision ntot31
      COMMON /HSSN31/ SIG31,SIG31E,T31GMA,T31MAX(NREG31),
     +                XX31(50,NDIM31),
     +                FFGO31,DNCG31,FFLO31,DNCL31,GOLD31,
     +                SI31,SI2N31,SWGT31,SCHI31,IT31,
     +                NM31(NREG31),NDO31,
     +                NTOT31,NCAL31,NCA131,NCA231,IBIM31,JCOR31,
     +                LGLO31,LLOC31
      PARAMETER(NDIM32=5,NBIN32=6,NREG32=7776)
C---                     NREG32=NBIN32**NDIM32
      LOGICAL LGLO32,LLOC32
      double precision ntot32
      COMMON /HSSN32/ SIG32,SIG32E,T32GMA,T32MAX(NREG32),
     +                XX32(50,NDIM32),
     +                FFGO32,DNCG32,FFLO32,DNCL32,GOLD32,
     +                SI32,SI2N32,SWGT32,SCHI32,IT32,
     +                NM32(NREG32),NDO32,
     +                NTOT32,NCAL32,NCA132,NCA232,IBIM32,JCOR32,
     +                LGLO32,LLOC32
      PARAMETER(NDIM33=5,NBIN33=6,NREG33=7776)
C---                     NREG33=NBIN33**NDIM33
      LOGICAL LGLO33,LLOC33
      double precision ntot33
      COMMON /HSSN33/ SIG33,SIG33E,T33GMA,T33MAX(NREG33),
     +                XX33(50,NDIM33),
     +                FFGO33,DNCG33,FFLO33,DNCL33,GOLD33,
     +                SI33,SI2N33,SWGT33,SCHI33,IT33,
     +                NM33(NREG33),NDO33,
     +                NTOT33,NCAL33,NCA133,NCA233,IBIM33,JCOR33,
     +                LGLO33,LLOC33
      PARAMETER(NDIM34=5,NBIN34=6,NREG34=7776)
C---                     NREG34=NBIN34**NDIM34
      LOGICAL LGLO34,LLOC34
      double precision ntot34
      COMMON /HSSN34/ SIG34,SIG34E,T34GMA,T34MAX(NREG34),
     +                XX34(50,NDIM34),
     +                FFGO34,DNCG34,FFLO34,DNCL34,GOLD34,
     +                SI34,SI2N34,SWGT34,SCHI34,IT34,
     +                NM34(NREG34),NDO34,
     +                NTOT34,NCAL34,NCA134,NCA234,IBIM34,JCOR34,
     +                LGLO34,LLOC34
      PARAMETER(NDM3CC=5)
      PARAMETER(NBN31C=6,NRG31C=7776)
C---                     NRG31C=NBN31C**NDM3CC
      LOGICAL LGL31C,LLC31C
      double precision ntt31c
      COMMON /HSSC31/ SIG31C,SG31EC,T31GMC,T31MXC(NRG31C),
     +                XX31C(50,5),
     +                FFG31C,DNG31C,FFL31C,DNL31C,GLD31C,
     +                SI31C,S2N31C,SWT31C,SCH31C,IT31C,
     +                NM31C(NRG31C),NDO31C,
     +                NTT31C,NCL31C,NC131C,NC231C,IBM31C,JCR31C,
     +                LGL31C,LLC31C
      PARAMETER(NBN32C=6,NRG32C=7776)
C---                     NRG32C=NBN32C**NDM3CC
      LOGICAL LGL32C,LLC32C
      double precision ntt32c
      COMMON /HSSC32/ SIG32C,SG32EC,T32GMC,T32MXC(NRG32C),
     +                XX32C(50,5),
     +                FFG32C,DNG32C,FFL32C,DNL32C,GLD32C,
     +                SI32C,S2N32C,SWT32C,SCH32C,IT32C,
     +                NM32C(NRG32C),NDO32C,
     +                NTT32C,NCL32C,NC132C,NC232C,IBM32C,JCR32C,
     +                LGL32C,LLC32C
      PARAMETER(NBN33C=6,NRG33C=7776)
C---                     NRG33C=NBN33C**NDM3CC
      LOGICAL LGL33C,LLC33C
      double precision ntt33c
      COMMON /HSSC33/ SIG33C,SG33EC,T33GMC,T33MXC(NRG33C),
     +                XX33C(50,5),
     +                FFG33C,DNG33C,FFL33C,DNL33C,GLD33C,
     +                SI33C,S2N33C,SWT33C,SCH33C,IT33C,
     +                NM33C(NRG33C),NDO33C,
     +                NTT33C,NCL33C,NC133C,NC233C,IBM33C,JCR33C,
     +                LGL33C,LLC33C
      PARAMETER(NDM3EL=4)
      PARAMETER(NBN31E=10,NRG31E=10000)
C---                     NRG31E=NBN31E**NDM3EL
      LOGICAL LGL31E,LLC31E
      double precision ntt31e
      COMMON /HSSE31/ SIG31L,SG31EE,T31GME,T31MXE(NRG31E),
     +                XX31E(50,5),
     +                FFG31E,DNG31E,FFL31E,DNL31E,GLD31E,
     +                SI31E,S2N31E,SWT31E,SCH31E,IT31E,
     +                NM31E(NRG31E),NDO31E,
     +                NTT31E,NCL31E,NC131E,NC231E,IBM31E,JCR31E,
     +                LGL31E,LLC31E
      PARAMETER(NBN32E=10,NRG32E=10000)
C---                     NRG32E=NBN32E**NDM3EL
      LOGICAL LGL32E,LLC32E
      double precision ntt32e
      COMMON /HSSE32/ SIG32L,SG32EE,T32GME,T32MXE(NRG32E),
     +                XX32E(50,5),
     +                FFG32E,DNG32E,FFL32E,DNL32E,GLD32E,
     +                SI32E,S2N32E,SWT32E,SCH32E,IT32E,
     +                NM32E(NRG32E),NDO32E,
     +                NTT32E,NCL32E,NC132E,NC232E,IBM32E,JCR32E,
     +                LGL32E,LLC32E
      PARAMETER(NBN33E=10,NRG33E=10000)
C---                     NRG33E=NBN33E**NDM3EL
      LOGICAL LGL33E,LLC33E
      double precision ntt33e
      COMMON /HSSE33/ SIG33L,SG33EE,T33GME,T33MXE(NRG33E),
     +                XX33E(50,5),
     +                FFG33E,DNG33E,FFL33E,DNL33E,GLD33E,
     +                SI33E,S2N33E,SWT33E,SCH33E,IT33E,
     +                NM33E(NRG33E),NDO33E,
     +                NTT33E,NCL33E,NC133E,NC233E,IBM33E,JCR33E,
     +                LGL33E,LLC33E
C-----------------------------
      CHARACTER CHNAME*45
      COMMON /HSNAMC/ CHNAME(20)
      COMMON /HSNUME/ SIGTOT,SIGTRR,SIGG(20),SIGGRR(20),NEVENT,NEVE(20)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD,IDPVR
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSELEP/ IDIPOL
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMAX,NPYMIN
      REAL            PYSTOP,PYSLAM
      COMMON /HSALFS/ PAR111,PAR112,PARL11,PARL19,MST111,MST115
      REAL            PAR111,PAR112,PARL11,PARL19
      INTEGER         MST111,MST115
      COMMON /INTINI/ IINI
cs..Sophia
      double precision WSOPHIA
      COMMON /SOPHCT/ WSOPHIA
C
C--- STANDARD KINEMATICS
      DATA EELE, EPRO  / 27.5D0, 920D0 /
      DATA SIGTOT,SIGTRR,SIGG,SIGGRR /42*0D0/
      DATA NEVENT,NEVE               /21*0/
      DATA LUNIN,LUNOUT,LUNTES,LUNRND,LUNDAT
     *    /    5,     6,     6,    10,    11/
      DATA NINP,NOUTP /   5,    6/
      DATA ISDINP,ISDOUT / 0, 0/
      DATA NPOVEG   / 10000 /
      DATA INT2, INT3   / 20*0/
      DATA ISAM2, ISAM3 / 20*0/
      DATA IOPLOT, IPRINT /  0, 0 /
      DATA IWEIGS / 0 /
      DATA PONULL,CHNULL,ICPASY / 1.0D0, 1.0D0, 0 /
      DATA IHSONL / 1 /
      DATA DELEPS, DELTA, EGMIN, IOPEGM   / 2D-2,  0D0,  0D0, 0 /
      DATA XIRDEL /2D-2/
      DATA XMIN, XMAX, Q2MIN, Q2MAX, YMIN, YMAX, WMIN, ICUT
     *    /1D-5,  1D0,   4D0,   1D8, 1D-2,  1D0,  5D0,    3 /
      DATA THEMIN,THEMAX,CTHMIN,CTHCON,PTMIN,PTXM0
     *    /0D0   ,180D0, 1D0,   1D15,  0D0,  0D0  /
      DATA TCUTQ,TCUTQS / 0.25, 0.25/
      DATA POLARI, HPOLAR, LLEPT, LQUA / 0D0, 0D0, -1, 0/
C
C---DEFINE GSW PARAMETERS
      DATA LPARIN / 2, 1, 3, 1, 0, 0, 2, 0, 0, 0, 0, 1/
      DATA ICODE /14400/, ILIB /2/, ILQMOD /1/, IDPVR /100/
      DATA NPYMIN /1/, NPYMAX/6/
      DATA PAR111 /0.2/, PAR112/0.25/, MST111/1/, MST115/0/
      DATA PARL11 /0.01/, PARL19/0.03/
      DATA LPAR / 1, 1, 3, 2, 0, 14400, 2, 0, 0, 0,
     *            1, 1, 0, 0, 0, 0, 0, 2, 0, 0/
      DATA IPDFOP / 1 /
      DATA IFLOPT,LQCD,LTM,LHT / 4*0 /
      DATA IDIPOL / 0 /
      DATA HNA,HNZ / 1D0,1D0 /
      DATA INUMOD / 0 /
      DATA IINI / 0 /
C
      DATA WSOPHIA /1.5D0/
C
      DATA  SIG2, SIG2E, T2GGMA, T2GMAX   , NM2,     NDO2
     *    / 0D0 ,   0D0,    0D0, NREG2N*0D0, NREG2N*0,   50/
      DATA FFGO2, DNCG2, FFLO2, DNCL2, GOLD2  /5*0D0/
      DATA LGLO2,LLOC2 /2*.FALSE./
C
      DATA  SIG2C, SIG2EC, T2GMAC, T2MAXC   , NM2C   ,  NDO2C
     *    /  0D0 ,    0D0,    0D0, NREG2C*0D0, NREG2C*0,     50/
      DATA FFGO2C, DNCG2C, FFLO2C, DNCL2C, GOLD2C  /5*0D0/
      DATA LGLO2C,LLOC2C /2*.FALSE./
C
      DATA  SIG2L, SIG2EE, T2GMAE, T2MAXE   , NM2E   ,  NDO2E
     *    /  0D0 ,    0D0,    0D0, NREG2E*0D0, NREG2E*0,     50/
      DATA FFGO2E, DNCG2E, FFLO2E, DNCL2E, GOLD2E  /5*0D0/
      DATA LGLO2E,LLOC2E /2*.FALSE./
C
      DATA  SIG31, SIG31E, T31GMA, T31MAX    , NM31    , NDO31
     *     / 0D0 ,    0D0,    0D0, NREG31*0D0, NREG31*0,    50/
      DATA FFGO31,DNCG31,FFLO31,DNCL31,GOLD31 /5*0D0/
      DATA LGLO31,LLOC31 /2*.FALSE./
C
      DATA  SIG32, SIG32E, T32GMA, T32MAX    , NM32    , NDO32
     *     / 0D0 ,    0D0,    0D0, NREG32*0D0, NREG32*0,    50/
      DATA FFGO32,DNCG32,FFLO32,DNCL32,GOLD32 /5*0D0/
      DATA LGLO32,LLOC32 /2*.FALSE./
C
      DATA  SIG33, SIG33E, T33GMA, T33MAX    , NM33    , NDO33
     *     / 0D0 ,    0D0,    0D0, NREG33*0D0, NREG33*0,    50/
      DATA FFGO33,DNCG33,FFLO33,DNCL33,GOLD33 /5*0D0/
      DATA LGLO33,LLOC33 /2*.FALSE./
C
      DATA  SIG34, SIG34E, T34GMA, T34MAX    , NM34    , NDO34
     *     / 0D0 ,    0D0,    0D0, NREG34*0D0, NREG34*0,    50/
      DATA FFGO34,DNCG34,FFLO34,DNCL34,GOLD34 /5*0D0/
      DATA LGLO34,LLOC34 /2*.FALSE./
C
      DATA  SIG31C, SG31EC, T31GMC, T31MXC    , NM31C   , NDO31C
     *     /  0D0 ,    0D0,    0D0, NRG31C*0D0, NRG31C*0,     50/
      DATA FFG31C,DNG31C,FFL31C,DNL31C,GLD31C /5*0D0/
      DATA LGL31C,LLC31C /2*.FALSE./
      DATA  SIG32C, SG32EC, T32GMC, T32MXC    , NM32C   , NDO32C
     *     /  0D0 ,    0D0,    0D0, NRG32C*0D0, NRG32C*0,     50/
      DATA FFG32C,DNG32C,FFL32C,DNL32C,GLD32C /5*0D0/
      DATA LGL32C,LLC32C /2*.FALSE./
      DATA  SIG33C, SG33EC, T33GMC, T33MXC    , NM33C   , NDO33C
     *     /  0D0 ,    0D0,    0D0, NRG33C*0D0, NRG33C*0,     50/
      DATA FFG33C,DNG33C,FFL33C,DNL33C,GLD33C /5*0D0/
      DATA LGL33C,LLC33C /2*.FALSE./
C
      DATA  SIG31L, SG31EE, T31GME, T31MXE    , NM31E   , NDO31E
     *     /  0D0 ,    0D0,    0D0, NRG31E*0D0, NRG31E*0,     50/
      DATA FFG31E,DNG31E,FFL31E,DNL31E,GLD31E /5*0D0/
      DATA LGL31E,LLC31E /2*.FALSE./
      DATA  SIG32L, SG32EE, T32GME, T32MXE    , NM32E   , NDO32E
     *     /  0D0 ,    0D0,    0D0, NRG32E*0D0, NRG32E*0,     50/
      DATA FFG32E,DNG32E,FFL32E,DNL32E,GLD32E /5*0D0/
      DATA LGL32E,LLC32E /2*.FALSE./
      DATA  SIG33L, SG33EE, T33GME, T33MXE    , NM33E   , NDO33E
     *     /  0D0 ,    0D0,    0D0, NRG33E*0D0, NRG33E*0,     50/
      DATA FFG33E,DNG33E,FFL33E,DNL33E,GLD33E /5*0D0/
      DATA LGL33E,LLC33E /2*.FALSE./
      DATA CHNAME /
     &   'NEUTRAL CURRENT / ELASTIC + SOFT&VIRTUAL     ',
     &   'CHARGED CURRENT / ELASTIC + SOFT&VIRTUAL     ',
     &   'ELASTIC EP / NON-RADIATIVE + SOFT&VIRTUAL    ',
     &   'UNDEFINED                                    ',
     &   'UNDEFINED                                    ',
     &   'NEUTRAL CURRENT / LEPT. INITIAL STATE RADIAT.',
     &   'NEUTRAL CURRENT / LEPT. FINAL STATE RADIAT.  ',
     &   'NEUTRAL CURRENT / LEPT. COMPTON CONTRIBUTION ',
     &   'NEUTRAL CURRENT / QUARKONIC RADIATION        ',
     &   'UNDEFINED                                    ',
     &   'UNDEFINED                                    ',
     &   'CHARGED CURRENT / LEPT. INITIAL STATE RADIAT.',
     &   'CHARGED CURRENT / QUARK. INITIAL STATE RAD.  ',
     &   'CHARGED CURRENT / QUARK. FINAL STATE RADIAT. ',
     &   'ELASTIC TAIL / INITIAL STATE RADIATION       ',
     &   'ELASTIC TAIL / FINAL STATE RADIATION         ',
     &   'ELASTIC TAIL / COMTPON PART                  ',
     &   'UNDEFINED                                    ',
     &   'UNDEFINED                                    ',
     &   'UNDEFINED                                    '/
C
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      Function NextUn()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE HSPRLG
C---INITIALIZATION OF KINEMATICS / SETTING OF PARAMETERS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSIRCX/ XIRDEL
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSNUME/ SIGTOT,SIGTRR,SIGG(20),SIGGRR(20),NEVENT,NEVE(20)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSISGM/ TCUTQ,TCUTQS
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMPLEX(kind=8) CMW2,CMZ2
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSCBMS/ CMW2,CMZ2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSELEP/ IDIPOL
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      REAL            PYSTOP,PYSLAM
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMAX,NPYMIN
      COMMON /HSLTYP/ LEPIN1

      WRITE(LUNOUT,901)
 901  FORMAT(////,
     1'**************************************************',
     2'*****************************'/)
      WRITE(LUNOUT,'(2(20X,A/))')
     &    'FINAL DEFINITION OF RUN PARAMETERS',
     &    '           FOR HERACLES           '
      WRITE(LUNOUT,902)
 902  FORMAT(
     1'**************************************************',
     2'*****************************')
C
C---ENERGIES AND MOMENTA IN THE HERA LAB SYSTEM
C---ELECTRON BEAM
c      IF(EELE.LE.1D0) EELE=3D1

c-HS-20.02.2014:
c allow any electron energy

      IF(POLARI.LT.-1D0) POLARI=-1D0
      IF(POLARI.GT.+1D0) POLARI=+1D0
      IF(LLEPT.NE.-1.AND.LLEPT.NE.1) LLEPT=-1
C
C---ELECTRO-WEAK PARAMETERS: OPTIONS FOR VIRTUAL&SOFT CONTRIBUTION
C---INCLUDING DEFINITION OF MASSES
      CALL HSSETP
C---TRANSFER TOP-MASS TO PYSTFU
      PYSTOP=REAL(MT)
C
C---PROTON BEAM
c      IF(EPRO.LT.1D0) EPRO=82D1
      IF (EPRO.LE.MPRO) EPRO=MPRO+0.00001D0

c-HS-20.02.2014:
c allow any proton energy

      PELE=DSQRT((EELE-MEI)*(EELE+MEI))
      PPRO=DSQRT((EPRO-MPRO)*(EPRO+MPRO))
      SP=2D0*(EELE*EPRO+PELE*PPRO)+MPRO2+MEI2
      GSP=SP-MEI2-MPRO2
      IF(HPOLAR.LT.-1D0.OR.HPOLAR.GT.1D0) HPOLAR=0D0
C
C---KINEMATICAL CUTS
C   DEFINITION OF THE KINEMATICAL CUTS FOR INTEGRATION
C   ICUT=1 : CUTS IN X AND LOWER LIMIT FOR Q2
C       =2 : CUTS IN X AND LOWER LIMITS FOR Q2 AND W
C       =3 : CUTS IN X, Y, AND LOWER LIMITS FOR Q2 AND W
      ICX01=0
      ICX02=0
      ICX03=0
      ICX04=0
      IF (XMIN.LE.0D0) THEN
        XMIN=MEI2/GSP
        ICX01=1
      ENDIF
      IF (XMAX.GE.1D0) THEN
        XMAX=1D0-MEI2/GSP
        ICX02=1
      ENDIF
      IF (YMIN.LE.0D0) THEN
        YMIN=MEI2/GSP
        ICX03=1
      ENDIF
      Q2MNY=MEI2*YMIN*YMIN/(1D0-YMIN)
      IF (Q2MIN.LE.Q2MNY) THEN
        Q2MIN=Q2MNY
        ICX04=1
      ENDIF

c...for version 4.6.11: allow ICUT=4 for elastic scattering 
c   with min/max scattering angle from input
c      IF(ICUT.LT.1.OR.ICUT.GT.3) ICUT=3
      ICUTO2=0
      ICUTO3=0
      ICUTO4=0
      IF(ICUT.EQ.1) THEN
        IF(Q2MIN.GE.XMIN*GSP) THEN
          XMIN=Q2MIN/GSP
          ICUTO2=1
        ENDIF
        IF (XMIN.GT.XMAX) ICUTO2=4
      ELSEIF(ICUT.EQ.2) THEN
        IF(WMIN.LT.MPRO) WMIN=MPRO
        Q2MIN1=XMIN*(WMIN*WMIN-MPRO2)/(1D0-XMIN)
        XMAX1=1D0-(WMIN*WMIN-MPRO2)/GSP
        IF (Q2MIN1.GT.Q2MIN) THEN
          Q2MIN=Q2MIN1
          ICUTO2=3
        ENDIF
        IF (XMAX1.LT.XMAX) THEN
          XMAX=XMAX1
          ICUTO2=2
        ENDIF
        IF (XMIN.GT.XMAX) ICUTO2=4
      ELSEIF(ICUT.EQ.3) THEN
        IF(WMIN.LT.MPRO) WMIN=MPRO
        YMIN1=Q2MIN/(XMAX*GSP)
        YMIN2=(WMIN*WMIN-MPRO2)/(1D0-XMIN)/GSP
        IF(YMIN1.GT.YMIN.OR.YMIN2.GT.YMIN) ICUTO2=1
        YMIN=MAX(YMIN,YMIN1,YMIN2)
        IF(YMAX.LT.YMIN) ICUTO2=4
        XMIN1=Q2MIN/(YMAX*GSP)
        IF(XMIN1.GT.XMIN.AND.ICUTO2.LE.1) ICUTO2=1
        XMIN=MAX(XMIN,XMIN1)
        XMAX1=1D0-(WMIN*WMIN-MPRO2)/YMAX/GSP
        IF (XMAX1.LT.XMAX) ICUTO2=2
        XMAX=MIN(XMAX,XMAX1)
        IF(XMAX.LT.XMIN) ICUTO2=4
        XMAX2=Q2MAX/YMIN/GSP
        IF (XMAX2.LT.XMAX) ICUTO2=2
        XMAX=MIN(XMAX,XMAX2)
        IF(XMAX.LT.XMIN) ICUTO2=4
        Q2MIN1=XMIN*YMIN*GSP
        Q2MAX1=XMAX*YMAX*GSP
        IF (Q2MIN1.GT.Q2MIN.OR.Q2MAX1.LT.Q2MAX) ICUTO2=3
        Q2MIN=MAX(Q2MIN,Q2MIN1)
        Q2MAX=MIN(Q2MAX,Q2MAX1)
C...CUT ON ELECTRON SCATTERING ANGLE
        IF (THEMIN.NE.0D0) THEN
          CTHMIN=DCOS(THEMIN/180d0*PI)
          CTHCON=SP*(1D0+CTHMIN)/4D0/EELE/EELE/(1D0-CTHMIN)
          YMINT1=1D0/(1D0+XMAX*CTHCON)
         ELSE
          CTHMIN=1D0
          CTHCON=1D15
          YMINT1=0D0
        ENDIF
        Q2MINT=XMIN*YMINT1*GSP
        IF (YMIN.LT.YMINT1) THEN
          YMIN=YMINT1
          ICUTO3=1
        ENDIF
        IF (Q2MIN.LE.Q2MINT) THEN
          Q2MIN=Q2MINT
          ICUTO3=1
        ENDIF
C...CUT ON ELECTRON TRANSVERSE MOMENTUM PTMIN
        PTM2=PTMIN*PTMIN
        PTXM0=4D0*PTM2/SP
        YMINP1=(1D0-DSQRT(1D0-PTXM0/XMAX))/2D0
        YMAXP1=(1D0+DSQRT(1D0-PTXM0/XMAX))/2D0
        Q2MPT1=YMINP1*XMAX*GSP
        Q2MPT2=YMAXP1*XMAX*GSP
        IF (XMIN.LT.PTXM0) THEN
          XMIN=PTXM0
          ICUTO4=1
        ENDIF
        XMIN4=0D0
        IF (YMIN.GT.0.5D0) XMIN4=PTM2/SP/YMIN/(1D0-YMIN)
        IF (YMAX.LT.0.5D0) XMIN4=PTM2/SP/YMAX/(1D0-YMAX)
        IF (XMIN.LT.XMIN4) THEN
          XMIN=XMIN4
          ICUTO4=1
        ENDIF
        XMIN5=0D0
        IF (Q2MIN.GT.2D0*PTM2) XMIN5=Q2MIN*Q2MIN/SP/(Q2MIN-PTM2)
        IF (Q2MAX.GT.Q2MPT2) THEN
          Q2MAX=Q2MPT2
          ICUTO4=1
        ENDIF
        IF (Q2MAX.LT.2D0*PTM2) XMIN5=Q2MAX*Q2MAX/SP/(Q2MAX-PTM2)
        IF (XMIN.LT.XMIN5) THEN
          XMIN=XMIN5
          ICUTO4=1
        ENDIF
        IF (YMIN.LT.YMINP1) THEN
          YMIN=YMINP1
          ICUTO4=1
        ENDIF
        IF (YMAX.GT.YMAXP1) THEN
          YMAX=YMAXP1
          ICUTO4=1
        ENDIF
        IF (Q2MIN.LE.Q2MPT1) THEN
          Q2MIN=Q2MPT1
          ICUTO4=1
        ENDIF
      ENDIF
      GMIN=-1D0/Q2MIN
C
C---LOWER LIMIT IN PHOTON ENERGY (OPTIONAL)
      IF (IOPEGM.GT.0) THEN
        DO 2 I=1,5
          INT2(I)=0
          ISAM2(I)=0
  2     CONTINUE
      ENDIF
C
C---PRINT KINEMATICS / BEAM PROPERTIES
C---ELECTRON BEAM
      WRITE(LUNOUT,'(///A/)')
     *    ' *****  PROPERTIES OF THE LEPTON BEAM  *****'
      WRITE(LUNOUT,'(10X,A,I3)')
     *      ' LEPTON TYPE =  ',LEPIN1
      WRITE(LUNOUT,'(10X,A,F8.1,A)')
     *      ' ENERGY OF INCIDENT LEPTON =  ',EELE,' GEV'
      WRITE(LUNOUT,'(10X,A,F8.1,A)')
     *      ' MOMENTUM OF INCIDENT LEPTON =',PELE,' GEV/C'
      WRITE(LUNOUT,'(10X,A,F11.8,A)')
     *      ' LEPTON MASS =',MEI,' GEV/C**2'
      WRITE(LUNOUT,'(10X,A,I3)')
     *      ' CHARGE OF INCIDENT LEPTON =',LLEPT
      WRITE(LUNOUT,'(10X,A,F8.4)')
     *      ' DEGREE OF LEPTON BEAM POLARIZATION =', POLARI
C
C---PROTON BEAM
      IF (INT(HNA).EQ.1.AND.INT(HNZ).EQ.1) THEN
        WRITE(LUNOUT,'(///A/)')
     *      ' *****  PROPERTIES OF THE PROTON BEAM  *****'
        WRITE(LUNOUT,'(10X,A,F10.3,A)')
     *      ' ENERGY OF INCIDENT PROTON =  ',EPRO,' GEV'
        WRITE(LUNOUT,'(10X,A,F10.3,A)')
     *      ' MOMENTUM OF INCIDENT PROTON =',PPRO,' GEV/C'
        WRITE(LUNOUT,'(10X,A,F8.4)')
     *      ' DEGREE OF PROTON BEAM POLARIZATION =', HPOLAR
        WRITE(LUNOUT,'(10X,A,F7.4,A)')
     *      ' PROTON MASS =',MPRO,' GEV/C**2'
        WRITE(LUNOUT,'(//10X,A,1PE12.5,A)')
     *      ' CMS ENERGY SQUARED  S =',SP,' GEV**2'
      ELSE
        WRITE(LUNOUT,'(///A/)')
     *      ' *****  PROPERTIES OF THE TARGET BEAM  *****'
        WRITE(LUNOUT,'(10X,A,F4.0,A,F4.0)')
     *      ' A-NUCLEUS = ',HNA,'   Z-NUCLEUS = ',HNZ
        WRITE(LUNOUT,'(10X,A,F10.3,A)')
     *      ' ENERGY PER INCIDENT NUCLEON =  ',EPRO,' GEV'
        WRITE(LUNOUT,'(10X,A,F10.3,A)')
     *      ' MOMENTUM PER INCIDENT NUCLEON =',PPRO,' GEV/C'
        WRITE(LUNOUT,'(10X,A,F8.4)')
     *      ' DEGREE OF NUCLEON BEAM POLARIZATION =', HPOLAR
        WRITE(LUNOUT,'(10X,A,F7.4,A)')
     *      ' NUCLEON MASS =',MPRO,' GEV/C**2'
        WRITE(LUNOUT,'(//10X,A,1PE12.5,A)')
     *      ' CMS ENERGY SQUARED  S =',SP,' GEV**2'
      ENDIF
C
C---KINEMATICAL CUTS FOR GENERATED EVENTS
      WRITE(LUNOUT,'(///A/)')
     *     ' *****  KINEMATICAL LIMITS FOR GENERATED EVENTS  *****'
      WRITE(LUNOUT,'(10X,A,1PE12.5,2X,A,1PE12.5)')
     *     ' XMIN=',XMIN,' XMAX=',XMAX
      WRITE(LUNOUT,'(10X,A,1PE12.5,A)')
     *     ' Q2MIN=',Q2MIN,' GEV**2, '
      WRITE(LUNOUT,'(10X,A,1PE12.5,A)')
     *     ' Q2MAX=',Q2MAX,' GEV**2'
      WRITE(LUNOUT,'(10X,A,1PE12.5,A,17X,A)')
     *     ' WMIN=', WMIN, ' GEV',' (ACTIVE ONLY FOR ICUT=2)'
      WRITE(LUNOUT,'(10X,A,1PE12.5,2X,A,1PE12.5,1X,A)')
     *     ' YMIN=',YMIN,' YMAX=',YMAX, ' (ACTIVE ONLY FOR ICUT=3)'
      WRITE(LUNOUT,'(10X,A,I2)')
     *     ' ICUT=',ICUT
C
      IF(ICUTO2.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' NOTE: XMIN AND/OR YMIN MODIFIED FOR CONSISTENCY WITH'
        WRITE(LUNOUT,'(10X,A)')
     &     ' MINIMUM ALLOWED Q**2 OR MINIMUM ALLOWED W'
      ELSEIF(ICUTO2.EQ.2) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' NOTE: XMAX MODIFIED FOR CONSISTENCY WITH'
        WRITE(LUNOUT,'(10X,A)')
     &     ' MIN/MAX ALLOWED Q**2 AND MAX/MIN ALLOWED Y'
      ELSEIF(ICUTO2.EQ.3) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' NOTE: Q2MIN/Q2MAX MODIFIED FOR CONSISTENCY WITH'
        WRITE(LUNOUT,'(10X,A)')
     &     ' MINIMUM ALLOWED W, X, Y'
      ELSEIF(ICUTO2.EQ.4) THEN
        WRITE(LUNOUT,'(/10X,A/10X,A/10X,A)')
     &     ' NOTE: INCONSISTENT KINEMATICAL LIMITS: ',
     &     '       XMIN > XMAX AND/OR YMIN > YMAX ',
     &     '       AFTER CONSISTENCY CHECK QITH Q**2_MIN AND W_MIN. ',
     &     '       EXECUTION STOPPED'
        STOP
      ENDIF
      IF (ICUTO3.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,2A)')
     &     ' NOTE: YMIN AND/OR Q2MIN MODIFIED FOR CONSISTENCY WITH',
     &     ' CUT ON ELECTRON SCATTERING ANGLE'
      ENDIF
      IF (ICUTO4.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,2A)')
     &     ' NOTE: YMIN/MAX, XMIN AND/OR Q2MIN/MAX MODIFIED FOR',
     &     ' CONSISTENCY WITH CUT ON ELECTRON TRANSVERSE MOMENTUM'
      ENDIF
      IF (ICX01.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' XMIN MODIFIED, VALUES <= 0 NOT ALLOWED '
      ENDIF
      IF (ICX02.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' XMAX MODIFIED, VALUES >= 1 NOT ALLOWED '
      ENDIF
      IF (ICX03.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' YMIN MODIFIED, VALUES <= 0 NOT ALLOWED '
      ENDIF
      IF (ICX04.EQ.1) THEN
        WRITE(LUNOUT,'(/10X,A)')
     &     ' Q2MIN MODIFIED, VALUES <= 0 NOT ALLOWED '
      ENDIF
C
C---LOWER LIMIT IN PHOTON ENERGY (OPTIONAL)
      IF (IOPEGM.GT.0) THEN
        WRITE(LUNOUT,'(/10X,A,1PE10.3,A)')
     &       ' MINIMUM PHOTON ENERGY REQUESTED: EGMIN = ',EGMIN,' GEV'
        WRITE(LUNOUT,'(10X,A)')
     &       ' (EVENT SAMPLING FOR HARD-PHOTON CONTRIBUTIONS ONLY)'
      ENDIF
C
C---PARAMETERS FOR GSW THEORY
      WRITE(LUNOUT,'(///A/)')
     *    ' *****  PARAMETERS FOR EL-WEAK THEORY  *****'
      WRITE(LUNOUT,'(10X,A,I2,10X,A/29X,A/)')
     *    ' IFIX =', LPAR(4), ' IFIX=1 : MW FIXED',
     *                        '     =2 : GF FIXED'
      WRITE(LUNOUT,212) MW,MZ,SW2
      WWIDTH=-DIMAG(CMW2)/DSQRT(DREAL(CMW2))
      ZWIDTH=-DIMAG(CMZ2)/DSQRT(DREAL(CMZ2))
      WRITE(LUNOUT,213) WWIDTH,ZWIDTH
      WRITE(LUNOUT,214) MU,MC,MD,MB,MS,MT
      WRITE(LUNOUT,215) MH
212   FORMAT(10X,' BOSON MASSES:   MW = ',F10.4,' GEV,',/
     F       10X,'                 MZ = ',F10.4,' GEV,   SW2 = ',F10.4)
213   FORMAT(10X,' BOSON WIDTHS:   GW = ',F10.4,' GEV,    GZ = ',F10.4)
214   FORMAT(10X,
     F       ' FERMION MASSES: MU = ',F10.4,' GEV,  MC = ',F10.4,' GEV',
     F /,10X,'                 MD = ',F10.4,' GEV,  MB = ',F10.4,' GEV',
     F /,10X,'                 MS = ',F10.4,' GEV,  MT = ',F10.4,' GEV')
215   FORMAT(10X,' HIGGS MASS:     MH = ',F10.4,' GEV')
C
C---PARTON DISTRIBUTION
      WRITE(LUNOUT,'(///A/A/)')
     *    ' *****  OPTIONS FOR PARTON DISTRIBUTIONS OR      *****',
     *    ' *****  STRUCTURE FUNCTIONS                      *****'
      CALL HSWPDF
C
C---ELASTIC SCATTERING
      IF (IDIPOL.NE.0) THEN
      WRITE(LUNOUT,'(///A/A/)')
     *    ' *****  ELASTIC SCATTERING INCLUDED               *****',
     *    ' *****  WITH DIPOLE FORM FACTOR FROM STEIN ET AL. *****'
      ENDIF
C
C---ANGULAR CUTS FOR QUARKONIC BREMSSTRAHLUNG
      IF(INT3(4).GE.1 .OR. ISAM3(4).GE.1 .OR.
     &   ((INT2(1).GE.1.OR.ISAM2(1).GE.1).AND.LPARIN(5).GE.1) ) THEN
        WRITE(LUNOUT,'(//A/)')
     &          ' *****  ANGULAR CUTS FOR QUARKONIC RADIATION  *****'
        WRITE(LUNOUT,'(10X,A,1PE10.3,A)')
     &    ' TCUTQ  =', TCUTQ, ' RAD',
     &    ' TCUTQS =', TCUTQS, ' RAD'
      ENDIF
C
C---OPTIONS FOR VIRTUAL&SOFT CONTRIBUTION
      IF(INT2(1).EQ.1.OR.INT2(2).EQ.1
     *   .OR.ISAM2(1).GE.1.OR.ISAM2(2).GE.1) THEN
        WRITE(LUNOUT,'(//A/)')
     *    ' *****  OPTIONS FOR VIRT&SOFT CONTRIBUTION         *****'
        WRITE(LUNOUT,216) LPAR
216     FORMAT(/,' PARAMETER LIST',/
     F,' *******************************************************',/
     F,' **      BORN CROSS SECTION:    LPAR( 1) = ',I6,'     **',/
     F,' **      1-LOOP CORRECTIONS:    LPAR( 2) = ',I6,'     **',/
     F,' **      HIGHER ORDERS:         LPAR( 3) = ',I6,'     **',/
     F,' **      MW OR GMU FIXED:       LPAR( 4) = ',I6,'     **',/
     F,' **      MASS PARAM FROM INPUT: LPAR( 5) = ',I6,'     **',/
     F,' **      STRUCTURE FUNCTIONS:   LPAR( 6) = ',I6,'     **',/
     F,' **      SIGMA-GAMMA:           LPAR( 7) = ',I6,'     **',/
     F,' **      SIGMA-GAMMA-Z:         LPAR( 8) = ',I6,'     **',/
     F,' **      SIGMA-Z:               LPAR( 9) = ',I6,'     **',/
     F,' **      SIGMA-W:               LPAR(10) = ',I6,'     **',/
     F,' **      QED CORRECTIONS:       LPAR(11) = ',I6,'     **',/
     F,' **      LEPTONIC QED:          LPAR(12) = ',I6,'     **',/
     F,' **      HADRONIC QED:          LPAR(13) = ',I6,'     **',/
     F,' **      LEPTON-QUARK-INTRF:    LPAR(14) = ',I6,'     **',/
     F,' **      WEAK CORRECTIONS:      LPAR(15) = ',I6,'     **',/
     F,' **      WEAK BOXES:            LPAR(16) = ',I6,'     **',/
     F,' **      GAMMA OR/AND Z         LPAR(17) = ',I6,'     **',/
     F,' **      NOT USED:              LPAR(18) = ',I6,'     **',/
     F,' **      NOT USED:              LPAR(19) = ',I6,'     **',/
     F,' **      NOT USED:              LPAR(20) = ',I6,'     **',/
     F,' *******************************************************'/)
      ENDIF
C
C---ELASTIC CONTRIBUTIONS COMPATIBLE WITH KINEMATIC CUTS?
      IF (XMAX.LT.(1D0-MEI2/GSP-XIRDEL)
     &   .AND.(INT2(3).NE.0.OR.ISAM2(3).NE.0))THEN
        INT2(3)=0
        ISAM2(3)=0
        WRITE(LUNOUT,'(/A,3(/A,10X,A))')
     &    ' ********** WARNING **********',
     &    ' ***',' VALUE OF XMAX DOES NOT ALLOW NON-RADIATIVE ',
     &    ' ***',' ELASTIC EP SCATTERING ',
     &    ' ***',' INT2(3) AND ISAM2(3) SET TO 0 '
      ENDIF
      IF (WMIN.GT.MPRO.AND.(INT2(3).NE.0.OR.ISAM2(3).NE.0.OR.
     &   INT3(10).NE.0.OR.INT3(11).NE.0.OR.INT3(12).NE.0.OR.
     &   ISAM3(10).NE.0.OR.ISAM3(11).NE.0.OR.ISAM3(12).NE.0)) THEN
        INT2(3)=0
        INT3(10)=0
        INT3(11)=0
        INT3(12)=0
        ISAM2(3)=0
        ISAM3(10)=0
        ISAM3(11)=0
        ISAM3(12)=0
        WRITE(LUNOUT,'(/A,3(/A,10X,A))')
     &    ' ********** WARNING **********',
     &    ' ***',' VALUE OF WMIN DOES NOT ALLOW ',
     &    ' ***',' ELASTIC AND QUASI-ELASTIC EP SCATTERING ',
     &    ' ***',' INT2(3), ISAM2(3), INT3(10-12) AND ISAM3(10-12) SET
     &TO 0 '
      ENDIF
C
C---OPTIONS FOR INTEGRATION / SAMPLING
      WRITE(LUNOUT,'(///A/)')
     *    ' *****  OPTIONS FOR INTEGRATION / SAMPLING  *****'
C---CHANNELS FOR CHARGED CURRENT NOT YET DEFINED
      IF (INT3(8).NE.0.OR.ISAM3(8).NE.0) THEN
        INT3(8)=0
        ISAM3(8)=0
        LPAR(13)=0
        WRITE(LUNOUT,'(10X,A,/,10X,A)')
     *  ' CHARGED CURRENT: QUARKONIC RADIATION NOT YET ACTIVATED ',
     *  ' ICC32, ISCC32 AND LPARIN(5) SET TO ZERO ' 
      ENDIF
      IF (INT3(9).NE.0.OR.ISAM3(9).NE.0) THEN
        INT3(9)=0
        ISAM3(9)=0
        LPAR(14)=0
        WRITE(LUNOUT,'(10X,A,/,10X,A)')
     *' CHARGED CURRENT: LEPTON-QUARK INTERFERENCE NOT YET ACTIVATED ',
     *' ICC33, ISCC33 AND LPARIN(6) SET TO ZERO ' 
      ENDIF

      DO 1 I=1,5
        IF(INT2(I).LT.0.OR.INT2(I).GT.1) INT2(I)=0
 1    CONTINUE
      WRITE(LUNOUT,'(10X,A,8I4)')
     *      ' INT2(I) ', INT2
      WRITE(LUNOUT,'(10X,A,15I4)')
     *      ' INT3(I) ', INT3
      WRITE(LUNOUT,'(10X,A,8I4)')
     *      ' ISAM2(I)',ISAM2
      WRITE(LUNOUT,'(10X,A,15I4)')
     *      ' ISAM3(I)', ISAM3
      IF(LPAR(14).GT.0) THEN
        ISM3TT=ISAM3(1)+ISAM3(2)+ISAM3(3)+ISAM3(4)
        IF ((ISM3TT.GT.0).AND.(ISAM3(1).EQ.0.OR.ISAM3(2).EQ.0.OR.
     *                         ISAM3(3).EQ.0.OR.ISAM3(4).EQ.0)   ) THEN
        WRITE(LUNOUT,'(/A,3(/A,10X,A))')
     &    ' ********** WARNING **********',
     &    ' ***',' LEPTON-QUARK INTERFERENCE CAN BE INCLUDED ONLY IF',
     &    ' ***',' AT THE SAME TIME ALL LEPTONIC AND QUARKONIC ',
     &    ' ***',' CHANNELS ARE REQUESTED,  LPARIN(6) SET TO 0 '
        LPARIN(6)=0
        LPAR(14)=0
        ENDIF
      ENDIF
C
C---NUMBER OF REQUESTED EVENTS
      IF (ISAM2(1).NE.0 .OR. ISAM2(2).NE.0.OR.ISAM2(3).NE.0
     &    .OR. ISAM3(1).GT.0 .OR. ISAM3(2).GT.0 .OR. ISAM3(3).GT.0
     &    .OR. ISAM3(4).GT.0 .OR. ISAM3(7).GT.0 .OR. ISAM3(9).GT.0
     &    .OR. ISAM3(10).GT.0 .OR. ISAM3(11).GT.0 .OR. ISAM3(12).GT.0
     &   ) THEN
        WRITE(LUNOUT,'(///A/)')
     *              ' *****  NUMBER OF EVENTS TO BE SAMPLED  *****'
        WRITE(LUNOUT,'(10X,A,I12)')
     *    ' NUMBER OF EVENTS REQUESTED  NEVENT =',NEVENT
      ENDIF
C---
      WRITE(LUNOUT,'(///)')
      RETURN
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE HSSETP
C---SETTING OF ELECTROWEAK PARAMETERS
C*****
C   MODIFIED VERSION
C   HJM 1/12/89
C   HS 7/7/91
C*****
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C
      COMPLEX(kind=8) HSSRWW,HSSRZZ
      COMPLEX(kind=8) CMW2,CMZ2
      EXTERNAL HSSRWW,HSSRZZ
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSLTYP/ LEPIN1
      COMMON /HSASY0/ PONULL,CHNULL,ICPASY
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSDELR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /HSCBMS/ CMW2,CMZ2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSKNCC/ SXNRCC,SX1NCC
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
C
C---DEFINE CONSTANTS
      PI=4D0*DATAN(1D0)
      ALPHA=1D0/137.0359895D0
      ALP1PI=ALPHA/PI
      ALP2PI=ALPHA/2D0/PI
      ALP4PI=ALPHA/4D0/PI
      E=DSQRT(4D0*PI*ALPHA)
      GF=1.166389D-5
      AGF0=PI*ALPHA/GF/DSQRT(2D0)
      SXNORM=PI*ALPHA*ALPHA/2D0*3.8938D5
      SX1NRM=ALPHA*ALPHA*ALPHA/16D0/PI*3.8938D5

C---DEFINE PARAMETERS OF THE ELECTROWEAK STANDARD MODEL
      ME=.51099906D-3
      MMY=.105658387D0
      MTAU=1.7841D0
      MU=.062D0
      MD=.083D0
      MS=.215D0
      MC=1.5D0
      MB=4.5D0
      IF (LPAR(5).NE.1) THEN
        MT=172.5D0
        MH=125.25D0
      ENDIF
      MH2=MH*MH
      ME2=ME*ME
      MMY2=MMY*MMY
      MTAU2=MTAU*MTAU
      MU2=MU*MU
      MD2=MD*MD
      MS2=MS*MS
      MC2=MC*MC
      MB2=MB*MB
      MT2=MT*MT
C
      IF (LPAR(5).NE.1) MZ=91.1876D0
      MZ2=MZ*MZ
      IF (LPAR(4).GT.1) THEN
        MW=HSPGFX()
      ELSE
        IF (LPAR(5).NE.1) MW=80.379D0
        MWKEEP=MW
        MW1=HSPGFX()
        MW=MWKEEP
      ENDIF
      MW2=MW*MW
      CW=MW/MZ
      CW2=CW*CW
      SW2=1D0-CW2
      SW=DSQRT(SW2)
      WWIDTH=DIMAG(HSSRWW(MW2))/MW
      ZWIDTH=DIMAG(HSSRZZ(MZ2))/MZ
      CMW2=MW*DCMPLX(MW,-WWIDTH)
      CMZ2=MZ*DCMPLX(MZ,-ZWIDTH)

C---NORMALIZATION OF NC AND CC CROSS SECTIONS
C---ON-MASS SHELL SCHEME
      IF (LPAR(4).EQ.1) THEN
        B=1D0/4D0/CW/SW
        SXNRCC=SXNORM/SW2/SW2
        SX1NCC=SX1NRM/SW2/SW2
        ELSE
C---MODIFIED ON-MASS SHELL SCHEME (NORMALIZATION TO G-MU)
        B=MZ/SQRT(AGF0)/4D0
        SXNRCC=SXNORM*MW2*MW2/AGF0/AGF0
        SX1NCC=SX1NRM*MW2*MW2/AGF0/AGF0
        IF (LPAR(4).EQ.3) B=B*SQRT(1D0-DELTAR)
      ENDIF
C---DEFINE FERMION GAUGE BOSON COUPLING CONSTANTS
      VAFI(2,1,1)=0D0
      VAFI(2,2,1)=0D0
      VAFI(2,3,1)=0D0
      VAFI(2,1,2)=-B
      VAFI(2,2,2)=B
      VAFI(2,3,2)=-B
      VAFI(1,1,1)=1D0
      VAFI(1,2,1)=-2D0/3D0
      VAFI(1,3,1)=1D0/3D0
      VAFI(1,1,2)=B*(4D0*SW2-1D0)
      VAFI(1,2,2)=B*(1D0-8D0*SW2/3D0)
      VAFI(1,3,2)=B*(4D0*SW2/3D0-1D0)

C----------------------------------------------------------------------
C---MASSES USED FOR HARD BREMSSTRAHLUNG KINEMATICS
      IF (ABS(LEPIN1).EQ.11) THEN
        MEI=ME
        MEF=ME
C...for muon scattering
       ELSEIF(ABS(LEPIN1).EQ.13) THEN
        MEI=MMY
        MEF=MMY
      ENDIF

C---QUARK MASSES USED AS REGULATORS IN THE LEPTONIC BREMSSTRAHLUNG
C---DO NOT USE FOR VERY SMALL X
      MPRO=938.28D-3
      MPRO2=MPRO*MPRO
      MQI=MU
      MQF=MU
      MEF2=MEF*MEF
      MQF2=MQF*MQF
      MEI2=MEI*MEI
      MQI2=MQI*MQI
C
      IGAMMA = 1
      IZ = 2
      IEL = 1
      IFU = 2
      IFD = 3
      INDV = 1
      INDA = 2
C
      DO 11 IF=IEL,IFD
        DO 12 IB1=IGAMMA,IZ
          DO 13 IB2=IGAMMA,IZ
          FLIND(INDV,IF,IB1,IB2)=
     *     2D0*(VAFI(INDV,IF,IB1)*VAFI(INDV,IF,IB2)
     *         +VAFI(INDA,IF,IB1)*VAFI(INDA,IF,IB2))
 13       CONTINUE
 12     CONTINUE
 11   CONTINUE
      DO 21 IF=IEL,IFD
        DO 22 IB1=IGAMMA,IZ
          DO 23 IB2=IGAMMA,IZ
          FLIND(INDA,IF,IB1,IB2)=
     *     2D0*(VAFI(INDV,IF,IB1)*VAFI(INDA,IF,IB2)
     *         +VAFI(INDA,IF,IB1)*VAFI(INDV,IF,IB2))
 23       CONTINUE
 22     CONTINUE
 21   CONTINUE
C
      DO 31 IVB1 = IGAMMA, IZ
       DO 32 IVB2 = IGAMMA, IZ
        DO 33 IFERM = IFU, IFD
        AFIJ(IFERM,IVB1,IVB2)=
     *  FLIND(INDV,IFERM,IVB1,IVB2)
     *  *(PONULL*FLIND(INDV,IEL,IVB1,IVB2) 
     *   - POLARI*FLIND(INDA,IEL,IVB1,IVB2))
        BFIJ(IFERM,IVB1,IVB2)=
     *  FLIND(INDA,IFERM,IVB1,IVB2)
     *  *(PONULL*FLIND(INDA,IEL,IVB1,IVB2) 
     *   - POLARI*FLIND(INDV,IEL,IVB1,IVB2))
 33     CONTINUE
 32    CONTINUE
 31   CONTINUE
C
      DO 41 IVB1 = IGAMMA, IZ
       DO 42 IVB2 = IGAMMA, IZ
        DO 43 IFERM = IFU, IFD
        AGIJ(IFERM,IVB1,IVB2)=
     *  FLIND(INDA,IFERM,IVB1,IVB2)
     *  *(PONULL*FLIND(INDV,IEL,IVB1,IVB2) 
     *    - POLARI*FLIND(INDA,IEL,IVB1,IVB2))
        BGIJ(IFERM,IVB1,IVB2)=
     *  FLIND(INDV,IFERM,IVB1,IVB2)
     *  *(PONULL*FLIND(INDA,IEL,IVB1,IVB2) 
     *    - POLARI*FLIND(INDV,IEL,IVB1,IVB2))
 43     CONTINUE
 42    CONTINUE
 41   CONTINUE
C
      RETURN
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      FUNCTION HSPGFX()
C
C   THIS PROGRAM DETERMINES MW FOR GIVEN MZ VIA MU-LIFETIME
C   ( W. HOLLIK, HAMBURG )
C   ( CHANGED AD 6.8.88)
C   ( CHANGED HS 14.5.91)
C
C   INPUT:   Z0 MASS = MZ, HIGGS MASS = MH, TOP MASS= MT
C   OUTPUT:  W MASS = MW
C   OUTPUT ON COMMONS: DELTA-R = DELTAR, WEAK MIXING ANGLE =SW2
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSRWW
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSDELR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      DIMENSION DR2(20)

C
      I=1
C---SW2 = SIN**2 THETA-W, START WITH APPROXIMATE VALUE
      SQ=DSQRT(1D0-4D0*AGF0/MZ2/(1D0-.07D0))
      SW2=(1D0-SQ)/2D0
      DR2(1)=0.07D0
      CW2=1D0-SW2
      SW=DSQRT(SW2)
      CW=DSQRT(CW2)
      MW=MZ*CW
      MW2=MW*MW

C---DELTAR MEANS THE QUANTITY 'DELTA R' IN MU LIFETIME FORMULA
51    CONTINUE
      LPR15K=LPAR(15)
      LPR7K=LPAR(7)
      LPAR(15)=1
      LPAR(7)=2
      DELTAR=DREAL(HSSRWW(0D0))/MW2
     *    +ALP4PI/SW2*(6D0+(3.5D0/SW2-2D0)*DLOG(CW2))
      LPAR(15)=LPR15K
      LPAR(7)=LPR7K
C---LEADING TWO-LOOP TERMS
      DRHOT=3D0*ALP4PI/4D0/SW2/CW2*MT2/MZ2
      DALPMZ=0.0602D0
      XGMT=GF*MT2/SQRT(2D0)/8D0/PI/PI
      ALAMB2=0.14D0**2D0
      NF=5
      ALPQCD=0D0
      DLRQCD=0D0
      IF (LPAR(3).GE.2) THEN
        ALPQCD=4D0*PI/(11D0-2D0/3D0*NF)/DLOG(    MT2/ALAMB2)
        DLRQCD=CW2/SW2*2D0*ALPQCD/PI*(PI*PI/3D0+1D0)*XGMT
      ENDIF
      DELTR2=0D0
      IF (LPAR(3).GE.1) THEN
        DELTR2=CW2/SW2*(CW2/SW2*DRHOT*DRHOT/(1D0-DALPMZ)
     *                  -3D0*XGMT*XGMT*(19D0-2D0*PI*PI) )
      ENDIF
      DELTAR=DELTAR+DELTR2+DLRQCD
      DRPIW2=DELTR2+DLRQCD
      SQ=DSQRT(1D0-4D0*AGF0/MZ2/(1D0-DELTAR))

C---THE CORRECTED VALUE FOR SIN**2 THETA-W
      SW2=(1D0-SQ)/2D0
      DR2(I+1)=DELTAR
      DS2=DABS(DR2(I+1)-DR2(I))
      CW2=1.D0-SW2
      CW=DSQRT(CW2)
      SW=DSQRT(SW2)
C---THE CORRECTED VALUE FOR THE W MASS
      MW=MZ*CW
      MW2=MW*MW
C---DELTA-RHO
      BTOP4=1D0/(1D0+DRHOT)/(1D0-3D0*XGMT*(1D0+XGMT*(19D0-2D0*PI*PI)))
      BTOP4=SQRT(BTOP4)

      IF(DS2.LT.1D-8) THEN
         HSPGFX=MW
         RETURN
      END IF
      I=I+1
      IF(I.LE.20) GOTO  51
      HSPGFX=MW
      WRITE(LUNOUT,1)
    1 FORMAT(' WARNING: CALCULATION OF MW IN HSPGFX DID NOT CONVERGE',/
     F      ,' AFTER 20 ITERATIONS. CHECK INPUT MASSES',/)
      RETURN
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE HSDELO(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      IF (IOPEGM.GT.0) THEN
        DELTA=EGMIN
        RETURN
      ELSE
        SIGMA=EPRO/EELE
        XMY=Y-(1D0-X*Y)*SIGMA
        XPY=Y+(1D0-X*Y)*SIGMA
        OMEGA=2D0*EELE*SIGMA*Y*(1D0-X)
     *          /(XPY+DSQRT(4D0*X*Y*SIGMA*(1D0-Y)+XMY*XMY))
        EQUA=X*SIGMA*EELE
        EES=EELE*(1D0-Y+X*Y*SIGMA)
        DELTA=DMIN1(EELE,EES,OMEGA)*DELEPS
      ENDIF
      END

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE HSEVTG
C---
C   EVENT SAMPLING
C   FRACTIONS FROM DIFFERENT CONTRIBUTIONS ACCORDING TO CROSS SECTIONS
C---
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      EXTERNAL HSNCG2, HSTSK1, HSTSK2, HSK1TS, HSK1K3
      EXTERNAL HSCCG2, HSCCKL, HSCCQI, HSCCQF
      EXTERNAL HSELG2, HSELK1, HSELK2, HSELCO
      PARAMETER (NCHN2=1,NCHC2=2,NCHE2=3)
      PARAMETER (NCHN31=6,NCHN32=7,NCHN33=8,NCHN34=9)
      PARAMETER (NCHC31=12,NCHC32=13,NCHC33=14)
      PARAMETER (NCHE31=15,NCHE32=16,NCHE33=16)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSNUME/ SIGTOT,SIGTRR,SIGG(20),SIGGRR(20),NEVENT,NEVE(20)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
C-----------------------------
      PARAMETER(NDIM2=2,NBIN2=50)
      PARAMETER(NREG2N=2500)
      LOGICAL LGLO2,LLOC2
      double precision ntot2
      COMMON /HSSNC2/ SIG2,SIG2E,T2GGMA,T2GMAX(NREG2N),
     +                XX2(50,2),
     +                FFGO2,DNCG2,FFLO2,DNCL2,GOLD2,
     +                NTOT2,NM2(NREG2N),NDO2,
     +                NCAL2,NCA12,NCA22,IBIM2,JCOR2,
     +                LGLO2,LLOC2
      PARAMETER(NREG2C=2500)
      LOGICAL LGLO2C,LLOC2C
      double precision ntot2c
      COMMON /HSSCC2/ SIG2C,SIG2EC,T2GMAC,T2MAXC(NREG2C),
     +                XX2C(50,2),
     +                FFGO2C,DNCG2C,FFLO2C,DNCL2C,GOLD2C,
     +                NTOT2C,NM2C(NREG2C),NDO2C,
     +                NCAL2C,NCA12C,NCA22C,IBIM2C,JCOR2C,
     +                LGLO2C,LLOC2C
      PARAMETER(NREG2E=50)
      PARAMETER(NDIM1=1)
      LOGICAL LGLO2E,LLOC2E
      double precision ntot2e
      COMMON /HSSEL2/ SIG2L,SIG2EE,T2GMAE,T2MAXE(NREG2E),
     +                XX2E(50,1),
     +                FFGO2E,DNCG2E,FFLO2E,DNCL2E,GOLD2E,
     +                NTOT2E,NM2E(NREG2E),NDO2E,
     +                NCAL2E,NCA12E,NCA22E,IBIM2E,JCOR2E,
     +                LGLO2E,LLOC2E
      PARAMETER(NDIM31=5,NBIN31=6,NREG31=7776)
      LOGICAL LGLO31,LLOC31
      double precision ntot31
      COMMON /HSSN31/ SIG31,SIG31E,T31GMA,T31MAX(NREG31),
     +                XX31(50,NDIM31),
     +                FFGO31,DNCG31,FFLO31,DNCL31,GOLD31,
     +                SI31,SI2N31,SWGT31,SCHI31,IT31,
     +                NM31(NREG31),NDO31,
     +                NTOT31,NCAL31,NCA131,NCA231,IBIM31,JCOR31,
     +                LGLO31,LLOC31
      PARAMETER(NDIM32=5,NBIN32=6,NREG32=7776)
      LOGICAL LGLO32,LLOC32
      double precision ntot32
      COMMON /HSSN32/ SIG32,SIG32E,T32GMA,T32MAX(NREG32),
     +                XX32(50,NDIM32),
     +                FFGO32,DNCG32,FFLO32,DNCL32,GOLD32,
     +                SI32,SI2N32,SWGT32,SCHI32,IT32,
     +                NM32(NREG32),NDO32,
     +                NTOT32,NCAL32,NCA132,NCA232,IBIM32,JCOR32,
     +                LGLO32,LLOC32
      PARAMETER(NDIM33=5,NBIN33=6,NREG33=7776)
      LOGICAL LGLO33,LLOC33
      double precision ntot33
      COMMON /HSSN33/ SIG33,SIG33E,T33GMA,T33MAX(NREG33),
     +                XX33(50,NDIM33),
     +                FFGO33,DNCG33,FFLO33,DNCL33,GOLD33,
     +                SI33,SI2N33,SWGT33,SCHI33,IT33,
     +                NM33(NREG33),NDO33,
     +                NTOT33,NCAL33,NCA133,NCA233,IBIM33,JCOR33,
     +                LGLO33,LLOC33
      PARAMETER(NDIM34=5,NBIN34=6,NREG34=7776)
      LOGICAL LGLO34,LLOC34
      double precision ntot34
      COMMON /HSSN34/ SIG34,SIG34E,T34GMA,T34MAX(NREG34),
     +                XX34(50,NDIM34),
     +                FFGO34,DNCG34,FFLO34,DNCL34,GOLD34,
     +                SI34,SI2N34,SWGT34,SCHI34,IT34,
     +                NM34(NREG34),NDO34,
     +                NTOT34,NCAL34,NCA134,NCA234,IBIM34,JCOR34,
     +                LGLO34,LLOC34
      PARAMETER(NDM3CC=5)
      PARAMETER(NBN31C=6,NRG31C=7776)
      LOGICAL LGL31C,LLC31C
      double precision ntt31c
      COMMON /HSSC31/ SIG31C,SG31EC,T31GMC,T31MXC(NRG31C),
     +                XX31C(50,5),
     +                FFG31C,DNG31C,FFL31C,DNL31C,GLD31C,
     +                SI31C,S2N31C,SWT31C,SCH31C,IT31C,
     +                NM31C(NRG31C),NDO31C,
     +                NTT31C,NCL31C,NC131C,NC231C,IBM31C,JCR31C,
     +                LGL31C,LLC31C
      PARAMETER(NBN32C=6,NRG32C=7776)
      LOGICAL LGL32C,LLC32C
      double precision ntt32c
      COMMON /HSSC32/ SIG32C,SG32EC,T32GMC,T32MXC(NRG32C),
     +                XX32C(50,5),
     +                FFG32C,DNG32C,FFL32C,DNL32C,GLD32C,
     +                SI32C,S2N32C,SWT32C,SCH32C,IT32C,
     +                NM32C(NRG32C),NDO32C,
     +                NTT32C,NCL32C,NC132C,NC232C,IBM32C,JCR32C,
     +                LGL32C,LLC32C
      PARAMETER(NBN33C=6,NRG33C=7776)
      LOGICAL LGL33C,LLC33C
      double precision ntt33c
      COMMON /HSSC33/ SIG33C,SG33EC,T33GMC,T33MXC(NRG33C),
     +                XX33C(50,5),
     +                FFG33C,DNG33C,FFL33C,DNL33C,GLD33C,
     +                SI33C,S2N33C,SWT33C,SCH33C,IT33C,
     +                NM33C(NRG33C),NDO33C,
     +                NTT33C,NCL33C,NC133C,NC233C,IBM33C,JCR33C,
     +                LGL33C,LLC33C
      PARAMETER(NDM3EL=4)
      PARAMETER(NBN31E=10,NRG31E=10000)
      LOGICAL LGL31E,LLC31E
      double precision ntt31e
      COMMON /HSSE31/ SIG31L,SG31EE,T31GME,T31MXE(NRG31E),
     +                XX31E(50,5),
     +                FFG31E,DNG31E,FFL31E,DNL31E,GLD31E,
     +                SI31E,S2N31E,SWT31E,SCH31E,IT31E,
     +                NM31E(NRG31E),NDO31E,
     +                NTT31E,NCL31E,NC131E,NC231E,IBM31E,JCR31E,
     +                LGL31E,LLC31E
      PARAMETER(NBN32E=10,NRG32E=10000)
      LOGICAL LGL32E,LLC32E
      double precision ntt32e
      COMMON /HSSE32/ SIG32L,SG32EE,T32GME,T32MXE(NRG32E),
     +                XX32E(50,5),
     +                FFG32E,DNG32E,FFL32E,DNL32E,GLD32E,
     +                SI32E,S2N32E,SWT32E,SCH32E,IT32E,
     +                NM32E(NRG32E),NDO32E,
     +                NTT32E,NCL32E,NC132E,NC232E,IBM32E,JCR32E,
     +                LGL32E,LLC32E
      PARAMETER(NBN33E=10,NRG33E=10000)
      LOGICAL LGL33E,LLC33E
      double precision ntt33e
      COMMON /HSSE33/ SIG33L,SG33EE,T33GME,T33MXE(NRG33E),
     +                XX33E(50,5),
     +                FFG33E,DNG33E,FFL33E,DNL33E,GLD33E,
     +                SI33E,S2N33E,SWT33E,SCH33E,IT33E,
     +                NM33E(NRG33E),NDO33E,
     +                NTT33E,NCL33E,NC133E,NC233E,IBM33E,JCR33E,
     +                LGL33E,LLC33E
      CHARACTER CHNAME*45
      COMMON /HSNAMC/ CHNAME(20)
      DIMENSION PROCON(20),ICONTI(20),EFFIC(20)
      PARAMETER(NEV2=1,NEV3NC=1,NEV3CC=1,NEVEL=1)
      LOGICAL LFIRST(20)
      DATA LFIRST /20*.TRUE./
C-----------------------------------------------------------------------
      WRITE(LUNOUT,'(//A/)')
     * ' CROSS SECTIONS ACTUALLY APPLIED FOR SAMPLING (IN NANOBARN): '
      WRITE(LUNOUT,'(A,1PE12.4,A,1PE11.4/)')
     *     ' TOTAL CROSS SECTION,  SIGTOT =              ',
     *       SIGTOT,' +/-',SIGTRR
      DO 110 I=1,20
        IF(SIGG(I).GT.0D0)
     *    WRITE(LUNOUT,'(A,1PE12.4,A,1PE11.4)')
     *          CHNAME(I),SIGG(I),' +/-', SIGGRR(I)
 110  CONTINUE
C-----------------------------------------------------------------------
C---CUMULATIVE PROBABILITIES FOR CHOSING THE ACTUAL CONTRIBUTION
      PROCON(1)=SIGG(1)/SIGTOT
      DO 101 I=2,20
        PROCON(I)=PROCON(I-1) + SIGG(I)/SIGTOT
  101 CONTINUE
      DO 102 I=1,5
        ICONTI(I)=ISAM2(I)
  102 CONTINUE
      DO 103 I=1,15
        ICONTI(I+5)=ISAM3(I)
  103 CONTINUE
C-----------------------------------------------------------------------
C---EVENT SAMPLING
      DO 1000 I=1,NEVENT
C---CHANNEL SELECTION
        RNC=HSRNDM(-1)
        DO 1001 NC=1,20
          IF(RNC.LE.PROCON(NC)) THEN
            NCA=NC
            GOTO 1002
          ENDIF
 1001   CONTINUE
        NCA=999
        WRITE(LUNOUT,'(A,I4/A)')
     &        ' HSEVTG: ERROR IN CHANNEL SELECTION - NCA=',NCA,
     &        '         EXECUTION STOPPED'
        STOP
C---GENERATION OF ONE SINGLE EVENT IN THE SELECTED CHANNEL
 1002   CONTINUE
        IF (NCA.EQ.1) THEN
           GOTO 1
        ELSEIF (NCA.EQ.2) THEN
           GOTO 2
        ELSEIF (NCA.EQ.3) THEN
           GOTO 3
        ELSEIF (NCA.EQ.4) THEN
           GOTO 4
        ELSEIF (NCA.EQ.5) THEN
           GOTO 5
        ELSEIF (NCA.EQ.6) THEN
           GOTO 6
        ELSEIF (NCA.EQ.7) THEN
           GOTO 7
        ELSEIF (NCA.EQ.8) THEN
           GOTO 8
        ELSEIF (NCA.EQ.9) THEN
           GOTO 9
        ELSEIF (NCA.EQ.10) THEN
           GOTO 10
        ELSEIF (NCA.EQ.11) THEN
           GOTO 11
        ELSEIF (NCA.EQ.12) THEN
           GOTO 12
        ELSEIF (NCA.EQ.13) THEN
           GOTO 13
        ELSEIF (NCA.EQ.14) THEN
           GOTO 14
        ELSEIF (NCA.EQ.15) THEN
           GOTO 15
        ELSEIF (NCA.EQ.16) THEN
           GOTO 16
        ELSEIF (NCA.EQ.17) THEN
           GOTO 17
        ELSEIF (NCA.EQ.18) THEN
           GOTO 18
        ELSEIF (NCA.EQ.19) THEN
           GOTO 19
        ELSEIF (NCA.EQ.20) THEN
           GOTO 20
        ENDIF
C
C---NON RADIATIVE EVENTS - NEUTRAL CURRENT / BORN + VIRTUAL&SOFT
  1     CONTINUE
          IF(IPRINT.GE.3)
     *      WRITE(LUNTES,'(//A/A,5X,2I3,I8,2I4)')
     *           ' CALL HSGENM ',' NCHN2,NDIM2,NEV2,NDO2,NBIN2',
     *                             NCHN2,NDIM2,NEV2,NDO2,NBIN2
          CALL HSGENM(HSNCG2,NCHN2,NDIM2,NEV2,ICONTI(1),T2GGMA,T2GMAX,
     &                GOLD2,FFGO2,FFLO2,DNCG2,DNCL2,LLOC2,LGLO2,
     &                NTOT2,NM2,NCAL2,NCA12,NCA22,IBIM2,JCOR2,
     &                XX2,NDO2,NBIN2,NREG2N)
          GOTO 1010
C---NON RADIATIVE EVENTS - CHARGED CURRENT / BORN + VIRTUAL&SOFT
  2     CONTINUE
          CALL HSGENM(HSCCG2,NCHC2,NDIM2,NEV2,ICONTI(2),T2GMAC,T2MAXC,
     &                GOLD2C,FFGO2C,FFLO2C,DNCG2C,DNCL2C,LLOC2C,LGLO2C,
     &                NTOT2C,NM2C,NCAL2C,NCA12C,NCA22C,IBIM2C,JCOR2C,
     &                XX2C,NDO2C,NBIN2,NREG2C)
          GOTO 1010
C---NON RADIATIVE EVENTS - ELASTIC EP / BORN + VIRTUAL&SOFT
  3     CONTINUE
          IF(IPRINT.GE.3)
     *      WRITE(LUNTES,'(//A/A,5X,2I3,I8,2I4)')
     *           ' CALL HSGENM ',' NCHE2,NDIM1,NEVEL,NDO2E,NBIN2',
     *                             NCHE2,NDIM1,NEVEL,NDO2E,NBIN2
          CALL HSGENM(HSELG2,NCHE2,NDIM1,NEVEL,ICONTI(3),T2GMAE,T2MAXE,
     &                GOLD2E,FFGO2E,FFLO2E,DNCG2E,DNCL2E,LLOC2E,LGLO2E,
     &                NTOT2E,NM2E,NCAL2E,NCA12E,NCA22E,IBIM2E,JCOR2E,
     &                XX2E,NDO2E,NBIN2,NREG2E)
          GOTO 1010
  4     CONTINUE
  5     CONTINUE
        WRITE(LUNOUT,'(A,I4/A)')
     &        ' HSEVTG: UNDEFINED CHANNEL SELECTED- NCA=',NCA,
     &        '         EXECUTION STOPPED'
        STOP
C---RADIATIVE EVENTS - NEUTRAL CURRENT / INITIAL STATE LEPTONIC RAD.
  6     CONTINUE
          IF(IPRINT.GE.3)
     *      WRITE(LUNTES,'(//A/A,5X,2I3,I8,2I4)')
     *           ' CALL HSGENM ',' NCHN31,NDIM31,NEV3NC,NDO31,NBIN31',
     *                             NCHN31,NDIM31,NEV3NC,NDO31,NBIN31
          CALL HSGENM(HSTSK1,NCHN31,NDIM31,NEV3NC,ICONTI(6),
     &                T31GMA,T31MAX,
     &                GOLD31,FFGO31,FFLO31,DNCG31,DNCL31,LLOC31,LGLO31,
     &                NTOT31,NM31,NCAL31,NCA131,NCA231,IBIM31,JCOR31,
     &                XX31,NDO31,NBIN31,NREG31)
          GOTO 1010
C---RADIATIVE EVENTS - NEUTRAL CURRENT / FINAL STATE LEPTONIC RAD.
  7     CONTINUE
          IF(IPRINT.GE.3)
     *      WRITE(LUNTES,'(//A/A,5X,2I3,I8,2I4)')
     *           ' CALL HSGENM ',' NCHN32,NDIM32,NEV3NC,NDO32,NBIN32',
     *                             NCHN32,NDIM32,NEV3NC,NDO32,NBIN32
          CALL HSGENM(HSTSK2,NCHN32,NDIM32,NEV3NC,ICONTI(7),
     &                T32GMA,T32MAX,
     &                GOLD32,FFGO32,FFLO32,DNCG32,DNCL32,LLOC32,LGLO32,
     &                NTOT32,NM32,NCAL32,NCA132,NCA232,IBIM32,JCOR32,
     &                XX32,NDO32,NBIN32,NREG32)
          GOTO 1010
C---RADIATIVE EVENTS - NEUTRAL CURRENT / COMPTON CONTRIBUTION
  8     CONTINUE
          IF(IPRINT.GE.1)
     *      WRITE(LUNTES,'(//A/A,5X,2I3,I8,2I4)')
     *           ' CALL HSGENM ',' NCHN33,NDIM33,NEV3NC,NDO33,NBIN33',
     *                             NCHN33,NDIM33,NEV3NC,NDO33,NBIN33
          CALL HSGENM(HSK1TS,NCHN33,NDIM33,NEV3NC,ICONTI(8),
     &                T33GMA,T33MAX,
     &                GOLD33,FFGO33,FFLO33,DNCG33,DNCL33,LLOC33,LGLO33,
     &                NTOT33,NM33,NCAL33,NCA133,NCA233,IBIM33,JCOR33,
     &                XX33,NDO33,NBIN33,NREG33)
          GOTO 1010
C---RADIATIVE EVENTS - NEUTRAL CURRENT / QUARKONIC RADIATION
  9     CONTINUE
          IF(IPRINT.GE.1)
     &      WRITE(LUNTES,'(//A/A,5X,2I3,I8,2I4)')
     &           ' CALL HSGENM ',' NCHN34,NDIM34,NEV3NC,NDO34,NBIN34',
     &                             NCHN34,NDIM34,NEV3NC,NDO34,NBIN34
          CALL HSGENM(HSK1K3,NCHN34,NDIM34,NEV3NC,ICONTI(9),
     &                T34GMA,T34MAX,
     &                GOLD34,FFGO34,FFLO34,DNCG34,DNCL34,LLOC34,LGLO34,
     &                NTOT34,NM34,NCAL34,NCA134,NCA234,IBIM34,JCOR34,
     &                XX34,NDO34,NBIN34,NREG34)
          GOTO 1010
 10     CONTINUE
 11     CONTINUE
        GOTO 20
C---RADIATIVE EVENTS - CHARGED CURRENT / LEPTONIC RADIATION
 12     CONTINUE
          CALL HSGENM(HSCCKL,NCHC31,NDM3CC,NEV3CC,ICONTI(12),T31GMC,
     &                                                          T31MXC,
     &                GLD31C,FFG31C,FFL31C,DNG31C,DNL31C,LLC31C,LGL31C,
     &                NTT31C,NM31C,NCL31C,NC131C,NC231C,IBM31C,JCR31C,
     &                XX31C,NDO31C,NBN31C,NRG31C)
          GOTO 1010
C---RADIATIVE EVENTS - CHARGED CURRENT / (KQ) CHANNEL
 13     CONTINUE
          CALL HSGENM(HSCCQI,NCHC32,NDM3CC,NEV3CC,ICONTI(13),T32GMC,
     &                                                          T32MXC,
     &                GLD32C,FFG32C,FFL32C,DNG32C,DNL32C,LLC32C,LGL32C,
     &                NTT32C,NM32C,NCL32C,NC132C,NC232C,IBM32C,JCR32C,
     &                XX32C,NDO32C,NBN32C,NRG32C)
          GOTO 1010
C---RADIATIVE EVENTS - CHARGED CURRENT / (KQS) CHANNEL
 14     CONTINUE
          CALL HSGENM(HSCCQF,NCHC33,NDM3CC,NEV3CC,ICONTI(14),T33GMC,
     &                                                          T33MXC,
     &                GLD33C,FFG33C,FFL33C,DNG33C,DNL33C,LLC33C,LGL33C,
     &                NTT33C,NM33C,NCL33C,NC133C,NC233C,IBM33C,JCR33C,
     &                XX33C,NDO33C,NBN33C,NRG33C)
          GOTO 1010
C---RADIATIVE EVENTS - ELASTIC TAIL / INITIAL STATE
 15     CONTINUE
          IF(IPRINT.GE.3)
     *      WRITE(LUNTES,'(//A/A,5X,2I3,I8,2I4)')
     *           ' CALL HSGENM ',' NCHE31,NDM3EL,NEVEL,NDO31E,NBN31E',
     *                             NCHE31,NDM3EL,NEVEL,NDO31E,NBN31E
          CALL HSGENM(HSELK1,NCHE31,NDM3EL,NEVEL,ICONTI(15),T31GME,
     &                                                          T31MXE,
     &                GLD31E,FFG31E,FFL31E,DNG31E,DNL31E,LLC31E,LGL31E,
     &                NTT31E,NM31E,NCL31E,NC131E,NC231E,IBM31E,JCR31E,
     &                XX31E,NDO31E,NBN31E,NRG31E)
          GOTO 1010
C---  RADIATIVE EVENTS - ELASTIC TAIL / FINAL STATE
 16     CONTINUE
          IF(IPRINT.GE.3)
     *      WRITE(LUNTES,'(//A/A,5X,2I3,I8,2I4)')
     *           ' CALL HSGENM ',' NCHE32,NDM3EL,NEVEL,NDO32E,NBN32E',
     *                             NCHE32,NDM3EL,NEVEL,NDO32E,NBN32E
          CALL HSGENM(HSELK2,NCHE32,NDM3EL,NEVEL,ICONTI(16),T32GME,
     &                                                          T32MXE,
     &                GLD32E,FFG32E,FFL32E,DNG32E,DNL32E,LLC32E,LGL32E,
     &                NTT32E,NM32E,NCL32E,NC132E,NC232E,IBM32E,JCR32E,
     &                XX32E,NDO32E,NBN32E,NRG32E)
          GOTO 1010
C---RADIATIVE EVENTS - ELASTIC TAIL / COMPTON
 17     CONTINUE
          IF(IPRINT.GE.3)
     *      WRITE(LUNTES,'(//A/A,5X,2I3,I8,2I4)')
     *           ' CALL HSGENM ',' NCHE33,NDM3EL,NEVEL,NDO33E,NBN33E',
     *                             NCHE33,NDM3EL,NEVEL,NDO33E,NBN33E
          CALL HSGENM(HSELCO,NCHE33,NDM3EL,NEVEL,ICONTI(17),T33GME,
     &                                                          T33MXE,
     &                GLD33E,FFG33E,FFL33E,DNG33E,DNL33E,LLC33E,LGL33E,
     &                NTT33E,NM33E,NCL33E,NC133E,NC233E,IBM33E,JCR33E,
     &                XX33E,NDO33E,NBN33E,NRG33E)
          GOTO 1010
 18     CONTINUE
 19     CONTINUE
 20     CONTINUE
        WRITE(LUNOUT,'(A,I4/A)')
     &        ' HSEVTG: UNDEFINED CHANNEL SELECTED- NCA=',NCA,
     &        '         EXECUTION STOPPED'
        STOP
C---
 1010   CONTINUE
        NEVE(NCA)=NEVE(NCA) + 1
        IF(LFIRST(NCA)) THEN
          LFIRST(NCA)=.FALSE.
          ICONTI(NCA)=2
        ENDIF
 1000 CONTINUE
C
C---Efficiencies
      IF (NCAL2.NE.0) EFFIC(1)=DFLOAT(NEVE(1))/DFLOAT(NCAL2)
      IF (NCAL2C.NE.0) EFFIC(2)=DFLOAT(NEVE(2))/DFLOAT(NCAL2C)
      IF (NCAL2E.NE.0) EFFIC(3)=DFLOAT(NEVE(3))/DFLOAT(NCAL2E)
      IF (NCAL31.NE.0) EFFIC(6)=DFLOAT(NEVE(6))/DFLOAT(NCAL31)
      IF (NCAL32.NE.0) EFFIC(7)=DFLOAT(NEVE(7))/DFLOAT(NCAL32)
      IF (NCAL33.NE.0) EFFIC(8)=DFLOAT(NEVE(8))/DFLOAT(NCAL33)
      IF (NCAL34.NE.0) EFFIC(9)=DFLOAT(NEVE(9))/DFLOAT(NCAL34)
      IF (NCL31C.NE.0) EFFIC(12)=DFLOAT(NEVE(12))/DFLOAT(NCL31C)
      IF (NCL32C.NE.0) EFFIC(13)=DFLOAT(NEVE(13))/DFLOAT(NCL32C)
      IF (NCL33C.NE.0) EFFIC(14)=DFLOAT(NEVE(14))/DFLOAT(NCL33C)
      IF (NCL31E.NE.0) EFFIC(15)=DFLOAT(NEVE(15))/DFLOAT(NCL31E)
      IF (NCL32E.NE.0) EFFIC(16)=DFLOAT(NEVE(16))/DFLOAT(NCL32E)
      IF (NCL33E.NE.0) EFFIC(17)=DFLOAT(NEVE(17))/DFLOAT(NCL33E)

      WRITE(LUNOUT,'(//2A/)') 'NUMBERS OF GENERATED EVENTS',
     &' AND EFFICIENCIES'
      WRITE(LUNOUT,'(A,I12)') 'TOTAL EVENT NUMBER', NEVENT
      DO 3911 I=1,20
        IF(NEVE(I).GT.0) WRITE(LUNOUT,'(A,I12,5X,1PE12.4)') 
     &  CHNAME(I),NEVE(I),EFFIC(I)
 3911 CONTINUE
C
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSDTIN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      PARAMETER(NDIM2=2,NBIN2=50)
      PARAMETER(NREG2N=2500)
C---                     NREG2N=NBIN2**NDIM2
      LOGICAL LGLO2,LLOC2
      double precision ntot2
      COMMON /HSSNC2/ SIG2,SIG2E,T2GGMA,T2GMAX(NREG2N),
     +                XX2(50,2),
     +                FFGO2,DNCG2,FFLO2,DNCL2,GOLD2,
     +                NTOT2,NM2(NREG2N),NDO2,
     +                NCAL2,NCA12,NCA22,IBIM2,JCOR2,
     +                LGLO2,LLOC2
      PARAMETER(NREG2C=2500)
C---                     NREG2C=NBIN2**NDIM2
      LOGICAL LGLO2C,LLOC2C
      double precision ntot2c
      COMMON /HSSCC2/ SIG2C,SIG2EC,T2GMAC,T2MAXC(NREG2C),
     +                XX2C(50,2),
     +                FFGO2C,DNCG2C,FFLO2C,DNCL2C,GOLD2C,
     +                NTOT2C,NM2C(NREG2C),NDO2C,
     +                NCAL2C,NCA12C,NCA22C,IBIM2C,JCOR2C,
     +                LGLO2C,LLOC2C
      PARAMETER(NREG2E=50)
      PARAMETER(NDIM1=1)
C---                     NREG2E=NBIN2**NDIM1
      LOGICAL LGLO2E,LLOC2E
      double precision ntot2e
      COMMON /HSSEL2/ SIG2L,SIG2EE,T2GMAE,T2MAXE(NREG2E),
     +                XX2E(50,1),
     +                FFGO2E,DNCG2E,FFLO2E,DNCL2E,GOLD2E,
     +                NTOT2E,NM2E(NREG2E),NDO2E,
     +                NCAL2E,NCA12E,NCA22E,IBIM2E,JCOR2E,
     +                LGLO2E,LLOC2E
C-----------------
      PARAMETER(NDIM31=5,NBIN31=6,NREG31=7776)
C---                     NREG31=NBIN31**NDIM31
      LOGICAL LGLO31,LLOC31
      double precision ntot31
      COMMON /HSSN31/ SIG31,SIG31E,T31GMA,T31MAX(NREG31),
     +                XX31(50,NDIM31),
     +                FFGO31,DNCG31,FFLO31,DNCL31,GOLD31,
     +                SI31,SI2N31,SWGT31,SCHI31,IT31,
     +                NM31(NREG31),NDO31,
     +                NTOT31,NCAL31,NCA131,NCA231,IBIM31,JCOR31,
     +                LGLO31,LLOC31
C
      PARAMETER(NDIM32=5,NBIN32=6,NREG32=7776)
C---                     NREG32=NBIN32**NDIM32
      LOGICAL LGLO32,LLOC32
      double precision ntot32
      COMMON /HSSN32/ SIG32,SIG32E,T32GMA,T32MAX(NREG32),
     +                XX32(50,NDIM32),
     +                FFGO32,DNCG32,FFLO32,DNCL32,GOLD32,
     +                SI32,SI2N32,SWGT32,SCHI32,IT32,
     +                NM32(NREG32),NDO32,
     +                NTOT32,NCAL32,NCA132,NCA232,IBIM32,JCOR32,
     +                LGLO32,LLOC32
C
      PARAMETER(NDIM33=5,NBIN33=6,NREG33=7776)
C---                     NREG33=NBIN33**NDIM33
      LOGICAL LGLO33,LLOC33
      double precision ntot33
      COMMON /HSSN33/ SIG33,SIG33E,T33GMA,T33MAX(NREG33),
     +                XX33(50,NDIM33),
     +                FFGO33,DNCG33,FFLO33,DNCL33,GOLD33,
     +                SI33,SI2N33,SWGT33,SCHI33,IT33,
     +                NM33(NREG33),NDO33,
     +                NTOT33,NCAL33,NCA133,NCA233,IBIM33,JCOR33,
     +                LGLO33,LLOC33
C
      PARAMETER(NDIM34=5,NBIN34=6,NREG34=7776)
C---                     NREG34=NBIN34**NDIM34
      LOGICAL LGLO34,LLOC34
      double precision ntot34
      COMMON /HSSN34/ SIG34,SIG34E,T34GMA,T34MAX(NREG34),
     +                XX34(50,NDIM34),
     +                FFGO34,DNCG34,FFLO34,DNCL34,GOLD34,
     +                SI34,SI2N34,SWGT34,SCHI34,IT34,
     +                NM34(NREG34),NDO34,
     +                NTOT34,NCAL34,NCA134,NCA234,IBIM34,JCOR34,
     +                LGLO34,LLOC34
C
      PARAMETER(NDM3CC=5)
      PARAMETER(NBN31C=6,NRG31C=7776)
C---                     NRG31C=NBN31C**NDM3CC
      LOGICAL LGL31C,LLC31C
      double precision ntt31c
      COMMON /HSSC31/ SIG31C,SG31EC,T31GMC,T31MXC(NRG31C),
     +                XX31C(50,5),
     +                FFG31C,DNG31C,FFL31C,DNL31C,GLD31C,
     +                SI31C,S2N31C,SWT31C,SCH31C,IT31C,
     +                NM31C(NRG31C),NDO31C,
     +                NTT31C,NCL31C,NC131C,NC231C,IBM31C,JCR31C,
     +                LGL31C,LLC31C
C
      PARAMETER(NBN32C=6,NRG32C=7776)
C---                     NRG32C=NBN32C**NDM3CC
      LOGICAL LGL32C,LLC32C
      double precision ntt32c
      COMMON /HSSC32/ SIG32C,SG32EC,T32GMC,T32MXC(NRG32C),
     +                XX32C(50,5),
     +                FFG32C,DNG32C,FFL32C,DNL32C,GLD32C,
     +                SI32C,S2N32C,SWT32C,SCH32C,IT32C,
     +                NM32C(NRG32C),NDO32C,
     +                NTT32C,NCL32C,NC132C,NC232C,IBM32C,JCR32C,
     +                LGL32C,LLC32C
      PARAMETER(NBN33C=6,NRG33C=7776)
C---                     NRG33C=NBN33C**NDM3CC
      LOGICAL LGL33C,LLC33C
      double precision ntt33c
      COMMON /HSSC33/ SIG33C,SG33EC,T33GMC,T33MXC(NRG33C),
     +                XX33C(50,5),
     +                FFG33C,DNG33C,FFL33C,DNL33C,GLD33C,
     +                SI33C,S2N33C,SWT33C,SCH33C,IT33C,
     +                NM33C(NRG33C),NDO33C,
     +                NTT33C,NCL33C,NC133C,NC233C,IBM33C,JCR33C,
     +                LGL33C,LLC33C
      PARAMETER(NDM3EL=4)
      PARAMETER(NBN31E=10,NRG31E=10000)
C---                     NRG31C=NBN31E**NDM3EL
      LOGICAL LGL31E,LLC31E
      double precision ntt31e
      COMMON /HSSE31/ SIG31L,SG31EE,T31GME,T31MXE(NRG31E),
     +                XX31E(50,5),
     +                FFG31E,DNG31E,FFL31E,DNL31E,GLD31E,
     +                SI31E,S2N31E,SWT31E,SCH31E,IT31E,
     +                NM31E(NRG31E),NDO31E,
     +                NTT31E,NCL31E,NC131E,NC231E,IBM31E,JCR31E,
     +                LGL31E,LLC31E
      PARAMETER(NBN32E=10,NRG32E=10000)
C---                     NRG32E=NBN32E**NDM3EL
      LOGICAL LGL32E,LLC32E
      double precision ntt32e
      COMMON /HSSE32/ SIG32L,SG32EE,T32GME,T32MXE(NRG32E),
     +                XX32E(50,5),
     +                FFG32E,DNG32E,FFL32E,DNL32E,GLD32E,
     +                SI32E,S2N32E,SWT32E,SCH32E,IT32E,
     +                NM32E(NRG32E),NDO32E,
     +                NTT32E,NCL32E,NC132E,NC232E,IBM32E,JCR32E,
     +                LGL32E,LLC32E
      PARAMETER(NBN33E=10,NRG33E=10000)
C---                     NRG33E=NBN33E**NDM3EL
      LOGICAL LGL33E,LLC33E
      double precision ntt33e
      COMMON /HSSE33/ SIG33L,SG33EE,T33GME,T33MXE(NRG33E),
     +                XX33E(50,5),
     +                FFG33E,DNG33E,FFL33E,DNL33E,GLD33E,
     +                SI33E,S2N33E,SWT33E,SCH33E,IT33E,
     +                NM33E(NRG33E),NDO33E,
     +                NTT33E,NCL33E,NC133E,NC233E,IBM33E,JCR33E,
     +                LGL33E,LLC33E
C
C-------------------------
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C
C---READ PARAMETER DEFINITIONS OF THE PREVIOUS RUN
C---CROSS CHECK WITH CURRENT ONES
C
      CALL HSTPAR
C
C---NEUTRAL CURRENT - NONRADIATIVE CHANNEL
C
      READ(LUNDAT,ERR=9900,END=9900,IOSTAT=IOS) SIG2,SIG2E,T2GGMA,NDO2
      READ(LUNDAT,ERR=9900,END=9900,IOSTAT=IOS) (T2GMAX(I),I=1,NREG2N)
      READ(LUNDAT,ERR=9900,END=9900,IOSTAT=IOS) (NM2(I),I=1,NREG2N)
      READ(LUNDAT,ERR=9900,END=9900,IOSTAT=IOS)
     &    ((XX2(I,II),I=1,NDO2),II=1,NDIM2)
      READ(LUNDAT,ERR=9900,END=9900,IOSTAT=IOS)
     &                FFGO2,DNCG2,FFLO2,DNCL2,GOLD2,
     &                NTOT2,NCAL2,NCA12,NCA22,IBIM2,JCOR2,
     &                LGLO2,LLOC2
C
 9900 CONTINUE
      IF(IOS.NE.0. OR. IPRINT.GE.2) THEN
        WRITE(LUNTES,'(//A/)') ' TEST PRINT HSDTIN'
        WRITE(LUNTES,'(/A)') ' SIG2, SIG2E, T2GGMA, NDO2'
        WRITE(LUNTES,11) SIG2,SIG2E,T2GGMA,NDO2
        IF(IPRINT.GE.3) THEN
          WRITE(LUNTES,'(/,A)') ' T2GMAX(NREG2N)'
          WRITE(LUNTES,12) (T2GMAX(I),I=1,NREG2N)
          WRITE(LUNTES,'(/,A)') ' NM2(NREG2N)'
          WRITE(LUNTES,13) (NM2(I),I=1,NREG2N)
          WRITE(LUNTES,'(/,A)') ' XX2(...)'
          WRITE(LUNTES,12) ((XX2(I,II),I=1,NDO2),II=1,NDIM2)
        ENDIF
C
        WRITE(LUNTES,'(/A/A/A)') ' FFGO2,DNCG2,FFLO2,DNCL2',
     &                   'NTOT2,NCAL2,NCA12,NCA22,IBIM2,JCOR2',
     &                   'LGLO2,LLOC2'
        WRITE(LUNTES,'(4(1PE13.4)/1PE25.16/5I10/2L4)')
     &                   FFGO2,DNCG2,FFLO2,DNCL2,
     &                   NTOT2,NCAL2,NCA12,NCA22,IBIM2,JCOR2,
     &                   LGLO2,LLOC2
        IF(IOS.NE.0) GOTO 9999
      ENDIF
C
C---NEUTRAL CURRENT - INITIAL STATE LEPTONIC RADIATION
C
      READ(LUNDAT,ERR=9901,END=9901,IOSTAT=IOS)
     &     SIG31,SIG31E,T31GMA,NDO31
      IF(IPRINT.GE.5) THEN
        WRITE(LUNTES,'(///A)') ' HSDTIN :  SIG31,SIG31E,T31GMA,NDO31'
        WRITE(LUNTES,*) SIG31,SIG31E,T31GMA,NDO31
      ENDIF
      READ(LUNDAT,ERR=9901,END=9901,IOSTAT=IOS) (T31MAX(I),I=1,NREG31)
      READ(LUNDAT,ERR=9901,END=9901,IOSTAT=IOS) (NM31(I),I=1,NREG31)
      READ(LUNDAT,ERR=9901,END=9901,IOSTAT=IOS)
     &     IT31,SI31,SI2N31,SWGT31,SCHI31
      READ(LUNDAT,ERR=9901,END=9901,IOSTAT=IOS)
     &     ((XX31(I,II),I=1,NDO31),II=1,NDIM31)
      READ(LUNDAT,ERR=9901,END=9901,IOSTAT=IOS)
     &                FFGO31,DNCG31,FFLO31,DNCL31,GOLD31,
     &                NTOT31,NCAL31,NCA131,NCA231,IBIM31,JCOR31,
     &                LGLO31,LLOC31
C
 9901 CONTINUE
C
      IF(IPRINT.GE.6 .OR. IOS.NE.0) THEN
        WRITE(LUNTES,'(/A/5X,5(1PD13.5))')
     &       ' HSDTIN :    FFGO31,DNCG31,FFLO31,DNCL31,GOLD31',
     &                    FFGO31,DNCG31,FFLO31,DNCL31,GOLD31
        WRITE(LUNTES,'(/A/5X,1PE25.16,I8,4I5)')
     &       ' HSDTIN :    NTOT31,NCAL31,NCA131,NCA231,IBIM31,JCOR31',
     &                    NTOT31,NCAL31,NCA131,NCA231,IBIM31,JCOR31
        WRITE(LUNTES,'(/A,2L4)')
     &       ' HSDTIN:     LGLO31,LLOC31',  LGLO31,LLOC31
        IF(IOS.NE.0) GOTO 9999
      ENDIF
C
C---NEUTRAL CURRENT - FINAL STATE LEPTONIC RADIATION
C
      READ(LUNDAT,ERR=9902,END=9902,IOSTAT=IOS)
     &     SIG32,SIG32E,T32GMA,NDO32
      READ(LUNDAT,ERR=9902,END=9902,IOSTAT=IOS) (T32MAX(I),I=1,NREG32)
      READ(LUNDAT,ERR=9902,END=9902,IOSTAT=IOS) (NM32(I),I=1,NREG32)
      READ(LUNDAT,ERR=9902,END=9902,IOSTAT=IOS)
     &     IT32,SI32,SI2N32,SWGT32,SCHI32
      READ(LUNDAT,ERR=9902,END=9902,IOSTAT=IOS)
     &    ((XX32(I,II),I=1,NDO32),II=1,NDIM32)
      READ(LUNDAT,ERR=9902,END=9902,IOSTAT=IOS)
     &     FFGO32,DNCG32,FFLO32,DNCL32,GOLD32,
     &     NTOT32,NCAL32,NCA132,NCA232,IBIM32,JCOR32,
     &     LGLO32,LLOC32
C
 9902 CONTINUE
C
      IF(IPRINT.GE.6 .OR. IOS.NE.0) THEN
        WRITE(LUNTES,'(/A/5X,5(1PD13.5))')
     &       ' HSDTIN :    FFGO32,DNCG32,FFLO32,DNCL32,GOLD32',
     &                    FFGO32,DNCG32,FFLO32,DNCL32,GOLD32
        WRITE(LUNTES,'(/A/5X,1PE25.16,I8,4I5)')
     &       ' HSDTIN :    NTOT32,NCAL32,NCA132,NCA232,IBIM32,JCOR32',
     &                    NTOT32,NCAL32,NCA132,NCA232,IBIM32,JCOR32
        WRITE(LUNTES,'(/A,2L4)')
     &       ' HSDTIN:     LGLO32,LLOC32',  LGLO32,LLOC32
        IF(IOS.NE.0) GOTO 9999
      ENDIF
C
C---NEUTRAL CURRENT - COMPTON CHANNEL
C
      READ(LUNDAT,ERR=9903,END=9903,IOSTAT=IOS)
     &     SIG33,SIG33E,T33GMA,NDO33
      IF(IPRINT.GE.5) THEN
        WRITE(LUNTES,'(///A)') ' HSDTIN :  SIG33,SIG33E,T33GMA,NDO33'
        WRITE(LUNTES,*) SIG33,SIG33E,T33GMA,NDO33
      ENDIF
      READ(LUNDAT,ERR=9903,END=9903,IOSTAT=IOS) (T33MAX(I),I=1,NREG33)
      READ(LUNDAT,ERR=9903,END=9903,IOSTAT=IOS) (NM33(I),I=1,NREG33)
      READ(LUNDAT,ERR=9903,END=9903,IOSTAT=IOS)
     &     IT33,SI33,SI2N33,SWGT33,SCHI33
      READ(LUNDAT,ERR=9903,END=9903,IOSTAT=IOS)
     &     ((XX33(I,II),I=1,NDO33),II=1,NDIM33)
      READ(LUNDAT,ERR=9903,END=9903,IOSTAT=IOS)
     &     FFGO33,DNCG33,FFLO33,DNCL33,GOLD33,
     &     NTOT33,NCAL33,NCA133,NCA233,IBIM33,JCOR33,
     &     LGLO33,LLOC33
C
 9903 CONTINUE
C
      IF(IPRINT.GE.6 .OR. IOS.NE.0) THEN
        WRITE(LUNTES,'(/A/5X,5(1PD13.5))')
     &       ' HSDTIN :    FFGO33,DNCG33,FFLO33,DNCL33,GOLD33',
     &                    FFGO33,DNCG33,FFLO33,DNCL33,GOLD33
        WRITE(LUNTES,'(/A/5X,1PE25.16,I8,4I5)')
     &       ' HSDTIN :    NTOT33,NCAL33,NCA133,NCA233,IBIM33,JCOR33',
     &                    NTOT33,NCAL33,NCA133,NCA233,IBIM33,JCOR33
        WRITE(LUNTES,'(/A,2L4)')
     &       ' HSDTIN:     LGLO33,LLOC33',  LGLO33,LLOC33
        IF(IOS.NE.0) GOTO 9999
      ENDIF
C
C---NEUTRAL CURRENT - QUARKONIC RADIATION
C
      READ(LUNDAT,ERR=9904,END=9904,IOSTAT=IOS)
     &     SIG34,SIG34E,T34GMA,NDO34
      IF(IPRINT.GE.5) THEN
        WRITE(LUNTES,'(///A)') ' HSDTIN :  SIG34,SIG34E,T34GMA,NDO34'
        WRITE(LUNTES,*) SIG34,SIG34E,T34GMA,NDO34
      ENDIF
      READ(LUNDAT,ERR=9904,END=9904,IOSTAT=IOS) (T34MAX(I),I=1,NREG34)
      READ(LUNDAT,ERR=9904,END=9904,IOSTAT=IOS) (NM34(I),I=1,NREG34)
      READ(LUNDAT,ERR=9904,END=9904,IOSTAT=IOS)
     &     IT34,SI34,SI2N34,SWGT34,SCHI34
      READ(LUNDAT,ERR=9904,END=9904,IOSTAT=IOS)
     &     ((XX34(I,II),I=1,NDO34),II=1,NDIM34)
      READ(LUNDAT,ERR=9904,END=9904,IOSTAT=IOS)
     &     FFGO34,DNCG34,FFLO34,DNCL34,GOLD34,
     &     NTOT34,NCAL34,NCA134,NCA234,IBIM34,JCOR34,
     &     LGLO34,LLOC34
C
 9904 CONTINUE
C
      IF(IPRINT.GE.6 .OR. IOS.NE.0) THEN
        WRITE(LUNTES,'(/A/5X,5(1PD13.5))')
     &       ' HSDTIN :    FFGO34,DNCG34,FFLO34,DNCL34,GOLD34',
     &                    FFGO34,DNCG34,FFLO34,DNCL34,GOLD34
        WRITE(LUNTES,'(/A/5X,1PE25.16,I8,4I5)')
     &       ' HSDTIN :    NTOT34,NCAL34,NCA134,NCA234,IBIM34,JCOR34',
     &                    NTOT34,NCAL34,NCA134,NCA234,IBIM34,JCOR34
        WRITE(LUNTES,'(/A,2L4)')
     &       ' HSDTIN:     LGLO34,LLOC34',  LGLO34,LLOC34
        IF(IOS.NE.0) GOTO 9999
      ENDIF
C
C---CHARGED CURRENT - NONRADIATIVE CHANNEL
C
      READ(LUNDAT,ERR=9905,END=9905,IOSTAT=IOS)SIG2C,SIG2EC,T2GMAC,NDO2C
      READ(LUNDAT,ERR=9905,END=9905,IOSTAT=IOS) (T2MAXC(I),I=1,NREG2C)
      READ(LUNDAT,ERR=9905,END=9905,IOSTAT=IOS) (NM2C(I),I=1,NREG2C)
      READ(LUNDAT,ERR=9905,END=9905,IOSTAT=IOS)
     &    ((XX2C(I,II),I=1,NDO2C),II=1,NDIM2)
      READ(LUNDAT,ERR=9905,END=9905,IOSTAT=IOS)
     &                FFGO2C,DNCG2C,FFLO2C,DNCL2C,GOLD2C,
     &                NTOT2C,NCAL2C,NCA12C,NCA22C,IBIM2C,JCOR2C,
     &                LGLO2C,LLOC2C
C
C---CHARGED CURRENT - (KP) CHANNEL
C
 9905 CONTINUE
      READ(LUNDAT,ERR=9906,END=9906,IOSTAT=IOS)
     &     SIG31C,SG31EC,T31GMC,NDO31C
      READ(LUNDAT,ERR=9906,END=9906,IOSTAT=IOS) (T31MXC(I),I=1,NRG31C)
      READ(LUNDAT,ERR=9906,END=9906,IOSTAT=IOS) (NM31C(I),I=1,NRG31C)
      READ(LUNDAT,ERR=9906,END=9906,IOSTAT=IOS)
     &     IT31C,SI31C,S2N31C,SWT31C,SCH31C
      READ(LUNDAT,ERR=9906,END=9906,IOSTAT=IOS)
     &     ((XX31C(I,II),I=1,NDO31C),II=1,NDM3CC)
      READ(LUNDAT,ERR=9906,END=9906,IOSTAT=IOS)
     &                FFG31C,DNG31C,FFL31C,DNL31C,GLD31C,
     &                NTT31C,NCL31C,NC131C,NC231C,IBM31C,JCR31C,
     &                LGL31C,LLC31C
C
C---CHARGED CURRENT - (KQ) CHANNEL
C
 9906 CONTINUE
      READ(LUNDAT,ERR=9907,END=9907,IOSTAT=IOS)
     &     SIG32C,SG32EC,T32GMC,NDO32C
      READ(LUNDAT,ERR=9907,END=9907,IOSTAT=IOS) (T32MXC(I),I=1,NRG32C)
      READ(LUNDAT,ERR=9907,END=9907,IOSTAT=IOS) (NM32C(I),I=1,NRG32C)
      READ(LUNDAT,ERR=9907,END=9907,IOSTAT=IOS)
     &     IT32C,SI32C,S2N32C,SWT32C,SCH32C
      READ(LUNDAT,ERR=9907,END=9907,IOSTAT=IOS)
     &     ((XX32C(I,II),I=1,NDO32C),II=1,NDM3CC)
      READ(LUNDAT,ERR=9907,END=9907,IOSTAT=IOS)
     &                FFG32C,DNG32C,FFL32C,DNL32C,GLD32C,
     &                NTT32C,NCL32C,NC132C,NC232C,IBM32C,JCR32C,
     &                LGL32C,LLC32C
C
C---CHARGED CURRENT - (KQS) CHANNEL
C
 9907 CONTINUE
      READ(LUNDAT,ERR=9908,END=9908,IOSTAT=IOS)
     &     SIG33C,SG33EC,T33GMC,NDO33C
      READ(LUNDAT,ERR=9908,END=9908,IOSTAT=IOS) (T33MXC(I),I=1,NRG33C)
      READ(LUNDAT,ERR=9908,END=9908,IOSTAT=IOS) (NM33C(I),I=1,NRG33C)
      READ(LUNDAT,ERR=9908,END=9908,IOSTAT=IOS)
     &     IT33C,SI33C,S2N33C,SWT33C,SCH33C
      READ(LUNDAT,ERR=9908,END=9908,IOSTAT=IOS)
     &     ((XX33C(I,II),I=1,NDO33C),II=1,NDM3CC)
      READ(LUNDAT,ERR=9908,END=9908,IOSTAT=IOS)
     &                FFG33C,DNG33C,FFL33C,DNL33C,GLD33C,
     &                NTT33C,NCL33C,NC133C,NC233C,IBM33C,JCR33C,
     &                LGL33C,LLC33C
C
C
C---ELASTIC NON-RADIATIVE CHANNEL
C
 9908 CONTINUE
      READ(LUNDAT,ERR=9909,END=9909,IOSTAT=IOS)SIG2L,SIG2EE,T2GMAE,NDO2E
      READ(LUNDAT,ERR=9909,END=9909,IOSTAT=IOS) (T2MAXE(I),I=1,NREG2E)
      READ(LUNDAT,ERR=9909,END=9909,IOSTAT=IOS) (NM2E(I),I=1,NREG2E)
      READ(LUNDAT,ERR=9909,END=9909,IOSTAT=IOS)
     &    (XX2E(I,1),I=1,NDO2E)
      READ(LUNDAT,ERR=9909,END=9909,IOSTAT=IOS)
     &                FFGO2E,DNCG2E,FFLO2E,DNCL2E,GOLD2E,
     &                NTOT2E,NCAL2E,NCA12E,NCA22E,IBIM2E,JCOR2E,
     &                LGLO2E,LLOC2E
C
C---ELASTIC TAIL (KP) CHANNEL
C
 9909 CONTINUE
      READ(LUNDAT,ERR=9910,END=9910,IOSTAT=IOS)
     &     SIG31L,SG31EE,T31GME,NDO31E
      READ(LUNDAT,ERR=9910,END=9910,IOSTAT=IOS) (T31MXE(I),I=1,NRG31E)
      READ(LUNDAT,ERR=9910,END=9910,IOSTAT=IOS) (NM31E(I),I=1,NRG31E)
      READ(LUNDAT,ERR=9910,END=9910,IOSTAT=IOS)
     &     IT31E,SI31E,S2N31E,SWT31E,SCH31E
      READ(LUNDAT,ERR=9910,END=9910,IOSTAT=IOS)
     &     ((XX31E(I,II),I=1,NDO31E),II=1,NDM3EL)
      READ(LUNDAT,ERR=9910,END=9910,IOSTAT=IOS)
     &                FFG31E,DNG31E,FFL31E,DNL31E,GLD31E,
     &                NTT31E,NCL31E,NC131E,NC231E,IBM31E,JCR31E,
     &                LGL31E,LLC31E
C
C---ELASTIC TAIL (KPS) CHANNEL
C
 9910 CONTINUE
      READ(LUNDAT,ERR=9911,END=9911,IOSTAT=IOS)
     &     SIG32L,SG32EE,T32GME,NDO32E
      READ(LUNDAT,ERR=9911,END=9911,IOSTAT=IOS) (T32MXE(I),I=1,NRG32E)
      READ(LUNDAT,ERR=9911,END=9911,IOSTAT=IOS) (NM32E(I),I=1,NRG32E)
      READ(LUNDAT,ERR=9911,END=9911,IOSTAT=IOS)
     &     IT32E,SI32E,S2N32E,SWT32E,SCH32E
      READ(LUNDAT,ERR=9911,END=9911,IOSTAT=IOS)
     &     ((XX32E(I,II),I=1,NDO32E),II=1,NDM3EL)
      READ(LUNDAT,ERR=9911,END=9911,IOSTAT=IOS)
     &                FFG32E,DNG32E,FFL32E,DNL32E,GLD32E,
     &                NTT32E,NCL32E,NC132E,NC232E,IBM32E,JCR32E,
     &                LGL32E,LLC32E
C
C---ELASTIC TAIL COMPTON PART
C
 9911 CONTINUE
      READ(LUNDAT,ERR=9912,END=9912,IOSTAT=IOS)
     &     SIG33L,SG33EE,T33GME,NDO33E
      READ(LUNDAT,ERR=9912,END=9912,IOSTAT=IOS) (T33MXE(I),I=1,NRG33E)
      READ(LUNDAT,ERR=9912,END=9912,IOSTAT=IOS) (NM33E(I),I=1,NRG33E)
      READ(LUNDAT,ERR=9912,END=9912,IOSTAT=IOS)
     &     IT33E,SI33E,S2N33E,SWT33E,SCH33E
      READ(LUNDAT,ERR=9912,END=9912,IOSTAT=IOS)
     &     ((XX33E(I,II),I=1,NDO33E),II=1,NDM3EL)
      READ(LUNDAT,ERR=9912,END=9912,IOSTAT=IOS)
     &                FFG33E,DNG33E,FFL33E,DNL33E,GLD33E,
     &                NTT33E,NCL33E,NC133E,NC233E,IBM33E,JCR33E,
     &                LGL33E,LLC33E
C
 9912 CONTINUE
      IF(IOS.NE.0) GOTO 9999
      RETURN
C------------------------------
C                    EXECUTION HALTED BECAUSE OF READING ERROR
 9999 CONTINUE
        WRITE(LUNOUT,'(A,I3/A,I3/A)')
     &          ' ERROR READING UNIT =',LUNDAT,
     &          ' IOSTAT =', IOS, ' EXECUTION HALTED IN HSDTIN'
      STOP
C
 11   FORMAT(3(1PD15.5),3I5)
 12   FORMAT(5(1PD15.5))
 13   FORMAT(10I8)
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSDOUT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     &                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      PARAMETER(NDIM2=2,NBIN2=50)
      PARAMETER(NREG2N=2500)
C---                     NREG2N=NBIN2**NDIM2
      LOGICAL LGLO2,LLOC2
      double precision ntot2
      COMMON /HSSNC2/ SIG2,SIG2E,T2GGMA,T2GMAX(NREG2N),
     +                XX2(50,2),
     +                FFGO2,DNCG2,FFLO2,DNCL2,GOLD2,
     +                NTOT2,NM2(NREG2N),NDO2,
     +                NCAL2,NCA12,NCA22,IBIM2,JCOR2,
     +                LGLO2,LLOC2
      PARAMETER(NREG2C=2500)
C---                     NREG2C=NBIN2**NDIM2
      LOGICAL LGLO2C,LLOC2C
      double precision ntot2c
      COMMON /HSSCC2/ SIG2C,SIG2EC,T2GMAC,T2MAXC(NREG2C),
     +                XX2C(50,2),
     +                FFGO2C,DNCG2C,FFLO2C,DNCL2C,GOLD2C,
     +                NTOT2C,NM2C(NREG2C),NDO2C,
     +                NCAL2C,NCA12C,NCA22C,IBIM2C,JCOR2C,
     +                LGLO2C,LLOC2C
      PARAMETER(NREG2E=50)
      PARAMETER(NDIM1=1)
C---                     NREG2E=NBIN2**NDIM1
      LOGICAL LGLO2E,LLOC2E
      double precision ntot2e
      COMMON /HSSEL2/ SIG2L,SIG2EE,T2GMAE,T2MAXE(NREG2E),
     +                XX2E(50,1),
     +                FFGO2E,DNCG2E,FFLO2E,DNCL2E,GOLD2E,
     +                NTOT2E,NM2E(NREG2E),NDO2E,
     +                NCAL2E,NCA12E,NCA22E,IBIM2E,JCOR2E,
     +                LGLO2E,LLOC2E
C-----------------
      PARAMETER(NDIM31=5,NBIN31=6,NREG31=7776)
C---                     NREG31=NBIN31**NDIM31
      LOGICAL LGLO31,LLOC31
      double precision ntot31
      COMMON /HSSN31/ SIG31,SIG31E,T31GMA,T31MAX(NREG31),
     +                XX31(50,NDIM31),
     +                FFGO31,DNCG31,FFLO31,DNCL31,GOLD31,
     +                SI31,SI2N31,SWGT31,SCHI31,IT31,
     +                NM31(NREG31),NDO31,
     +                NTOT31,NCAL31,NCA131,NCA231,IBIM31,JCOR31,
     +                LGLO31,LLOC31
      PARAMETER(NDIM32=5,NBIN32=6,NREG32=7776)
C---                     NREG32=NBIN32**NDIM32
      LOGICAL LGLO32,LLOC32
      double precision ntot32
      COMMON /HSSN32/ SIG32,SIG32E,T32GMA,T32MAX(NREG32),
     +                XX32(50,NDIM32),
     +                FFGO32,DNCG32,FFLO32,DNCL32,GOLD32,
     +                SI32,SI2N32,SWGT32,SCHI32,IT32,
     +                NM32(NREG32),NDO32,
     +                NTOT32,NCAL32,NCA132,NCA232,IBIM32,JCOR32,
     +                LGLO32,LLOC32
      PARAMETER(NDIM33=5,NBIN33=6,NREG33=7776)
C---                     NREG33=NBIN33**NDIM33
      LOGICAL LGLO33,LLOC33
      double precision ntot33
      COMMON /HSSN33/ SIG33,SIG33E,T33GMA,T33MAX(NREG33),
     +                XX33(50,NDIM33),
     +                FFGO33,DNCG33,FFLO33,DNCL33,GOLD33,
     +                SI33,SI2N33,SWGT33,SCHI33,IT33,
     +                NM33(NREG33),NDO33,
     +                NTOT33,NCAL33,NCA133,NCA233,IBIM33,JCOR33,
     +                LGLO33,LLOC33
      PARAMETER(NDIM34=5,NBIN34=6,NREG34=7776)
C---                     NREG34=NBIN34**NDIM34
      LOGICAL LGLO34,LLOC34
      double precision ntot34
      COMMON /HSSN34/ SIG34,SIG34E,T34GMA,T34MAX(NREG34),
     +                XX34(50,NDIM34),
     +                FFGO34,DNCG34,FFLO34,DNCL34,GOLD34,
     +                SI34,SI2N34,SWGT34,SCHI34,IT34,
     +                NM34(NREG34),NDO34,
     +                NTOT34,NCAL34,NCA134,NCA234,IBIM34,JCOR34,
     +                LGLO34,LLOC34
C
      PARAMETER(NDM3CC=5)
      PARAMETER(NBN31C=6,NRG31C=7776)
C---                     NRG31C=NBN31C**NDM3CC
      LOGICAL LGL31C,LLC31C
      double precision ntt31c
      COMMON /HSSC31/ SIG31C,SG31EC,T31GMC,T31MXC(NRG31C),
     +                XX31C(50,5),
     +                FFG31C,DNG31C,FFL31C,DNL31C,GLD31C,
     +                SI31C,S2N31C,SWT31C,SCH31C,IT31C,
     +                NM31C(NRG31C),NDO31C,
     +                NTT31C,NCL31C,NC131C,NC231C,IBM31C,JCR31C,
     +                LGL31C,LLC31C
C
      PARAMETER(NBN32C=6,NRG32C=7776)
C---                     NRG32C=NBN32C**NDM3CC
      LOGICAL LGL32C,LLC32C
      double precision ntt32c
      COMMON /HSSC32/ SIG32C,SG32EC,T32GMC,T32MXC(NRG32C),
     +                XX32C(50,5),
     +                FFG32C,DNG32C,FFL32C,DNL32C,GLD32C,
     +                SI32C,S2N32C,SWT32C,SCH32C,IT32C,
     +                NM32C(NRG32C),NDO32C,
     +                NTT32C,NCL32C,NC132C,NC232C,IBM32C,JCR32C,
     +                LGL32C,LLC32C
      PARAMETER(NBN33C=6,NRG33C=7776)
C---                     NRG33C=NBN33C**NDM3CC
      LOGICAL LGL33C,LLC33C
      double precision ntt33c
      COMMON /HSSC33/ SIG33C,SG33EC,T33GMC,T33MXC(NRG33C),
     +                XX33C(50,5),
     +                FFG33C,DNG33C,FFL33C,DNL33C,GLD33C,
     +                SI33C,S2N33C,SWT33C,SCH33C,IT33C,
     +                NM33C(NRG33C),NDO33C,
     +                NTT33C,NCL33C,NC133C,NC233C,IBM33C,JCR33C,
     +                LGL33C,LLC33C
      PARAMETER(NDM3EL=4)
      PARAMETER(NBN31E=10,NRG31E=10000)
C---                     NRG31C=NBN31E**NDM3EL
      LOGICAL LGL31E,LLC31E
      double precision ntt31e
      COMMON /HSSE31/ SIG31L,SG31EE,T31GME,T31MXE(NRG31E),
     +                XX31E(50,5),
     +                FFG31E,DNG31E,FFL31E,DNL31E,GLD31E,
     +                SI31E,S2N31E,SWT31E,SCH31E,IT31E,
     +                NM31E(NRG31E),NDO31E,
     +                NTT31E,NCL31E,NC131E,NC231E,IBM31E,JCR31E,
     +                LGL31E,LLC31E
      PARAMETER(NBN32E=10,NRG32E=10000)
C---                     NRG32E=NBN32E**NDM3EL
      LOGICAL LGL32E,LLC32E
      double precision ntt32e
      COMMON /HSSE32/ SIG32L,SG32EE,T32GME,T32MXE(NRG32E),
     +                XX32E(50,5),
     +                FFG32E,DNG32E,FFL32E,DNL32E,GLD32E,
     +                SI32E,S2N32E,SWT32E,SCH32E,IT32E,
     +                NM32E(NRG32E),NDO32E,
     +                NTT32E,NCL32E,NC132E,NC232E,IBM32E,JCR32E,
     +                LGL32E,LLC32E
      PARAMETER(NBN33E=10,NRG33E=10000)
C---                     NRG33E=NBN33E**NDM3EL
      LOGICAL LGL33E,LLC33E
      double precision ntt33e
      COMMON /HSSE33/ SIG33L,SG33EE,T33GME,T33MXE(NRG33E),
     +                XX33E(50,5),
     +                FFG33E,DNG33E,FFL33E,DNL33E,GLD33E,
     +                SI33E,S2N33E,SWT33E,SCH33E,IT33E,
     +                NM33E(NRG33E),NDO33E,
     +                NTT33E,NCL33E,NC133E,NC233E,IBM33E,JCR33E,
     +                LGL33E,LLC33E
C---------------------------
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C-----------------------------------------------------------------------
C
C---WRITE PARAMETER DEFINITIONS OF THE CURRENT RUN
      CALL HSWRPA
C
C---WRITE INFORMATION FOR SAMPLING -------------------------------------
C---NEUTRAL CURRENT - NONRADIATIVE CHANNEL
      IF(IPRINT.GE.2) THEN
        WRITE(LUNTES,'(//2A,I3)')
     &       ' *** WRITE SAMPLING INFORMATION FOR NEUTRAL CURRENT/',
     &       'NON-RADIATIVE CHANNEL ONTO UNIT', LUNDAT
      ENDIF
      ISET=0
      CALL HSWRSA(ISET,NREG2N,NDIM2,
     &            SIG2,SIG2E,T2GGMA,NDO2,T2GMAX,NM2,XX2,
     &            FFGO2,DNCG2,FFLO2,DNCL2,GOLD2,
     &            NTOT2,NCAL2,NCA12,NCA22,IBIM2,JCOR2,
     &            LGLO2,LLOC2,
     &            IT31,SI31,SI2N31,SWGT31,SCHI31)
C***    NOTE: LAST 6 VARIABLES NOT USED IN THIS CASE!!!
C
C---NEUTRAL CURRENT - INITIAL STATE LEPTONIC RADIATION
      ISET=1
      CALL HSWRSA(ISET,NREG31,NDIM31,
     &            SIG31,SIG31E,T31GMA,NDO31,T31MAX,NM31,XX31,
     &            FFGO31,DNCG31,FFLO31,DNCL31,GOLD31,
     &            NTOT31,NCAL31,NCA131,NCA231,IBIM31,JCOR31,
     &            LGLO31,LLOC31,
     &            IT31,SI31,SI2N31,SWGT31,SCHI31)
C
C---NEUTRAL CURRENT - FINAL STATE LEPTONIC RADIATION
      ISET=1
      CALL HSWRSA(ISET,NREG32,NDIM32,
     &            SIG32,SIG32E,T32GMA,NDO32,T32MAX,NM32,XX32,
     &            FFGO32,DNCG32,FFLO32,DNCL32,GOLD32,
     &            NTOT32,NCAL32,NCA132,NCA232,IBIM32,JCOR32,
     &            LGLO32,LLOC32,
     &            IT32,SI32,SI2N32,SWGT32,SCHI32)
C
C---NEUTRAL CURRENT - COMPTON CHANNEL
      ISET=1
      CALL HSWRSA(ISET,NREG33,NDIM33,
     &            SIG33,SIG33E,T33GMA,NDO33,T33MAX,NM33,XX33,
     &            FFGO33,DNCG33,FFLO33,DNCL33,GOLD33,
     &            NTOT33,NCAL33,NCA133,NCA233,IBIM33,JCOR33,
     &            LGLO33,LLOC33,
     &            IT33,SI33,SI2N33,SWGT33,SCHI33)
C
C---NEUTRAL CURRENT - QUARKONIC RADIATION
      ISET=1
      CALL HSWRSA(ISET,NREG34,NDIM34,
     &            SIG34,SIG34E,T34GMA,NDO34,T34MAX,NM34,XX34,
     &            FFGO34,DNCG34,FFLO34,DNCL34,GOLD34,
     &            NTOT34,NCAL34,NCA134,NCA234,IBIM34,JCOR34,
     &            LGLO34,LLOC34,
     &            IT34,SI34,SI2N34,SWGT34,SCHI34)
C
C---CHARGED CURRENT - NONRADIATIVE CHANNEL
      IF(IPRINT.GE.2) THEN
        WRITE(LUNTES,'(//2A,I3)')
     &       ' *** WRITE SAMPLING INFORMATION FOR CHARGED CURRENT/',
     &       'NON-RADIATIVE CHANNEL ONTO UNIT', LUNDAT
      ENDIF
      ISET=0
      CALL HSWRSA(ISET,NREG2C,NDIM2,
     &            SIG2C,SIG2EC,T2GMAC,NDO2C,T2MAXC,NM2C,XX2C,
     &            FFGO2C,DNCG2C,FFLO2C,DNCL2C,GOLD2C,
     &            NTOT2C,NCAL2C,NCA12C,NCA22C,IBIM2C,JCOR2C,
     &            LGLO2C,LLOC2C,
     &            IT31C,SI31C,S2N31C,SWT31C,SCH31C)
C
C---CHARGED CURRENT - (KP) CHANNEL
      ISET=1
      CALL HSWRSA(ISET,NRG31C,NDM3CC,
     &            SIG31C,SG31EC,T31GMC,NDO31C,T31MXC,NM31C,XX31C,
     &            FFG31C,DNG31C,FFL31C,DNL31C,GLD31C,
     &            NTT31C,NCL31C,NC131C,NC231C,IBM31C,JCR31C,
     &            LGL31C,LLC31C,
     &            IT31C,SI31C,S2N31C,SWT31C,SCH31C)
C
C---CHARGED CURRENT - (KQ) CHANNEL
      ISET=1
      CALL HSWRSA(ISET,NRG32C,NDM3CC,
     &            SIG32C,SG32EC,T32GMC,NDO32C,T32MXC,NM32C,XX32C,
     &            FFG32C,DNG32C,FFL32C,DNL32C,GLD32C,
     &            NTT32C,NCL32C,NC132C,NC232C,IBM32C,JCR32C,
     &            LGL32C,LLC32C,
     &            IT32C,SI32C,S2N32C,SWT32C,SCH32C)
C
C---CHARGED CURRENT - (KQS) CHANNEL
      ISET=1
      CALL HSWRSA(ISET,NRG33C,NDM3CC,
     &            SIG33C,SG33EC,T33GMC,NDO33C,T33MXC,NM33C,XX33C,
     &            FFG33C,DNG33C,FFL33C,DNL33C,GLD33C,
     &            NTT33C,NCL33C,NC133C,NC233C,IBM33C,JCR33C,
     &            LGL33C,LLC33C,
     &            IT33C,SI33C,S2N33C,SWT33C,SCH33C)
C
C---ELASTIC EP - NONRADIATIVE CHANNEL
      IF(IPRINT.GE.2) THEN
        WRITE(LUNTES,'(//2A,I3)')
     &       ' *** WRITE SAMPLING INFORMATION FOR ELASTIC EP /',
     &       'NON-RADIATIVE CHANNEL ONTO UNIT', LUNDAT
      ENDIF
      ISET=0
      CALL HSWRSA(ISET,NREG2E,NDIM1,
     &            SIG2L,SIG2EE,T2GMAE,NDO2E,T2MAXE,NM2E,XX2E,
     &            FFGO2E,DNCG2E,FFLO2E,DNCL2E,GOLD2E,
     &            NTOT2E,NCAL2E,NCA12E,NCA22E,IBIM2E,JCOR2E,
     &            LGLO2E,LLOC2E,
     &            IT31E,SI31E,S2N31E,SWT31E,SCH31E)
C
C---ELASTIC TAIL (KP) CHANNEL
      ISET=1
      CALL HSWRSA(ISET,NRG31E,NDM3EL,
     &            SIG31L,SG31EE,T31GME,NDO31E,T31MXE,NM31E,XX31E,
     &            FFG31E,DNG31E,FFL31E,DNL31E,GLD31E,
     &            NTT31E,NCL31E,NC131E,NC231E,IBM31E,JCR31E,
     &            LGL31E,LLC31E,
     &            IT31E,SI31E,S2N31E,SWT31E,SCH31E)
C
C---ELASTIC TAIL (KPS) CHANNEL
      ISET=1
      CALL HSWRSA(ISET,NRG32E,NDM3EL,
     &            SIG32L,SG32EE,T32GME,NDO32E,T32MXE,NM32E,XX32E,
     &            FFG32E,DNG32E,FFL32E,DNL32E,GLD32E,
     &            NTT32E,NCL32E,NC132E,NC232E,IBM32E,JCR32E,
     &            LGL32E,LLC32E,
     &            IT32E,SI32E,S2N32E,SWT32E,SCH32E)
C
C---ELASTIC TAIL COMPTON PART
      ISET=1
      CALL HSWRSA(ISET,NRG33E,NDM3EL,
     &            SIG33L,SG33EE,T33GME,NDO33E,T33MXE,NM33E,XX33E,
     &            FFG33E,DNG33E,FFL33E,DNL33E,GLD33E,
     &            NTT33E,NCL33E,NC133E,NC233E,IBM33E,JCR33E,
     &            LGL33E,LLC33E,
     &            IT33E,SI33E,S2N33E,SWT33E,SCH33E)
C
      WRITE(LUNOUT,'(/2A,I3)')
     &       ' *** SAMPLING INFORMATION FOR ALL CHANNELS',
     &       ' WRITTEN ONTO UNIT LUNDAT=',LUNDAT
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSTPAR
C---
C   TEST CURRENTLY DEFINED RUN PARAMETERS
C   AGAINST THE INPUT READ FROM UNIT LUNDAT
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C
      PARAMETER(TESACC=1D-10)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD,IDPVR
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSISGM/ TCUTQ,TCUTQS
      DIMENSION LPARIT(12)
C
      READ(LUNDAT,ERR=9901,END=9901,IOSTAT=IOS)
     &    TXMIN,TXMAX,TQ2MIN,TYMIN,TYMAX,TWMIN,ICUTT
      READ(LUNDAT,ERR=9901,END=9901,IOSTAT=IOS) TTCUTQ,TTCTQS
      READ(LUNDAT,ERR=9901,END=9901,IOSTAT=IOS)
     &      (LPARIT(I),I=1,12)
      READ(LUNDAT,ERR=9901,END=9901,IOSTAT=IOS)
     &    ICODET,ILIBT,ILQMODT
      READ(LUNDAT,ERR=9901,END=9901,IOSTAT=IOS)
     &    TEELE,TEPRO,TPOLAR,LLEPTT
      READ(LUNDAT,ERR=9901,END=9901,IOSTAT=IOS)  TEGMIN
C
      IF((ABS(TXMIN-XMIN).GT.TESACC)
     &    .OR. (ABS(TXMIN-XMIN).GT.TESACC)
     &    .OR. (ABS(TXMAX-XMAX).GT.TESACC)
     &    .OR. (ABS(TQ2MIN-Q2MIN).GT.TESACC)
     &    .OR. (ABS(TYMIN-YMIN).GT.TESACC)
     &    .OR. (ABS(TYMAX-YMAX).GT.TESACC)
     &    .OR. (ICUTT.NE.ICUT))                 THEN
        WRITE(LUNOUT,'(/2A,I3//A/)')
     &      ' CURRENT KINEMATICAL CUTS INCONSISTENT WITH THE DATA SET',
     &      ' TO BE READ FROM UNIT',LUNDAT,
     &      ' CUTS DEFINED IN THE DATA SET:'
        WRITE(LUNOUT,'(10X,A,1PE10.3,5X,A,1PE10.3)')
     *     ' XMIN=',TXMIN,' XMAX=',TXMAX
        WRITE(LUNOUT,'(10X,A,1PE10.3)')
     *     ' Q2MIN=',TQ2MIN
        WRITE(LUNOUT,'(10X,A,1PE10.3,26X,A)')
     *     ' WMIN=', TWMIN, ' (ACTIVE ONLY FOR ICUT>1)'
        WRITE(LUNOUT,'(10X,A,1PE10.3,5X,A,1PE10.3,5X,A)')
     *     ' YMIN=',TYMIN,' YMAX=',TYMAX, ' (ACTIVE ONLY FOR ICUT=3)'
        WRITE(LUNOUT,'(10X,A,I2)')
     *     ' ICUT=',ICUTT
C
        WRITE(LUNOUT,'(/A/)')  ' CURRENTLY DEFINED PARAMETERS:'
        WRITE(LUNOUT,'(10X,A,1PE10.3,5X,A,1PE10.3)')
     *     ' XMIN=',XMIN,' XMAX=',XMAX
        WRITE(LUNOUT,'(10X,A,1PE10.3)')
     *     ' Q2MIN=',Q2MIN
        WRITE(LUNOUT,'(10X,A,1PE10.3,26X,A)')
     *     ' WMIN=', WMIN, ' (ACTIVE ONLY FOR ICUT>1)'
        WRITE(LUNOUT,'(10X,A,1PE10.3,5X,A,1PE10.3,5X,A)')
     *     ' YMIN=',YMIN,' YMAX=',YMAX, ' (ACTIVE ONLY FOR ICUT=3)'
        WRITE(LUNOUT,'(10X,A,I2)')
     *     ' ICUT=',ICUT
        WRITE(LUNOUT,100)
        STOP
      ENDIF
C
      IF(INT3(4).GT.1 .OR. ISAM3(4).GE.1) THEN
        IF(ABS(TTCUTQ-TCUTQ).GT.TESACC
     &      .OR. ABS(TTCTQS-TCUTQS).GT.TESACC) THEN
          WRITE(LUNOUT,'(/2A,I3//A/)')
     &      ' CURRENT ANGULAR CUTS FOR QUARKONIC BREMSSTRAHLUNG',
     &      ' INCONSISTENT WITH THE DATA SET TO BE READ FROM UNIT',
     &        LUNDAT,
     &      ' CUTS DEFINED IN THE DATA SET:'
          WRITE(LUNOUT,'(10X,A,1PE10.3,5X,A,1PE10.3)')
     &     ' TCUTQ=',TTCUTQ,' TCUTQS=',TTCTQS
          WRITE(LUNOUT,'(/A/)')  ' CURRENTLY DEFINED PARAMETERS:'
          WRITE(LUNOUT,'(10X,A,1PE10.3,5X,A,1PE10.3)')
     &     ' TCUTQ=',TCUTQ,' TCUTQS=',TCUTQS
          WRITE(LUNOUT,100)
          STOP
        ENDIF
      ENDIF
C
      DO 11 I=1,12
        IF(LPARIN(I).NE.LPARIT(I)) GOTO 12
  11  CONTINUE
      GOTO 13
  12  CONTINUE
      WRITE(LUNOUT,'(//2A,I3//A,I3,A)')
     &  ' CURRENT GSW-PARAMETERS INCONSISTENT WITH DATA SET ',
     &  ' TO BE READ FROM UNIT ',LUNDAT,
     &  ' PARAMETERS FROM DATA SET ON UNIT ',LUNDAT, ' :'
      WRITE(LUNOUT,'(12I3)') LPARIT
      WRITE(LUNOUT,'(/A/)')  ' CURRENTLY DEFINED PARAMETERS:'
      WRITE(LUNOUT,'(12I3)') LPARIN
      WRITE(LUNOUT,100)
      STOP
C
  13  CONTINUE
C
      IF(ICODET.NE.ICODE
     &  .OR.ILIBT.NE.ILIB
     &  .OR.ILQMODT.NE.ILQMOD) THEN
        WRITE(LUNOUT,'(//2A,I3//A,I3,A)')
     &    ' OTHER PARTON DISTRIBUTIONS USED FOR THE DATA SET ',
     &    ' TO BE READ FROM UNIT',LUNDAT,
     &    ' PARAMETER FROM DATA SET ON UNIT',LUNDAT,' : '
        WRITE(LUNOUT,'(10X,A,I7,A,I3,A,I3)')
     &      ' ICODE = ',ICODET,' ILIB = ',ILIBT,' ILQMOD = ',ILQMODT
        WRITE(LUNOUT,'(/A)')
     &    ' CURRENTLY DEFINED PARAMETERS:'
        WRITE(LUNOUT,'(10X,A,I7,A,I3,A,I3)')
     &      ' ICODE = ',ICODE,' ILIB = ',ILIB,' ILQMOD = ',ILQMOD
C
        WRITE(LUNOUT,100)
        STOP
      ENDIF
C
      IF((ABS(TEELE-EELE).GT.TESACC)
     &    .OR. (ABS(TEPRO-EPRO).GT.TESACC)
     &    .OR. (ABS(TPOLAR-POLARI).GT.TESACC)
     &    .OR. (LLEPTT.NE.LLEPT))                 THEN
        WRITE(LUNOUT,'(/2A,I3//A,I3/)')
     &      ' CURRENT BEAM PARAMETERS INCONSISTENT WITH THE DATA SET',
     &      ' TO BE READ FROM UNIT',LUNDAT,
     &      ' PARAMETERS FROM DSTS SET ON UNIT',LUNDAT,' :'
        WRITE(LUNOUT,'(10X,A,F8.1,A)')
     &      ' ENERGY OF INCIDENT ELECTRON =',TEELE,' GEV'
        WRITE(LUNOUT,'(10X,A,I3)')
     &                  ' CHARGE OF INCIDENT ELECTRON =',LLEPTT
        WRITE(LUNOUT,'(10X,A,F8.4)')
     &      ' DEGREE OF BEAM POLARIZATION =', TPOLAR
        WRITE(LUNOUT,'(10X,A,F8.1,A)')
     &      ' ENERGY OF INCIDENT PROTON =',TEPRO,' GEV'
C
        WRITE(LUNOUT,'(/A/)')  ' CURRENTLY DEFINED PARAMETERS:'
        WRITE(LUNOUT,'(10X,A,F8.1,A)')
     &      ' ENERGY OF INCIDENT ELECTRON =',EELE,' GEV'
        WRITE(LUNOUT,'(10X,A,I3)')
     &                  ' CHARGE OF INCIDENT ELECTRON =',LLEPT
        WRITE(LUNOUT,'(10X,A,F8.4)')
     &      ' DEGREE OF BEAM POLARIZATION =', POLARI
        WRITE(LUNOUT,'(10X,A,F8.1,A)')
     *      ' ENERGY OF INCIDENT PROTON =',EPRO,' GEV'
C
        WRITE(LUNOUT,100)
        STOP
      ENDIF
C
      IF(ABS(TEGMIN-EGMIN).GT.TESACC) THEN
        WRITE(LUNOUT,'(//2A,I3//A,I3,A,1PE11.3,A)')
     &      ' CURRENT CUT ON PHOTON ENERGY INCONSISTENT WITH DATA SET ',
     &       ' TO BE READ FROM UNIT',LUNDAT,
     &       ' CUT FROM DATA SET ON UNIT',LUNDAT,' :  EGMIN=',
     &       TEGMIN,' GEV'
        WRITE(LUNOUT,'(/A,1PE11.3,A)')
     &       ' CURRENTLY DEFINED PARAMETER: EGMIN=',EGMIN,' GEV'
        WRITE(LUNOUT,100)
        STOP
      ENDIF
C
      RETURN
C
 9901 CONTINUE
      WRITE(LUNOUT,'(/A,I3/A,I3/A)')
     &    ' ***  ERROR IN HSTPAR READING DATA FROM UNIT ',LUNDAT,
     &    ' ***  IOS=',IOS,
     &    ' ***  EXECUTION STOPPED'
      WRITE(LUNOUT,'(/A/5X,A/5X,6(1PE12.4),I3/5X,A,I7,A,I3,A,I3/
     &               5X,A/5X,3(1PE12.4),I3/5X,A,1PE12.4)')
     &    ' ***  ACTUAL VALUES OF PARAMETERS TO BE READ:',
     &    ' TXMIN,TXMAX,TQ2MIN,TYMIN,TYMAX,TWMIN,ICUTT',
     &      TXMIN,TXMAX,TQ2MIN,TYMIN,TYMAX,TWMIN,ICUTT,
     &    ' ICODET=',ICODET,' ILIBT=', ILIBT,' ILQMODT=',ILQMODT,
     &    ' TEELE,TEPRO,TPOLAR,LLEPTT',
     &      TEELE,TEPRO,TPOLAR,LLEPTT,
     &    ' TEGMIN=', TEGMIN
      STOP
C
 100  FORMAT(/' ***  EXECUTION STOPPED IN SUBROUTINE HSTPAR')
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSWRPA
C---
C   WRITE PARAMETER DEFINITIONS OF THE CURRENT RUN
C   FOR ALL CHANNELS TO UNIT LUNDAT
C---
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     &                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD,IDPVR
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSISGM/ TCUTQ,TCUTQS
C----------------
      WRITE(LUNDAT) XMIN,XMAX,Q2MIN,YMIN,YMAX,WMIN,ICUT
      WRITE(LUNDAT) TCUTQ,TCUTQS
      WRITE(LUNDAT) (LPARIN(I),I=1,12)
      WRITE(LUNDAT) ICODE,ILIB,ILQMOD
      WRITE(LUNDAT) EELE,EPRO,POLARI,LLEPT
      WRITE(LUNDAT) EGMIN
C
      IF(IPRINT.GE.2) THEN
        WRITE(LUNTES,'(//A)') ' *** TEST PRINT HSWRPA'
        WRITE(LUNTES,'(A/A,I3)')
     &      ' *** PARAMETERS OF THE ACTUAL RUN',
     &      ' *** WRITTEN ONTO UNIT', LUNDAT
        WRITE(LUNTES,'(/A)') ' XMIN,XMAX,Q2MIN,YMIN,YMAX,WMIN,ICUT'
        WRITE(LUNTES,'(6(1PE12.3),I4)')
     &               XMIN,XMAX,Q2MIN,YMIN,YMAX,WMIN,ICUT
        WRITE(LUNTES,'(/A)') ' TCUTQ,TCUTQS'
        WRITE(LUNTES,'(2(1PE12.3))') TCUTQ,TCUTQS
        WRITE(LUNTES,'(/A)') ' (LPARIN(I),I=1,12)'
        WRITE(LUNTES,'(12I2)') (LPARIN(I),I=1,12)
        WRITE(LUNTES,'(/A,I3,A,I8,A,I3)')
     &    ' ICODE = ',ICODE,'ILIB = ',ILIB,' ILQMOD = ',ILQMOD
        WRITE(LUNTES,'(/A)') ' EELE,EPRO,POLARI,LLEPT'
        WRITE(LUNTES,'(3(1PE12.3),I5)') EELE,EPRO,POLARI,LLEPT
        WRITE(LUNTES,'(A,1PE12.3)') ' EGMIN=',EGMIN
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSWRSA(ISET,NREGX,NDIMX,
     &                  SIG,SIGE,TGMAX,NDOX,TLMAX,NM,XX,
     &                  FFGOX,DNCGX,FFLOX,DNCLX,GOLDX,
     &                  NTOTX,NCALX,NCA1X,NCA2X,IBIMX,JCORX,
     &                  LGLOX,LLOCX,
     &                  IT,SI,SI2,SWGT,SCHI)
C---
C   WRITE INFORMATION FOR SAMPLING FOR ONE CONTRIBUTION TO UNIT LUNDAT
C   ISET.GT.0 : RADIATIVE CHANNELS (INTEGRATION INFORMATION FROM VEGAS)
C---
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision ntotx
      LOGICAL LGLOX,LLOCX
C
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     &                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      DIMENSION TLMAX(NREGX),NM(NREGX),XX(NDOX,NDIMX)
C---------------------------------
      WRITE(LUNDAT) SIG,SIGE,TGMAX,NDOX
      IF(IPRINT.GE.6) THEN
        WRITE(LUNTES,'(/A,3I5)') ' HSWRSA:  ISET, NREGX, NDIMX',
     &                                      ISET, NREGX, NDIMX
        WRITE(LUNTES,'(/A)') ' HSWRSA:  SIG,SIGE,TGMAX,NDOX'
        WRITE(LUNTES,*) SIG,SIGE,TGMAX,NDOX
      ENDIF
      WRITE(LUNDAT) (TLMAX(I),I=1,NREGX)
      WRITE(LUNDAT) (NM(I),I=1,NREGX)
C
      IF(ISET.GT.0) WRITE(LUNDAT) IT,SI,SI2,SWGT,SCHI
      WRITE(LUNDAT) ((XX(I,II),I=1,NDOX),II=1,NDIMX)
      WRITE(LUNDAT) FFGOX,DNCGX,FFLOX,DNCLX,GOLDX,
     &              NTOTX,NCALX,NCA1X,NCA2X,IBIMX,JCORX,
     &              LGLOX,LLOCX
      IF(IPRINT.GE.6) THEN
        WRITE(LUNTES,'(/A/5X,5(1PD13.5))')
     &       ' HSWRSA :    FFGOX,DNCGX,FFLOX,DNCLX,GOLDX',
     &                    FFGOX,DNCGX,FFLOX,DNCLX,GOLDX
        WRITE(LUNTES,'(/A/5X,1PE25.16,I8,4I5)')
     &       ' HSWRSA :    NTOTX,NCALX,NCA1X,NCA2X,IBIMX,JCORX',
     &                    NTOTX,NCALX,NCA1X,NCA2X,IBIMX,JCORX
        WRITE(LUNTES,'(/A,2L4)')
     &       ' HSWRSA:     LGLOX,LLOCX',  LGLOX,LLOCX
      ENDIF
      IF(IPRINT.GE.6) THEN
         WRITE(LUNTES,'(/A)') ' HSWRSA: XX(NDOX,NDIMX)'
         WRITE(LUNTES,2) ((XX(I,J),I=1,NDOX),J=1,NDIMX)
         WRITE(LUNTES,'(/A)') ' HSWRSA: TLMAX(NREGX)'
         WRITE(LUNTES,2) (TLMAX(I),I=1,NREGX)
         WRITE(LUNTES,'(/,A)') ' HSWRSA: NM2(NDOX,NDIMX)'
         WRITE(LUNTES,3) (NM(I),I=1,NREGX)
      ENDIF
      RETURN
 2    FORMAT(5(1PD15.5))
 3    FORMAT(10I8)
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSINIT(IFUN,EPSO,NBIN2,NDO2,SIG2,SIG2E,XX2)
C
C   LAST CHANGE  18/04/90  HJM
C   LAST CHANGE  24/07/90  HJM
C   LAST CHANGE  13/12/90  HJM         EXTENSION FOR CHARGED CURRENT
C   Last change  22/01/97  HS          NC / CC corrected: FUN -> IFUN
C******************
C   INITIALIZATION FOR EP EVENT GENERATION
C   2 --> 2 PROCESS WITH SOFT AND VIRTUAL CORRECTIONS FROM H.SP.
C           (CROSS SECTION IN NANOBARN)
C
C   --> CONTRIBUTION SIG2 TO THE OVERALL CROSS SECTION
C
C   --> TABLE OF MAXIMUM FUNCTION VALUES IN THE NREG2 (X,G) INTERVALS
C       NUMBERING PROCEDURE AS IN HSESTM, HSGENM..
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      INTEGER MINPTS,MAXPTS
      PARAMETER (NXINT=100)
      PARAMETER (LENWRK=2000)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSINTL/ XL,XU
C
      DIMENSION XX2(50,2)
      DIMENSION BLOW(2),BUP(2)
      DIMENSION DFXT(NXINT),XXINT(NXINT)
      DATA BLOW/0D0,0D0/, BUP/1D0,1D0/
      DATA MINPTS,MAXPTS / 0, 1000/
C
C   PARAMETERS
C
      ACCURA=EPSO
      NDO2=NBIN2
C
C   INTEGRATION
C
      XL=XMIN
      DLOGX=LOG(XMAX/XMIN)/FLOAT(NXINT)
      SIG2=0D0
      SIG2E=0D0
      IF(IPRINT.GT.1)
     *  WRITE(LUNTES,'(A)') ' IX, XL, XU,  DFXT(IX), EPSF, EPSO'
C
      DO 1 IX=1,NXINT
        XU=XMIN*EXP(FLOAT(IX)*DLOGX)
        IF(XU.GT.XMAX) THEN
          XU=XMAX
        ENDIF
        XXINT(IX)=XL
        IFAIL=1
        NDIMEN=2
C       WRITE(LUNTES,'(/A/I2,4F5.1,2I6,1PE10.1,I6)')
C    &    ' NDIMEN,BLOW(2),BUP(2),MINPTS,MAXPTS,ACCURA,LENWRK',
C    &     NDIMEN,BLOW,BUP,MINPTS,MAXPTS,ACCURA,LENWRK
        CALL DX1FCF(BLOW,BUP,IFUN,ACCURA,ACCFIN,RESULT,IFAIL)
        DFXT(IX)=RESULT
        IF(IFAIL.NE.0.OR.IPRINT.GT.1)
     &    WRITE(LUNTES,'(A,I5,/,3(1PD15.6),A,I5)')
     &      ' D01FCF DID NOT MEET REQUIRED ACCURACY IN BIN ',IX,
     &                 XL,XU,RESULT,' IFAIL = ',IFAIL
        XL=XU
        SIG2=SIG2 + DFXT(IX)
        SIG2E=SIG2E + ACCFIN*RESULT
 1    CONTINUE
C
      IF(IPRINT.GT.1) THEN
        WRITE(LUNOUT,'(///A,5X,1PE12.4,A,1PE12.4,A)')
     *        ' CROSS SECTION VALUE SIG2 (WITH ERROR ESTIMATE):',
     *        SIG2, ' +/- ', SIG2E, '  NB'
        WRITE(LUNTES,'(A,1PD10.1)') ' RELATIVE ACCURACY REQUIRED:',
     &        ACCURA
      ENDIF

      if (ioplot.lt.0) return
C
      DO 2 IX=1,NXINT
        DFXT(IX)=DFXT(IX)/SIG2
 2    CONTINUE
C
      IF(IPRINT.GT.1) THEN
         WRITE(LUNTES,'(6D15.5)') (DFXT(I),I=1,NXINT)
      ENDIF
C
C   INVERSION OF THE X-DISTRIBUTION
C
      IF(IPRINT.GT.2)
     *  WRITE(LUNTES,'(2A)') ' JX,IX, G,DG, FI, DFXT(IX), DX,',
     *                   ' XMIN, XX2(JX)'
      G=0D0
      DG=1D0/DFLOAT(NDO2)
      IX=1
      FI=DFXT(1)
      DO 3 JX=1,NDO2-1
        G=G + DG
   4    CONTINUE
        IF(G.GT.FI) THEN
          IX=IX + 1
          FI=FI + DFXT(IX)
          GOTO 4
        ENDIF
        IF (IX.LT.100) THEN
          XXINT1=XXINT(IX+1)
          ELSE
          XXINT1=XMAX
        ENDIF
        DX=XXINT1-XXINT(IX)
        XX2(JX,1)=XXINT1-(FI-G)*DX/DFXT(IX)
        IF(IPRINT.GT.2) WRITE(LUNTES,'(2I5/7(1PD13.5))')
     *        JX,IX,G,DG,FI,DFXT(IX),DX,XMIN,XX2(JX,1)
 3    CONTINUE
      XX2(NDO2,1)=XMAX
C***     HJM 24/07/90                RESCALING OF X-INTERVALS,
C                                    SINCE HSTRIT ASSUMES INTERVAL (0,1)
      DX=XMAX-XMIN
      DO 6 I=1,NDO2
        XX2(I,1)=(XX2(I,1)-XMIN)/DX
  6   CONTINUE
C
      DO 5 I=1,NDO2
        XX2(I,2)=DFLOAT(I)*DG
  5   CONTINUE
C
      IF(IPRINT.GT.1) THEN
        WRITE(LUNTES,'(A,/,4(5(1PD15.5)/))')
     &                     ' XX2(I,1)', (XX2(IX,1),IX=1,NDO2)
        WRITE(LUNTES,'(A,/,4(5(1PD15.5)/))')
     &                     ' XX2(I,2)', (XX2(IG,2),IG=1,NDO2)
      ENDIF
      IF(IPRINT.GT.1) WRITE(LUNTES,'(A)') ' HSINIT FINISHED'
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSGLOW(X)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
C
      GSP=SP-MEI2-MPRO2
      IF(ICUT.EQ.1) THEN
C                                   CUT IN EXTERNALLY DEFINED Q**2(MIN)
        HSGLOW=-1D0/Q2MIN
      ELSEIF(ICUT.EQ.2) THEN
C                                   CUT IN W / MAXIMUM Q**2
        QQ2MIN=X*(WMIN**2-MPRO2)/(1D0 - X)
        HSGLOW=-1D0/MAX(Q2MIN,QQ2MIN)
      ELSEIF(ICUT.EQ.3) THEN
C                                   CUT IN Y AND W
        QQ2MIN=X*(WMIN**2-MPRO2)/(1D0 - X)
        QQQ2MN=X*YMIN*GSP
C                                   CUT ON ELECTRON SCATTERING ANGLE
C                                   (MASSES NEGLECTED)
        YMINTH=1D0/(1D0+X*CTHCON)
        QT2MIN=X*YMINTH*GSP
C                                   CUT ON ELECTRON TRANSVERSE MOMENTUM
C                                   (MASSES NEGLECTED)
        QP2MIN=X*GSP/2D0*(1D0-DSQRT(1D0-PTXM0/X))
        HSGLOW=-1D0/MAX(Q2MIN,QQ2MIN,QQQ2MN,QT2MIN,QP2MIN)
      ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSGLOW'
        STOP
      ENDIF
      GMIN=HSGLOW
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSGUPP(X)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSPCUT/ PTMIN,PTXM0
C
      GSP=SP-MEI2-MPRO2
      YMAXX=X*(1D0-4D0*MEI2*MPRO2/GSP/GSP)
     *        /(X*(1D0+X*MPRO2/GSP)+MEI2/GSP)
      IF(ICUT.EQ.1) THEN
C                                   CUT IN EXTERNALLY DEFINED Q**2(MIN)
        GMAX=-1D0/(X*YMAXX*GSP)
      ELSEIF(ICUT.EQ.2) THEN
C                                   CUT IN W / MAXIMUM Q**2
        GMAX=-1D0/(X*YMAXX*GSP)
      ELSEIF(ICUT.EQ.3) THEN
C                                   CUT IN Y
        Q2MAX1=X*MIN(YMAX,YMAXX)*GSP
C                                   CUT ON ELECTRON TRANSVERSE MOMENTUM
C                                   (MASSES NEGLECTED)
        QP2MAX=X*GSP/2D0*(1D0+DSQRT(1D0-PTXM0/X))
        GMAX=-1D0/MIN(Q2MAX1,QP2MAX,Q2MAX)
      ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSGUPP'
        STOP
      ENDIF
C
      HSGUPP=GMAX
      IF(HSGUPP.LT.GMIN) THEN
        HSGUPP=GMIN
      ENDIF
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSNCG1(NDIMEN,ARGUM)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      DOUBLE PRECISION HSNCG1
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSINTL/ XL,XU
      DIMENSION ARGUM(NDIMEN)
C
      DX=XU-XL
      X=XL+ARGUM(1)*DX
      GL=HSGLOW(X)
      GU=HSGUPP(X)
      DG=DMAX1(GU-GL,0D0)
      G=GL+ARGUM(2)*DG
      Q2=-1D0/G
      IF(IPRINT.GT.20) WRITE(LUNTES,'(A,4D15.6)')
     &                      ' HSNCG1: X, G, Q2',X,G,Q2
      HSNCG1=Q2**2*HSNC22(X,Q2)*DX*DG
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSNCG2(XARG)
C***
C   LAST CHANGE 24/07/90 BY HJM
C***
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      DOUBLE PRECISION HSNCG2
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      DIMENSION XARG(2)
C
      DX=XMAX-XMIN
      X=XMIN+XARG(1)*DX
      Z=XARG(2)
      GL=HSGLOW(X)
      GU=HSGUPP(X)
      DG=DMAX1(GU-GL,0D0)
      G=GL+Z*DG
      Q2=-1D0/G
      IF(IPRINT.GT.20) WRITE(LUNTES,'(A,4D15.6)')
     &                      ' HSNCG2: X, Z, G, Q2',X,Z,G,Q2
      HSNCG2=Q2**2*HSNC22(X,Q2)*DG*DX
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSNC22(X,Q2)
C
C   D2SIG/DX*DQ2 WITH COMPLETE 1-LOOP SOFT AND VIRTUAL CORRECTIONS
C   FROM H. SPIESBERGER
C
C        NOTE: OUTPUT FROM H.S.  DSIG / DX*DY
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
C---------------------------------------------------------------------
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKPXY/ XX,Y
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      DATA NEVERR /0/
      LOGICAL LERR
      DATA LERR /.TRUE./
C
      XX=X
      IF(IPRINT.GT.20)
     &  WRITE(LUNTES,'(A/3(1PD13.5),2(0PF8.3),2I3)')
     *         ' HSNC22: SP, X, Q2, POLARI, HPOLAR, LLEPT, LQUA',
     *                   SP, X, Q2, POLARI, HPOLAR, LLEPT, LQUA
      Y=Q2/X/(SP-MEI2-MPRO2)
      HSNC22=HSSGNC(X,Y,LLEPT,POLARI,HPOLAR,LQUA)/X/SP
C
      IF(HSNC22.LT.0D0.and.IOPLOT.ge.0) THEN
        HSNC22=0D0
        NEVERR=NEVERR+1
        IF (NEVERR.LT.10) THEN
         WRITE(LUNTES,'(A,/,4(1PD13.5),2I5,2(0PF8.3)/A/2(6(1PD13.5)/))')
     +     ' HSNC22: X, Y, Q2, HSNC22, LLEPT, LQUA, POLARI, HPOLAR',
     +               X, Y, Q2, HSNC22, LLEPT, LQUA, POLARI, HPOLAR,
     +     ' HSPDFQ: QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT',
     +               QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
         WRITE(LUNTES,'(2A,/,2(6(1PD13.5)/))')
     +     ' HSDPDF: DQU,DQBU,DQD,DQBD,DQS,DQBS,',
     +                                  'DQC,DQBC,DQB,DQBB,DQT,DQBT',
     +               DQU,DQBU,DQD,DQBD,DQS,DQBS,
     +                                   DQC,DQBC,DQB,DQBB,DQT,DQBT
        ELSEIF (LERR) THEN
         LERR=.FALSE.
         WRITE(LUNTES,'(A,I3,A)')
     &     ' ERROR HSNC22 < 0 HAS OCCURED ',NEVERR,
     &     ' TIMES, NO FURTHER WARNINGS ARE PRINTED'
        ELSE
        ENDIF
      ENDIF
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE VEGAS(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH,
     &                 NDO,IT,SI,SI2,SWGT,SCHI,XI)
C
C  SUBROUTINE PERFORMS N-DIMENSIONAL MONTE CARLO INTEG*N
C    - BY G.P.  LEPAGE  SEPT 1976/ (REV) APR 1978
C
C    -FTN5 VERSION 21-8-1984
C    -HBOOK/HPLOT INTERFACE 6-1-1985
C
C    - LAST MODIFICATIONS HJM 7/12/89
C
C  AUTHOR                                       : G. P. LEPAGE
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL FXN
C*    REAL ZZF1,ZZW
C
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      PARAMETER(NDIMX=5)
C
      DIMENSION XIN(50),R(50),DX(10),IA(10),KG(10),DT(10)
      DIMENSION XL(NDIMX),XU(NDIMX),QRAN(NDIMX),X(NDIMX)
      COMMON /VGASIO/ NINP,NOUTP
C*    COMMON /VGB2/   NDO,IT,SI,SI2,SWGT,SCHI,XI(50,NDIMX),SCALLS,
C*   *            D(50,NDIMX),DI(50,NDIMX)
      DIMENSION XI(50,NDIMX),D(50,NDIMX),DI(50,NDIMX)
C***
      COMMON /VGRES/  S1,S2,S3,S4
C*    REAL S1,S2,S3,S4
C
      DATA XL,XU/NDIMX*0D0,NDIMX*1D0/
      DATA NDMX/50/,ALPH/1.5D0/,ONE/1D0/,MDS/1/
C
C*    CALL VGDAT
C
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)' VEGAS CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
      IF(NPRN.GT.0)THEN
         IPR=0
      ELSE
         IPR=1
      ENDIF
      NDO=1
      DO 1 J=1,NDIM
         XI(1,J)=ONE
1     CONTINUE
C
      ENTRY VEGAS1(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH,
     &                 NDO,IT,SI,SI2,SWGT,SCHI,XI)
C...INITIALIZES CUMMULATIVE VARIABLES,BUT NOT GRID
C*    CALL VGDAT
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)' VEGAS1 CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
C*    IF(MOD(IGRAPH,100) .GT. 0 )THEN
C*       IF(MOD(IGRAPH,10) .GT. 0)THEN
C*          NOW=MOD(IGRAPH,100)/10
C*       ELSE
C*          NOW=MOD(IGRAPH,10)
C*       ENDIF
C*       ZZF1  = SNGL(F1)
C*       ZZW   = SNGL(W)
C*       CALL INPLOT(NOW,ZZF1,ZZW)
C*       CALL INPLOT(NOW,F1,W)
C*    ENDIF
C*    IF(IGRAPH .GE. 100) THEN
C*       NOW=MOD(IGRAPH,1000)/100
C*       CALL HVBOOK(NOW,F1,W)
C*    ENDIF
C
      IT=0
      SI=0.
      SI2=SI
      SWGT=SI
      SCHI=SI
C*    SCALLS=SI
C
      ENTRY VEGAS2(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH,
     &                 NDO,IT,SI,SI2,SWGT,SCHI,XI)
C...NO INITIALIZATION
C*    CALL VGDAT
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)' VEGAS2 CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
C***                                 HJM 14/12/89
      ITMXX=ITMX + IT
C***
      ND=NDMX
      NG=1
      IF(MDS.NE.0)THEN
         NG=INT((NCALL/2.)**(1./NDIM))
         MDS=1
         IF((2*NG-NDMX).GE.0)THEN
            MDS=-1
            NPG=NG/NDMX+1
            ND=NG/NPG
            NG=NPG*ND
         ENDIF
      ENDIF
C
      K=NG**NDIM
      NPG=NCALL/K
      IF(NPG.LT.2) NPG=2
      CALLS=NPG*K
      DXG=ONE/NG
      DV2G=DXG**(2*NDIM)/NPG/NPG/(NPG-ONE)
      XND=ND
      NDM=ND-1
      DXG=DXG*XND
      XJAC=ONE
      DO 3 J=1,NDIM
         DX(J)=XU(J)-XL(J)
         XJAC=XJAC*DX(J)
3     CONTINUE
C
C  REBIN PRESERVING BIN DENSITY
C
      IF(ND.NE.NDO)THEN
         RC=NDO/XND
         DO 7 J=1,NDIM
            K=0
            XN=0D0
            DR=XN
            I=K
4           K=K+1
            DR=DR+ONE
            XO=XN
            XN=XI(K,J)
5           IF(RC.GT.DR) GO TO 4
            I=I+1
            DR=DR-RC
            XIN(I)=XN-(XN-XO)*DR
            IF(I.LT.NDM) GO TO 5
            DO 6  I=1,NDM
               XI(I,J)=XIN(I)
6           CONTINUE
            XI(ND,J)=ONE
7        CONTINUE
         NDO=ND
      ENDIF
C
       IF(NPRN.NE.0.AND.NPRN.NE.10)WRITE(NOUTP,200)NDIM,CALLS,IT,ITMX
     * ,ACC,MDS,ND
       IF(NPRN.EQ.10)WRITE(NOUTP,290)NDIM,CALLS,ITMX,ACC,MDS,ND
C
      ENTRY VEGAS3(FXN,ACC,NDIM,NCALL,ITMX,NPRN,IGRAPH,
     &                 NDO,IT,SI,SI2,SWGT,SCHI,XI)
C     - MAIN INTEGRATION LOOP
      IF(ITMX.LE.0)THEN
         WRITE(NOUTP,199)' VEGAS3 CALLED WITH AT MAX LESS EQUAL ZERO'//
     +                   ' ITERATIONS. NO EXECUTION.'
         RETURN
      ENDIF
9     CONTINUE
      IT=IT+1
      TI=0.
      TSI=TI
C*    IF(MOD(IGRAPH,100) .GT. 0 )THEN
C*       ZZF1  = SNGL(F1)
C*       ZZW   = SNGL(W)
C*       CALL REPLOT(NOW,ZZF1,ZZW)
C*       CALL REPLOT(NOW,F1,W)
C*    ENDIF
C*    IF(IGRAPH .GE. 100) THEN
C*       CALL HVRSET(NOW,F1,W)
C*    ENDIF
C
      DO 101 J=1,NDIM
         KG(J)=1
         DO 102 I=1,ND
            D(I,J)=TI
            DI(I,J)=TI
 102     CONTINUE
 101  CONTINUE
C
11    FB=0D0
      F2B=FB
      K=0
C
12    CONTINUE
      K=K+1
      DO 121 J=1,NDIM
C*       QRAN(J)=VGRAN(0.0)
         QRAN(J)=HSRNDM(-1)
121   CONTINUE
      WGT=XJAC
      DO 15 J=1,NDIM
         XN=(KG(J)-QRAN(J))*DXG+ONE
         IA(J)=INT(XN)
         IAJ=IA(J)
         IAJ1=IAJ-1
         IF(IAJ.LE.1)THEN
            XO=XI(IAJ,J)
            RC=(XN-IAJ)*XO
         ELSE
            XO=XI(IAJ,J)-XI(IAJ1,J)
            RC=XI(IAJ1,J)+(XN-IAJ)*XO
         ENDIF
         X(J)=XL(J)+RC*DX(J)
         WGT=WGT*XO*XND
15    CONTINUE
C
      F=FXN(X)*WGT
      F1=F/CALLS
      W=WGT/CALLS
C*    IF(MOD(IGRAPH,100) .GT. 0 )THEN
C*       ZZF1  = SNGL(F1)
C*       ZZW   = SNGL(W)
C*       CALL XPLOT(NOW,ZZF1,ZZW)
C*       CALL XPLOT(NOW,F1,W)
C*    ENDIF
C*    IF(IGRAPH .GE. 100) THEN
C*       CALL HVFILL(NOW,F1,W)
C*    ENDIF
C
      F2=F*F
      FB=FB+F
      F2B=F2B+F2
      DO 16 J=1,NDIM
         IAJ=IA(J)
         DI(IAJ,J)=DI(IAJ,J)+F/CALLS
         IF(MDS.GE.0)  D(IAJ,J)=D(IAJ,J)+F2
16    CONTINUE
      IF(K.LT.NPG) GO TO 12
C
      F2B=F2B*NPG
      F2B=DSQRT(F2B)
      F2B=DABS((F2B-FB)*(F2B+FB))
      TI=TI+FB
      TSI=TSI+F2B
      IF(MDS.LT.0)THEN
         DO 17 J=1,NDIM
            IAJ=IA(J)
            D(IAJ,J)=D(IAJ,J)+F2B
17       CONTINUE
      ENDIF
      K=NDIM
19    KG(K)=MOD(KG(K),NG)+1
      IF(KG(K).NE.1) GO TO 11
      K=K-1
      IF(K.GT.0) GO TO 19
C
C  FINAL RESULTS FOR THIS ITERATION
C
      TI=TI/CALLS
      TSI=TSI*DV2G
      TI2=TI*TI
      IF(TSI .EQ. 0D0)THEN
         WGT = 0D0
      ELSE
         WGT=TI2/TSI
      ENDIF
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SCHI=SCHI+TI2*WGT
      IF(SWGT .EQ. 0D0)THEN
         AVGI=TI
      ELSE
         AVGI=SI/SWGT
      ENDIF
      IF(SI2 .EQ. 0D0)THEN
         SD=TSI
      ELSE
         SD=SWGT*IT/SI2
      ENDIF
C*    SCALLS=SCALLS+CALLS
      CHI2A=0D0
      IF(IT.GT.1)CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-1)
      IF(SD .NE. 0D0)THEN
         SD=ONE/SD
         SD=SQRT(SD)
      ELSE
         SD=TSI
      ENDIF
      IF(NPRN.NE.0)THEN
         TSI=SQRT(TSI)
C        IF(NPRN.NE.10)WRITE(NOUTP,201)IPR,IT,TI,TSI,AVGI,SD,CHI2A
         IF(NPRN.NE.10)WRITE(NOUTP,201)    IT,TI,TSI,AVGI,SD,CHI2A
         IF(NPRN.EQ.10)WRITE(NOUTP,203)IT,TI,TSI,AVGI,SD,CHI2A
         IF(NPRN.LT.0)THEN
            DO 20 J=1,NDIM
               WRITE(NOUTP,202)J
               WRITE(NOUTP,204)(XI(I,J),DI(I,J),D(I,J),I=1,ND)
20         CONTINUE
         ENDIF
      ENDIF
C
C   REFINE GRID
C
      IF(SD .NE. 0D0)THEN
         REL = DABS(SD/AVGI)
      ELSE
         REL = 0D0
      ENDIF
      IF(REL.LE.DABS(ACC).OR.IT.GE.ITMX)NOW=2
      S1=AVGI
      S2=SD
      S3=TI
      S4=TSI
C*    IGRPH=MOD(IGRAPH,100)
C*    IF(IGRPH .GT. 0 .AND. IGRPH .LT. 10)THEN
C*       ZZF1  = SNGL(F1)
C*       ZZW   = SNGL(W)
C*       CALL PLOTIT(NOW,ZZF1,ZZW)
C*       CALL PLOTIT(NOW,F1,W)
C*    ELSE IF(IGRPH .GE. 10 )THEN
C*       ZZF1  = SNGL(F1)
C*       ZZW   = SNGL(W)
C*       CALL PLOTTA(NOW,ZZF1,ZZW)
C*       CALL PLOTTA(NOW,F1,W)
C*    ENDIF
C*    IF(IGRAPH .GE. 100) THEN
C*       CALL HVEDIT(NOW,F1,W)
C*    ENDIF
C
      DO 23 J=1,NDIM
         XO=D(1,J)
         XN=D(2,J)
         D(1,J)=(XO+XN)/2D0
         DT(J)=D(1,J)
         DO 22 I=2,NDM
            D(I,J)=XO+XN
            XO=XN
            XN=D(I+1,J)
            D(I,J)=(D(I,J)+XN)/3D0
            DT(J)=DT(J)+D(I,J)
22       CONTINUE
         D(ND,J)=(XN+XO)/2D0
         DT(J)=DT(J)+D(ND,J)
23    CONTINUE
C
      DO 28 J=1,NDIM
         RC=0D0
         DO 24 I=1,ND
            R(I)=0D0
            IF(D(I,J).GT.0D0)THEN
               XO=DT(J)/D(I,J)
               R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
            ENDIF
            RC=RC+R(I)
24       CONTINUE
         RC=RC/XND
         K=0
         XN=0D0
         DR=XN
         I=K
25       K=K+1
         DR=DR+R(K)
         XO=XN
         XN=XI(K,J)
26       IF(RC.GT.DR) GO TO 25
         I=I+1
         DR=DR-RC
         IF(DR .EQ. 0D0)THEN
            XIN(I)=XN
         ELSE
            XIN(I)=XN-(XN-XO)*DR/R(K)
         ENDIF
         IF(I.LT.NDM) GO TO 26
         DO 27 I=1,NDM
            XI(I,J)=XIN(I)
27       CONTINUE
         XI(ND,J)=ONE
28    CONTINUE
C
C***  IF(IT.LT.ITMX.AND.DABS(ACC).LT.REL)GO TO 9        HJM 14/12/89
      IF(IT.LT.ITMXX.AND.DABS(ACC).LT.REL)GO TO 9
C
      S1=AVGI
      S2=SD
      S3=CHI2A
      RETURN
C
199   FORMAT(A)
200   FORMAT(' INPUT PARAMETERS FOR VEGAS   NDIM=',I3,
     +       '  NCALL=',F8.0/28X,'  IT=',I5,'  ITMX =',I5/28X,
     +       '  ACC=',1PD11.3 / 28X,'  MDS=',I3,'   ND=',I4//)
290   FORMAT(' VEGAS  NDIM=',I3,'  NCALL=',F8.0,'  ITMX =',I5,
     +       '  ACC=',1PD11.3,'  MDS=',I3,'   ND=',I4)
C201   FORMAT(/I1,'INTEGRATION BY VEGAS'/' ITERATION NO',I3,
201   FORMAT(' INTEGRATION BY VEGAS'/' ITERATION NO',I3,
     +       '.   INTEGRAL =',1PD16.8/20X,'STD DEV  =',1PD12.4/
     +       ' ACCUMULATED RESULTS.   INTEGRAL =',1PD16.8/
     +   24X,'STD DEV  =',1PD12.4 / 24X,'CHI**2 PER ITN   =',1PD12.4)
202   FORMAT(' DATA FOR AXIS',I2 / 7X,'X',7X,'  DELT I  ',
     +       2X,' CONVCE    ',11X,'X',7X,'  DELT I  ',2X,' CONVCE     '
     +       ,11X,'X',7X,'  DELT I  ',2X,' CONVCE     '/)
204   FORMAT(1X,3D12.4,5X,3D12.4,5X,3D12.4)
203   FORMAT(1X,I3,D20.8,D12.4,D20.8,D12.4,D12.4)
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSTRIT(F,X,NDIM,NCONT,XI,NDO)
C
C  AUTHOR      : J. VERMASEREN
C                MODIFIED 2/88 BY HJM
C                LAST CHANGE 23/07/90  BY HJM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION HSTRIT
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      DIMENSION X(NDIM),XI(50,NDIM)
      DIMENSION Z(10),R(20),NCALL(20)
      DATA NCALL /20*0/
      EXTERNAL F
C
      IF(NCALL(NCONT).EQ.0) THEN
        NCALL(NCONT)=1
        R(NCONT)=NDO**NDIM
      ENDIF
C
      W=R(NCONT)
      DO 4 I=1,NDIM
         XX=X(I)*NDO
         J=INT(XX)
         JJ=J+1
         Y=XX-J
         IF(J.LE.0)THEN
            DD=XI(1,I)
         ELSE
            DD=XI(JJ,I)-XI(J,I)
         ENDIF
         Z(I)=XI(JJ,I)-DD*(1.-Y)
         W=W*DD
4     CONTINUE
C
      FZ=F(Z)
      HSTRIT=W*FZ
C
      IF(HSTRIT.LE.0D0.AND.IPRINT.GE.5) THEN
        WRITE(LUNOUT,'(A,5(1PE12.4)/15X,5(1PE12.4))')
     *        ' HSTRIT - Z(I):',Z
        WRITE(LUNOUT,'(A,1PE12.4)')  ' HSTRIT -  W  :',W
        WRITE(LUNOUT,'(A,1PE12.4)')  ' HSTRIT - F(Z):',FZ
        WRITE(LUNOUT,'(A/50(10(1PE12.4)/))')  ' HSTRIT - XI  :',XI
      ENDIF
C
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSESTM(F,NCONT,NDIM,NPOIN,NDO,MBIN,FFMAX,FMAX,XI,
     &                  IBIMAX,NREGX)
C
C  AUTHOR      : J. VERMASEREN
C                MODIFIED 2/88 HJM
C                MODIFIED 2/97 HS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL F
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      DIMENSION X(10),N(10)
      DIMENSION FMAX(NREGX),XI(NDO,NDIM)

      IF (IOPLOT.LT.0) RETURN

      FFMAX=0D0
      DO 5 J=1,NREGX
        FMAX(J)=0D0
 5    CONTINUE
C
      SUM=0D0
      SUM2=0D0
      SUM2P=0D0
      IF(IPRINT.GE.1) WRITE(LUNTES,200) MBIN,NREGX,NPOIN

C...DETERMINATION OF GLOBAL/LOCAL MAXIMA
      DO 1 J=1,NREGX
         JJ=J-1
         DO 2 K=1,NDIM
            JJJ=JJ/MBIN
            N(K)=JJ-JJJ*MBIN
            JJ=JJJ
2        CONTINUE
         FSUM=0D0
         FSUM2=0D0
         DO 3 M=1,NPOIN
            DO 4 K=1,NDIM
               X(K)=(HSRNDM(-1)+N(K))/MBIN
4           CONTINUE
            Z=HSTRIT(F,X,NDIM,NCONT,XI,NDO)
            IF(Z.GT.FMAX(J)) FMAX(J)=Z
            FSUM=FSUM + Z
            FSUM2=FSUM2 + Z*Z
3        CONTINUE
         IF(FMAX(J).GT.FFMAX) THEN
           FFMAX=FMAX(J)
           IBIMAX=J
         ENDIF
C
         AV=FSUM/FLOAT(NPOIN)
         AV2=FSUM2/FLOAT(NPOIN)
         SIG2=AV2-AV*AV
         SIG=SQRT(DABS(SIG2))
         SUM=SUM+AV
         SUM2=SUM2+AV2
         SUM2P=SUM2P+SIG2
         EFF=10000.
         IF(FMAX(J).NE.0) EFF=FMAX(J)/AV
         IF(IPRINT.GE.5) WRITE(LUNTES,100) J,AV,SIG,FMAX(J),EFF,
     +                                 (N(KJ),KJ=1,NDIM)
1     CONTINUE
      SUM=SUM/NREGX
      SUM2=SUM2/NREGX
      SUM2P=SUM2P/NREGX
      SIG=SQRT(SUM2-SUM*SUM)
      SIGP=SQRT(SUM2P)
      EFF1=0D0
      DO 6 J=1,NREGX
         EFF1=EFF1+FMAX(J)
6     CONTINUE
      EFF1=EFF1/(NREGX*SUM)
      EFF2=FFMAX/SUM
      IF(IPRINT.GE.1) THEN
        WRITE(LUNTES,101)SUM,SIG,SIGP,FFMAX,EFF1,EFF2
        WRITE (LUNTES,'(/A,1PE13.5,5X,A,I5)')
     &        ' GLOBAL MAXIMUM FFMAX=',FFMAX,' IN BIN IBIMAX=',IBIMAX
      ENDIF
C
100   FORMAT(I6,3X,G13.6,G12.4,G13.6,3X,F8.2,3X,10I2)
101   FORMAT(' THE AVERAGE FUNCTION VALUE =',G14.6/
     +       ' THE OVERALL STD DEV        =',G14.4/
     +       ' THE AVERAGE STD DEV        =',G14.4/
     +       ' THE MAXIMUM FUNCTION VALUE =',G14.6/
     +       ' THE AVERAGE INEFFICIENCY   =',G14.3/
     +       ' THE OVERALL INEFFICIENCY   =',G14.3/)
200   FORMAT(' SUBROUTINE HSESTM USES A',I3,'**NDIM DIVISION'/
     + ' THIS RESULTS IN ',I7,' CUBES'/
     + ' THE PROGRAM PUT ',I5,' POINTS IN EACH CUBE TO FIND',
     + ' STARTING VALUES FOR THE MAXIMA'//)
C
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSGENM(F,NCONT,NDIM,NEVENT,ICONTI,FFMAX,FMAX,GLAST,
     &                  FFMOLD,FMLOLD,DNCGLO,DNCLOC,LLOCAL,LGLOB,
     &                  NTOT,NM,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,
     &                  XI,NDO,MBIN,NREG)
C
C  ORIGINAL AUTHOR      : J. VERMASEREN
C                MODIFIED 2/88 HJM
C                LAST CHANGE   7/12/89  HJM
C                LAST CHANGE  19/03/90  HJM
C                MAJOR MODIFICATION  12/04/90  HJM
C                LAST CHANGE  23/07/90  HJM
C                MODIFIED 12/02/97 HS
C
C---------------------------------------------------------------------
C
C   VARIABLES:
C   **********
C
C   F :          ACTUAL DISTRIBUTION FUNCTION FOR SAMPLING
C
C   NCONT :      ACTUAL NUMBER OF THE CORRESPONDING CONTRIBUTION
C                TO THE MATRIX ELEMENT
C
C   NDIM  :      DIMENSION OF THE COORDINATE VECTOR
C
C   NEVENT:      NUMBER OF EVENTS TO BE SAMPLED
C
C   ICONTI = 1   NEW SEQUENCE OF EVENTS TO BE SAMPLED,
C                I.E. INFORMATION FROM POTENTIAL PREVIOUS SAMPLING RUNS
C                     IS TO BE DISCARDED
C          = 2   EXTENSION OF AN EXISTING EVENT SAMPLE
C
C   FFMAX :      ESTIMATED GLOBAL MAXIMUM
C
C   FMAX  :      ESTIMATES FOR LOCAL MAXIMA (FIELD)
C
C   GLAST :      LAST CALCULATED FUNCTION VALUE
C                (NEEDED FOR CONTINUED SAMPLING, ICONTI=1)
C
C   NTOT  :      EFFECTIVE TOTAL NUMBER OF TRIALS FOR SAMPLING
C
C   NM(I) :      NUMBER OF TRIALS IN REGION I
C
C   NCALL :      TOTAL NUMBER OF ACTUAL FUNCTION EVALUATIONS
C
C   MCALL1 :     NUMBER OF FUNCTION EVALUATIONS DURING CORRECTION
C                PROCEDURE FOR WRONG LOCAL MAXIMUM
C
C   MCALL2 :     NUMBER OF FUNCTION EVALUATIONS DURING CORRECTION
C                PROCEDURE FOR WRONG GLOBAL MAXIMUM
C
C   IBIMAX :     NUMBER OF THE BIN CONTAINING THE CURRENT GLOBAL MAXIMUM
C
C   FFMOLD :     GLOBAL MAXIMUM FROM LAST SAMPLING RUN
C                (NOT YET CORRECTED EVEN IF NECESSARY!)
C
C   FMLOLD :     LOCAL MAXIMUM FOR THE BIN SELECTED FOR THE LAST EVENT
C                IN THE PREVIOUS SAMPLING RUN
C                (NOT YET CORRECTED EVEN IF NECESSARY!)
C
C   DNCGLO :     NUMBER OF EVENTS TO BE ADDITIONALLY SAMPLED
C                TO CONTINUE THE CORRECTION PROCEDURE
C                FOR WRONG ESTIMATION OF GLOBAL MAXIMUM
C                (INFORMATION FROM POTENTIAL PREVIOUS RUN)
C
C   DNCLOC :     NUMBER OF EVENTS TO BE ADDITIONALLY SAMPLED
C                IN REGION **JCOR** TO CONTINUE THE CORRECTION PROCEDURE
C                FOR WRONG ESTIMATION OF LOCAL MAXIMUM
C                (INFORMATION FROM POTENTIAL PREVIOUS RUN)
C
C   LGLOB :      (LOGICAL) IF TRUE CORRECTION FOR WRONG ESTIMATE OF
C                GLOBAL MAXIMUM NECESSARY
C
C   LLOCAL :     (LOGICAL) IF TRUE CORRECTION FOR WRONG ESTIMATE OF
C                LOCAL MAXIMUM NECESSARY
C
C   XI    :      ACTUAL SUBDIVISION OF NDIM-DIMENSIONAL COORDINATE SPACE
C                AS PROVIDED BY VEGAS
C
C   NDO   :      NUMBER OF INTERVALS PER COORDINATE AXIS FROM VEGAS
C                FROM VEGAS-SUBDIVISION
C
C   MBIN  :      NUMBER OF BINS PER AXIS FOR SUBVOLUMES USED IN EVENT
C                SAMPLING TO INCREASE THE SAMPLING EFFICIENCY
C
C---------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL F
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      PARAMETER (NDIMX=10)
      DIMENSION X(NDIMX),N(NDIMX),FMAX(NREG),NM(NREG),XI(NDO,NDIM)
      double precision ntot,ntold
      LOGICAL LGLOB,LLOCAL
C----------------------------------------------------------------------
C
C  COUNTERS:
C  *********
C
C  NTOT  :   CURRENT TOTAL NUMBER OF TRIALS INCLUDING THE ONES FROM
C            EARLIER RUNS FOR CONTINUED SAMPLING
C  NN    :   CURRENT NUMBER OF TRIALS IN THE ACTUAL RUN
C  NCALL :   TOTAL NUMBER OF ACTUAL FUNCTION CALLS
C  MCALL1:   TOTAL NUMBER OF FUNCTION CALLS DURING CORRECTION PROCEDURE
C            FOR UNDERESTIMATED LOCAL MAXIMA
C  MCALL2:   TOTAL NUMBER OF FUNCTION CALLS DURING CORRECTION PROCEDURE
C            FOR UNDERESTIMATED GLOBAL MAXIMUM
C  NEV   :   CURRENT NUMBER OF ACCEPTED EVENTS IN THE ACTUAL RUN
C  MEV1  :   CURRENT NUMBER OF EVENTS ACCEPTED IN THE ACTUAL RUN
C            DURING CORRECTION PROCEDURE FOR UNDERESTIMATED LOCAL MAXIMA
C  MEV2  :   CURRENT NUMBER OF EVENTS ACCEPTED IN THE ACTUAL RUN
C            DURING CORRECTION PROCEDURE FOR UNDERESTIMATED GLOBAL MAX.
C
C-----------------------------------------------------------------------

      IF(IPRINT.GE.2) THEN
        WRITE(LUNTES,'(//A/A,5X,2I3,I8,I3,2(1PE13.6))')
     +     ' ENTRY HSGENM ',' NCONT,NDIM,NEVENT,ICONTI,FFMAX, GLAST',
     +                        NCONT,NDIM,NEVENT,ICONTI,FFMAX,GLAST
        WRITE(LUNTES,'(5X,A/5X,1PE25.16,4I8,4I5)')
     &     ' NTOT,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,NDO,MBIN,NREG :',
     &       NTOT,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,NDO,MBIN,NREG
        WRITE(LUNTES,'(5X,A/5X,4(1PD12.5),2L4)')
     &     ' FFMOLD, FMLOLD, DNCGLO, DNCLOC, LLOCAL, LGLOB ',
     &       FFMOLD,FMLOLD,DNCGLO,DNCLOC,LLOCAL,LGLOB
        WRITE(LUNTES,'(A/50(10(1PE12.4)/))')  ' HSGENM - XI  :',XI
      ENDIF

C...Basic initialization
      NN=0
      NEV=0
      MEV1=0
      MEV2=0
      AMI=1D0/FLOAT(MBIN)

      IF(ICONTI.EQ.1) THEN
C...New start of event sampling discarding potentially existing 
C...information on trial numbers and corrections from previous runs
        LGLOB=.FALSE.
        LLOCAL=.FALSE.
        DNCGLO=0D0
        DNCLOC=0D0
        JCOR=0
        NTOT=0D0
        NCALL=0
        MCALL1=0
        MCALL2=0
        DO 1 I=1,NREG
          NM(I)=0
    1   CONTINUE

      ELSEIF(ICONTI.EQ.2) THEN
C...Restart for continued sampling
        G=GLAST
        GLOMAX=FFMAX
C...Reset parameter for random evaluation of coordinate vector, 
C...needed in case of continued correction
         JJ=JCOR - 1
         DO 3 K=1,NDIM
            JJJ=JJ/MBIN
            N(K)=JJ-JJJ*MBIN
            JJ=JJJ
  3      CONTINUE
C...Test whether correction of local maximum has to be continued
        IF(LLOCAL) GOTO 290
C...Test whether correction of global maximum has to be continued
        IF(LGLOB) GOTO 390
C...Test whether last function value exceeds local/global maxima
        J=JCOR
        GOTO 190

      ELSE
        WRITE(LUNOUT,'(/A,I5/A,I5/A)')
     &       ' *** HSGENM: WRONG INPUT FOR ICONTI, ICONTI=',ICONTI,
     &       ' ***         CALLED WITH NCONT=',NCONT,
     &       ' *** EXECUTION STOPPED'
        stop
      ENDIF

  100 CONTINUE
C...Entry after correctino procedures for local/global maxima
      LLOCAL=.FALSE.
      LGLOB=.FALSE.
      DNCGLO=0D0
      DNCLOC=0D0

  110 CONTINUE
C...Unrestricted event sampling
      NTOT=NTOT + 1
      NN=NN + 1
      J=INT(HSRNDM(-1)*FLOAT(NREG) + 1D0)
      Y=HSRNDM(-1)*FFMAX
      NM(J)=NM(J)+1
C...Rejection on the basis of estimated local maxima
      IF(Y.GT.FMAX(J)) GOTO 110
      JJ=J-1
      DO 10 K=1,NDIM
        JJJ=JJ/MBIN
        N(K)=JJ-JJJ*MBIN
        X(K)=(HSRNDM(-1)+N(K))*AMI
        JJ=JJJ
 10   CONTINUE
      G=HSTRIT(F,X,NDIM,NCONT,XI,NDO)
      NCALL=NCALL + 1
C...Rejection on the basis of actual function value
      IF (Y.GE.G) GOTO 110

C...Event accepted
      CALL HSACPT(NCONT)
      NEV=NEV + 1
      IF(NEV.GE.NEVENT) THEN
C...Prepare potential correction for further runs if the last
C...calculated function value is larger than the local/global maximum
        JCOR=J
        GOTO 400
      ENDIF

C*********************************
C...Entry for restart with continued sampling. No further correction
C...from previous runs necessary: test the last function value
C...evaluated in the prevsious runs.
 190  CONTINUE
      IF (G.LE.FMAX(J)) GOTO 110
      JCOR=J

      IF (JCOR.EQ.IBIMAX) THEN
C...Correction of global maximum in bin IBIMAX
        LGLOB=.TRUE.
        GLOMAX=G
        FFMOLD=FFMAX
        GOTO 300
      ENDIF

C...Correct for underestimation of local maximum in bin J.
C...Correction for global maximum later at 300
      LLOCAL=.TRUE.
      FMLOLD=FMAX(J)
      GLOMAX=FFMAX

 200  CONTINUE
C...Entry for correction loop with continued correction of maximum
      IF(G.GT.GLOMAX) THEN
C...New local maximum larger than current global one
        GLOMAX=G
        IF(.NOT.LGLOB) THEN
C...LGLOB=.TRUE. possible at this point for repeated correction in 
C...the same bin
          LGLOB=.TRUE.
          FFMOLD=FFMAX
          DNCLOC=DNCLOC + (NM(JCOR)-1)*(FFMOLD-FMAX(JCOR))/FFMOLD
        ENDIF
      ELSE
C...New value smaller than global maximum
        DNCLOC=DNCLOC + (NM(JCOR)-1)*(G-FMAX(JCOR))/FFMAX
      ENDIF

      IF(IPRINT.GE.4) THEN
        WRITE(LUNTES,'(/A,I2/A,I5,A,2(1PD14.5))')
     &       ' HSGENM / NCONT=', NCONT,
     &       ' CORRECTION OF LOCAL MAXIMUM IN BIN', JCOR,
     &       ' / OLD-NEW :',FMAX(JCOR),G
      ENDIF
      FMAX(JCOR)=G

 210  CONTINUE
C...Entry for correction loop without continued correction of the 
C...local maximum
      IF(DNCLOC.LT.1D0) THEN
        IF(HSRNDM(-1).GT.DNCLOC) THEN
          DNCLOC=0D0
          LLOCAL=.FALSE.
          GOTO 300
        ENDIF
        DNCLOC=1D0
      ENDIF
      DNCLOC=DNCLOC - 1D0
      DO 211 K=1,NDIM
        X(K)=(HSRNDM(-1)+N(K))*AMI
 211  CONTINUE
      G=HSTRIT(F,X,NDIM,NCONT,XI,NDO)
      MCALL1=MCALL1 + 1
      NCALL=NCALL+1
C...Events accepted only if function value larger than old maximum, 
C...i.e., falling into the bin where events might have been discarded
C...because of wrong maximum
      IF (G.LT.FMLOLD) GOTO 210
      IF (G.LT.(HSRNDM(-1)*FMAX(JCOR))) GOTO 290
      CALL HSACPT(NCONT)
      MEV1=MEV1 + 1
      NEV=NEV + 1
      IF(NEV.GE.NEVENT) THEN
C...Remember potentially outstanding correction for wrong estimation
C...of global maximum
        IF (LGLOB) THEN
          DNCGLO=DNCGLO + NTOT*(GLOMAX-FFMAX)/FFMOLD
          FFMAX=GLOMAX
          IF(DNCGLO.LE.0D0) THEN
            WRITE(LUNTES,'(2A/A,3(1PE13.5))')
     &         ' HSGENM - ERROR AT EXIT: WRONG DEFINITION OF DNCGLO',
     &         ' (BEFORE 290 CONTINUE)', ' G, GLOMAX, FFMAX',
     &          G, GLOMAX, FFMAX
            WRITE(LUNTES,'(A,4I5,2L4/2(A,5(1PE20.12)/))')
     &       ' JCOR, J, IBIMAX, ICONTI, LGLOB, LLOCAL',
     &         JCOR, J, IBIMAX, ICONTI, LGLOB, LLOCAL,
     &       ' FFMAX, GLOMAX, FFMOLD, G, GLAST',
     &         FFMAX, GLOMAX, FFMOLD, G,GLAST,
     &       ' FMAX(JCOR), DNCGLO, DNCLOC, FMLOLD',
     &         FMAX(JCOR), DNCGLO, DNCLOC, FMLOLD
            WRITE(LUNTES,'(5X,A/5X,1PE25.16,I20,3I10,4I5)')
     &       ' NTOT,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,NDO,MBIN,NREG :',
     &         NTOT,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,NDO,MBIN,NREG
          ENDIF
        ENDIF
        IF (DNCLOC.LE.0D0) LLOCAL=.FALSE.
        GOTO 400
      ENDIF

C****************************
C...Entry for restart with continued sampling / correction for wrong
C...local maximum to be continued: test the function value from the 
C...previous run
 290  CONTINUE
      IF (G.LE.FMAX(JCOR)) THEN
        GOTO 210
      ELSE
C...Further correction of estimated local maximum. New local maximum
C...could become global one
        GOTO 200
      ENDIF

 300  CONTINUE
      IF (.NOT.LGLOB) GOTO 100
C...Correction for underestimation of global maximum
      IBIMAX=JCOR
      NTOLD=NTOT

 320  CONTINUE
C...Entry if repeated correction of global maximum is necessary
      DNCGLO=DNCGLO + NTOLD*(GLOMAX-FFMAX)/FFMOLD
      IF(IPRINT.GE.3) THEN
        WRITE(LUNTES,'(/A,I2/A,I5,A,2(1PD14.5))')
     &       ' HSGENM / NCONT=', NCONT,
     &       ' CORRECTION OF GLOBAL MAXIMUM IN BIN', JCOR,
     &       ' / OLD-NEW :',FFMAX,GLOMAX
      ENDIF
      IF(GLOMAX.LT.FFMAX.OR.DNCGLO.LE.0D0) THEN
        WRITE(LUNOUT,'(/A,I2,A/A,4I5,2L4/2(A,5(1PE20.12)/))')
     &       ' HSGENM / NCONT=', NCONT,
     &       ' : WRONG CORRECTION OF GLOBAL MAXIMUM',
     &       ' JCOR, J, IBIMAX, ICONTI, LGLOB, LLOCAL',
     &         JCOR, J, IBIMAX, ICONTI, LGLOB, LLOCAL,
     &       ' FFMAX, GLOMAX, FFMOLD, G, GLAST',
     &         FFMAX, GLOMAX, FFMOLD, G,GLAST,
     &       ' FMAX(JCOR), DNCGLO, DNCLOC, FMLOLD',
     &         FMAX(JCOR), DNCGLO, DNCLOC, FMLOLD
            WRITE(LUNOUT,'(5X,A/5X,2(1PE25.16),3I10,4I5)')
     &       ' NTOT,ntold,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,NDO,MBIN :',
     &         NTOT,ntold,NCALL,MCALL1,MCALL2,IBIMAX,JCOR,NDO,MBIN
            WRITE(LUNOUT,'(/A)')
     &       ' **** EXECUTION STOPPED '
            STOP
      ENDIF
      FFMAX=GLOMAX
      FMAX(JCOR)=GLOMAX

 310  CONTINUE
      IF(DNCGLO.LT.1D0) THEN
        IF(HSRNDM(-1).GT.DNCGLO) THEN
          DNCGLO=0D0
          LGLOB=.FALSE.
          GOTO 100
        ENDIF
        DNCGLO=1D0
      ENDIF
      DNCGLO=DNCGLO - 1D0
      J=INT(HSRNDM(-1)*FLOAT(NREG) + 1D0)
      NM(J)=NM(J)+1
      NTOT=NTOT + 1
      NN=NN + 1
C...Accept additional events only in the bin containing the 
C...underestimated global maximum
      IF (J.NE.JCOR) GOTO 310
      DO 311 K=1,NDIM
        X(K)=(HSRNDM(-1)+N(K))*AMI
 311  CONTINUE
      G=HSTRIT(F,X,NDIM,NCONT,XI,NDO)
      MCALL2=MCALL2 + 1
      NCALL=NCALL + 1

C...Event accepted only if function value larger than old maximum to 
C...correct for too high relative normalization of all other bins.
      IF(G.LE.FFMOLD) GOTO 310
C...Standard rejection procedure (step 2) for events in the required
C...region
      IF (G.GE.(HSRNDM(-1)*FFMAX)) THEN
        CALL HSACPT(NCONT)
        MEV2=MEV2 + 1
        NEV=NEV + 1
        IF(NEV.GE.NEVENT) THEN
          IF (DNCGLO.LE.0D0) LGLOB=.FALSE.
C...Remember potentially outstanding correction for wrong estimation 
C...of the global maximum
          GOTO 400
        ENDIF
      ENDIF

C****************************
C...Entry for restart with continued sampling / correction of wrong
C...global maximum to be continued: Test the function value of the 
C...previous run
 390  CONTINUE
      IF (G.GT.FFMAX) THEN
        GLOMAX=G
        GOTO 320
      ELSE
        GOTO 310
      ENDIF

 400  CONTINUE
      GLAST=G

      IF(IPRINT.GE.3) THEN
        WRITE(LUNTES,
     &         '(/A/5X,A/5X,I3,4(1PD12.5)/5X,A/5X,2(1PD12.5),2L4)')
     &  ' *** TEST PRINT HSGENM / EXIT',
     &     ' NCONT, FFMAX, FFMOLD, FMLOLD, GLAST ',
     &       NCONT,FFMAX, FFMOLD,FMLOLD, GLAST,
     &     ' DNCGLO, DNCLOC, LLOCAL, LGLOB ',
     &       DNCGLO, DNCLOC, LLOCAL, LGLOB
        WRITE(LUNTES,'(5X,A/5X,1PE25.16,5I8,2I5)')
     &     ' NTOT, NCALL, MCALL1, MCALL2, MEV1, MEV2, IBIMAX, JCOR',
     &       NTOT,NCALL,MCALL1,MCALL2,MEV1,MEV2,IBIMAX,JCOR
        WRITE(LUNTES,'(/A)') ' *****  ARRAY NM(IBIN)'
        WRITE(LUNTES,'(10I8)')  (NM(I),I=1,NREG)
        WRITE(LUNTES,'(/A)') ' *****  LOCAL MAXIMA FMAX(I) '
        WRITE(LUNTES,'(3(I5,1PE15.6))')  (I,FMAX(I),I=1,NREG)
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C#######################################################################
C
C    Function/routine names changed (29.01.98, HS)
C
C    THIS IS A PACKAGE CONTAINING A RANDOM NUMBER GENERATOR AND
C    SERVICE ROUTINES.
C    THE ALGORITHM IS FROM
C      'TOWARD A UNVERSAL RANDOM NUMBER GENERATOR'
C      G.MARSAGLIA, A.ZAMAN ;  FLORIDA STATE UNIV. PREPRINT FSU-SCRI-87-
C    IMPLEMENTATION BY K. HAHN  DEC. 88,
C    THIS GENERATOR SHOULD NOT DEPEND ON THE HARD WARE ( IF A REAL HAS
C    AT LEAST 24 SIGNIFICANT BITS IN INTERNAL REPRESENTATION ),
C    THE PERIOD IS ABOUT 2**144,
C    TIME FOR ONE CALL AT IBM-XT IS ABOUT 0.7 MILLISECONDS,
C    THE PACKAGE CONTAINS
C      FUNCTION HSRNDM(I)                     : GENERATOR
C      SUBROUTINE HSRNST(NA1,NA2,NA3,NB4)   : INITIALIZATION
C      SUBROUTINE HSRNIN(U,C,CD,CM,I,J)     : PUT SEED TO GENERATOR
C      SUBROUTINE HSRNOU(U,C,CD,CM,I,J)     : TAKE SEED FROM GENERATOR
C      SUBROUTINE HSRNTE(IO)                : TEST OF GENERATOR
C---
C    FUNCTION HSRNDM(I)
C       GIVES UNIFORMLY DISTRIBUTED RANDOM NUMBERS  IN c0..1)
C       I  - DUMMY VARIABLE, NOT USED
C    SUBROUTINE HSRNST(NA1,NA2,NA3,NB1)
C       INITIALIZES THE GENERATOR, MUST BE CALLED BEFORE USING HSRNDM
C       NA1,NA2,NA3,NB1  - VALUES FOR INITIALIZING THE GENERATOR
C                          NA? MUST BE IN 1..178 AND NOT ALL 1
C                          12,34,56  ARE THE STANDARD VALUES
C                          NB1 MUST BE IN 1..168
C                          78  IS THE STANDARD VALUE
C    SUBROUTINE HSRNIN(U,C,CD,CM,I,J)
C       PUTS SEED TO GENERATOR ( BRINGS GENERATOR IN THE SAME STATUS
C       AS AFTER THE LAST HSRNOU CALL )
C       U(97),C,CD,CM,I,J  - SEED VALUES AS TAKEN FROM HSRNOU
C    SUBROUTINE HSRNOU(U,C,CD,CM,I,J)
C       TAKES SEED FROM GENERATOR
C       U(97),C,CD,CM,I,J  - SEED VALUES
C    SUBROUTINE HSRNTE(IO)
C       TEST OF THE GENERATOR
C       IO     - DEFINES OUTPUT
C                  = 0  OUTPUT ONLY IF AN ERROR IS DETECTED
C                  = 1  OUTPUT INDEPENDEND ON AN ERROR
C       HSRNTE USES HSRNIN AND HSRNOU TO BRING GENERATOR TO SAME STATUS
C       AS BEFORE CALL OF HSRNTE
C#######################################################################
C===========================================================
      FUNCTION HSRNDM(IDUMMY)
C     REAL*4 U,C,CD,CM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      real r
      COMMON /HSRNDC/ U(97),C,CD,CM,I,J
      COMMON /HSRDIO/ ISDINP,ISDOUT
C
      if (isdinp.lt.0) then
        call ranlux(r,1)
        hsrndm=dble(r)
        return
      endif
      
      HSRNDM = U(I)-U(J)
      IF ( HSRNDM.LT.0D0 ) HSRNDM = HSRNDM+1D0
      U(I) = HSRNDM
      I    = I-1
      IF ( I.EQ.0 ) I = 97
      J    = J-1
      IF ( J.EQ.0 ) J = 97
      C    = C-CD
      IF ( C.LT.0.0 ) C = C+CM
      HSRNDM = HSRNDM-C
      IF ( HSRNDM.LT.0D0 ) HSRNDM = HSRNDM+1D0
      RETURN
      END
C===========================================================
      SUBROUTINE HSRNST(NA1,NA2,NA3,NB1)
C     REAL*4 U,C,CD,CM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /HSRNDC/ U(97),C,CD,CM,I,J
C
      MA1 = NA1
      MA2 = NA2
      MA3 = NA3
      MB1 = NB1
      I   = 97
      J   = 33
      DO 20 II2 = 1,97
        S = 0D0
        T = 0.5D0
        DO 10 II1 = 1,24
          MAT  = MOD(MOD(MA1*MA2,179)*MA3,179)
          MA1  = MA2
          MA2  = MA3
          MA3  = MAT
          MB1  = MOD(53*MB1+1,169)
          IF ( MOD(MB1*MAT,64).GE.32 ) S = S+T
          T    = 0.5D0*T
 10     CONTINUE
        U(II2) = S
 20     CONTINUE
      C  =   362436.D0/16777216.D0
      CD =  7654321.D0/16777216.D0
      CM = 16777213.D0/16777216.D0
      RETURN
      END
C==========================================================
      SUBROUTINE HSRNIN(UIN,CIN,CDIN,CMIN,IIN,JIN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     REAL*4 UIN(97),CIN,CDIN,CMIN
C     REAL*4 U,C,CD,CM
C
      COMMON /HSRNDC/ U(97),C,CD,CM,I,J
      DIMENSION UIN(97)
C
      DO 10 KKK = 1,97
        U(KKK) = UIN(KKK)
 10   CONTINUE
      C  = CIN
      CD = CDIN
      CM = CMIN
      I  = IIN
      J  = JIN
      RETURN
      END
C==========================================================
      SUBROUTINE HSRNOU(UOUT,COUT,CDOUT,CMOUT,IOUT,JOUT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     REAL*4 UOUT(97),COUT,CDOUT,CMOUT
C     REAL*4 U,C,CD,CM
      COMMON /HSRNDC/ U(97),C,CD,CM,I,J
      DIMENSION UOUT(97)
C
      DO 10 KKK = 1,97
        UOUT(KKK) = U(KKK)
 10   CONTINUE
      COUT  = C
      CDOUT = CD
      CMOUT = CM
      IOUT  = I
      JOUT  = J
      RETURN
      END
C==========================================================
      SUBROUTINE HSRNTE(IO)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
C
C     REAL*4 UU(97)
C     REAL*4 U(6),X(6),D(6)
C
      DIMENSION UU(97)
      DIMENSION U(6),X(6),D(6)
      DATA U / 6533892.0D0 , 14220222.0D0 ,  7275067.0D0 ,
     &         6172232.0D0 ,  8354498.0D0 , 10633180.0D0 /
C
      CALL HSRNOU(UU,CC,CCD,CCM,II,JJ)
      CALL HSRNST(12,34,56,78)
      DO 10 II1 = 1,20000
        XX      = HSRNDM(-1)
 10   CONTINUE
      SD        = 0.0D0
      DO 20 II2 = 1,6
        X(II2)  = 4096.D0*(4096.D0*HSRNDM(-1))
        D(II2)  = X(II2)-U(II2)
        SD      = SD+D(II2)
 20   CONTINUE
      CALL HSRNIN(UU,CC,CCD,CCM,II,JJ)
      IF ( IO.EQ.1 .OR. SD.NE.0D0 )
     &   WRITE(LUNTES,50) (U(I),X(I),D(I),I=1,6)
      RETURN
50    FORMAT('  === TEST OF THE RANDOM-GENERATOR ===',/,
     &       '    EXPECTED VALUE    CALCULATED VALUE     DIFFERENCE',/,
     &       6(F17.1,F20.1,F15.3,/),
     &       '  === END OF TEST ;',
     &       '  GENERATOR HAS THE SAME STATUS AS BEFORE CALLING HSRNTE')
C==END OF RANDOM GENERATOR PACKAGE==========================
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE SFECFE(SFE,CFE)
      ENTRY         COSI(SFE,CFE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C********************************************************************
C
C     SUBROUTINE OF FLUKA TO GIVE SINE AND COSINE OF AN
C     RANDOM ANGLE UNIFORMLY DISTRIBUTED BETWEEN 0 AND 2*PI
C********************************************************************
C
 10   X=2.0*HSRNDM(-1)-1.0
      Y=HSRNDM(-1)
      X2=X*X
      Y2=Y*Y
      IF (X2+Y2.GT.1.0) GO TO 10
      CFE=(X2-Y2)/(X2+Y2)
      SFE=2.0*X*Y/(X2+Y2)
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSACPT(ICONT)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSKPXY/ XX,Y
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSLTYP/ LEPIN1
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSCHNN/ ICHNN
      COMMON /HSONLY/ IHSONL
      PARAMETER (NMXHEP=2000)
      COMMON /HEPEVT/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &                JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),
     &                PHEP(5,NMXHEP),VHKK(4,NMXHEP)
      COMMON /SOPHCT/ WSOPHIA
      COMMON /SPPASS/ NSOPH,NSPOUT,NFAILP,NSPACC
      COMMON /HSPASS/ NREJCW
      LOGICAL IELAST
      LOGICAL IDJSPH
      DIMENSION IQFLAV(-6:6)
      DIMENSION IQFLCC(-6:6)
      DATA IQFLAV /-6,-5,-4,-3,-1,-2,0,2,1,3,4,5,6/
      DATA IQFLCC /-5,-6,-3,-4,-2,-1,0,1,2,4,3,6,5/
      DATA NCEVE /0/

      IELAST=.FALSE.
      IDJSPH=.FALSE.
      ICHNN=ICONT
      IF (ICONT.EQ.3.OR.ICONT.EQ.15.OR.ICONT.EQ.16.OR.ICONT.EQ.17)
     &  IELAST=.TRUE.
      X=XX

C---------------------------------------      COUNT ACCEPTED EVENTS
      NCEVE=NCEVE + 1
      NEVHEP=NCEVE

C----------------------------------------     HADRONIC MASS
      GSP=SP-MEI2-MPRO2
      IF (ICONT.GT.5) THEN
        OMEGA=(PH*DKQ+PQH*DKP)/(PH*EQH+PQH*EH)
        DCTHGA=(DKP/OMEGA-MEF2/2D0/EH)/PH
        CTHGA=1D0-DCTHGA
        DKPRO=OMEGA*(EPRO+CTHGA*PPRO)
      ENDIF
      IF (ICONT.LE.2) THEN
        W2=Y*(1D0-X)*GSP+MPRO2
      ELSEIF (IELAST) THEN
        W2=MPRO2
      ELSE
        W2=Y*(1D0-X)*GS+MPRO2-2D0*(DKP-DKPS+DKPRO)
      ENDIF
      W=DSQRT(W2)
      IF (W.LE.WSOPHIA.AND.
     *    (ICONT.LE.10.OR.ICONT.GE.15)) IDJSPH=.TRUE.
      IFAIL=0
      IF ((.NOT.IDJSPH).AND.IHSONL.EQ.0) CALL HSWCUT(W,1,12,IFAIL)
      IF (IFAIL.NE.0) THEN 
        NREJCW=NREJCW+1
        RETURN
      ENDIF 

C----------------------------------------     SAMPLE QUARK TYPE
      IF (IPDFOP.GE.1.AND.(.NOT.IELAST).AND.(.NOT.IDJSPH)) THEN
  100   CALL HSFLAV(SNGL(W2),IQF,IQFR)
C...Check if energy in jet system is enough for fragmentation.
        IFAIL=0
        IF (IHSONL.EQ.0) CALL HSWCUT(W,IQF,IQFR,IFAIL)
        IF (IFAIL.NE.0) GOTO 100
        ELSE
        IQF=0
        IQFR=0
      ENDIF

C...test print
      IF(NEVHEP.LE.10.AND.IPRINT.GE.3) THEN
        WRITE(LUNTES,'(///A/5X,2I5,1PE15.4)')
     *        ' HSACPT : NEVHEP, IQF, CQP(12)',
     *                   NEVHEP, IQF
        WRITE(LUNTES,'(2(5X,6(1PE15.4)/))') CQP
      ENDIF

C------------------------------         RECONSTRUCT KINEMATICS
      IF(ICONT.LE.3) THEN
C------------------------------         NON-RADIATIVE CHANNELS
        Q2=X*Y*GSP
C -----------------------------         DEFINE  COMMON /HSIKP/
        S=X*SP
        T=-Q2
        U=-X*SP+Q2
        SS=S
        TS=T
        US=U
        DKP=0D0
        DKPS=0D0
        DKQ=0D0
        DKQS=0D0
C
        AMFEL=(SP-MPRO2-MEI2)/2D0-Q2/2D0/X
        BMFEL=Q2/2D0+MEI2
        EFEL=(BMFEL*PPRO+AMFEL*PELE)/(PELE*EPRO+PPRO*EELE)
        PFEL=SQRT((EFEL-MEF)*(EFEL+MEF))
        CTFEL=EELE/PELE*EFEL/PFEL - BMFEL/PELE/PFEL
        STFEL2=1D0-CTFEL*CTFEL
        IF (STFEL2.GT.0D0) THEN
          STFEL=DSQRT(STFEL2)
          ELSE
          STFEL=0D0
        ENDIF
        CALL SFECFE(SPHIEL,CPHIEL)
        EFQU=EELE+X*EPRO-EFEL
        PFQU=EFQU
        CTFQU=(PELE-X*PPRO-PFEL*CTFEL)/PFQU
        STFQU=SQRT((1D0+CTFQU)*(1D0-CTFQU))
C
C------------------------------        FILL STANDARD COMMON
        NHEP=4
C------------------------------        FINAL ELECTRON : IHEP=1
        IHEP=1
        ISTHEP(IHEP)=1
        IF (ICONT.EQ.1.OR.ICONT.EQ.3) THEN
          IDHEP(IHEP)=LEPIN1
          ELSE
          IF (ABS(LEPIN1).EQ.11) IDHEP(IHEP)=-LLEPT*12
          IF (ABS(LEPIN1).EQ.13) IDHEP(IHEP)=-LLEPT*14
        ENDIF
        PHEP(1,IHEP)=CPHIEL*STFEL*PFEL
        PHEP(2,IHEP)=SPHIEL*STFEL*PFEL
        PHEP(3,IHEP)=CTFEL*PFEL
        PHEP(4,IHEP)=EFEL
        PHEP(5,IHEP)=MEF
        JMOHEP(1,IHEP)=0
        JMOHEP(2,IHEP)=0
        JDAHEP(1,IHEP)=0
        JDAHEP(2,IHEP)=0
C-------------------------------       HADRONIC FINAL STATE
C-------------------------------       FOR (QUASI)ELASTIC EP:
C                                      IT IS A PROTON
        IF (IELAST) THEN
          IHEP=2
          ISTHEP(IHEP)=1
          IDHEP(IHEP)=2212
          PHEP(1,IHEP)=-CPHIEL*STFQU*PFQU
          PHEP(2,IHEP)=-SPHIEL*STFQU*PFQU
          PHEP(3,IHEP)=CTFQU*PFQU
          PHEP(4,IHEP)=EFQU
          PHEP(5,IHEP)=0D0
          JMOHEP(1,IHEP)=0
          JMOHEP(2,IHEP)=0
          JDAHEP(1,IHEP)=0
          JDAHEP(2,IHEP)=0
C-------------------------------       IF PARTON MODEL IS USED:
C                                      SCATTERED QUARK
        ELSEIF (IPDFOP.GE.1) THEN
          IHEP=2
          ISTHEP(IHEP)=1
          IF (ICONT.EQ.1) THEN
            IDHEP(IHEP)=IQFLAV(IQF)
            ELSE
            IDHEP(IHEP)=IQFLCC(IQF)
          ENDIF
          PHEP(1,IHEP)=-CPHIEL*STFQU*PFQU
          PHEP(2,IHEP)=-SPHIEL*STFQU*PFQU
          PHEP(3,IHEP)=CTFQU*PFQU
          PHEP(4,IHEP)=EFQU
          PHEP(5,IHEP)=0D0
          JMOHEP(1,IHEP)=0
          JMOHEP(2,IHEP)=0
          JDAHEP(1,IHEP)=0
          JDAHEP(2,IHEP)=0
C-------------------------------       SPECTATOR QUARK SYSTEM
          IHEP=4
          ISTHEP(IHEP)=1
          IDHEP(IHEP)=90
          PHEP(1,IHEP)=          -(PHEP(1,1)+PHEP(1,2))
          PHEP(2,IHEP)=          -(PHEP(2,1)+PHEP(2,2))
          PHEP(3,IHEP)=-PPRO+PELE-(PHEP(3,1)+PHEP(3,2))
          PHEP(4,IHEP)= EPRO+EELE-(PHEP(4,1)+PHEP(4,2))
          PHEP(5,IHEP)=0D0
        ELSE
C-------------------------------    HADRONIC FINAL STATE : IHEP=2
C                                   FOR THE STRUCTURE FUNCTIONS OPTION
          IHEP=2
          ISTHEP(IHEP)=1
          IDHEP(IHEP)=93
          PHEP(1,IHEP)=-CPHIEL*STFQU*PFQU
          PHEP(2,IHEP)=-SPHIEL*STFQU*PFQU
          PHEP(3,IHEP)=PELE-PPRO-PHEP(3,1)
          PHEP(4,IHEP)=EELE+PPRO-PHEP(4,1)
          PHEP(5,IHEP)=W
          JMOHEP(1,IHEP)=0
          JMOHEP(2,IHEP)=0
          JDAHEP(1,IHEP)=0
          JDAHEP(2,IHEP)=0
C-------------------------------       PROTON REMNANT NOT SEPARATED
          IHEP=4
          ISTHEP(4)=0
          DO 1 IAX=1,5
            PHEP(IAX,IHEP)=0D0
 1        CONTINUE
        ENDIF
C-------------------------------       NO PHOTON : IHEP=3
        IHEP=3
        ISTHEP(IHEP)=0
        DO 2 IAX=1,5
          PHEP(IAX,IHEP)=0D0
 2      CONTINUE
C
C------------------------------------------------------------------------------
C------------------------------   ***  HARD PHOTON CONTRIBUTIONS
      ELSEIF(ICONT.GT.5) THEN
        Q2=X*Y*GSP
        NHEP=4
C-------------------------------       FINAL ELECTRON : IHEP=1
        IHEP=1
        ISTHEP(IHEP)=1
        IF (ICONT.LE.10.OR.ICONT.GE.15) THEN
          IDHEP(IHEP)=LEPIN1
          ELSE
          IF (ABS(LEPIN1).EQ.11) IDHEP(IHEP)=-LLEPT*12
          IF (ABS(LEPIN1).EQ.13) IDHEP(IHEP)=-LLEPT*14
        ENDIF
        PHEP(4,IHEP)=ESH
        PHEP(3,IHEP)=COSEH*PSH
        PHEP(1,IHEP)=SINEH*PSH
        PHEP(2,IHEP)=0D0
        PHEP(5,IHEP)=MEF
C-------------------------------       FINAL STATE PHOTON : IHEP=3
        IHEP=3
        ISTHEP(IHEP)=1
        IDHEP(IHEP)=22
C---ALREADY DONE:
C       OMEGA=(PH*DKQ + PQH*DKP) / (PH*EQH + PQH*EH)
C       DCTHGA=(DKP/OMEGA-MEF2/2D0/EH)/PH
C       CTHGA=1D0-DCTHGA
        STHGA=DSQRT(DCTHGA*(2D0-DCTHGA))
        PPP2=PSH*SINEH
        STCPHG=(ESH - DKPS/OMEGA - PSH*CTHGA*COSEH) / PPP2
        PHEP(4,IHEP)=OMEGA
        PHEP(3,IHEP)=OMEGA*CTHGA
        PHEP(1,IHEP)=OMEGA*STCPHG
        AKY2=PHEP(4,IHEP)**2-PHEP(1,IHEP)**2-PHEP(3,IHEP)**2
        IF (AKY2.LT.0D0) THEN
          PHEP(2,IHEP)=0D0
          ELSE
          PHEP(2,IHEP)=DSQRT(AKY2)*DSIGN(1D0,0.5D0-HSRNDM(-1))
        ENDIF
        PHEP(5,IHEP)=0D0
C-------------------------------       HADRONIC FINAL STATE
        IF (IELAST) THEN
C-------------------------------
C                                      FOR (QUASI)ELASTIC EP: PROTON
          IHEP=2
          ISTHEP(IHEP)=1
          IDHEP(IHEP)=2212
          ETOT=EELE+EQH
          PZTOT=PELE-PQH
          PHEP(4,IHEP)=ETOT-PHEP(4,1)-PHEP(4,3)
          PHEP(1,IHEP)=-PHEP(1,1)-PHEP(1,3)
          PHEP(2,IHEP)=-PHEP(2,1)-PHEP(2,3)
          PHEP(3,IHEP)=PZTOT-PHEP(3,1)-PHEP(3,3)
          PHEP(5,IHEP)=0D0
        ELSEIF (IPDFOP.GE.1) THEN
C-------------------------------
C                                      FOR PARTON MODEL:
C                                      SCATTERED QUARK : IHEP=2
C                                      USE E-P CONSERVATION WITHOUT
C                                      QUARK PT
          IHEP=2
          ISTHEP(IHEP)=1
          IF (ICONT.LE.10) THEN
            IDHEP(IHEP)=IQFLAV(IQF)
            ELSE
            IDHEP(IHEP)=IQFLCC(IQF)
          ENDIF
          ETOT=EELE + EQH
          PZTOT=PELE - PQH
          PHEP(4,IHEP)=ETOT - PHEP(4,1) - PHEP(4,3)
          PHEP(1,IHEP)= - PHEP(1,1) - PHEP(1,3)
          PHEP(2,IHEP)= - PHEP(2,1) - PHEP(2,3)
          PHEP(3,IHEP)=PZTOT - PHEP(3,1) - PHEP(3,3)
          PHEP(5,IHEP)=0D0
C-------------------------------       SPECTATOR QUARK SYSTEM
          IHEP=4
          ISTHEP(IHEP)=1
          IDHEP(IHEP)=90
          PHEP(1,IHEP)=-(PHEP(1,1)+PHEP(1,2)+PHEP(1,3))
          PHEP(2,IHEP)=-(PHEP(2,1)+PHEP(2,2)+PHEP(2,3))
          PHEP(3,IHEP)=-PPRO+PELE
     &                 -(PHEP(3,1)+PHEP(3,2)+PHEP(3,3))
          PHEP(4,IHEP)= EPRO+EELE
     &                 -(PHEP(4,1)+PHEP(4,2)+PHEP(4,3))
          PHEP(5,IHEP)=0D0
        ELSE
C-------------------------------    HADRONIC FINAL STATE : IHEP=2
C                                   FOR STRUCTURE FUNCTIONS OPTION
          IHEP=2
          ISTHEP(IHEP)=1
          IDHEP(IHEP)=93
          ETOT=EELE+EPRO
          PZTOT=PELE-PPRO
          PHEP(4,IHEP)=ETOT-PHEP(4,1)-PHEP(4,3)
          PHEP(1,IHEP)=-PHEP(1,1)-PHEP(1,3)
          PHEP(2,IHEP)=-PHEP(2,1)-PHEP(2,3)
          PHEP(3,IHEP)=PZTOT-PHEP(3,1)-PHEP(3,3)
          PHEP(5,IHEP)=W
C-------------------------------       PROTON REMNANT NOT SEPARATED
          IHEP=4
          ISTHEP(4)=0
          DO 3 IAX=1,5
            PHEP(IAX,IHEP)=0D0
 3        CONTINUE
        ENDIF
C-------------------------------       ROTATE WHOLE EVENT
        PHIEVE=HSRNDM(-1)*2D0*PI
        CROT=DCOS(PHIEVE)
        SROT=DSIN(PHIEVE)
        DO 4 INHEP=1,4
          PNHEPX= CROT*PHEP(1,INHEP)+SROT*PHEP(2,INHEP)
          PNHEPY=-SROT*PHEP(1,INHEP)+CROT*PHEP(2,INHEP)
          PHEP(1,INHEP)=PNHEPX
          PHEP(2,INHEP)=PNHEPY
 4      CONTINUE
C
      ELSE
        STOP
      ENDIF
C
      IF (IELAST) CALL DJGELV
      IF (IHSONL.EQ.0.AND..NOT.IELAST) THEN
        CALL DJGVAR(ICHNN,X,Y,Q2)
        NSPACC=0
        IF (IDJSPH) THEN
          CALL DJGSPH(13,IMODE)
        ELSE
          CALL DJGEVT
        ENDIF
      ENDIF
C-------------------------
C---Define beam particles for /HEPEVT/
chs...ISTHEP according to H1 standard
      NHEP = NHEP + 2
C---incoming electron:
      IHEP=NHEP-1
      ISTHEP(IHEP)=201
      IDHEP(IHEP)=LEPIN1
      PHEP(5,IHEP)=0D0
      PHEP(4,IHEP)=EELE
      PHEP(1,IHEP)=0D0
      PHEP(2,IHEP)=0D0
      PHEP(3,IHEP)=PELE
      VHKK(4,IHEP)=0D0
      VHKK(1,IHEP)=0D0
      VHKK(2,IHEP)=0D0
      VHKK(3,IHEP)=0D0
      JMOHEP(1,IHEP)=0
      JMOHEP(2,IHEP)=0
      JDAHEP(1,IHEP)=0
      JDAHEP(2,IHEP)=0
C incoming nucleon:
      IHEP = NHEP
      ISTHEP(IHEP)=201
C incoming proton:
      IF (INT(HNA).EQ.1.AND.INT(HNZ).EQ.1) THEN
        IDHEP(IHEP)=2212
C incoming neutron:
      ELSEIF (INT(HNA).EQ.1.AND.INT(HNZ).EQ.0) THEN 
        IDHEP(IHEP)=2112
C other nuclei
      ELSE
        IDHEP(IHEP)=1000000000+INT(HNZ)*10000+INT(HNA)*10
      ENDIF
      PHEP(5,IHEP)=0D0
      PHEP(4,IHEP)=EPRO
      PHEP(1,IHEP)=0D0
      PHEP(2,IHEP)=0D0
      PHEP(3,IHEP)=-PPRO
      VHKK(4,IHEP)=0D0
      VHKK(1,IHEP)=0D0
      VHKK(2,IHEP)=0D0
      VHKK(3,IHEP)=0D0
      JMOHEP(1,IHEP)=0
      JMOHEP(2,IHEP)=0
      JDAHEP(1,IHEP)=0
      JDAHEP(2,IHEP)=0
C
      CALL HSUSER(2,X,Y,Q2)
C
      RETURN
      END
C
C***********************************************************************
C
C***********************************************************************
C
C   GIVE EXTERNAL WEIGHT TO EVENT
C
      SUBROUTINE HSWGTX(X,Y,IACPT)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSWGTC/ IWEIGS
      LOGICAL LFIRST
      DATA LFIRST/.TRUE./
 
      IF (LFIRST) THEN
c        Q20=100D0
        Q20=2D0
        LFIRST=.FALSE.
      ENDIF
 
      IACPT=1
      WI=1D0
      IF (IWEIGS.EQ.1) THEN
        WI=X
      ELSEIF (IWEIGS.EQ.2) THEN
        Q2=X*Y*SP
        IF (Q2.LT.Q20) THEN
          WI=Q2/Q20
          ELSE
          WI=1D0
        ENDIF
      ELSE
        RETURN
      ENDIF
      IF (HSRNDM(-1).LE.WI) RETURN
      IACPT=0
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFLAV(W2,IFL,IFLR)

C...Choose flavour of struck quark and the
C...corresponding flavour of the target remnant jet.

      COMMON /HSCUMS/  CQP(12)
      DOUBLE PRECISION CQP,HSRNDM,R

C...May want to check hadronic mass before choosing flavor? 
      IF (W2.LE.0.) THEN
         CONTINUE
      ENDIF
      
      NFL=6
      R=HSRNDM(-1)*CQP(12)
      DO 30 I=1,2*NFL
      IFL=I
      IF(R.LE.CQP(I)) GOTO 40
   30 CONTINUE
   40 CONTINUE
      IF(MOD(IFL,2).EQ.1) THEN
        IFL=(IFL+1)/2
      ELSE
        IFL=-IFL/2
      ENDIF

      IFLR=-IFL

      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   VERSION 4.2 FOR MODEL-INEDEPENDENT CALCULATION OF LEPTONIC
C   QED-CORRECTIONS IN DEEP INELASTIC ELECTRON PROTON SCATTERING
C   AND FOR VERY SMALL X
C   AND INCLUDING QUARKONIC RADIATION
C
C   MAJOR CHANGES (SMALL X AND STRFCT) 24.FEB.91: AK
C   FOR SOFT+VIRTUAL CORRECTIONS CHANGED 13.3.91: HS
C   MODIFICATIONS FOR LOW Q2 INCLUDED    28.8.92: HS
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   MAXIMAL PHOTON ENERGY AS FUNCTION OF XS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSOMAX(X,Y,XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
C
      CALL HSFLAB(XS)
      STU=(XS-X)*Y*GS
      STUM=(XS-X)*Y*GS+XS*XS*MPRO2+MEI2+MEF2
      HO1 = STU/2D0/STUM*( EH + XS*EPRO - ESH
     *                     + DSQRT(PSH*PSH*SINEH*SINEH
     *                             + (PH - XS*PPRO - PSH*COSEH)
     *                              *(PH - XS*PPRO - PSH*COSEH)))
      HO2 = STU*STU/4D0/STUM/HO1
      HSOMAX = DMAX1(HO1,HO2)
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   MAXIMAL PHOTON ENERGY AS FUNCTION OF XS (FOR CC)
C   FROM CK0MAX (epcctot.f, 28/08/98)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCKMX(X,Y,XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
C
      S=XS*SP
      T=-X*Y*SP
      U=-XS*(1D0-Y)*SP
      SIGMA=(S-MEI2-MQI2)/2D0
      TAU=-(T-MEI2)/2D0
      QRHO=-(U-MQI2)/2D0
      QH1=-(SIGMA*PELE+EELE*DSQRT(SIGMA*SIGMA-MEI2*MQI2))
      QPQH=-(SIGMA*SIGMA-EELE*EELE*MQI2)/QH1
      QEQH=DSQRT(QPQH*QPQH+MQI2)
      EESH=(QPQH*TAU+PELE*QRHO)/(PELE*QEQH+QPQH*EELE)
      PPSH=EESH
      ECOS=(QRHO*EELE-TAU*QEQH)/PPSH/(PELE*QEQH+QPQH*EELE)
      ESIN=DSQRT(1D0-ECOS*ECOS)
      SMF=MEI2+MQI2+MQF2
      STU=S+T+U-SMF
      STUM=STU+MQF2
      CK0MAX=STU/2D0/STUM*(EELE+QEQH-EESH
     *  +DSQRT(PPSH*PPSH*ESIN*ESIN
     *         +(PELE-QPQH-PPSH*ECOS)*(PELE-QPQH-PPSH*ECOS)))
      HSCKMX=CK0MAX
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   MINIMAL PHOTON ENERGY AS FUNCTION OF XS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSOMIN(X,Y,XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
C
      CALL HSFLAB(XS)
      STU = (XS-X)*Y*GS
      STUM = (XS-X)*Y*GS + XS*XS*MPRO2 + MEI2 + MEF2
      HO1 = STU/2D0/STUM*( EH + XS*EPRO - ESH
     *                     + DSQRT(PSH*PSH*SINEH*SINEH
     *                             + (PH - XS*PPRO - PSH*COSEH)
     *                              *(PH - XS*PPRO - PSH*COSEH)))
      HO2 = STU*STU/4D0/STUM/HO1
      HSOMIN = DMIN1(HO1,HO2)
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   XSMIN FROM OMEGA-MAX(XS)=DELTA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSXSMN(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C
      XS1 = X
      XS2 = 1D0
      DO 90, N = 1,70,1
        XS3 = (XS1+XS2)/2D0
        OM3 = HSOMAX(X,Y,XS3)
        IF (OM3.LT.DELTA) THEN
          XS1 = XS3
        ELSE
          XS2 = XS3
        ENDIF
90    CONTINUE
      HSXSMN = XS3

C...CHECK ON ES > ME
      IF (X.LT.1D-6) THEN
      BXS=GU*(TP+MEI2)/(GU*GU-4D0*MEI2*MPRO2)
      DXS=1D0-(GU*GU-4D0*MEI2*MPRO2)/GU/GU
     *        *(TP*(TP-2D0*MEI2)-7D0*MEI2*MEI2)/(TP+MEI2)/(TP+MEI2)
      IF (DXS.LT.0D0) RETURN
      XESMIN=-BXS*(1D0+DSQRT(DXS))
      XESMAX=-BXS*(1D0-DSQRT(DXS))
      HSXSMN=DMAX1(HSXSMN,XESMIN,XESMAX)
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   XSMIN FROM OMEGA-MAX(XS)=DELTA (FOR CC)
C   FROM CXSMIN (epcctot.f, 28/08/98)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCXSM(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C
      XS1=X+(MEI2+MQI2+MQF2)/SP/Y
      XSMM=XS1
      XS2=1D0
      DO 90, N=1,70,1
        XS3=(XS1+XS2)/2D0
        OM3=HSCKMX(X,Y,XS3)
        IF (OM3.LT.DELTA) THEN
          XS1=XS3
        ELSE
          XS2=XS3
        ENDIF
90    CONTINUE
      HSCXSM=DMAX1(XS3,XSMM)
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   XSCUT FROM OMEGA-MIN(XS)=DELTA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSXSCT(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C
      XS1 = X
      XS2 = 1D0
      DO 100, N = 1,70,1
        XS3 = (XS1+XS2)/2D0
        OM3 = HSOMIN(X,Y,XS3)
        IF (OM3.LT.DELTA) THEN
          XS1 = XS3
        ELSE
          XS2 = XS3
        ENDIF
100   CONTINUE
      HSXSCT = XS3
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   PARAMETERS OF THE HERA-LAB-SYSTEM REAL*8
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFLAB(XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
C
      EH = EELE
      PH = PELE
      EQH = XS*EPRO
      PQH = XS*PPRO
      ESH = -(PPRO*(TP-MEI2-MEF2)+PH*(UP-MPRO2-MEF2))
     &              /2D0/(EH*PPRO+PH*EPRO)
      PSH = DSQRT((ESH+MEF)*(ESH-MEF))
      COSEH = (TP-2D0*MEF2+2D0*EH*ESH)/2D0/PH/PSH
      SINH2=1D0-COSEH*COSEH
      IF(SINH2.LE.0D0) THEN
        SINEH=0D0
        ELSE
        SINEH=DSQRT(SINH2)
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   PARAMETERS OF THE HERA-LAB-SYSTEM REAL*8 (FOR CC)
C   FROM CCLABEA (epcctot.f, 28/08/98)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSCLAB(X,Y,XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
C
      EH=EELE
      PH=PELE
      S=SP*XS
      CHI=(S-MEI2-MQI2)/2D0
      QH1=-(CHI*PELE+EELE*DSQRT(CHI*CHI-MEI2*MQI2))
      PQH=-(CHI*CHI-EELE*EELE*MQI2)/QH1
      EQH=DSQRT(PQH*PQH+MQI2)
      GS=SP-MEI2-MPRO2
      GU=-(1D0-Y)*GS
      TP=-X*Y*GS
      UP=GU+MPRO2
      ESH=-(PPRO*(TP-MEI2)+PH*(UP-MPRO2))/2D0/(EH*PPRO+PH*EPRO)
      PSH=ESH
      COSEH=(TP+2D0*EH*ESH)/2D0/PH/PSH
      SINH2=1D0-COSEH*COSEH
      IF(SINH2.LE.0D0) THEN
        SINEH=0D0
        ELSE
        SINEH=DSQRT(SINH2)
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   PARAMETERS OF THE PSEUDO-CMS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFCMS(X,Y,XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSCMS1/ COSQ,SINQ
      COMMON /HSPSPC/ IPHSPC
C
      IPHSPC=0
      MXS = XS*MPRO
      MXS2 = MXS*MXS
      ETOT = DSQRT((XS-X)*Y*GS+MXS2+MEI2+MEF2)
      EPTN = (Y*GS+2D0*XS*MPRO2+(MEI2+MEF2)/XS)/2D0/ETOT
      PPTN = DSQRT((EPTN+MPRO)*(EPTN-MPRO))
      EQ = XS*EPTN
      PQ = XS*PPTN
      EEL = ((XS-X*Y)*GS+MEI2)/2D0/ETOT
      IF (EEL.LE.MEI) THEN
        IPHSPC=1
        EEL=MEI*(1D0+1D-12)
      ENDIF
      PEL = DSQRT((EEL+MEI)*(EEL-MEI))
      ES = ((XS-Y*(XS-X))*GS-MEF2)/2D0/ETOT
      IF (ES.LE.MEF) THEN
        IPHSPC=1
        ES=MEF*(1D0+1D-12)
      ENDIF
      PS = DSQRT((ES+MEF)*(ES-MEF))
      OMEGA=(XS-X)*Y*GS/2D0/ETOT
      COSE=(-X*Y*GS+2D0*EEL*ES-MEI2-MEF2)/2D0/PEL/PS
      IF (COSE.GE.1D0.OR.COSE.LE.-1D0) THEN
        COSE=DSIGN(1D0,COSE)
        SINE=0D0
        ELSE
        SINE=DSQRT(1D0-COSE*COSE)
      ENDIF
      U = -XS*(1D0-Y)*GS + XS*XS*MPRO2 + MEF2
      COSQ = (U + 2D0*ES*EQ-MEF2-MXS2)/2D0/PS/PQ
      IF (COSQ.GE.1D0.OR.COSQ.LE.-1D0) THEN
        COSQ=DSIGN(1D0,COSQ)
        SINQ=0D0
        ELSE
        SINQ=DSQRT(1D0-COSQ*COSQ)
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   PARAMETERS OF THE PSEUDO-CMS (FOR CC)
C   FROM CCMSEA (epcctot.f, 28/08/98)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSCCMS(X,Y,XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSCMS1/ COSQ,SINQ
      COMMON /HSPSPC/ IPHSPC
C
      IPHSPC=0
      EQ=(XS*Y*SP-MEI2)/2D0/DSQRT(Y*(XS-X)*SP-MEI2-MQI2)
      PQ=DSQRT(EQ*EQ-MQI2)
      EEL=((XS-X*Y)*SP-MQI2)/(XS*Y*SP-MEI2)*EQ
      IF (EEL.LE.MEI) THEN
        IPHSPC=1
        EEL=MEI*(1D0+1D-12)
      ENDIF
      PEL=DSQRT((EEL+MEI)*(EEL-MEI))
      ES=((XS-Y*(XS-X))*SP+MEI2+MQI2)/(XS*Y*SP-MEI2)*EQ
      IF (ES.LE.0D0) THEN
        IPHSPC=1
        ES=0D0
        RETURN
      ENDIF
      PS=ES
      COSE=(-X*Y*SP+2D0*EEL*ES-MEI2)/2D0/PEL/PS
      IF (COSE.GE.1D0.OR.COSE.LE.-1D0) THEN
        COSE=DSIGN(1D0,COSE)
        SINE=0D0
        ELSE
        SINE=DSQRT(1D0-COSE*COSE)
      ENDIF
      COSQ=(2D0*ES*(ES-EEL)+X*Y*SP+MEI2)/2D0/PQ/PS
      IF (COSQ.GE.1D0.OR.COSQ.LE.-1D0) THEN
        COSQ=DSIGN(1D0,COSQ)
        SINQ=0D0
        ELSE
        SINQ=DSQRT(1D0-COSQ*COSQ)
      ENDIF
      SMF=MEI2+MQI2+MQF2
      OMEGA=((XS-X)*Y*SP-SMF)/2D0/(EEL+EQ-ES)
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INVARIANTS CALCULATED IN KP - CMS - SYSTEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFIV1(X,Y,XS,A1,TTS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
C
      CALL HSFIVM(X,Y,XS)
      DKP = A1/2D0
      DKPS = (TTS-T+A1)/2D0
      DKQS =(XS-X)*Y*GS/2D0
      DKQ = DKQS +(TTS-T)/2D0
      SS = S - 2D0*(DKP+DKQ)
      US = U + 2D0*(DKPS-DKQ)
      GX = -TTS/XS
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INVARIANTS CALCULATED IN KPS - CMS - SYSTEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFIV2(X,Y,XS,A2,TTS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
C
      CALL HSFIVM(X,Y,XS)
      DKP = (T+A2-TTS)/2D0
      DKPS = A2/2D0
      DKQS =(XS-X)*Y*GS/2D0
      DKQ = DKQS +(TTS-T)/2D0
      SS = S - 2D0*(DKP+DKQ)
      US = U + 2D0*(DKPS-DKQ)
      GX = -TTS/XS
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INVARIANTS CALCULATED IN KQ - CMS - SYSTEM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFIV3(X,Y,XS,A1,A3)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGIKP/ GS,GU,GX,TP,UP

      CALL HSFIVM(X,Y,XS)
      DKP = A1/2D0
      DKQ = A3/2D0
      DKQS =(XS-X)*Y*GS/2D0
      DKPS = DKP+DKQ-DKQS
      TS = T + A3 - 2D0*DKQS
      SS = S - A1 - A3
      US = U + A1 - 2D0*DKQS
      GX = -TS/XS
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   MANDELSTAM INVARIANTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFIVM(X,Y,XS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
C
      S = XS*GS + XS*XS*MPRO2 + MEI2
      T = -X*Y*GS
      U = -XS*(1D0-Y)*GS + XS*XS*MPRO2 + MEF2
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   CAPITAL INVARIANTS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFIVC(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
C
      GS = SP-MEI2-MPRO2
      TP = -X*Y*GS
      GU = -(1D0-Y)*GS
      UP = GU+MPRO2+MEF2
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND OF KP TERM (INITIAL STATE RADIATION)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSTSK1(X)
C
C  X(1) -->  XX
C  X(2) -->  G=-1/Q**2   -->  Q**2=-1/G-->  Y
C  X(3) -->  LOG(XS-XX)
C  X(4) -->  LOG(2*K.P)
C  X(5) -->  TS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSKPXY/ XX,Y
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSXSLM/ XSMIN,XSCUT
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSASY0/ PONULL,CHNULL,ICPASY
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSFIJK/ F1(2,2),F2(2,2),F3(2,2)
      COMMON /HSGIJK/ G1(2,2),G3(2,2),G4(2,2)
      COMMON /HSPSPC/ IPHSPC
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSWGTC/ IWEIGS
      DIMENSION X(5),DBOS(2),DSBOS(2),RUNALP(2)
      COMPLEX(kind=8) HSSRGG,CG
C
C---X-VALUE
      XMINL=DLOG(XMIN)
      XMAXL=DLOG(XMAX)
      XXL=XMINL+(XMAXL-XMINL)*X(1)
      XX=DEXP(XXL)
C---Y-VALUE
      GS=SP-MEI2-MPRO2
      YMAXX=XX*(1D0-4D0*MEI2*MPRO2/GS/GS)/(XX*(1D0+XX*MPRO2/GS)+MEI2/GS)
      IF(ICUT.LT.3) THEN
C---CUT IN EXTERNALLY DEFINED Q**2(MIN)
        GMAX=-1D0/(XX*YMAXX*GS)
        ELSEIF(ICUT.EQ.3) THEN
C---CUT IN Y
        QQ2MIN=XX*YMIN*GS
C---CUT ON ELECTRON SCATTERING ANGLE (MASSES NEGLECTED)
        YMINTH=1D0/(1D0+XX*CTHCON)
        QT2MIN=XX*YMINTH*GS
C---CUT ON ELECTRON TRANSVERSE MOMENTUM (MASSES NEGLECTED)
        QP2MIN=XX*GS/2D0*(1D0-DSQRT(1D0-PTXM0/XX))
        QP2MAX=XX*GS/2D0*(1D0+DSQRT(1D0-PTXM0/XX))
        YP2MAX=QP2MAX/GS/XX
        GMIN=-1D0/MAX(Q2MIN,QQ2MIN,QT2MIN,QP2MIN)
        GMAX=-1D0/(XX*MIN(YMAX,YMAXX,YP2MAX)*GS)
        GMAX=MIN(GMAX,-1D0/Q2MAX)
        ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSTSK1'
        STOP
      ENDIF
      DG2=DMAX1(GMAX-GMIN,0D0)
C---CUT IN W LATER
      GACT=GMIN+X(2)*DG2
      Y=-1D0/(XX*GACT*GS)
C
C---EXTERNAL WEIGHT
      IACPT=1
      IF (IWEIGS.GT.0) THEN
        CALL HSWGTX(XX,Y,IACPT)
        IF (IACPT.EQ.0) THEN
          HSTSK1=0D0
          RETURN
        ENDIF
      ENDIF
C
      CALL HSFIVC(XX,Y)
      CALL HSDELO(XX,Y)
      XSMIN=HSXSMN(XX,Y)
C     XSCUT=HSXSCT(XX,Y)
C
      IF(IPRINT.GT.30) THEN
        WRITE(LUNTES,230)SP,XX,Y,XSMIN,XSCUT,DELEPS,DELTA
230     FORMAT(' ***************************************************',/
     F        ,' SP = ',1PD12.3,/
     F        ,' X = ',D12.6,'   Y = ',D12.6,/
     F        ,' XSMIN = ',D17.11,'   XSCUT = ',D17.11,/
     F        ,' DELEPS = ',D12.6,'   DELTA = ',D14.8,/
     F       ,' ***************************************************',//)
      ENDIF
C
      XSMAX = 1D0
      XSMINI = XSMIN
      IF(XSMAX.LE.XSMINI) THEN
        HSTSK1=0D0
        RETURN
      ENDIF
C--- SUBSTITUTION UV=LN(XS-X)
      UVMIN = DLOG(XSMINI-XX)
      UVMAX = DLOG(XSMAX-XX)
      UV = UVMIN + (UVMAX-UVMIN)*X(3)
      XS = XX+DEXP(UV)
      CALL HSFCMS(XX,Y,XS)
      IF(IPHSPC.EQ.1) THEN
        HSTSK1=0D0
        RETURN
      ENDIF
      CALL HSFLAB(XS)
      CALL HSLZK1(ZMIN,ZMAX)
      IF(IPHSPC.EQ.1) THEN
        HSTSK1=0D0
        RETURN
      ENDIF
C---SUBSTITUTION V = LN(A1)
      A1MAX = 2D0*OMEGA*(EEL-PEL*ZMIN)
      IF ((ZMAX.GE.0.9999D0).AND.(EEL/MEI.GT.1D3)) THEN
        A1MIN = 2D0*OMEGA*MEI2/2D0/EEL
      ELSE
        A1MIN = 2D0*OMEGA*(EEL-PEL*ZMAX)
      ENDIF
      VMAX = DLOG(A1MAX)
      VMIN = DLOG(A1MIN)
      V= VMIN + (VMAX-VMIN)*X(4)
      A1 = DEXP(V)
C---NO SUBSTITUTION FOR TS
      CALL HSLTS1(A1,XX,Y,XS,TSMIN,TSMAX,TSM,TSP,CFKP)
      IF(IPHSPC.EQ.1)THEN
       HSTSK1=0D0
       RETURN
      ENDIF
      TS=TSMIN + (TSMAX-TSMIN)*X(5)
      SQGRM2=-CFKP*(TS-TSM)*(TS-TSP)
      IF (SQGRM2.LE.0D0) THEN
        HSTSK1=0D0
        RETURN
      ENDIF
      SQGRAM = DSQRT(DABS(SQGRM2))
      CALL HSFIV1(XX,Y,XS,A1,TS)
C
      A2=2D0*DKPS
      DBOS(1)=1D0
      DBOS(2)=T/(T-MZ2)
      DSBOS(1)=1D0
      RUNALP(1)=1D0
      IF (LPAR(3).GE.3) THEN
        CG=HSSRGG(TS)
        DSBOS(1)=1D0/(1D0+DREAL(CG)/TS)
        RUNALP(1)=DSBOS(1)
      ENDIF
      DSBOS(2)=TS/(TS-MZ2)
      RUNALP(2)=1D0
      IF(LPAR(14).EQ.0)THEN
        SU1=0D0
        SU2=0D0
        PROPI1=0D0
      ELSE
        SU1 = S*S+SS*SS+U*U+US*US
        SU2 = (S-U)*(S+U) + (SS-US)*(SS+US)
        PROPI1= - 1D0/8D0/XS
     &       * (SS/DKPS/DKQS + S/DKP/DKQ + U/DKPS/DKQ + US/DKP/DKQS)
      ENDIF

      R1=4D0*GX*(
     &     -(T+6D0*MEI2)/(A1-TS)/(A1-TS)/A1
     &     +1D0/(A1-TS)/A1
     &     + 2D0/(A1+A2)/A1
     &     -8D0*MEI2*MEI2/(A1*A2+TS*TS)/(A1+A2)/A1
     &     +4D0*MEI2*MEI2/(A1*A1+TS*TS)/A1/A1
     &     -2D0*MEI2/(A1-TS)/A1/A1 )
C
      FAC1 = -(GS*GS + GU*GU - GX*(GS+GU) -4D0*MEI2*MPRO2)
      FAC2 = -GS*GX + MPRO2*T
      FAC3 = (-2D0*GS*GU + GX*(GS+GU-GX))*2D0*MEI2
      FAC4 = GU*(GX-GU)*2D0*MEI2
      R2=4D0*(
     &     -FAC1/(A1+A2-TS)*(1D0/(A1+A2)+1D0/(A1-TS))/A1
     &     +FAC2/(A1-TS)/(A1-TS)/A1
     &     +FAC3/(A1*A2+TS*TS)/(A1+A2)/A1
     &     -2D0*MPRO2/(A1+A2)/A1
     &     +2D0*MEI2*MPRO2/(A1-TS)/(A1-TS)/A1
     &     +2D0*MEI2*MPRO2/(A1-TS)/A1/A1
     &     -MPRO2/(A1-TS)/A1
     &     +FAC4/(A1*A1+TS*TS)/A1/A1  )
      FK1=-2D0*MEI2*(GS-GU)
      R3= -4D0*(
     &     (GS-GU)/(A1+A2)/A1
     &     +FK1/(A1+A2-TS)*(-1D0/(A1+A2)-1D0/(A1-TS))/A1
     &     -(GX/2D0-GS)/(A1-TS)/A1
     &     +GX*T/2D0/(A1-TS)/(A1-TS)/A1
     &     -(GX-2D0*GU)*MEI2/(A1-TS)/(A1-TS)/A1
     &     -(GX-2D0*GU)*MEI2/(A1-TS)/A1/A1)
      DO 20 IFL=1,12
        CQP(IFL)=0D0
 20   CONTINUE

      SUMME=0D0
      IF (IPDFOP.EQ.0) THEN
C..................................STRUCTURE FUNCTIONS
         CALL HSSTRF(XS,-TS)
         DO 11 IB1=1,2
          DO 12 IB2=1,2
             CQP(12)=CQP(12)
     &         +CHNULL*F1(IB1,IB2)*R1*RUNALP(IB1)*RUNALP(IB2)
     &         +CHNULL*F2(IB1,IB2)*R2*RUNALP(IB1)*RUNALP(IB2)
     &         +F3(IB1,IB2)*DFLOAT(LLEPT)*R3*RUNALP(IB1)*RUNALP(IB2)
     &       -HPOLAR*(
     &         +G1(IB1,IB2)*DFLOAT(LLEPT)*R3*2D0*RUNALP(IB1)*RUNALP(IB2)
     &         +CHNULL*G4(IB1,IB2)*R2*RUNALP(IB1)*RUNALP(IB2)
     &         +CHNULL*G3(IB1,IB2)*R1*RUNALP(IB1)*RUNALP(IB2))
 12       CONTINUE
 11      CONTINUE
        SUMME=CQP(12)
      ELSEIF (IPDFOP.EQ.1) THEN
C..................................PARTON DENSITIES
        CALL HSPVER(XS,-TS)
        IF (HPOLAR.NE.0D0) THEN
          CALL HSDPVR(XS,-TS)
        ENDIF
        CALL HSIROT
        AUP=0D0
        AUM=0D0
        ADP=0D0
        ADM=0D0
        GUP=0D0
        GUM=0D0
        GDP=0D0
        GDM=0D0
        FUPI = 0D0
        FDOI = 0D0
        FUPIB =0D0
        FDOIB =0D0
        GUPI = 0D0
        GDOI = 0D0
        GUPIB =0D0
        GDOIB =0D0
        IF(LPAR(14).NE.0)THEN
          SU3 = (S*S+U*U)/2D0
          SU4 = (S-U)*(S+U)/2D0
          SU5 = (SS*SS+US*US)/2D0
          SU6 = (SS-US)*(SS+US)/2D0
          PPHA1= -T/DKQ/DKQS
          PPHA2= -4D0*MQF2/DKQS/DKQS
          PPHA3= -4D0*MQI2/DKQ/DKQ
          RALEP =
     &      GX/2D0*(
     &        +(T+4D0*MEI2)*(1D0/A2/TS/TS-1D0/A1/TS/TS)
     &        +2D0/TS/TS - 1D0/A1/TS + 1D0/A2/TS + 2D0/A1/A2
     &        -8D0*MEI2*MEI2/A1/A2/TS/TS
     &        +4D0*MEI2*MEI2/TS/TS*(1D0/A1/A1+1D0/A2/A2)
     &        +2D0*MEI2/TS*(1D0/A1/A1+1D0/A2/A2)  )
     &     +XS*(
     &        +FAC1/TS/A1/A2
     &        +FAC2/A1/TS/TS
     &        +(GU*GX - MPRO2*T)/A2/TS/TS
     &        +FAC3/A1/A2/TS/TS
     &        +FAC4/A1/A1/TS/TS
     &        +GS*(GX-GS)*2D0*MEI2/A2/A2/TS/TS
     &        -2D0*MEI2*MPRO2/TS*(1D0/A1/A1+1D0/A2/A2)
     &        -2D0*MPRO2*(1D0/TS/TS+1D0/A1/A2)
     &        +MPRO2/TS*(1D0/A1-1D0/A2) )
          RBLEP =
     &       -((GS-GU)/A1/A2
     &        -2D0*MEI2*(GS-GU)/A1/A2/TS
     &        +(GX/2D0-GS)/A1/TS
     &        +(GX/2D0-GU)/A2/TS
     &        +GX*T/2D0/A1/TS/TS
     &        +GX*T/2D0/A2/TS/TS
     &        +(GX-2D0*GU)*MEI2/A1/A1/TS
     &        +(2D0*GS-GX)*MEI2/A2/A2/TS)
          PLEP=0D0
          PHAD=0D0
          PFRCH=0D0
          GLEP=0D0
          GHAD=0D0
          GFRCH=0D0
          DO 511 IB1=1,2
            DO 512 IB2=1,2
              PLEP=(AFIJ(2,IB1,IB2)*RALEP
     &               +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*RBLEP)
     &             *DSBOS(IB1)*DSBOS(IB2) + PLEP
              PHAD=((AFIJ(2,IB1,IB2)*SU1
     &              -BFIJ(2,IB1,IB2)*SU2*LLEPT)*PPHA1
     &             +(AFIJ(2,IB1,IB2)*SU3
     &              -BFIJ(2,IB1,IB2)*SU4*LLEPT)*PPHA2
     &             +(AFIJ(2,IB1,IB2)*SU5
     &              -BFIJ(2,IB1,IB2)*SU6*LLEPT)*PPHA3)
     &            *DBOS(IB1)*DBOS(IB2)/T/T*4D0/9D0/XS/8D0 + PHAD
              PFRCH=(AFIJ(2,IB1,IB2)*(R1/8D0+R2/4D0*XS)
     &              +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3/4D0)
     &             *DSBOS(IB1)*DSBOS(IB2) + PFRCH
              GLEP=(AGIJ(2,IB1,IB2)*RALEP
     &              +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*RBLEP)
     &             *DSBOS(IB1)*DSBOS(IB2) + GLEP
              GHAD=((AGIJ(2,IB1,IB2)*SU1
     &              -BGIJ(2,IB1,IB2)*SU2*LLEPT)*PPHA1
     &             +(AGIJ(2,IB1,IB2)*SU3
     &              -BGIJ(2,IB1,IB2)*SU4*LLEPT)*PPHA2
     &             +(AGIJ(2,IB1,IB2)*SU5
     &              -BGIJ(2,IB1,IB2)*SU6*LLEPT)*PPHA3)
     &            *DBOS(IB1)*DBOS(IB2)/T/T*4D0/9D0/XS/8D0 + GHAD
              GFRCH=(AGIJ(2,IB1,IB2)*(R1/8D0+R2/4D0*XS)
     &              +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*R3/4D0)
     &             *DSBOS(IB1)*DSBOS(IB2) + GFRCH
 512        CONTINUE
 511      CONTINUE
          SUTOT=PHAD+PLEP
          FRCH=PFRCH/SUTOT
          SGTOT=GHAD+GLEP
          GRCH=GFRCH/SGTOT
        ELSE
          FRCH=0D0
          GRCH=0D0
        ENDIF
        DO 501 IB1=1,2
         DO 502 IB2=1,2
          IF(LPAR(12).NE.0)THEN
             AUP=(CHNULL*AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUP
             AUM=(CHNULL*AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &            -DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUM
             ADP=(CHNULL*AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BFIJ(3,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADP
             ADM=(CHNULL*AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &            -DFLOAT(LLEPT)*BFIJ(3,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADM
             GUP=(CHNULL*AGIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + GUP
             GUM=(-CHNULL*AGIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + GUM
             GDP=(CHNULL*AGIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BGIJ(3,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + GDP
             GDM=(-CHNULL*AGIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BGIJ(3,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + GDM
          ENDIF
         IF(LPAR(14).NE.0)THEN
          PROPI2 = DBOS(IB1)*DSBOS(IB2)/T/TS
          AB12UP = AFIJ(2,IB1,IB2)*SU1*LLEPT-CHNULL*BFIJ(2,IB1,IB2)*SU2
          AB12UM = AFIJ(2,IB1,IB2)*SU1*LLEPT+CHNULL*BFIJ(2,IB1,IB2)*SU2
          AB12DP = AFIJ(3,IB1,IB2)*SU1*LLEPT-CHNULL*BFIJ(3,IB1,IB2)*SU2
          AB12DM = AFIJ(3,IB1,IB2)*SU1*LLEPT+CHNULL*BFIJ(3,IB1,IB2)*SU2
          GB12UP = AGIJ(2,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(2,IB1,IB2)*SU2
          GB12UM =-AGIJ(2,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(2,IB1,IB2)*SU2
          GB12DP = AGIJ(3,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(3,IB1,IB2)*SU2
          GB12DM =-AGIJ(3,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(3,IB1,IB2)*SU2
          FUPI = AB12UP*PROPI1 * PROPI2     + FUPI
          FDOI = AB12DP*PROPI1 * PROPI2     + FDOI
          FUPIB =AB12UM*PROPI1 * PROPI2     + FUPIB
          FDOIB =AB12DM*PROPI1 * PROPI2     + FDOIB
          GUPI = GB12UP*PROPI1 * PROPI2     + GUPI
          GDOI = GB12DP*PROPI1 * PROPI2     + GDOI
          GUPIB =GB12UM*PROPI1 * PROPI2     + GUPIB
          GDOIB =GB12DM*PROPI1 * PROPI2     + GDOIB
         ENDIF
 502     CONTINUE
 501    CONTINUE
        FUPI=FUPI*FRCH
        FDOI=FDOI*FRCH
        FUPIB=FUPIB*FRCH
        FDOIB=FDOIB*FRCH
        GUPI=GUPI*GRCH
        GDOI=GDOI*GRCH
        GUPIB=GUPIB*GRCH
        GDOIB=GDOIB*GRCH

        CQP(1) =  QU *(AUP + FUPI  *  2D0/3D0 )
     &  - HPOLAR*DQU *(GUP + GUPI  *  2D0/3D0 )
        CQP(2) =  QBU*(AUM + FUPIB *(-2D0/3D0)) 
     &  - HPOLAR*DQBU*(GUM + GUPIB *(-2D0/3D0)) + CQP(1)
        CQP(3) =  QD *(ADP + FDOI  *(-1D0/3D0))
     &  - HPOLAR*DQD *(GDP + GDOI  *(-1D0/3D0)) + CQP(2)
        CQP(4) =  QBD*(ADM + FDOIB *  1D0/3D0 )
     &  - HPOLAR*DQBD*(GDM + GDOIB *  1D0/3D0 ) + CQP(3)
        CQP(5) =  QS *(ADP + FDOI  *(-1D0/3D0))
     &  - HPOLAR*DQS *(GDP + GDOI  *(-1D0/3D0)) + CQP(4)
        CQP(6) =  QBS*(ADM + FDOIB *  1D0/3D0 )
     &  - HPOLAR*DQBS*(GDM + GDOIB *  1D0/3D0 ) + CQP(5)
        CQP(7) =  QC *(AUP + FUPI  *  2D0/3D0 )
     &  - HPOLAR*DQC *(GUP + GUPI  *  2D0/3D0 ) + CQP(6)
        CQP(8) =  QBC*(AUM + FUPIB *(-2D0/3D0))
     &  - HPOLAR*DQBC*(GUM + GUPIB *(-2D0/3D0)) + CQP(7)
        CQP(9) =  QB *(ADP + FDOI  *(-1D0/3D0))
     &  - HPOLAR*DQB *(GDP + GDOI  *(-1D0/3D0)) + CQP(8)
        CQP(10) = QBB*(ADM + FDOIB *  1D0/3D0 )
     &  - HPOLAR*DQBB*(GDM + GDOIB *  1D0/3D0 ) + CQP(9)
        CQP(11) = QT *(AUP + FUPI  *  2D0/3D0 )
     &  - HPOLAR*DQT *(GUP + GUPI  *  2D0/3D0 ) + CQP(10)
        CQP(12) = QBT*(AUM + FUPIB *(-2D0/3D0))
     &  - HPOLAR*DQBT*(GUM + GUPIB *(-2D0/3D0)) + CQP(11)
        SUMME=CQP(12)
      ELSEIF (IPDFOP.GE.2) THEN
C..................................PARTON DISTRIBUTION FUNCTIONS
C                                  INCLUDING F_L
        CALL HSSTRF(XS,-TS)
        AUP=0D0
        AUM=0D0
        ADP=0D0
        ADM=0D0
        DO 521 IB1=1,2
         DO 522 IB2=1,2
            AUP=(CHNULL*AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &           +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3)
     &         /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUP
            AUM=(CHNULL*AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &           -DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3)
     &         /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUM
            ADP=(CHNULL*AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &           +DFLOAT(LLEPT)*BFIJ(3,IB1,IB2)*R3)
     &         /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADP
            ADM=(CHNULL*AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &           -DFLOAT(LLEPT)*BFIJ(3,IB1,IB2)*R3)
     &         /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADM
 522     CONTINUE
 521    CONTINUE
        CQP(1) =  QU *AUP
        CQP(2) =  QBU*AUM + CQP(1)
        CQP(3) =  QD *ADP + CQP(2)
        CQP(4) =  QBD*ADM + CQP(3)
        CQP(5) =  QS *ADP + CQP(4)
        CQP(6) =  QBS*ADM + CQP(5)
        CQP(7) =  QC *AUP + CQP(6)
        CQP(8) =  QBC*AUM + CQP(7)
        CQP(9) =  QB *ADP + CQP(8)
        CQP(10) = QBB*ADM + CQP(9)
        CQP(11) = QT *AUP + CQP(10)
        CQP(12) = QBT*AUM + CQP(11)
        SUMME=0D0
        DO 531 IB1=1,2
         DO 532 IB2=1,2
            SUMME=SUMME
     &      +F1(IB1,IB2)*CHNULL*R1*RUNALP(IB1)*RUNALP(IB2)
     &      +F2(IB1,IB2)*CHNULL*R2*RUNALP(IB1)*RUNALP(IB2)
     &      +F3(IB1,IB2)*DFLOAT(LLEPT)*R3*RUNALP(IB1)*RUNALP(IB2)
 532     CONTINUE
 531    CONTINUE
        RNORM=SUMME/CQP(12)
        DO 54 IF=1,12
          CQP(IF)=CQP(IF)*RNORM
 54     CONTINUE
      ENDIF

      HSTSK1 = SUMME*Y*2D0*SX1NRM/SQGRAM
     *      * (UVMAX-UVMIN)*(XS-XX)*(VMAX-VMIN)*A1*(TSMAX-TSMIN)
     *      * (XMAXL-XMINL)* XX * DG2/(XX*GS*GACT**2)
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR COS(THETA(K,P)) FOR KP-PARAMETRIZATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSLZK1(ZMIN,ZMAX)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSPSPC/ IPHSPC
C
      IPHSPC=0
      DEPS = DELTA/OMEGA*( PQH*EH + PH*EQH)
      SIGM = PQH*EEL + PH*EQ
      TAU  = PH*( PS*COSE - PEL ) + PQH*PEL
      F2U = - (TAU + SIGM - DEPS)
      F2O =   (TAU - SIGM + DEPS)
C
      IF ( (F2U.LT.0D0).AND.(F2O.LT.0D0) ) THEN
         ZMIN = -1D0
         ZMAX =  1D0
         RETURN
      ENDIF
C
      EZ = F2U*F2O
      IF ( EZ.LE.0D0 ) THEN
         PPP4 = PH*PH*PS*PS*SINE*SINE
         DZ = 4D0*PPP4*( PPP4 - EZ )
         AZ = TAU*TAU + PPP4
         BZ = - 2D0*TAU*(SIGM - DEPS)
         CZ = (SIGM - DEPS)*(SIGM - DEPS) - PPP4
         IF (TAU.GT.0D0) THEN
            ZMIN = -1D0
            IF (BZ.GT.0D0) THEN
               ZMAX = 2D0*CZ/(-BZ - DSQRT(DZ))
            ELSE
               ZMAX = (- BZ + DSQRT(DZ))/2D0/AZ
            ENDIF
            RETURN
         ELSE
            ZMAX = 1D0
            IF (BZ.GT.0D0) THEN
               ZMIN = (- BZ - DSQRT(DZ))/2D0/AZ
            ELSE
               ZMIN = 2D0*CZ/(-BZ + DSQRT(DZ))
            ENDIF
            RETURN
         ENDIF
      ENDIF
C
      IF ( (F2U.GT.0D0).AND.(F2O.GT.0D0) ) THEN
         PPP4 = PH*PH*PS*PS*SINE*SINE
         DZ = 4D0*PPP4*( PPP4 - EZ )
         IF (DZ.LT.0D0) THEN
            IPHSPC=1
            RETURN
         ELSE
         AZ = TAU*TAU + PPP4
         BZ = - 2D0*TAU*(SIGM - DEPS)
         CZ = (SIGM - DEPS)*(SIGM - DEPS) - PPP4
            IF (BZ.GT.0D0) THEN
               ZMIN = (- BZ - DSQRT(DZ))/2D0/AZ
               ZMAX = 2D0*CZ/(-BZ - DSQRT(DZ))
               if (zmin.ge.zmax.or.zmin.ge.1d0) iphspc=1
               RETURN
            ELSE
               ZMAX = (- BZ + DSQRT(DZ))/2D0/AZ
               ZMIN = 2D0*CZ/(-BZ + DSQRT(DZ))
               if (zmin.ge.zmax.or.zmin.ge.1d0) iphspc=1
               RETURN
            ENDIF
         ENDIF
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR TS FOR KP - PART
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSLTS1(A1,XX,Y,XS,TSMIN,TSMAX,TSM,TSP,CFKP)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSPSPC/ IPHSPC
C
      CALL HSFIVM(XX,Y,XS)
C---TS LIMIT FROM OMH > DELTA : TS > TGR
      TGR=(2D0*DELTA*(PH*EPRO+PPRO*EH)-PPRO*A1)*XS/PH
     &    -(XS-XX)*Y*GS + TP
C---CUT IN HADRONIC MASS
      IF (ICUT.GT.1.AND.WMIN.GT.MPRO) THEN
        W2CUT=WMIN*WMIN
        TWCUT=-XS*(W2CUT-MPRO2)/(1D0-XS)
      ELSE
        TWCUT=0D0
      ENDIF
C
      MXS=XS*MPRO
      MXS2=MXS*MXS
      A=((S+T-3D0*MEF2-MXS2)*(S+T-3D0*MEF2-MXS2)-4D0*MEF2*U)/16D0
      B=(-2D0*MEF2*(S+U-2D0*MEF2-2D0*MXS2)
     &            *(S+U-2D0*MEF2-2D0*MXS2)
     &   +(U+T+S-2D0*MEF2-2D0*MXS2)
     &     *(A1*S-2D0*T*MEF2-(A1+T)*(MEF2+MXS2))
     &   +T*(U*(S+T-A1)+A1*(MEF2-MXS2)
     &       -(MEF2-MXS2)*(MEF2-MXS2)+2D0*T*MEF2))/8D0
      C=((A1*(S+U-2D0*MEF2-2D0*MXS2)-T*(U-MEF2-MXS2))
     &  *(A1*(S+U-2D0*MEF2-2D0*MXS2)-T*(U-MEF2-MXS2))
     &  +4D0*T*T*MXS2*(A1-MEF2)-4D0*T*A1*A1*MXS2)/16D0
C
      CFKP=A
      DISK=1D0/16D0*(MEF2*(S+T+U-2D0*MEF2-2D0*MXS2)
     &                   *(S+T+U-2D0*MEF2-2D0*MXS2)
     &              +A1*(S+T+U-2D0*MEF2-2D0*MXS2)*(A1-S-T+MEF2+MXS2)
     &              +A1*A1*MXS2)
     &      *(MEF2*(S+U-2D0*MEF2-2D0*MXS2)
     &            *(S+U-2D0*MEF2-2D0*MXS2)
     &        +T*(MEF2+MXS2)*(S+U-MEF2-MXS2)
     &        -S*U*T + T*MXS2*(T-4D0*MEF2))
      IF (DISK.LE.0D0) THEN
        DISK = 0D0
      ENDIF
C
      IF (B.GE.0D0) THEN
        T1 = (-B-DSQRT(DISK))/2D0/A
        T2 = C/A/T1
      ELSE
        T2 = (-B+DSQRT(DISK))/2D0/A
        T1 = C/A/T2
      ENDIF
      TSM = DMIN1(T1,T2)
      TSP = DMAX1(T1,T2)
      TSMIN = DMAX1(TSM,TGR)
      TSMAX = DMIN1(TSP,TWCUT)
      IF (TSMIN.GE.TSMAX) THEN
        IPHSPC=1
      ELSE
        IPHSPC=0
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND OF KPS TERM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSTSK2(X)
C
C  X(1) -->  XX
C  X(2) -->  G=-1/Q**2   -->  Q**2=-1/G-->  Y
C  X(3) -->  LOG(XS-XX)
C  X(4) -->  LOG(2*K.PS)
C  X(5) -->  TS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSKPXY/ XX,Y
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSXSLM/ XSMIN,XSCUT
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSASY0/ PONULL,CHNULL,ICPASY
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSFIJK/ F1(2,2),F2(2,2),F3(2,2)
      COMMON /HSGIJK/ G1(2,2),G3(2,2),G4(2,2)
      COMMON /HSPSPC/ IPHSPC
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSWGTC/ IWEIGS
      DIMENSION X(5),DBOS(2),DSBOS(2),RUNALP(2)
      COMPLEX(kind=8) HSSRGG,CG
C
C---X-VALUE
      XMINL=DLOG(XMIN)
      XMAXL=DLOG(XMAX)
      XXL=XMINL+(XMAXL-XMINL)*X(1)
      XX=DEXP(XXL)
C---Y-VALUE
      GS=SP-MEI2-MPRO2
      YMAXX=XX*(1D0-4D0*MEI2*MPRO2/GS/GS)/(XX*(1D0+XX*MPRO2/GS)+MEI2/GS)
      IF(ICUT.LT.3) THEN
C---CUT IN EXTERNALLY DEFINED Q**2(MIN)
        GMAX=-1D0/(XX*YMAXX*GS)
        ELSEIF(ICUT.EQ.3) THEN
C---CUT IN Y
        QQ2MIN=XX*YMIN*GS
C---CUT ON ELECTRON SCATTERING ANGLE (MASSES NEGLECTED)
        YMINTH=1D0/(1D0+XX*CTHCON)
        QT2MIN=XX*YMINTH*GS
C---CUT ON ELECTRON TRANSVERSE MOMENTUM (MASSES NEGLECTED)
        QP2MIN=XX*GS/2D0*(1D0-DSQRT(1D0-PTXM0/XX))
        QP2MAX=XX*GS/2D0*(1D0+DSQRT(1D0-PTXM0/XX))
        YP2MAX=QP2MAX/GS/XX
        GMIN=-1D0/MAX(Q2MIN,QQ2MIN,QT2MIN,QP2MIN)
        GMAX=-1D0/(XX*MIN(YMAX,YMAXX,YP2MAX)*GS)
        GMAX=MIN(GMAX,-1D0/Q2MAX)
        ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSTSK2'
        STOP
      ENDIF
      DG2=DMAX1(GMAX-GMIN,0D0)
C---CUT IN W LATER
      GACT=GMIN+X(2)*DG2
      Y=-1D0/(XX*GACT*GS)
C
C---EXTERNAL WEIGHT
      IACPT=1
      IF (IWEIGS.GT.0) THEN
        CALL HSWGTX(XX,Y,IACPT)
        IF (IACPT.EQ.0) THEN
          HSTSK2=0D0
          RETURN
        ENDIF
      ENDIF
C
      CALL HSFIVC(XX,Y)
      CALL HSDELO(XX,Y)
      XSMIN=HSXSMN(XX,Y)
C     XSCUT=HSXSCT(XX,Y)
C
      IF(IPRINT.GT.30) THEN
        WRITE(LUNTES,230)SP,XX,Y,XSMIN,XSCUT,DELEPS,DELTA
230     FORMAT(' ***************************************************',/
     F        ,' SP = ',1PD12.3,/
     F        ,' X = ',D12.6,'   Y = ',D12.6,/
     F        ,' XSMIN = ',D17.11,'   XSCUT = ',D17.11,/
     F        ,' DELEPS = ',D12.6,'   DELTA = ',D14.8,/
     F       ,' ***************************************************',//)
      ENDIF
C
      XSMAX=1D0
      XSMINI=XSMIN
      IF(XSMAX.LE.XSMINI) THEN
        HSTSK2=0D0
        RETURN
      ENDIF
C--- SUBSTITUTION UV=LN(XS-X)
      UVMIN=DLOG(XSMINI-XX)
      UVMAX=DLOG(XSMAX-XX)
      UV=UVMIN+(UVMAX-UVMIN)*X(3)
      XS=XX+DEXP(UV)
      CALL HSFCMS(XX,Y,XS)
      IF(IPHSPC.EQ.1) THEN
        HSTSK2=0D0
        RETURN
      ENDIF
      CALL HSFLAB(XS)
      CALL HSLZK2(ZMIN,ZMAX)
      IF(IPHSPC.EQ.1) THEN
        HSTSK2=0D0
        RETURN
      ENDIF
C---SUBSTITUTION V = LN(A2)
      A2MAX=2D0*OMEGA*(ES-PS*ZMIN)
      IF ((ZMAX.GE.0.9999D0).AND.(ES/MEF.GT.1D3)) THEN
        A2MIN=2D0*OMEGA*MEF2/2D0/ES
      ELSE
        A2MIN = 2D0*OMEGA*(ES-PS*ZMAX)
      ENDIF
      VMAX=DLOG(A2MAX)
      VMIN=DLOG(A2MIN)
      V=VMIN+(VMAX-VMIN)*X(4)
      A2=DEXP(V)
C---NO SUBSTITUTION FOR = TS
      CALL HSLTS2(A2,XX,Y,XS,TSMIN,TSMAX,TSM,TSP,CFKPS)
      IF(IPHSPC.EQ.1)THEN
       HSTSK2=0D0
       RETURN
      ENDIF
      TS=TSMIN+(TSMAX-TSMIN)*X(5)
      SQGRM2=-CFKPS*(TS-TSM)*(TS-TSP)
      IF (SQGRM2.LE.0D0) THEN
        HSTSK2=0D0
        RETURN
      ENDIF
      SQGRAM=DSQRT(DABS(SQGRM2))
      CALL HSFIV2(XX,Y,XS,A2,TS)
C
      A1=2D0*DKP
      DBOS(1)=1D0
      DBOS(2)=T/(T-MZ2)
      DSBOS(1)=1D0
      RUNALP(1)=1D0
      IF (LPAR(3).GE.3) THEN
        CG=HSSRGG(TS)
        DSBOS(1)=1D0/(1D0+DREAL(CG)/TS)
        RUNALP(1)=DSBOS(1)
      ENDIF
      DSBOS(2)=TS/(TS-MZ2)
      RUNALP(2)=1D0
      IF(LPAR(14).EQ.0)THEN
        SU1=0D0
        SU2=0D0
        PROPI1=0D0
      ELSE
        SU1 = S*S+SS*SS+U*U+US*US
        SU2 = (S-U)*(S+U) + (SS-US)*(SS+US)
        PROPI1= - 1D0/8D0/XS
     &       * (SS/DKPS/DKQS + S/DKP/DKQ + U/DKPS/DKQ + US/DKP/DKQS)
      ENDIF

      R1=4D0*GX*(
     &      (T+2D0*MEI2)/(A2-TS)/(A2-TS)/A2
     &     -1D0/(A2-TS)/A2
     &     + 2D0/(A1+A2)/A2
     &     -8D0*MEI2*MEI2/(A1*A2+TS*TS)/(A1+A2)/A2
     &     +4D0*MEI2*MEI2/(A2*A2+TS*TS)/A2/A2
     &     -2D0*MEI2/(A2-TS)/A2/A2 )
      FAC1 =-(GS*GS + GU*GU - GX*(GS+GU) -4D0*MEI2*MPRO2)
      FAC2 = GU*GX - MPRO2*T
      FAC3 = (-2D0*GS*GU + GX*(GS+GU-GX))*2D0*MEI2
      FAC4 = GS*(GX-GS)*2D0*MEI2
      R2=4D0*(
     &     -FAC1/(A1+A2-TS)*(1D0/(A1+A2)+1D0/(A2-TS))/A2
     &     +FAC2/(A2-TS)/(A2-TS)/A2
     &     +FAC3/(A1*A2+TS*TS)/(A1+A2)/A2
     &     -2D0*MPRO2/(A1+A2)/A2
     &     +2D0*MEI2*MPRO2/(A2-TS)/(A2-TS)/A2
     &     +2D0*MEI2*MPRO2/(A2-TS)/A2/A2
     &     +MPRO2/(A2-TS)/A2
     &     +FAC4/(A2*A2+TS*TS)/A2/A2 )
      FK1=-2D0*MEI2*(GS-GU)
      R3= -4D0*(
     &     (GS-GU)/(A1+A2)/A2
     &     +FK1/(A1+A2-TS)*(-1D0/(A1+A2)-1D0/(A2-TS))/A2
     &     -(GX/2D0-GU)/(A2-TS)/A2
     &     +GX*T/2D0/(A2-TS)/(A2-TS)/A2
     &     -(2D0*GS-GX)*MEI2/(A2-TS)/(A2-TS)/A2
     &     -(2D0*GS-GX)*MEI2/(A2-TS)/A2/A2  )
      DO 20 IFL=1,12
        CQP(IFL)=0D0
 20   CONTINUE

      SUMME=0D0
      IF (IPDFOP.EQ.0) THEN
C..................................STRUCTURE FUNCTIONS
         CALL HSSTRF(XS,-TS)
         DO 11 IB1=1,2
          DO 12 IB2=1,2
             CQP(12)=CQP(12)
     &      +CHNULL*F1(IB1,IB2)*R1*RUNALP(IB1)*RUNALP(IB2)
     &      +CHNULL*F2(IB1,IB2)*R2*RUNALP(IB1)*RUNALP(IB2)
     &      +F3(IB1,IB2)*DFLOAT(LLEPT)*R3*RUNALP(IB1)*RUNALP(IB2)
     &     -HPOLAR*(
     &      +G1(IB1,IB2)*DFLOAT(LLEPT)*R3*2D0*RUNALP(IB1)*RUNALP(IB2)
     &      +CHNULL*G4(IB1,IB2)*R2*RUNALP(IB1)*RUNALP(IB2)
     &      +CHNULL*G3(IB1,IB2)*R1*RUNALP(IB1)*RUNALP(IB2))
 12       CONTINUE
 11      CONTINUE
        SUMME = CQP(12)
      ELSEIF (IPDFOP.EQ.1) THEN
C..................................PARTON DENSITIES
        CALL HSPVER(XS,-TS)
        IF (HPOLAR.NE.0D0) THEN
          CALL HSDPVR(XS,-TS)
        ENDIF
        CALL HSIROT
        AUP=0D0
        AUM=0D0
        ADP=0D0
        ADM=0D0
        GUP=0D0
        GUM=0D0
        GDP=0D0
        GDM=0D0
        FUPI = 0D0
        FDOI = 0D0
        FUPIB =0D0
        FDOIB =0D0
        GUPI = 0D0
        GDOI = 0D0
        GUPIB =0D0
        GDOIB =0D0
        IF(LPAR(14).NE.0)THEN
          SU3 = (S*S+U*U)/2D0
          SU4 = (S-U)*(S+U)/2D0
          SU5 = (SS*SS+US*US)/2D0
          SU6 = (SS-US)*(SS+US)/2D0
          PPHA1= -T/DKQ/DKQS
          PPHA2= -4D0*MQF2/DKQS/DKQS
          PPHA3= -4D0*MQI2/DKQ/DKQ
          RALEP =
     &      GX/2D0*(
     &        +(T+4D0*MEI2)*(1D0/A2/TS/TS-1D0/A1/TS/TS)
     &        +2D0/TS/TS - 1D0/A1/TS + 1D0/A2/TS + 2D0/A1/A2
     &        -8D0*MEI2*MEI2/A1/A2/TS/TS
     &        +4D0*MEI2*MEI2/TS/TS*(1D0/A1/A1+1D0/A2/A2)
     &        +2D0*MEI2/TS*(1D0/A1/A1+1D0/A2/A2)  )
     &     +XS*(
     &        +FAC1/TS/A1/A2
     &        +(-GS*GX + MPRO2*T)/A1/TS/TS
     &        +FAC2/A2/TS/TS
     &        +FAC3/A1/A2/TS/TS
     &        +GU*(GX-GU)*2D0*MEI2/A1/A1/TS/TS
     &        +FAC4/A2/A2/TS/TS
     &        -2D0*MEI2*MPRO2/TS*(1D0/A1/A1+1D0/A2/A2)
     &        -2D0*MPRO2*(1D0/TS/TS+1D0/A1/A2)
     &        +MPRO2/TS*(1D0/A1-1D0/A2) )
          RBLEP =
     &       -((GS-GU)/A1/A2
     &       -2D0*MEI2*(GS-GU)/A1/A2/TS
     &       +(GX/2D0-GS)/A1/TS
     &       +(GX/2D0-GU)/A2/TS
     &       +GX*T/2D0/A1/TS/TS
     &       +GX*T/2D0/A2/TS/TS
     &       +(GX-2D0*GU)*MEI2/A1/A1/TS
     &       +(2D0*GS-GX)*MEI2/A2/A2/TS)
          PLEP=0D0
          PHAD=0D0
          PFRCH=0D0
          GLEP=0D0
          GHAD=0D0
          GFRCH=0D0
          DO 511 IB1=1,2
            DO 512 IB2=1,2
              PLEP=(AFIJ(2,IB1,IB2)*RALEP
     &               +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*RBLEP)
     &             *DSBOS(IB1)*DSBOS(IB2) + PLEP
              PHAD =((AFIJ(2,IB1,IB2)*SU1
     &               -BFIJ(2,IB1,IB2)*SU2*LLEPT)*PPHA1
     &              +(AFIJ(2,IB1,IB2)*SU3
     &               -BFIJ(2,IB1,IB2)*SU4*LLEPT)*PPHA2
     &              +(AFIJ(2,IB1,IB2)*SU5
     &               -BFIJ(2,IB1,IB2)*SU6*LLEPT)*PPHA3)
     &             *DBOS(IB1)*DBOS(IB2)/T/T*4D0/9D0/XS/8D0 + PHAD
              PFRCH=(AFIJ(2,IB1,IB2)*(R1/8D0+R2/4D0*XS)
     &              +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3/4D0)
     &             *DSBOS(IB1)*DSBOS(IB2) + PFRCH
              GLEP=(AGIJ(2,IB1,IB2)*RALEP
     &              +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*RBLEP)
     &             *DSBOS(IB1)*DSBOS(IB2) + GLEP
              GHAD =((AGIJ(2,IB1,IB2)*SU1
     &               -BGIJ(2,IB1,IB2)*SU2*LLEPT)*PPHA1
     &              +(AGIJ(2,IB1,IB2)*SU3
     &               -BGIJ(2,IB1,IB2)*SU4*LLEPT)*PPHA2
     &              +(AGIJ(2,IB1,IB2)*SU5
     &               -BGIJ(2,IB1,IB2)*SU6*LLEPT)*PPHA3)
     &             *DBOS(IB1)*DBOS(IB2)/T/T*4D0/9D0/XS/8D0 + GHAD
              GFRCH=(AGIJ(2,IB1,IB2)*(R1/8D0+R2/4D0*XS)
     &              +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*R3/4D0)
     &             *DSBOS(IB1)*DSBOS(IB2) + GFRCH
 512        CONTINUE
 511      CONTINUE
          SUTOT=PHAD+PLEP
          FRCH=PFRCH/SUTOT
          SGTOT=GHAD+GLEP
          GRCH=GFRCH/SGTOT
        ELSE
          FRCH=0D0
          GRCH=0D0
        ENDIF
        DO 501 IB1=1,2
         DO 502 IB2=1,2
          IF(LPAR(12).NE.0)THEN
             AUP=(CHNULL*AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUP
             AUM=(CHNULL*AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &            -DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUM
             ADP=(CHNULL*AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BFIJ(3,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADP
             ADM=(CHNULL*AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &            -DFLOAT(LLEPT)*BFIJ(3,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADM
             GUP=(CHNULL*AGIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + GUP
             GUM=(-CHNULL*AGIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + GUM
             GDP=(CHNULL*AGIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BGIJ(3,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + GDP
             GDM=(-CHNULL*AGIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BGIJ(3,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + GDM
          ENDIF
         IF(LPAR(14).NE.0)THEN
          PROPI2 = DBOS(IB1)*DSBOS(IB2)/T/TS
          AB12UP = AFIJ(2,IB1,IB2)*SU1*LLEPT-CHNULL*BFIJ(2,IB1,IB2)*SU2
          AB12UM = AFIJ(2,IB1,IB2)*SU1*LLEPT+CHNULL*BFIJ(2,IB1,IB2)*SU2
          AB12DP = AFIJ(3,IB1,IB2)*SU1*LLEPT-CHNULL*BFIJ(3,IB1,IB2)*SU2
          AB12DM = AFIJ(3,IB1,IB2)*SU1*LLEPT+CHNULL*BFIJ(3,IB1,IB2)*SU2
          GB12UP = AGIJ(2,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(2,IB1,IB2)*SU2
          GB12UM =-AGIJ(2,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(2,IB1,IB2)*SU2
          GB12DP = AGIJ(3,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(3,IB1,IB2)*SU2
          GB12DM =-AGIJ(3,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(3,IB1,IB2)*SU2
          FUPI = AB12UP*PROPI1 * PROPI2      + FUPI
          FDOI = AB12DP*PROPI1 * PROPI2      + FDOI
          FUPIB =AB12UM*PROPI1 * PROPI2      + FUPIB
          FDOIB =AB12DM*PROPI1 * PROPI2      + FDOIB
          GUPI = GB12UP*PROPI1 * PROPI2      + GUPI
          GDOI = GB12DP*PROPI1 * PROPI2      + GDOI
          GUPIB =GB12UM*PROPI1 * PROPI2      + GUPIB
          GDOIB =GB12DM*PROPI1 * PROPI2      + GDOIB
         ENDIF
 502     CONTINUE
 501    CONTINUE
        FUPI=FUPI*FRCH
        FDOI=FDOI*FRCH
        FUPIB=FUPIB*FRCH
        FDOIB=FDOIB*FRCH
        GUPI=GUPI*GRCH
        GDOI=GDOI*GRCH
        GUPIB=GUPIB*GRCH
        GDOIB=GDOIB*GRCH
        
        CQP(1) =  QU *(AUP + FUPI  *  2D0/3D0 )
     &  - HPOLAR*DQU *(GUP + GUPI  *  2D0/3D0 )
        CQP(2) =  QBU*(AUM + FUPIB *(-2D0/3D0)) 
     &  - HPOLAR*DQBU*(GUM + GUPIB *(-2D0/3D0)) + CQP(1)
        CQP(3) =  QD *(ADP + FDOI  *(-1D0/3D0))
     &  - HPOLAR*DQD *(GDP + GDOI  *(-1D0/3D0)) + CQP(2)
        CQP(4) =  QBD*(ADM + FDOIB *  1D0/3D0 )
     &  - HPOLAR*DQBD*(GDM + GDOIB *  1D0/3D0 ) + CQP(3)
        CQP(5) =  QS *(ADP + FDOI  *(-1D0/3D0))
     &  - HPOLAR*DQS *(GDP + GDOI  *(-1D0/3D0)) + CQP(4)
        CQP(6) =  QBS*(ADM + FDOIB *  1D0/3D0 )
     &  - HPOLAR*DQBS*(GDM + GDOIB *  1D0/3D0 ) + CQP(5)
        CQP(7) =  QC *(AUP + FUPI  *  2D0/3D0 )
     &  - HPOLAR*DQC *(GUP + GUPI  *  2D0/3D0 ) + CQP(6)
        CQP(8) =  QBC*(AUM + FUPIB *(-2D0/3D0))
     &  - HPOLAR*DQBC*(GUM + GUPIB *(-2D0/3D0)) + CQP(7)
        CQP(9) =  QB *(ADP + FDOI  *(-1D0/3D0))
     &  - HPOLAR*DQB *(GDP + GDOI  *(-1D0/3D0)) + CQP(8)
        CQP(10) = QBB*(ADM + FDOIB *  1D0/3D0 )
     &  - HPOLAR*DQBB*(GDM + GDOIB *  1D0/3D0 ) + CQP(9)
        CQP(11) = QT *(AUP + FUPI  *  2D0/3D0 )
     &  - HPOLAR*DQT *(GUP + GUPI  *  2D0/3D0 ) + CQP(10)
        CQP(12) = QBT*(AUM + FUPIB *(-2D0/3D0))
     &  - HPOLAR*DQBT*(GUM + GUPIB *(-2D0/3D0)) + CQP(11)
        SUMME=CQP(12)
      ELSEIF (IPDFOP.GE.2) THEN
C..................................PARTON DISTRIBUTION FUNCTIONS
C                                  INCLUDING F_L
        CALL HSSTRF(XS,-TS)
        AUP=0D0
        AUM=0D0
        ADP=0D0
        ADM=0D0
        DO 521 IB1=1,2
         DO 522 IB2=1,2
            AUP=(CHNULL*AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &           +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3)
     &         /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUP
            AUM=(CHNULL*AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &           -DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3)
     &         /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUM
            ADP=(CHNULL*AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &           +DFLOAT(LLEPT)*BFIJ(3,IB1,IB2)*R3)
     &         /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADP
            ADM=(CHNULL*AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &           -DFLOAT(LLEPT)*BFIJ(3,IB1,IB2)*R3)
     &         /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADM
 522     CONTINUE
 521    CONTINUE
        CQP(1) =  QU *AUP
        CQP(2) =  QBU*AUM + CQP(1)
        CQP(3) =  QD *ADP + CQP(2)
        CQP(4) =  QBD*ADM + CQP(3)
        CQP(5) =  QS *ADP + CQP(4)
        CQP(6) =  QBS*ADM + CQP(5)
        CQP(7) =  QC *AUP + CQP(6)
        CQP(8) =  QBC*AUM + CQP(7)
        CQP(9) =  QB *ADP + CQP(8)
        CQP(10) = QBB*ADM + CQP(9)
        CQP(11) = QT *AUP + CQP(10)
        CQP(12) = QBT*AUM + CQP(11)
        SUMME=0D0
        DO 531 IB1=1,2
         DO 532 IB2=1,2
            SUMME=SUMME
     &     +F1(IB1,IB2)*CHNULL*R1*RUNALP(IB1)*RUNALP(IB2)
     &     +F2(IB1,IB2)*CHNULL*R2*RUNALP(IB1)*RUNALP(IB2)
     &     +F3(IB1,IB2)*DFLOAT(LLEPT)*R3*RUNALP(IB1)*RUNALP(IB2)
 532     CONTINUE
 531    CONTINUE
        RNORM=SUMME/CQP(12)
        DO 54 IF=1,12
          CQP(IF)=CQP(IF)*RNORM
 54     CONTINUE
      ENDIF

      HSTSK2 = SUMME*Y*2D0*SX1NRM/SQGRAM
     *      * (UVMAX-UVMIN)*(XS-XX)*(VMAX-VMIN)*A2*(TSMAX-TSMIN)
     *      * (XMAXL-XMINL) * XX * DG2/(XX*GS*GACT**2)
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR COS(THETA(K,PS)) FOR KPS-PARAMETRIZATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSLZK2(ZMIN,ZMAX)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSPSPC/ IPHSPC
C
      IPHSPC=0
      DEPS = DELTA/OMEGA*( PQH*EH + PH*EQH)
      SIGM = PQH*EEL + PH*EQ
      TAU  = PQH*PEL*COSE - PH*(PEL*COSE-PS)
      F2U = - (TAU + SIGM - DEPS)
      F2O =   (TAU - SIGM + DEPS)
C
      IF ( (F2U.LT.0D0).AND.(F2O.LT.0D0) ) THEN
         ZMIN = -1D0
         ZMAX =  1D0
         RETURN
      ENDIF
C
      EZ = F2U*F2O
      IF ( EZ.LE.0D0 ) THEN
         PPP4 = PEL*PEL*SINE*SINE*(PQH-PH)*(PQH-PH)
         DZ = 4D0*PPP4*( PPP4 - EZ )
         AZ = TAU*TAU + PPP4
         BZ = - 2D0*TAU*(SIGM - DEPS)
         CZ = (SIGM - DEPS)*(SIGM - DEPS) - PPP4
         IF (TAU.GT.0D0) THEN
            ZMIN = -1D0
            IF (BZ.GT.0D0) THEN
               ZMAX = 2D0*CZ/(-BZ - DSQRT(DZ))
            ELSE
               ZMAX = (- BZ + DSQRT(DZ))/2D0/AZ
            ENDIF
            RETURN
         ELSE
            ZMAX = 1D0
            IF (BZ.GT.0D0) THEN
               ZMIN = (- BZ - DSQRT(DZ))/2D0/AZ
            ELSE
               ZMIN = 2D0*CZ/(-BZ + DSQRT(DZ))
            ENDIF
            RETURN
         ENDIF
      ENDIF
C
      IF ( (F2U.GT.0D0).AND.(F2O.GT.0D0) ) THEN
         PPP4 = PEL*PEL*SINE*SINE*(PQH-PH)*(PQH-PH)
         DZ = 4D0*PPP4*( PPP4 - EZ )
         IF (DZ.LT.0D0) THEN
            IPHSPC=1
            RETURN
         ELSE
         AZ = TAU*TAU + PPP4
         BZ = - 2D0*TAU*(SIGM - DEPS)
         CZ = (SIGM - DEPS)*(SIGM - DEPS) - PPP4
            IF (BZ.GT.0D0) THEN
               ZMIN = (- BZ - DSQRT(DZ))/2D0/AZ
               ZMAX = 2D0*CZ/(-BZ - DSQRT(DZ))
               RETURN
            ELSE
               ZMAX = (- BZ + DSQRT(DZ))/2D0/AZ
               ZMIN = 2D0*CZ/(-BZ + DSQRT(DZ))
               RETURN
            ENDIF
         ENDIF
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR TS FOR KPS - PART
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSLTS2(A2,XX,Y,XS,TSMIN,TSMAX,TSM,TSP,CFKPS)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSPSPC/ IPHSPC
C
      CALL HSFIVM(XX,Y,XS)
C---TS LIMIT FROM OMH > DELTA
      TGR=2D0*XS*DELTA*(PH*EPRO+PPRO*EH)
     &   -XS*PPRO*(TP+A2)-PH*((XS-XX)*Y*GS-TP)
C---CUT IN HADRONIC MASS
      IF (ICUT.GT.1.AND.WMIN.GT.MPRO) THEN
        W2CUT=WMIN*WMIN
        TWCUT=-XS*(W2CUT-MPRO2)/(1D0-XS)
      ELSE
        TWCUT=0D0
      ENDIF
C
      MXS=XS*MPRO
      MXS2=MXS*MXS
      A=((U+T-3D0*MEF2-MXS2)*(U+T-3D0*MEF2-MXS2)-4D0*MEF2*S)/16D0
      B=(-2D0*MEF2*(U+S-2D0*MEF2-2D0*MXS2)*(U+S-2D0*MEF2-2D0*MXS2)
     &   +(U+S+T-2D0*MEF2-2D0*MXS2)
     &     *((A2-T)*(MEF2+MXS2)-2D0*T*MEF2-U*A2)
     &   +T*(S*(U+T+A2)-(MEF2-MXS2)*(MEF2-MXS2)
     &       -A2*(MEF2-MXS2)+2D0*T*MEF2))/8D0
      C=((A2*(U+S-2D0*MEF2-2D0*MXS2)+T*(S-MEF2-MXS2))
     &  *(A2*(U+S-2D0*MEF2-2D0*MXS2)+T*(S-MEF2-MXS2))
     &  -4D0*T*T*MXS2*(A2+MEF2) -4D0*T*A2*A2*MXS2)/16D0
      DISK=1D0/16D0*(MEF2*(S+U+T-2D0*MEF2-2D0*MXS2)
     &                   *(S+U+T-2D0*MEF2-2D0*MXS2)
     &               +A2*(S+U+T-2D0*MEF2-2D0*MXS2)*(A2+U+T-MEF2-MXS2)
     &               +A2*A2*MXS2)
     &       *(MEF2*(S+U-2D0*MEF2-2D0*MXS2)
     &             *(S+U-2D0*MEF2-2D0*MXS2)
     &         +T*(MEF2+MXS2)*(S+U-MEF2-MXS2)
     &         -U*T*S+T*MXS2*(T-4D0*MEF2))
      CFKPS=A
      IF (DISK.LE.0D0) THEN
        DISK=0D0
      ENDIF
C
      IF (B.GE.0D0) THEN
        T1=(-B-DSQRT(DISK))/2D0/A
        T2=C/A/T1
      ELSE
        T2=(-B+DSQRT(DISK))/2D0/A
        T1=C/A/T2
      ENDIF
      TSM=DMIN1(T1,T2)
      TSP=DMAX1(T1,T2)
C
      TSMIN=TSM
      TSMAX=TSP
c      write(*,*) ' in HSLTS2: '
c      write(*,*) ' tsmin,tsmax = ',tsmin,tsmax
C
      IF (PH.GT.(XS*PPRO)) THEN
        TGR=TGR/(PH-XS*PPRO)
        TSMIN=DMAX1(TSMIN,TGR)
      ENDIF
      IF (PH.LT.(XS*PPRO)) THEN
        TGR=TGR/(PH-XS*PPRO)
        TSMAX=DMIN1(TGR,TSMAX)
      ENDIF
      TSMAX=DMIN1(TSMAX,TWCUT)
C
      IF (TSMIN.GE.TSMAX
     .    .OR.TSMAX.GE.0D0) THEN
        IPHSPC=1
      ELSE
        IPHSPC=0
      ENDIF
c      write(*,*) ' 2: tsmin,tsmax = ',tsmin,tsmax
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND OF TS TERM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSK1TS(X)
C
C  X(1) -->  XX
C  X(2) -->  G=-1/Q**2   -->  Q**2=-1/G-->  Y
C  X(3) -->  LOG(XS-XX)
C  X(4) -->  LOG(-TS)
C  X(5) -->  2*K.P
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSKPXY/ XX,Y
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSXSLM/ XSMIN,XSCUT
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSASY0/ PONULL,CHNULL,ICPASY
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSFIJK/ F1(2,2),F2(2,2),F3(2,2)
      COMMON /HSGIJK/ G1(2,2),G3(2,2),G4(2,2)
      COMMON /HSPSPC/ IPHSPC
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSWGTC/ IWEIGS
      DIMENSION X(5),DBOS(2),DSBOS(2),RUNALP(2)
      COMPLEX(kind=8) HSSRGG,CG
C
C---X-VALUE
      XX=XMIN + (XMAX-XMIN)*X(1)
C---Y-VALUE
      GS=SP-MEI2-MPRO2
      YMAXX=XX*(1D0-4D0*MEI2*MPRO2/GS/GS)/(XX*(1D0+XX*MPRO2/GS)+MEI2/GS)
      IF(ICUT.LT.3) THEN
C---CUT IN EXTERNALLY DEFINED Q**2(MIN)
        GMAX=-1D0/(XX*YMAXX*GS)
        ELSEIF(ICUT.EQ.3) THEN
C---CUT IN Y
        QQ2MIN=XX*YMIN*GS
C---CUT ON ELECTRON SCATTERING ANGLE (MASSES NEGLECTED)
        YMINTH=1D0/(1D0+XX*CTHCON)
        QT2MIN=XX*YMINTH*GS
C---CUT ON ELECTRON TRANSVERSE MOMENTUM (MASSES NEGLECTED)
        QP2MIN=XX*GS/2D0*(1D0-DSQRT(1D0-PTXM0/XX))
        QP2MAX=XX*GS/2D0*(1D0+DSQRT(1D0-PTXM0/XX))
        YP2MAX=QP2MAX/GS/XX
        GMIN=-1D0/MAX(Q2MIN,QQ2MIN,QT2MIN,QP2MIN)
        GMAX=-1D0/(XX*MIN(YMAX,YMAXX,YP2MAX)*GS)
        GMAX=MIN(GMAX,-1D0/Q2MAX)
        ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSK1TS'
        STOP
      ENDIF
      DG2=DMAX1(GMAX-GMIN,0D0)
CUT IN W LATER
      GACT=GMIN+X(2)*DG2
      Y=-1D0/(XX*GACT*GS)
C
C---EXTERNAL WEIGHT
      IACPT=1
      IF (IWEIGS.GT.0) THEN
        CALL HSWGTX(XX,Y,IACPT)
        IF (IACPT.EQ.0) THEN
          HSK1TS=0D0
          RETURN
        ENDIF
      ENDIF
C
      CALL HSFIVC(XX,Y)
      CALL HSDELO(XX,Y)
      XSMIN = HSXSMN(XX,Y)
C     XSCUT = HSXSCT(XX,Y)
C
      IF(IPRINT.GT.30) THEN
        WRITE(LUNTES,230)SP,XX,Y,XSMIN,XSCUT,DELEPS,DELTA
230     FORMAT(' ***************************************************',/
     F        ,' SP = ',1PD12.3,/
     F        ,' X = ',D12.6,'   Y = ',D12.6,/
     F        ,' XSMIN = ',D17.11,'   XSCUT = ',D17.11,/
     F        ,' DELEPS = ',D12.6,'   DELTA = ',D14.8,/
     F       ,' ***************************************************',//)
      ENDIF
C
      XSMAX = 1D0
      XSMINI = XSMIN
      IF((XSMAX-XSMINI).LE.1D-14) THEN
        HSK1TS=0D0
        RETURN
      ENDIF
C
C--- SUBSTITUTION UV=LN(XS-X)
      UVMIN = DLOG(XSMINI-XX)
      UVMAX = DLOG(XSMAX-XX)
      UV = UVMIN + (UVMAX-UVMIN)*X(3)
      XS = XX+DEXP(UV)
      CALL HSFCMS(XX,Y,XS)
      IF(IPHSPC.EQ.1) THEN
        HSK1TS=0D0
        RETURN
      ENDIF
      CALL HSFLAB(XS)
      CALL HSLZTS(ZMIN,ZMAX)
      IF(IPHSPC.EQ.1) THEN
        HSK1TS=0D0
        RETURN
      ENDIF
      TSMIN=TP+2D0*OMEGA*(ES-EEL-PQ*ZMAX)
      IF(ZMIN.LT.-0.9999D0)THEN
        TSMAX=-XX*XX*XS*MPRO2/(XS-XX+XS*XS*MPRO2/Y/GS)
      ELSE
        TSMAX=TP+2D0*OMEGA*(ES-EEL-PQ*ZMIN)
      ENDIF
C---CUT ON HADRONIC MASS
      IF (ICUT.GT.1) THEN
        W2CUT=WMIN*WMIN
        TWCUT=-XS*(W2CUT-MPRO2)/(1D0-XS)
        TSMAX=DMIN1(TSMAX,TWCUT)
      ENDIF
C---SUBSTITUTION R = LN(-TS)
      RMAX=DLOG(-TSMIN)
      RMIN=DLOG(-TSMAX)
      R=RMIN+(RMAX-RMIN)*X(4)
      TS=-DEXP(R)
C---NO SUBSTITUTION FOR A1
      CALL HSL1TS(TS,XX,Y,XS,A1MIN,A1MAX,A1M,A1P,CFKP)
      IF(IPHSPC.EQ.1)THEN
       HSK1TS=0D0
       RETURN
      ENDIF
      A1=A1MIN+(A1MAX-A1MIN)*X(5)
      SQGRM2=-CFKP*(A1-A1M)*(A1-A1P)
      IF (SQGRM2.LE.0D0) THEN
        HSK1TS=0D0
        RETURN
      ENDIF
      SQGRAM=DSQRT(DABS(SQGRM2))
      CALL HSFIV1(XX,Y,XS,A1,TS)
C
      A2=2D0*DKPS
      DBOS(1)=1D0
      DBOS(2)=T/(T-MZ2)
      DSBOS(1)=1D0
      RUNALP(1)=1D0
      IF (LPAR(3).GE.3) THEN
        CG=HSSRGG(TS)
        DSBOS(1)=1D0/(1D0+DREAL(CG)/TS)
        RUNALP(1)=DSBOS(1)
      ENDIF
      DSBOS(2)=TS/(TS-MZ2)
      RUNALP(2)=1D0
      IF(LPAR(14).EQ.0)THEN
        SU1=0D0
        SU2=0D0
        PROPI1=0D0
      ELSE
        SU1=S*S+SS*SS+U*U+US*US
        SU2=(S-U)*(S+U)+(SS-US)*(SS+US)
        PROPI1=-1D0/8D0/XS
     &       * (SS/DKPS/DKQS + S/DKP/DKQ + U/DKPS/DKQ + US/DKP/DKQS)
      ENDIF

      R1=4D0*GX*(
     &     +(T+6D0*MEI2)/(A1-TS)/(A1-TS)/TS
     &     -(T+4D0*MEI2)/(A1-TS)/TS/TS
     &     -(T+2D0*MEI2)/(A2-TS)/(A2-TS)/TS
     &     -(T+4D0*MEI2)/(A2-TS)/TS/TS
     &     +2D0/TS/TS
     &     -1D0/(A1-TS)/TS
     &     +1D0/(A2-TS)/TS
     &     -8D0*MEI2*MEI2/(A1*A2+TS*TS)/TS/TS
     &     +4D0*MEI2*MEI2/(A1*A1+TS*TS)/TS/TS
     &     +4D0*MEI2*MEI2/(A2*A2+TS*TS)/TS/TS )
 
      FAC1 =-(GS*GS + GU*GU - GX*(GS+GU) -4D0*MEI2*MPRO2)
      FAC2 = -GS*GX + MPRO2*T
      FAC3 = GU*GX - MPRO2*T
      FAC4 = (-2D0*GS*GU + GX*(GS+GU-GX))*2D0*MEI2
      FAC5 = GU*(GX-GU)*2D0*MEI2
      FAC6 = GS*(GX-GS)*2D0*MEI2
      R2=4D0*(
     &      FAC1/(A1+A2-TS)*(1D0/(A1-TS)+1D0/(A2-TS))/TS
     &     -FAC2/(A1-TS)/(A1-TS)/TS
     &     +FAC2/(A1-TS)/TS/TS
     &     -FAC3/(A2-TS)/(A2-TS)/TS
     &     +FAC3/(A2-TS)/TS/TS
     &     +FAC4/(A1*A2+TS*TS)/TS/TS
     &     -2D0*MEI2*MPRO2/(A1-TS)/(A1-TS)/TS
     &     -2D0*MEI2*MPRO2/(A2-TS)/(A2-TS)/TS
     &     -2D0*MPRO2/TS/TS
     &     +MPRO2/(A1-TS)/TS
     &     -MPRO2/(A2-TS)/TS
     &     +FAC5/(A1*A1+TS*TS)/TS/TS
     &     +FAC6/(A2*A2+TS*TS)/TS/TS)
      FK1=-2D0*MEI2*(GS-GU)
      R3= -4D0*(
     &     +FK1/(A1+A2-TS)*(1D0/(A1-TS)+1D0/(A2-TS))/TS
     &     +(GX/2D0-GS)/(A1-TS)/TS
     &     +(GX/2D0-GU)/(A2-TS)/TS
     &     -GX*T/2D0/(A1-TS)/(A1-TS)/TS
     &     +GX*T/2D0/(A1-TS)/TS/TS
     &     -GX*T/2D0/(A2-TS)/(A2-TS)/TS
     &     +GX*T/2D0/(A2-TS)/TS/TS
     &     +(GX-2D0*GU)*MEI2/(A1-TS)/(A1-TS)/TS
     &     +(2D0*GS-GX)*MEI2/(A2-TS)/(A2-TS)/TS )
      DO 20 IFL=1,12
        CQP(IFL)=0D0
 20   CONTINUE

      SUMME=0D0
      IF (IPDFOP.EQ.0) THEN
C..................................STRUCTURE FUNCTIONS
         CALL HSSTRF(XS,-TS)
         DO 10 IB1=1,2
          DO 11 IB2=1,2
             CQP(12)=CQP(12)
     &         +CHNULL*F1(IB1,IB2)*R1*RUNALP(IB1)*RUNALP(IB2)
     &         +CHNULL*F2(IB1,IB2)*R2*RUNALP(IB1)*RUNALP(IB2)
     &         +F3(IB1,IB2)*DFLOAT(LLEPT)*R3*RUNALP(IB1)*RUNALP(IB2)
     &       -HPOLAR*(
     &         +G1(IB1,IB2)*DFLOAT(LLEPT)*R3*2D0*RUNALP(IB1)*RUNALP(IB2)
     &         +CHNULL*G4(IB1,IB2)*R2*RUNALP(IB1)*RUNALP(IB2)
     &         +CHNULL*G3(IB1,IB2)*R1*RUNALP(IB1)*RUNALP(IB2))
 11     CONTINUE
 10     CONTINUE
        SUMME = CQP(12)
      ELSEIF (IPDFOP.EQ.1) THEN
C..................................PARTON DENSITIES
        CALL HSPVER(XS,-TS)
        IF (HPOLAR.NE.0D0) THEN
          CALL HSDPVR(XS,-TS)
        ENDIF
        CALL HSIROT
        AUP=0D0
        AUM=0D0
        ADP=0D0
        ADM=0D0
        GUP=0D0
        GUM=0D0
        GDP=0D0
        GDM=0D0
        FUPI = 0D0
        FDOI = 0D0
        FUPIB =0D0
        FDOIB =0D0
        GUPI = 0D0
        GDOI = 0D0
        GUPIB =0D0
        GDOIB =0D0
        IF(LPAR(14).NE.0)THEN
          SU3 = (S*S+U*U)/2D0
          SU4 = (S-U)*(S+U)/2D0
          SU5 = (SS*SS+US*US)/2D0
          SU6 = (SS-US)*(SS+US)/2D0
          PPHA1= -T/DKQ/DKQS
          PPHA2= -4D0*MQF2/DKQS/DKQS
          PPHA3= -4D0*MQI2/DKQ/DKQ
          RALEP =
     &      GX/2D0*(
     &        +(T+4D0*MEI2)*(1D0/A2/TS/TS-1D0/A1/TS/TS)
     &        +2D0/TS/TS - 1D0/A1/TS + 1D0/A2/TS + 2D0/A1/A2
     &        -8D0*MEI2*MEI2/A1/A2/TS/TS
     &        +4D0*MEI2*MEI2/TS/TS*(1D0/A1/A1+1D0/A2/A2)
     &        +2D0*MEI2/TS*(1D0/A1/A1+1D0/A2/A2)  )
     &     +XS*(
     &        FAC1/TS/A1/A2
     &       +FAC2/A1/TS/TS
     &       +FAC3/A2/TS/TS
     &       +FAC4/A1/A2/TS/TS
     &       +FAC5/A1/A1/TS/TS
     &       +FAC6/A2/A2/TS/TS
     &       -2D0*MEI2*MPRO2/TS*(1D0/A1/A1+1D0/A2/A2)
     &       -2D0*MPRO2*(1D0/TS/TS+1D0/A1/A2)
     &       +MPRO2/TS*(1D0/A1-1D0/A2) )
          RBLEP =
     &      -((GS-GU)/A1/A2
     &       -2D0*MEI2*(GS-GU)/A1/A2/TS
     &       +(GX/2D0-GS)/A1/TS
     &       +(GX/2D0-GU)/A2/TS
     &       +GX*T/2D0/A1/TS/TS
     &       +GX*T/2D0/A2/TS/TS
     &       +(GX-2D0*GU)*MEI2/A1/A1/TS
     &       +(2D0*GS-GX)*MEI2/A2/A2/TS)
          PLEP=0D0
          PHAD=0D0
          PFRCH=0D0
          PFRCH0=0D0
          PFRCHC=0D0
          GLEP=0D0
          GHAD=0D0
          GFRCH=0D0
          GFRCH0=0D0
          GFRCHC=0D0
          DO 511 IB1=1,2
            DO 512 IB2=1,2
              PLEP=(AFIJ(2,IB1,IB2)*RALEP
     &               +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*RBLEP)
     &            *DSBOS(IB1)*DSBOS(IB2) + PLEP
              PHAD =((AFIJ(2,IB1,IB2)*SU1
     &               -BFIJ(2,IB1,IB2)*SU2*LLEPT)*PPHA1
     &              +(AFIJ(2,IB1,IB2)*SU3
     &               -BFIJ(2,IB1,IB2)*SU4*LLEPT)*PPHA2
     &              +(AFIJ(2,IB1,IB2)*SU5
     &               -BFIJ(2,IB1,IB2)*SU6*LLEPT)*PPHA3)
     &             *DBOS(IB1)*DBOS(IB2)/T/T*4D0/9D0/XS/8D0 + PHAD
              PFRCH=(AFIJ(2,IB1,IB2)*(R1/8D0+R2/4D0*XS)
     &              +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3/4D0)
     &             *DSBOS(IB1)*DSBOS(IB2) + PFRCH
              GLEP=(AGIJ(2,IB1,IB2)*RALEP
     &              +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*RBLEP)
     &            *DSBOS(IB1)*DSBOS(IB2) + GLEP
              GHAD =((AGIJ(2,IB1,IB2)*SU1
     &               -BGIJ(2,IB1,IB2)*SU2*LLEPT)*PPHA1
     &              +(AGIJ(2,IB1,IB2)*SU3
     &               -BGIJ(2,IB1,IB2)*SU4*LLEPT)*PPHA2
     &              +(AGIJ(2,IB1,IB2)*SU5
     &               -BGIJ(2,IB1,IB2)*SU6*LLEPT)*PPHA3)
     &             *DBOS(IB1)*DBOS(IB2)/T/T*4D0/9D0/XS/8D0 + GHAD
              GFRCH=(AGIJ(2,IB1,IB2)*(R1/8D0+R2/4D0*XS)
     &              +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*R3/4D0)
     &             *DSBOS(IB1)*DSBOS(IB2) + GFRCH
 512        CONTINUE
 511      CONTINUE
          SUTOT=PHAD+PLEP
          FRCH=PFRCH/SUTOT
          SGTOT=GHAD+GLEP
          GRCH=GFRCH/SGTOT
        ELSE
          FRCH=0D0
          GRCH=0D0
        ENDIF
        DO 501 IB1=1,2
         DO 502 IB2=1,2
          IF(LPAR(12).NE.0)THEN
             AUP=(CHNULL*AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUP
             AUM=(CHNULL*AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &            -DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUM
             ADP=(CHNULL*AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BFIJ(3,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADP
             ADM=(CHNULL*AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &            -DFLOAT(LLEPT)*BFIJ(3,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADM
             GUP=(CHNULL*AGIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + GUP
             GUM=(-CHNULL*AGIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + GUM
             GDP=(CHNULL*AGIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BGIJ(3,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + GDP
             GDM=(-CHNULL*AGIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &            +DFLOAT(LLEPT)*BGIJ(3,IB1,IB2)*R3)
     &          /4D0*DSBOS(IB1)*DSBOS(IB2)  + GDM
          ENDIF
         IF(LPAR(14).NE.0)THEN
          PROPI2 = DBOS(IB1)*DSBOS(IB2)/T/TS
          AB12UP = AFIJ(2,IB1,IB2)*SU1*LLEPT-CHNULL*BFIJ(2,IB1,IB2)*SU2
          AB12UM = AFIJ(2,IB1,IB2)*SU1*LLEPT+CHNULL*BFIJ(2,IB1,IB2)*SU2
          AB12DP = AFIJ(3,IB1,IB2)*SU1*LLEPT-CHNULL*BFIJ(3,IB1,IB2)*SU2
          AB12DM = AFIJ(3,IB1,IB2)*SU1*LLEPT+CHNULL*BFIJ(3,IB1,IB2)*SU2
          GB12UP = AGIJ(2,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(2,IB1,IB2)*SU2
          GB12UM =-AGIJ(2,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(2,IB1,IB2)*SU2
          GB12DP = AGIJ(3,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(3,IB1,IB2)*SU2
          GB12DM =-AGIJ(3,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(3,IB1,IB2)*SU2
          FUPI = AB12UP*PROPI1 * PROPI2     + FUPI
          FDOI = AB12DP*PROPI1 * PROPI2     + FDOI
          FUPIB =AB12UM*PROPI1 * PROPI2     + FUPIB
          FDOIB =AB12DM*PROPI1 * PROPI2     + FDOIB
          GUPI = GB12UP*PROPI1 * PROPI2     + GUPI
          GDOI = GB12DP*PROPI1 * PROPI2     + GDOI
          GUPIB =GB12UM*PROPI1 * PROPI2     + GUPIB
          GDOIB =GB12DM*PROPI1 * PROPI2     + GDOIB
         ENDIF
 502     CONTINUE
 501    CONTINUE
        FUPI=FUPI*FRCH
        FDOI=FDOI*FRCH
        FUPIB=FUPIB*FRCH
        FDOIB=FDOIB*FRCH
        GUPI=GUPI*GRCH
        GDOI=GDOI*GRCH
        GUPIB=GUPIB*GRCH
        GDOIB=GDOIB*GRCH

        CQP(1) =  QU *(AUP + FUPI  *  2D0/3D0 )
     &  - HPOLAR*DQU *(GUP + GUPI  *  2D0/3D0 )
        CQP(2) =  QBU*(AUM + FUPIB *(-2D0/3D0)) 
     &  - HPOLAR*DQBU*(GUM + GUPIB *(-2D0/3D0)) + CQP(1)
        CQP(3) =  QD *(ADP + FDOI  *(-1D0/3D0))
     &  - HPOLAR*DQD *(GDP + GDOI  *(-1D0/3D0)) + CQP(2)
        CQP(4) =  QBD*(ADM + FDOIB *  1D0/3D0 )
     &  - HPOLAR*DQBD*(GDM + GDOIB *  1D0/3D0 ) + CQP(3)
        CQP(5) =  QS *(ADP + FDOI  *(-1D0/3D0))
     &  - HPOLAR*DQS *(GDP + GDOI  *(-1D0/3D0)) + CQP(4)
        CQP(6) =  QBS*(ADM + FDOIB *  1D0/3D0 )
     &  - HPOLAR*DQBS*(GDM + GDOIB *  1D0/3D0 ) + CQP(5)
        CQP(7) =  QC *(AUP + FUPI  *  2D0/3D0 )
     &  - HPOLAR*DQC *(GUP + GUPI  *  2D0/3D0 ) + CQP(6)
        CQP(8) =  QBC*(AUM + FUPIB *(-2D0/3D0))
     &  - HPOLAR*DQBC*(GUM + GUPIB *(-2D0/3D0)) + CQP(7)
        CQP(9) =  QB *(ADP + FDOI  *(-1D0/3D0))
     &  - HPOLAR*DQB *(GDP + GDOI  *(-1D0/3D0)) + CQP(8)
        CQP(10) = QBB*(ADM + FDOIB *  1D0/3D0 )
     &  - HPOLAR*DQBB*(GDM + GDOIB *  1D0/3D0 ) + CQP(9)
        CQP(11) = QT *(AUP + FUPI  *  2D0/3D0 )
     &  - HPOLAR*DQT *(GUP + GUPI  *  2D0/3D0 ) + CQP(10)
        CQP(12) = QBT*(AUM + FUPIB *(-2D0/3D0))
     &  - HPOLAR*DQBT*(GUM + GUPIB *(-2D0/3D0)) + CQP(11)
        SUMME = CQP(12)
      ELSEIF (IPDFOP.GE.2) THEN
C..................................PARTON DISTRIBUTION FUNCTIONS
C                                  INCLUDING F_L
        CALL HSSTRF(XS,-TS)
        AUP=0D0
        AUM=0D0
        ADP=0D0
        ADM=0D0
        DO 521 IB1=1,2
         DO 522 IB2=1,2
            AUP=(CHNULL*AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &           +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3)
     &         /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUP
            AUM=(CHNULL*AFIJ(2,IB1,IB2)*(R1/2D0+R2*XS)
     &           -DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*R3)
     &         /4D0*DSBOS(IB1)*DSBOS(IB2)  + AUM
            ADP=(CHNULL*AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &           +DFLOAT(LLEPT)*BFIJ(3,IB1,IB2)*R3)
     &         /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADP
            ADM=(CHNULL*AFIJ(3,IB1,IB2)*(R1/2D0+R2*XS)
     &           -DFLOAT(LLEPT)*BFIJ(3,IB1,IB2)*R3)
     &         /4D0*DSBOS(IB1)*DSBOS(IB2)  + ADM
 522    CONTINUE
 521  CONTINUE
        CQP(1) =  QU *AUP
        CQP(2) =  QBU*AUM + CQP(1)
        CQP(3) =  QD *ADP + CQP(2)
        CQP(4) =  QBD*ADM + CQP(3)
        CQP(5) =  QS *ADP + CQP(4)
        CQP(6) =  QBS*ADM + CQP(5)
        CQP(7) =  QC *AUP + CQP(6)
        CQP(8) =  QBC*AUM + CQP(7)
        CQP(9) =  QB *ADP + CQP(8)
        CQP(10) = QBB*ADM + CQP(9)
        CQP(11) = QT *AUP + CQP(10)
        CQP(12) = QBT*AUM + CQP(11)
        SUMME=0D0
        DO 531 IB1=1,2
         DO 532 IB2=1,2
            SUMME=SUMME
     &      +F1(IB1,IB2)*CHNULL*R1*RUNALP(IB1)*RUNALP(IB2)
     &      +F2(IB1,IB2)*CHNULL*R2*RUNALP(IB1)*RUNALP(IB2)
     &      +F3(IB1,IB2)*DFLOAT(LLEPT)*R3*RUNALP(IB1)*RUNALP(IB2)
 532     CONTINUE
 531    CONTINUE
        RNORM=SUMME/CQP(12)
        DO 54 IF=1,12
          CQP(IF)=CQP(IF)*RNORM
 54     CONTINUE
      ENDIF

      HSK1TS = SUMME*Y*2D0*SX1NRM/SQGRAM
     *      * (UVMAX-UVMIN)*(XS-XX)*(A1MAX-A1MIN)*(RMAX-RMIN)*(-TS)
     *      * (XMAX-XMIN) * DG2/(XX*GS*GACT**2)
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR COS(THETA(K,P)) FOR KPS-KP PARAMETRIZATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSLZTS(ZMIN,ZMAX)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSPSPC/ IPHSPC
C
      IPHSPC=0
      CALP = (PS-PEL*COSE)/PQ
      IF (DABS(CALP).GT.1D0) THEN
        IPHSPC=1
        RETURN
      ENDIF
      ATAU = DELTA/OMEGA*( PQH*EH + PH*EQH)
      ALAM = PQH*EEL + PH*EQ
      AMU = PQH*PS*DSQRT((1D0-CALP)*(1D0+CALP))
      AXI = -PH*PQ+PQH*(PQ-PS*CALP)
C-----
      IF (AMU.GT.0D0) THEN
C
      F2U = (ATAU + AXI - ALAM)/AMU
      F2O = (ATAU - AXI - ALAM)/AMU
C
      IF ( (F2U.LT.0D0).AND.(F2O.LT.0D0) ) THEN
         ZMIN = -1D0
         ZMAX =  1D0
         RETURN
      ENDIF
C
      EZ = F2U*F2O
      IF ( EZ.LE.0D0 ) THEN
         AZ = AXI*AXI + AMU*AMU
         BZ = 2D0*AXI*(ALAM - ATAU)
         CZ = (ALAM-ATAU)*(ALAM-ATAU) - AMU*AMU
         DZ = BZ*BZ - 4D0*AZ*CZ
         IF(DZ.LT.0D0)DZ=0D0
         IF ((AXI/AMU).LT.0D0) THEN
            ZMIN = -1D0
            IF (BZ.GT.0D0) THEN
               ZMAX = 2D0*CZ/(-BZ - DSQRT(DZ))
            ELSE
               ZMAX = (- BZ + DSQRT(DZ))/2D0/AZ
            ENDIF
            RETURN
         ELSE
            ZMAX = 1D0
            IF (BZ.GT.0D0) THEN
               ZMIN = (- BZ - DSQRT(DZ))/2D0/AZ
            ELSE
               ZMIN = 2D0*CZ/(-BZ + DSQRT(DZ))
            ENDIF
            RETURN
         ENDIF
      ENDIF
C
      IF ( (F2U.GT.0D0).AND.(F2O.GT.0D0) ) THEN
         AZ = AXI*AXI + AMU*AMU
         BZ = 2D0*AXI*(ALAM - ATAU)
         CZ = (ALAM-ATAU)*(ALAM-ATAU) - AMU*AMU
         DZ = BZ*BZ - 4D0*AZ*CZ
         IF (DZ.LT.0D0) THEN
            ZMIN = 0D0
            ZMAX = 0D0
            RETURN
         ELSE
            IF (BZ.GT.0D0) THEN
               ZMIN = (- BZ - DSQRT(DZ))/2D0/AZ
               ZMAX = 2D0*CZ/(-BZ - DSQRT(DZ))
               RETURN
            ELSE
               ZMAX = (- BZ + DSQRT(DZ))/2D0/AZ
               ZMIN = 2D0*CZ/(-BZ + DSQRT(DZ))
               RETURN
            ENDIF
         ENDIF
      ENDIF
      ENDIF
C
C-----
      IF (AMU.EQ.0D0) THEN
        IF(AXI.GT.0D0) THEN
           ZMAX = 1D0
           ZMIN = DMAX1(-1D0,(ATAU-ALAM)/AXI)
           IF (ZMIN.GT.1D0) ZMIN = 1D0
           RETURN
         ELSEIF(AXI.LT.0D0) THEN
           ZMIN =-1D0
           ZMAX = DMIN1(1D0,(ATAU-ALAM)/AXI)
           IF (ZMAX.LT.-1D0) ZMAX =-1D0
           RETURN
         ELSE
           IF (ALAM.GE.ATAU) THEN
              ZMIN = -1D0
              ZMAX = 1D0
              RETURN
            ELSE
              IPHSPC=1
            ENDIF
         ENDIF
      ENDIF
      IF (AMU.LT.0D0)THEN
          IPHSPC=1
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR A1 AS FUNCTION OF TS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSL1TS(TTS,XX,Y,XS,A1MIN,A1MAX,A1M,A1P,CFKP)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSPSPC/ IPHSPC
C
      CALL HSFIVM(XX,Y,XS)
C----- A1 LIMIT FROM OMH > DELTA: A1 > A1GR
      A1GR = ( 2D0*DELTA*(PH*EPRO+PPRO*EH)
     &         - PH/XS*((XS-XX)*Y*GS + TTS - TP) )/PPRO
C
      MXS = XS*MPRO
      MXS2 = MXS*MXS
      A=((S+U-2D0*MEF2-2D0*MXS2)
     &  *(S+U-2D0*MEF2-2D0*MXS2)-4D0*T*MXS2)/16D0
      B=((U+T+S-2D0*MEF2-2D0*MXS2)
     &    *(S*TTS-U*T+(T-TTS)*(MEF2+MXS2))
     &   +T*(T-TTS)*(U-MEF2+MXS2))/8D0
      C=((TTS*(S+T)+U*T-(T+TTS)*(MEF2+MXS2))
     &  *(TTS*(S+T)+U*T-(T+TTS)*(MEF2+MXS2))
     &  -4D0*MEF2*(TTS*(U+S-2D0*MEF2-2D0*MXS2)
     &               *(U+T+S+TTS-2D0*MEF2-2D0*MXS2)
     &             +MXS2*(TTS-T)*(TTS-T)+TTS*TTS*T))/16D0
      DISK=1D0/16D0*(MEF2*(S+U-2D0*MEF2-2D0*MXS2)
     &                   *(S+U-2D0*MEF2-2D0*MXS2)
     &              +T*(MEF2+MXS2)*(S+U-MEF2-MXS2)
     &              -U*T*S + T*MXS2*(T-4D0*MEF2))
     &       *(TTS*(S+U+TTS-2D0*MEF2-2D0*MXS2)*(S+U+T-2D0*MEF2-2D0*MXS2)
     &        +MXS2*(TTS-T)*(TTS-T))
      CFKP = A
      IF (DISK.LE.0D0) THEN
        DISK = 0D0
      ENDIF
C
      IF (B.GE.0D0) THEN
        AH1 = (-B-DSQRT(DISK))/2D0/A
        AH2 = C/A/AH1
      ELSE
        AH2 = (-B+DSQRT(DISK))/2D0/A
        AH1 = C/A/AH2
      ENDIF
      A1M = DMIN1(AH1,AH2)
      A1P = DMAX1(AH1,AH2)
      A1MIN = DMAX1(A1M,A1GR)
      A1MAX = A1P
      IF (A1MIN.GE.A1MAX) THEN
        IPHSPC=1
      ELSE
        IPHSPC=0
      ENDIF
C
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND OF KQ TERM (QUARKONIC RADIATION)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSK1K3(X)
C
C  X(1) -->  XX
C  X(2) -->  G=-1/Q**2   -->  Q**2=-1/G-->  Y
C  X(3) -->  LOG(XS-XX)
C  X(4) -->  LOG(2*K.Q)
C  X(5) -->  2*K.P
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSKPXY/ XX,Y
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSXSLM/ XSMIN,XSCUT
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSASY0/ PONULL,CHNULL,ICPASY
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSFIJK/ F1(2,2),F2(2,2),F3(2,2)
      COMMON /HSPSPC/ IPHSPC
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSWGTC/ IWEIGS
      COMMON /HSISGM/ TCUTQ,TCUTQS
      DIMENSION X(5),DBOS(2),DSBOS(2)


C---X-VALUE
      XX=XMIN + (XMAX-XMIN)*X(1)
C---Y-VALUE
      GS=SP-MEI2-MPRO2
      YMAXX=XX*(1D0-4D0*MEI2*MPRO2/GS/GS)/(XX*(1D0+XX*MPRO2/GS)+MEI2/GS)
      IF(ICUT.LT.3) THEN
C---CUT IN EXTERNALLY DEFINED Q**2(MIN)
        GMAX=-1D0/(XX*YMAXX*GS)
        ELSEIF(ICUT.EQ.3) THEN
C---CUT IN Y
        QQ2MIN=XX*YMIN*GS
C---CUT ON ELECTRON SCATTERING ANGLE (MASSES NEGLECTED)
        YMINTH=1D0/(1D0+XX*CTHCON)
        QT2MIN=XX*YMINTH*GS
C---CUT ON ELECTRON TRANSVERSE MOMENTUM (MASSES NEGLECTED)
        QP2MIN=XX*GS/2D0*(1D0-DSQRT(1D0-PTXM0/XX))
        QP2MAX=XX*GS/2D0*(1D0+DSQRT(1D0-PTXM0/XX))
        YP2MAX=QP2MAX/GS/XX
        GMIN=-1D0/MAX(Q2MIN,QQ2MIN,QT2MIN,QP2MIN)
        GMAX=-1D0/(XX*MIN(YMAX,YMAXX,YP2MAX)*GS)
        GMAX=MIN(GMAX,-1D0/Q2MAX)
        ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSK1K3'
        STOP
      ENDIF
      DG2=DMAX1(GMAX-GMIN,0D0)
      GACT=GMIN+X(2)*DG2
      Y=-1D0/(XX*GACT*GS)
C
C---EXTERNAL WEIGHT
      IACPT=1
      IF (IWEIGS.GT.0) THEN
        CALL HSWGTX(XX,Y,IACPT)
        IF (IACPT.EQ.0) THEN
          HSK1K3=0D0
          RETURN
        ENDIF
      ENDIF
C
      CALL HSFIVC(XX,Y)
      CALL HSDELO(XX,Y)
      XSMIN = HSXSMN(XX,Y)
C     XSCUT = HSXSCT(XX,Y)

      XSMAX = 1D0
      XSMINI = XSMIN
      IF((XSMAX-XSMINI).LE.1D-14) THEN
        HSK1K3=0D0
        RETURN
      ENDIF

C---SUBSTITUTION UV=LN(XS-X)
      UVMIN = DLOG(XSMINI-XX)
      UVMAX = DLOG(XSMAX-XX)
      UV = UVMIN + (UVMAX-UVMIN)*X(3)
      XS = XX+DEXP(UV)
      CALL HSFCMS(XX,Y,XS)
      IF(IPHSPC.EQ.1) THEN
        HSK1K3=0D0
        RETURN
      ENDIF
      CALL HSFLAB(XS)
      CALL HSLZK3(ZMIN,ZMAX)
      IF(IPHSPC.EQ.1) THEN
        HSK1K3=0D0
        RETURN
      ENDIF

C---SUBSTITUTION V = LN(A3)
      A3MAX = 2D0*OMEGA*(EQ-PQ*ZMIN)
      IF(ZMAX.GE.0.99D0) THEN
        A3MIN = MQI2*OMEGA/EQ
      ELSE
        A3MIN = 2D0*OMEGA*(EQ-PQ*ZMAX)
      ENDIF

C---CUT FROM EXCLUSION OF INITIAL STATE COLLINEAR PHOTONS
      IF(TCUTQ.NE.0D0)A3MIN=MAX(A3MIN,2D0*DELTA*(EQH-PQH*DCOS(TCUTQ)))
C---CUT IN HADRONIC MASS
      IF (ICUT.GT.1) THEN
        W2CUT=WMIN*WMIN
        AWCUT=XS*Y*GS-XS*(W2CUT-MPRO2)/(1D0-XS)
        A3MAX=MIN(A3MAX,AWCUT)
      ENDIF

      IF (A3MAX.LE.0D0.OR.A3MIN.LE.0D0) THEN
        HSK1K3=0D0
        RETURN
      ENDIF
      VMAX = DLOG(A3MAX)
      VMIN = DLOG(A3MIN)
      V= VMIN + (VMAX-VMIN)*X(4)
      A3 = DEXP(V)

C---NO SUBSTITUTION FOR A1
      CALL HSL1K3(A3,XX,Y,XS,A1MIN,A1MAX,A1M,A1P,CFKP)
      IF(IPHSPC.EQ.1) THEN
        HSK1K3=0D0
        RETURN
      ENDIF
      A1=A1MIN + (A1MAX-A1MIN)*X(5)
      SQGRM2=-CFKP*(A1-A1M)*(A1-A1P)
      SQGRAM = DSQRT(DABS(SQGRM2))
      CALL HSFIV3(XX,Y,XS,A1,A3)
      CALL HSPVER(XS,-T)
      IF (HPOLAR.NE.0D0) THEN
        CALL HSDPVR(XS,-T)
      ENDIF
      CALL HSIROT

      OMH = (PH*DKQ+PQH*DKP)/(PH*EQH+PQH*EH)
C---SUBTRACTION OF INITIAL STATE MASS SINGULARITY
      DCKQH = (DKQ/OMH-EQH+PQH)/PQH
      CKQH = 1D0-DCKQH
      IF(DABS(CKQH).GT.1D0)CKQH=DSIGN(1D0,CKQH)
      THKQ = DACOS(CKQH)
      IF(THKQ.LT.TCUTQ)THEN
        HSK1K3=0D0
        RETURN
      ENDIF
C---SUBTRACTION OF FINAL STATE MASS SINGULARITY
      EQSH = EH+EQH-OMH-ESH
      PQSH = DSQRT(DABS((EQSH+MQF)*(EQSH-MQF)))
      DCKQSH = (DKQS/OMH-EQSH+PQSH)/PQSH
      CKQSH = 1D0-DCKQSH
      IF(DABS(CKQSH).GT.1D0)CKQSH=DSIGN(1D0,CKQSH)
      THKQS = DACOS(CKQSH)
      IF(THKQS.LT.TCUTQS)THEN
        HSK1K3=0D0
        RETURN
      ENDIF

      A2=2D0*DKPS
      DBOS(1)=1D0/T
      DBOS(2)=1D0/(T-MZ2)
      DSBOS(1)=1D0/TS
      DSBOS(2)=1D0/(TS-MZ2)
      SU1 = S*S+SS*SS+U*U+US*US
      SU2 = (S-U)*(S+U) + (SS-US)*(SS+US)
      SU3 = (S*S+U*U)/2D0
      SU4 = (S-U)*(S+U)/2D0
      SU5 = (SS*SS+US*US)/2D0
      SU6 = (SS-US)*(SS+US)/2D0
      PROP1= -T/DKQ/DKQS
      PROP2= -4D0*MQF2/DKQS/DKQS
      PROP3= -4D0*MQI2/DKQ/DKQ
      IF(LPAR(14).EQ.0)THEN
        PROPI1=0D0
      ELSE
        PROPI1= 
     &       - (SS/DKPS/DKQS + S/DKP/DKQ + U/DKPS/DKQ + US/DKP/DKQS)
      ENDIF
      DO 1 IFL=1,12
        CQP(IFL) = 0D0
 1    CONTINUE
      FUP = 0D0
      FDO = 0D0
      FUPB =0D0
      FDOB =0D0
      GUP = 0D0
      GDO = 0D0
      GUPB =0D0
      GDOB =0D0
      FUPI = 0D0
      FDOI = 0D0
      FUPIB =0D0
      FDOIB =0D0
      GUPI = 0D0
      GDOI = 0D0
      GUPIB =0D0
      GDOIB =0D0
      IF (LPAR(14).NE.0) THEN
        RALEP = 8D0*XS*(
     &     GX/2D0*(
     &      +(T+4D0*MEI2)*(1D0/A2/TS/TS-1D0/A1/TS/TS)
     &      +2D0/TS/TS - 1D0/A1/TS + 1D0/A2/TS + 2D0/A1/A2
     &      -8D0*MEI2*MEI2/A1/A2/TS/TS
     &      +4D0*MEI2*MEI2/TS/TS*(1D0/A1/A1+1D0/A2/A2)
     &      +2D0*MEI2/TS*(1D0/A1/A1+1D0/A2/A2)  )
     &     +XS*(
     &      -(GS*GS + GU*GU - GX*(GS+GU) -4D0*MEI2*MPRO2)/TS/A1/A2
     &      +(-GS*GX + MPRO2*T)/A1/TS/TS
     &      +(GU*GX - MPRO2*T)/A2/TS/TS
     &      +(-2D0*GS*GU + GX*(GS+GU-GX))*2D0*MEI2/A1/A2/TS/TS
     &      +GU*(GX-GU)*2D0*MEI2/A1/A1/TS/TS
     &      +GS*(GX-GS)*2D0*MEI2/A2/A2/TS/TS
     &      -2D0*MEI2*MPRO2/TS*(1D0/A1/A1+1D0/A2/A2)
     &      -2D0*MPRO2*(1D0/TS/TS+1D0/A1/A2)
     &      +MPRO2/TS*(1D0/A1-1D0/A2) ))
        RBLEP = -8D0*XS*(
     &     (GS-GU)/A1/A2
     &     -2D0*MEI2*(GS-GU)/A1/A2/TS
     &     +(GX/2D0-GS)/A1/TS
     &     +(GX/2D0-GU)/A2/TS
     &     +GX*T/2D0/A1/TS/TS
     &     +GX*T/2D0/A2/TS/TS
     &     +(GX-2D0*GU)*MEI2/A1/A1/TS
     &     +(2D0*GS-GX)*MEI2/A2/A2/TS)
        PLEP=0D0
        PHAD=0D0
        GLEP=0D0
        GHAD=0D0
        DO 511 IB1=1,2
         DO 512 IB2=1,2
         PLEP=(AFIJ(2,IB1,IB2)*RALEP
     &          +DFLOAT(LLEPT)*BFIJ(2,IB1,IB2)*RBLEP)
     &                *DSBOS(IB1)*DSBOS(IB2)*TS*TS  + PLEP
         PHAD =((AFIJ(2,IB1,IB2)*SU1 - BFIJ(2,IB1,IB2)*SU2*LLEPT)*PROP1
     &         +(AFIJ(2,IB1,IB2)*SU3 - BFIJ(2,IB1,IB2)*SU4*LLEPT)*PROP2
     &         +(AFIJ(2,IB1,IB2)*SU5 - BFIJ(2,IB1,IB2)*SU6*LLEPT)*PROP3)
     &         *DBOS(IB1)*DBOS(IB2)*4D0/9D0 + PHAD
         GLEP=(AGIJ(2,IB1,IB2)*RALEP
     &         +DFLOAT(LLEPT)*BGIJ(2,IB1,IB2)*RBLEP)
     &                *DSBOS(IB1)*DSBOS(IB2)*TS*TS  + GLEP
         GHAD =((AGIJ(2,IB1,IB2)*SU1 - BGIJ(2,IB1,IB2)*SU2*LLEPT)*PROP1
     &         +(AGIJ(2,IB1,IB2)*SU3 - BGIJ(2,IB1,IB2)*SU4*LLEPT)*PROP2
     &         +(AGIJ(2,IB1,IB2)*SU5 - BGIJ(2,IB1,IB2)*SU6*LLEPT)*PROP3)
     &         *DBOS(IB1)*DBOS(IB2)*4D0/9D0 + GHAD
 512     CONTINUE
 511    CONTINUE
        SUTOT=PHAD+PLEP
        FRCH=PHAD/SUTOT
        SGTOT=GHAD+GLEP
        GRCH=GHAD/SGTOT
      ELSE
        FRCH=0D0
        GRCH=0D0
      ENDIF
      DO 10 IB1 = 1,2
       DO 20 IB2 = 1,2
        PROP = DBOS(IB1)*DBOS(IB2)
        AB12UP=CHNULL*AFIJ(2,IB1,IB2)*SU1 - BFIJ(2,IB1,IB2)*SU2*LLEPT
        AB12UM=CHNULL*AFIJ(2,IB1,IB2)*SU1 + BFIJ(2,IB1,IB2)*SU2*LLEPT
        AB12DP=CHNULL*AFIJ(3,IB1,IB2)*SU1 - BFIJ(3,IB1,IB2)*SU2*LLEPT
        AB12DM=CHNULL*AFIJ(3,IB1,IB2)*SU1 + BFIJ(3,IB1,IB2)*SU2*LLEPT
        AB34UP=CHNULL*AFIJ(2,IB1,IB2)*SU3 - BFIJ(2,IB1,IB2)*SU4*LLEPT
        AB34UM=CHNULL*AFIJ(2,IB1,IB2)*SU3 + BFIJ(2,IB1,IB2)*SU4*LLEPT
        AB34DP=CHNULL*AFIJ(3,IB1,IB2)*SU3 - BFIJ(3,IB1,IB2)*SU4*LLEPT
        AB34DM=CHNULL*AFIJ(3,IB1,IB2)*SU3 + BFIJ(3,IB1,IB2)*SU4*LLEPT
        AB56UP=CHNULL*AFIJ(2,IB1,IB2)*SU5 - BFIJ(2,IB1,IB2)*SU6*LLEPT
        AB56UM=CHNULL*AFIJ(2,IB1,IB2)*SU5 + BFIJ(2,IB1,IB2)*SU6*LLEPT
        AB56DP=CHNULL*AFIJ(3,IB1,IB2)*SU5 - BFIJ(3,IB1,IB2)*SU6*LLEPT
        AB56DM=CHNULL*AFIJ(3,IB1,IB2)*SU5 + BFIJ(3,IB1,IB2)*SU6*LLEPT
        GB12UP=CHNULL*AGIJ(2,IB1,IB2)*SU1 - BGIJ(2,IB1,IB2)*SU2*LLEPT
        GB12UM=-CHNULL*AGIJ(2,IB1,IB2)*SU1 - BGIJ(2,IB1,IB2)*SU2*LLEPT
        GB12DP=CHNULL*AGIJ(3,IB1,IB2)*SU1 - BGIJ(3,IB1,IB2)*SU2*LLEPT
        GB12DM=-CHNULL*AGIJ(3,IB1,IB2)*SU1 - BGIJ(3,IB1,IB2)*SU2*LLEPT
        GB34UP=CHNULL*AGIJ(2,IB1,IB2)*SU3 - BGIJ(2,IB1,IB2)*SU4*LLEPT
        GB34UM=-CHNULL*AGIJ(2,IB1,IB2)*SU3 - BGIJ(2,IB1,IB2)*SU4*LLEPT
        GB34DP=CHNULL*AGIJ(3,IB1,IB2)*SU3 - BGIJ(3,IB1,IB2)*SU4*LLEPT
        GB34DM=-CHNULL*AGIJ(3,IB1,IB2)*SU3 - BGIJ(3,IB1,IB2)*SU4*LLEPT
        GB56UP=CHNULL*AGIJ(2,IB1,IB2)*SU5 - BGIJ(2,IB1,IB2)*SU6*LLEPT
        GB56UM=-CHNULL*AGIJ(2,IB1,IB2)*SU5 - BGIJ(2,IB1,IB2)*SU6*LLEPT
        GB56DP=CHNULL*AGIJ(3,IB1,IB2)*SU5 - BGIJ(3,IB1,IB2)*SU6*LLEPT
        GB56DM=-CHNULL*AGIJ(3,IB1,IB2)*SU5 - BGIJ(3,IB1,IB2)*SU6*LLEPT
       IF(LPAR(13).NE.0)THEN
        FUP = (AB12UP*PROP1+AB34UP*PROP2+AB56UP*PROP3) * PROP  + FUP
        FDO = (AB12DP*PROP1+AB34DP*PROP2+AB56DP*PROP3) * PROP  + FDO
        FUPB =(AB12UM*PROP1+AB34UM*PROP2+AB56UM*PROP3) * PROP  + FUPB
        FDOB =(AB12DM*PROP1+AB34DM*PROP2+AB56DM*PROP3) * PROP  + FDOB
        GUP = (GB12UP*PROP1+GB34UP*PROP2+GB56UP*PROP3) * PROP  + GUP
        GDO = (GB12DP*PROP1+GB34DP*PROP2+GB56DP*PROP3) * PROP  + GDO
        GUPB =(GB12UM*PROP1+GB34UM*PROP2+GB56UM*PROP3) * PROP  + GUPB
        GDOB =(GB12DM*PROP1+GB34DM*PROP2+GB56DM*PROP3) * PROP  + GDOB
       ENDIF
        IF(LPAR(14).NE.0)THEN
          PROPI2 = DBOS(IB1)*DSBOS(IB2)
          AB12UP = AFIJ(2,IB1,IB2)*SU1*LLEPT-CHNULL*BFIJ(2,IB1,IB2)*SU2
          AB12UM = AFIJ(2,IB1,IB2)*SU1*LLEPT+CHNULL*BFIJ(2,IB1,IB2)*SU2
          AB12DP = AFIJ(3,IB1,IB2)*SU1*LLEPT-CHNULL*BFIJ(3,IB1,IB2)*SU2
          AB12DM = AFIJ(3,IB1,IB2)*SU1*LLEPT+CHNULL*BFIJ(3,IB1,IB2)*SU2
          GB12UP = AGIJ(2,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(2,IB1,IB2)*SU2
          GB12UM =-AGIJ(2,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(2,IB1,IB2)*SU2
          GB12DP = AGIJ(3,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(3,IB1,IB2)*SU2
          GB12DM =-AGIJ(3,IB1,IB2)*SU1*LLEPT-CHNULL*BGIJ(3,IB1,IB2)*SU2
          FUPI = AB12UP*PROPI1 * PROPI2     + FUPI
          FDOI = AB12DP*PROPI1 * PROPI2     + FDOI
          FUPIB =AB12UM*PROPI1 * PROPI2     + FUPIB
          FDOIB =AB12DM*PROPI1 * PROPI2     + FDOIB
          GUPI = GB12UP*PROPI1 * PROPI2     + GUPI
          GDOI = GB12DP*PROPI1 * PROPI2     + GDOI
          GUPIB =GB12UM*PROPI1 * PROPI2     + GUPIB
          GDOIB =GB12DM*PROPI1 * PROPI2     + GDOIB
        ENDIF
20     CONTINUE
10    CONTINUE
      FUPI=FUPI*FRCH
      FDOI=FDOI*FRCH
      FUPIB=FUPIB*FRCH
      FDOIB=FDOIB*FRCH
      GUPI=GUPI*GRCH
      GDOI=GDOI*GRCH
      GUPIB=GUPIB*GRCH
      GDOIB=GDOIB*GRCH

      CQP(1) =  QU *(FUP  *4D0/9D0 + FUPI  *  2D0/3D0 )
     &  - HPOLAR*DQU *(GUP *4D0/9D0 + GUPI *  2D0/3D0 )
      CQP(2) =  QBU*(FUPB *4D0/9D0 + FUPIB *(-2D0/3D0))
     &  - HPOLAR*DQBU*(GUPB *4D0/9D0 + GUPIB *(-2D0/3D0)) + CQP(1)
      CQP(3) =  QD *(FDO  *1D0/9D0 + FDOI  *(-1D0/3D0))
     &  - HPOLAR*DQD *(GDO  *1D0/9D0 + GDOI  *(-1D0/3D0)) + CQP(2)
      CQP(4) =  QBD*(FDOB *1D0/9D0 + FDOIB *  1D0/3D0 )
     &  - HPOLAR*DQBD*(GDOB *1D0/9D0 + GDOIB *  1D0/3D0 ) + CQP(3)
      CQP(5) =  QS *(FDO  *1D0/9D0 + FDOI  *(-1D0/3D0))
     &  - HPOLAR*DQS *(GDO  *1D0/9D0 + GDOI  *(-1D0/3D0)) + CQP(4)
      CQP(6) =  QBS*(FDOB *1D0/9D0 + FDOIB *  1D0/3D0 )
     &  - HPOLAR*DQBS*(GDOB *1D0/9D0 + GDOIB *  1D0/3D0 ) + CQP(5)
      CQP(7) =  QC *(FUP  *4D0/9D0 + FUPI  *  2D0/3D0 )
     &  - HPOLAR*DQC *(GUP  *4D0/9D0 + GUPI  *  2D0/3D0 ) + CQP(6)
      CQP(8) =  QBC*(FUPB *4D0/9D0 + FUPIB *(-2D0/3D0))
     &  - HPOLAR*DQBC*(GUPB *4D0/9D0 + GUPIB *(-2D0/3D0)) + CQP(7)
      CQP(9) =  QB *(FDO  *1D0/9D0 + FDOI  *(-1D0/3D0))
     &  - HPOLAR*DQB *(GDO  *1D0/9D0 + GDOI  *(-1D0/3D0)) + CQP(8)
      CQP(10) =  QBB*(FDOB *1D0/9D0 + FDOIB*  1D0/3D0 )
     &  - HPOLAR*DQBB*(GDOB *1D0/9D0 + GDOIB*  1D0/3D0 ) + CQP(9)
      CQP(11) =  QT *(FUP  *4D0/9D0 + FUPI *  2D0/3D0 )
     &  - HPOLAR*DQT *(GUP  *4D0/9D0 + GUPI *  2D0/3D0 ) + CQP(10)
      CQP(12) =  QBT*(FUPB *4D0/9D0 + FUPIB*(-2D0/3D0))
     &  - HPOLAR*DQBT*(GUPB *4D0/9D0 + GUPIB*(-2D0/3D0)) + CQP(11)
      SUMME = CQP(12)
      HSK1K3 = SUMME*Y*SX1NRM/4D0/SQGRAM/XS
     *      * (UVMAX-UVMIN)*(XS-XX)*(VMAX-VMIN)*A3*(A1MAX-A1MIN)
     *      * (XMAX-XMIN) * DG2/(XX*GS*GACT**2)
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR COS(THETA(K,Q)) FOR KQ-PARAMETRIZATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSLZK3(ZMIN,ZMAX)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSCMS1/ COSQ,SINQ
      COMMON /HSPSPC/ IPHSPC
C
      IPHSPC=0
      DEPS = DELTA/OMEGA*( PQH*EH + PH*EQH)
      ALAM = PQH*EEL + PH*EQ
      TAU  = PQH*( PS*COSQ - PQ ) + PH*PQ
      F2U = -TAU - ALAM + DEPS
      F2O =  TAU - ALAM + DEPS
C
      IF ( (F2U.LT.0D0).AND.(F2O.LT.0D0) ) THEN
         ZMIN = -1D0
         ZMAX =  1D0
         RETURN
      ENDIF
C
      EZ = F2U*F2O
      IF ( EZ.LE.0D0 ) THEN
         PPP4 = PQH*PQH*PS*PS*SINQ*SINQ
         DZ = 4D0*PPP4*( PPP4 - EZ )
         AZ = TAU*TAU + PPP4
         BZ = 2D0*TAU*(DEPS - ALAM)
         CZ = (DEPS-ALAM)*(DEPS-ALAM) - PPP4
         IF (TAU.GT.0D0) THEN
            ZMIN = -1D0
            IF (BZ.GT.0D0) THEN
               ZMAX = 2D0*CZ/(-BZ - DSQRT(DZ))
            ELSE
               ZMAX = (- BZ + DSQRT(DZ))/2D0/AZ
            ENDIF
            RETURN
         ELSE
            ZMAX = 1D0
            IF (BZ.GT.0D0) THEN
               ZMIN = (- BZ - DSQRT(DZ))/2D0/AZ
            ELSE
               ZMIN = 2D0*CZ/(-BZ + DSQRT(DZ))
            ENDIF
            RETURN
         ENDIF
      ENDIF
C
      IF ( (F2U.GT.0D0).AND.(F2O.GT.0D0) ) THEN
         PPP4 = PQH*PQH*PS*PS*SINQ*SINQ
         DZ = 4D0*PPP4*( PPP4 - EZ )
         IF (DZ.LT.0D0) THEN
            IPHSPC=1
            RETURN
         ELSE
         AZ = TAU*TAU + PPP4
         BZ =  2D0*TAU*(DEPS - ALAM)
         CZ = (DEPS - ALAM)*(DEPS - ALAM) - PPP4
            IF (BZ.GT.0D0) THEN
               ZMIN = (- BZ - DSQRT(DZ))/2D0/AZ
               ZMAX = 2D0*CZ/(-BZ - DSQRT(DZ))
               RETURN
            ELSE
               ZMAX = (- BZ + DSQRT(DZ))/2D0/AZ
               ZMIN = 2D0*CZ/(-BZ + DSQRT(DZ))
               RETURN
            ENDIF
         ENDIF
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR A1 AS FUNCTION OF A3
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSL1K3(A3,XX,Y,XS,A1MIN,A1MAX,A1M,A1P,CFKP)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSPSPC/ IPHSPC
C
      IPHSPC=0
      CALL HSFIVM(XX,Y,XS)

C---A1 LIMIT FROM OMH > DELTA: A1 > A1GR
      A1GR = (2D0*DELTA*(PH*EQH+PQH*EH) - PH*A3)/PQH

      A = ((S+U-2D0*MEF2-2D0*MQF2)*(S+U-2D0*MEF2-2D0*MQF2) - 4D0*T*MQF2)
     &    / 16D0
      B =( (U+S+T-A3-2D0*MEF2-2D0*MQF2)
     &   *(2D0*T*MQF2 - (S+U-2D0*MEF2-2D0*MQF2)
     &                 *(S-MEF2-MQF2) )
     &    + T*A3*(S-U) )/8D0
      C=( ((U+S+T-A3-2D0*MEF2-2D0*MQF2)*(S-MEF2-MQF2)-A3*T)
     &   *((U+S+T-A3-2D0*MEF2-2D0*MQF2)*(S-MEF2-MQF2)-A3*T)
     &  -4D0*MEF2*( MQF2*(S+U+T-2D0*MEF2-2D0*MQF2)
     &                  *(S+U+T-2D0*MEF2-2D0*MQF2)
     &       -A3*(S+U-A3-2D0*MEF2-MQF2)*(S+U+T-2D0*MEF2-MQF2)
     &       -A3*MQF2*(T-MQF2)) )/16D0
      DISK=( MEF2*(U+S-2D0*MEF2-2D0*MQF2)
     &           *(U+S-2D0*MEF2-2D0*MQF2)
     &         +T*(U+S-MEF2-MQF2)*(MEF2+MQF2)
     &         -U*T*S + T*MQF2*(T-4D0*MEF2) )
     &   * ( MQF2*(S+U+T-2D0*MEF2-2D0*MQF2)
     &           *(S+U+T-2D0*MEF2-2D0*MQF2)
     &        -A3*(S+U+T-2D0*MEF2-2D0*MQF2)*(S+U-A3-2D0*MEF2)
     &        + A3*A3*MQF2 )/16D0

      CFKP = A
      IF (DISK.LE.0D0) THEN
        DISK = 0D0
      ENDIF
C
      IF (B.GE.0D0) THEN
        AH1 = (-B-DSQRT(DISK))/2D0/A
        AH2 = C/A/AH1
      ELSE
        AH2 = (-B+DSQRT(DISK))/2D0/A
        AH1 = C/A/AH2
      ENDIF
      A1M = DMIN1(AH1,AH2)
      A1P = DMAX1(AH1,AH2)
      A1MIN = DMAX1(A1M,A1GR)
      A1MAX = A1P
      IF (A1MIN.GE.A1MAX) THEN
        IPHSPC=1
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   CROSS SECTION FOR DEEP INELASTIC LEPTON PROTON SCATTERING:
C   MODIFIED FOR USE OF STRUCTURE FUNCTIONS AS INPUT: 13.3.91: HS
CCCCCCC 22. 9. 86   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC H. SPIESBERGER CCC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSSGNC(X,Y,LL,POL,HPOL,LQ)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSASY0/ PONULL,CHNULL,ICPASY
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSWGTC/ IWEIGS
      COMMON /HSFIJK/ F1(2,2),F2(2,2),F3(2,2)
      COMMON /HSGIJK/ G1(2,2),G3(2,2),G4(2,2)
C
C---EXTERNAL WEIGHT
      IACPT=1
      IF (IWEIGS.GT.0) THEN
        CALL HSWGTX(X,Y,IACPT)
        IF (IACPT.EQ.0) THEN
          HSSGNC=0D0
          RETURN
        ENDIF
      ENDIF

      SHAT=SP-MEI2-MPRO2
      T = -X*Y*SHAT
      SPN=SXNORM*SP*X
      CALL HSDELO(X,Y)
      IF (IPDFOP.EQ.1) THEN
        CALL HSPVER(X,-T)
        IF (HPOL.NE.0D0) THEN
          CALL HSDPVR(X,-T)
        ENDIF
        CALL HSIROT
        ELSE
        IF (LPAR(2).EQ.0) THEN
          CALL HSSTRF(X,-T)
          ELSE
          CALL HSSTR1(X,-T)
        ENDIF
      ENDIF
      DO 10 I=1,12
        CQP(I)=0D0
 10   CONTINUE
      HSSGNC=0D0
      SPU=1D0+(1D0-Y)*(1D0-Y)
      SMU=1D0-(1D0-Y)*(1D0-Y)
C
C...BORN CROSS SECTIONS
      IF (LPAR(1).EQ.1) THEN
        IF (LQ.EQ.1) THEN
         CALL HSSAB0(X,Y,POL,HPOL, 1,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU*CHNULL
     *                -LL*(PONULL*B3F+POL*A3F)*SMU)*QU
     *              +HPOL*((PONULL*A3F+POL*B3F)*SPU*CHNULL
     *                -LL*(PONULL*B1F+POL*A1F)*SMU)*DQU)
        ENDIF
        IF (LQ.EQ.2) THEN
         CALL HSSAB0(X,Y,POL,HPOL, 2,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU*CHNULL
     *                -LL*(PONULL*B3F+POL*A3F)*SMU)*QD
     *              +HPOL*((PONULL*A3F+POL*B3F)*SPU*CHNULL
     *                -LL*(PONULL*B1F+POL*A1F)*SMU)*DQD)
        ENDIF
        IF (LQ.EQ.3) THEN
         CALL HSSAB0(X,Y,POL,HPOL, 2,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU*CHNULL
     *                -LL*(PONULL*B3F+POL*A3F)*SMU)*QS
     *              +HPOL*((PONULL*A3F+POL*B3F)*SPU*CHNULL
     *                -LL*(PONULL*B1F+POL*A1F)*SMU)*DQS)
        ENDIF
        IF (LQ.EQ.4) THEN
         CALL HSSAB0(X,Y,POL,HPOL, 1,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU*CHNULL
     *                -LL*(PONULL*B3F+POL*A3F)*SMU)*QC
     *              +HPOL*((PONULL*A3F+POL*B3F)*SPU*CHNULL
     *                -LL*(PONULL*B1F+POL*A1F)*SMU)*DQC)
        ENDIF
        IF (LQ.EQ.5) THEN
         CALL HSSAB0(X,Y,POL,HPOL, 2,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU*CHNULL
     *                -LL*(PONULL*B3F+POL*A3F)*SMU)*QB
     *              +HPOL*((PONULL*A3F+POL*B3F)*SPU*CHNULL
     *                -LL*(PONULL*B1F+POL*A1F)*SMU)*DQB)
        ENDIF
        IF (LQ.EQ.6) THEN
         CALL HSSAB0(X,Y,POL,HPOL, 1,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU*CHNULL
     *                -LL*(PONULL*B3F+POL*A3F)*SMU)*QT
     *              +HPOL*((PONULL*A3F+POL*B3F)*SPU*CHNULL
     *                -LL*(PONULL*B1F+POL*A1F)*SMU)*DQT)
        ENDIF
        IF (LQ.EQ.-1) THEN
         CALL HSSAB0(X,Y,POL,HPOL, 1,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU*CHNULL
     *                +LL*(PONULL*B3F+POL*A3F)*SMU)*QBU
     *              -HPOL*((PONULL*A3F+POL*B3F)*SPU*CHNULL
     *                +LL*(PONULL*B1F+POL*A1F)*SMU)*DQBU)
        ENDIF
        IF (LQ.EQ.-2) THEN
         CALL HSSAB0(X,Y,POL,HPOL, 2,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU*CHNULL
     *                +LL*(PONULL*B3F+POL*A3F)*SMU)*QBD
     *              -HPOL*((PONULL*A3F+POL*B3F)*SPU*CHNULL
     *                +LL*(PONULL*B1F+POL*A1F)*SMU)*DQBD)
        ENDIF
        IF (LQ.EQ.-3) THEN
         CALL HSSAB0(X,Y,POL,HPOL, 2,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU*CHNULL
     *                +LL*(PONULL*B3F+POL*A3F)*SMU)*QBS
     *              -HPOL*((PONULL*A3F+POL*B3F)*SPU*CHNULL
     *                +LL*(PONULL*B1F+POL*A1F)*SMU)*DQBS)
        ENDIF
        IF (LQ.EQ.-4) THEN
         CALL HSSAB0(X,Y,POL,HPOL, 1,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU*CHNULL
     *                +LL*(PONULL*B3F+POL*A3F)*SMU)*QBC
     *              -HPOL*((PONULL*A3F+POL*B3F)*SPU*CHNULL
     *                +LL*(PONULL*B1F+POL*A1F)*SMU)*DQBC)
        ENDIF
        IF (LQ.EQ.-5) THEN
         CALL HSSAB0(X,Y,POL,HPOL, 2,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU*CHNULL
     *                +LL*(PONULL*B3F+POL*A3F)*SMU)*QBB
     *              -HPOL*((PONULL*A3F+POL*B3F)*SPU*CHNULL
     *                +LL*(PONULL*B1F+POL*A1F)*SMU)*DQBB)
        ENDIF
        IF (LQ.EQ.-6) THEN
         CALL HSSAB0(X,Y,POL,HPOL, 1,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU*CHNULL
     *                +LL*(PONULL*B3F+POL*A3F)*SMU)*QBT
     *              -HPOL*((PONULL*A3F+POL*B3F)*SPU*CHNULL
     *                +LL*(PONULL*B1F+POL*A1F)*SMU)*DQBT)
        ENDIF
C
        IF (LQ.EQ.0.AND.IPDFOP.EQ.1) THEN
          CALL HSSAB0(X,Y,POL,HPOL,1,A1U,A3U,B1U,B3U)
          CALL HSSAB0(X,Y,POL,HPOL,2,A1D,A3D,B1D,B3D)
          F1U=CHNULL*(PONULL*A1U+POL*B1U)*SPU
          F3U=-LL*(PONULL*B3U+POL*A3U)*SMU
          F1D=CHNULL*(PONULL*A1D+POL*B1D)*SPU
          F3D=-LL*(PONULL*B3D+POL*A3D)*SMU
          G3U=CHNULL*(PONULL*A3U+POL*B3U)*SPU
          G1U=-LL*(PONULL*B1U+POL*A1U)*SMU
          G3D=CHNULL*(PONULL*A3D+POL*B3D)*SPU
          G1D=-LL*(PONULL*B1D+POL*A1D)*SMU
          CQP(1)=         (F1U+F3U)*QU +HPOL*(G3U+G1U)*DQU
          CQP(2) =CQP(1) +(F1U-F3U)*QBU-HPOL*(G3U-G1U)*DQBU
          CQP(3) =CQP(2) +(F1D+F3D)*QD +HPOL*(G3D+G1D)*DQD
          CQP(4) =CQP(3) +(F1D-F3D)*QBD-HPOL*(G3D-G1D)*DQBD
          CQP(5) =CQP(4) +(F1D+F3D)*QS +HPOL*(G3D+G1D)*DQS
          CQP(6) =CQP(5) +(F1D-F3D)*QBS-HPOL*(G3D-G1D)*DQBS
          CQP(7) =CQP(6) +(F1U+F3U)*QC +HPOL*(G3U+G1U)*DQC
          CQP(8) =CQP(7) +(F1U-F3U)*QBC-HPOL*(G3U-G1U)*DQBC
          CQP(9) =CQP(8) +(F1D+F3D)*QB +HPOL*(G3D+G1D)*DQB
          CQP(10)=CQP(9) +(F1D-F3D)*QBB-HPOL*(G3D-G1D)*DQBB
          CQP(11)=CQP(10)+(F1U+F3U)*QT +HPOL*(G3U+G1U)*DQT
          CQP(12)=CQP(11)+(F1U-F3U)*QBT-HPOL*(G3U-G1U)*DQBT
          HSSGNC=SPN*CQP(12)
          DO 11 I=1,12
            CQP(I)=CQP(I)*SPN
 11       CONTINUE
        ELSEIF (LQ.EQ.0.AND.IPDFOP.EQ.0) THEN
          R1=CHNULL*Y*Y*(1D0+2D0*MEI2/T)
          R2=CHNULL*(1D0-Y+MPRO2*T/SHAT/SHAT)/X
          R3=-LL*Y*(1D0-Y/2D0)
          CQP(12)=0D0
          DO 121 IB1=1,2
          DO 122 IB2=1,2
          CQP(12)=CQP(12)+F1(IB1,IB2)*R1+F2(IB1,IB2)*R2+F3(IB1,IB2)*R3
     *        -HPOL*(G1(IB1,IB2)*R3*2D0+G3(IB1,IB2)*R1+G4(IB1,IB2)*R2)
 122      CONTINUE
 121      CONTINUE
          HSSGNC=8D0*SPN*CQP(12)/T/T
        ELSEIF (LQ.EQ.0.AND.IPDFOP.GE.2) THEN
          R1=CHNULL*Y*Y*(1D0+2D0*MEI2/T)
          R2=CHNULL*(1D0-Y+MPRO2*T/SHAT/SHAT)/X
          R3=-LL*Y*(1D0-Y/2D0)
          SUMME=0D0
          DO 131 IB1=1,2
          DO 132 IB2=1,2
          SUMME=SUMME+F1(IB1,IB2)*R1+F2(IB1,IB2)*R2+F3(IB1,IB2)*R3
     *       -HPOL*(G1(IB1,IB2)*R3*2D0+G3(IB1,IB2)*R1+G4(IB1,IB2)*R2)
 132      CONTINUE
 131      CONTINUE
          HSSGNC=8D0*SPN*SUMME/T/T
          CALL HSSAB0(X,Y,POL,HPOL,1,A1U,A3U,B1U,B3U)
          CALL HSSAB0(X,Y,POL,HPOL,2,A1D,A3D,B1D,B3D)
          F1U=CHNULL*(PONULL*A1U+POL*B1U)*SPU
          F3U=-LL*(PONULL*B3U+POL*A3U)*SMU
          F1D=CHNULL*(PONULL*A1D+POL*B1D)*SPU
          F3D=-LL*(PONULL*B3D+POL*A3D)*SMU
          G3U=CHNULL*(PONULL*A3U+POL*B3U)*SPU
          G1U=-LL*(PONULL*B1U+POL*A1U)*SMU
          G3D=CHNULL*(PONULL*A3D+POL*B3D)*SPU
          G1D=-LL*(PONULL*B1D+POL*A1D)*SMU
          CQP(1)=         (F1U+F3U)*QU +HPOL*(G3U+G1U)*DQU
          CQP(2) =CQP(1) +(F1U-F3U)*QBU-HPOL*(G3U-G1U)*DQBU
          CQP(3) =CQP(2) +(F1D+F3D)*QD +HPOL*(G3D+G1D)*DQD
          CQP(4) =CQP(3) +(F1D-F3D)*QBD-HPOL*(G3D-G1D)*DQBD
          CQP(5) =CQP(4) +(F1D+F3D)*QS +HPOL*(G3D+G1D)*DQS
          CQP(6) =CQP(5) +(F1D-F3D)*QBS-HPOL*(G3D-G1D)*DQBS
          CQP(7) =CQP(6) +(F1U+F3U)*QC +HPOL*(G3U+G1U)*DQC
          CQP(8) =CQP(7) +(F1U-F3U)*QBC-HPOL*(G3U-G1U)*DQBC
          CQP(9) =CQP(8) +(F1D+F3D)*QB +HPOL*(G3D+G1D)*DQB
          CQP(10)=CQP(9) +(F1D-F3D)*QBB-HPOL*(G3D-G1D)*DQBB
          CQP(11)=CQP(10)+(F1U+F3U)*QT +HPOL*(G3U+G1U)*DQT
          CQP(12)=CQP(11)+(F1U-F3U)*QBT-HPOL*(G3U-G1U)*DQBT
          RNORM=8D0*SUMME/T/T/CQP(12)*SPN
          DO 14 I=1,12
            CQP(I)=CQP(I)*RNORM
 14       CONTINUE
        ENDIF
      ENDIF
C
C...1-LOOP CORRECTIONS
      IF (LPAR(2).EQ.1) THEN
        IF (LQ.EQ.1) THEN
         CALL HSSAB1(X,Y,LL,POL,HPOL, 1,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU
     *                -LL*(PONULL*B3F+POL*A3F)*SMU)*QU
     *              +HPOL*((PONULL*A3F+POL*B3F)*SPU
     *                -LL*(PONULL*B1F+POL*A1F)*SMU)*DQU)
        ENDIF
        IF (LQ.EQ.2) THEN
         CALL HSSAB1(X,Y,LL,POL,HPOL, 2,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU
     *                -LL*(PONULL*B3F+POL*A3F)*SMU)*QD
     *              +HPOL*((PONULL*A3F+POL*B3F)*SPU
     *                -LL*(PONULL*B1F+POL*A1F)*SMU)*DQD)
        ENDIF
        IF (LQ.EQ.3) THEN
         CALL HSSAB1(X,Y,LL,POL,HPOL, 2,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU
     *                -LL*(PONULL*B3F+POL*A3F)*SMU)*QS
     *              +HPOL*((PONULL*A3F+POL*B3F)*SPU
     *                -LL*(PONULL*B1F+POL*A1F)*SMU)*DQS)
        ENDIF
        IF (LQ.EQ.4) THEN
         CALL HSSAB1(X,Y,LL,POL,HPOL, 1,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU
     *                -LL*(PONULL*B3F+POL*A3F)*SMU)*QC
     *              +HPOL*((PONULL*A3F+POL*B3F)*SPU
     *                -LL*(PONULL*B1F+POL*A1F)*SMU)*DQC)
        ENDIF
        IF (LQ.EQ.5) THEN
         CALL HSSAB1(X,Y,LL,POL,HPOL, 2,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU
     *                -LL*(PONULL*B3F+POL*A3F)*SMU)*QB
     *              +HPOL*((PONULL*A3F+POL*B3F)*SPU
     *                -LL*(PONULL*B1F+POL*A1F)*SMU)*DQB)
        ENDIF
        IF (LQ.EQ.6) THEN
         CALL HSSAB1(X,Y,LL,POL,HPOL, 1,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU
     *                -LL*(PONULL*B3F+POL*A3F)*SMU)*QT
     *              +HPOL*((PONULL*A3F+POL*B3F)*SPU
     *                -LL*(PONULL*B1F+POL*A1F)*SMU)*DQT)
        ENDIF
        IF (LQ.EQ.-1) THEN
         CALL HSSAB1(X,Y,-LL,POL,HPOL, 1,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU
     *                +LL*(PONULL*B3F+POL*A3F)*SMU)*QBU
     *              -HPOL*((PONULL*A3F+POL*B3F)*SPU
     *                +LL*(PONULL*B1F+POL*A1F)*SMU)*DQBU)
        ENDIF
        IF (LQ.EQ.-2) THEN
         CALL HSSAB1(X,Y,-LL,POL,HPOL, 2,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU
     *                +LL*(PONULL*B3F+POL*A3F)*SMU)*QBD
     *              -HPOL*((PONULL*A3F+POL*B3F)*SPU
     *                +LL*(PONULL*B1F+POL*A1F)*SMU)*DQBD)
        ENDIF
        IF (LQ.EQ.-3) THEN
         CALL HSSAB1(X,Y,-LL,POL,HPOL, 2,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU
     *                +LL*(PONULL*B3F+POL*A3F)*SMU)*QBS
     *              -HPOL*((PONULL*A3F+POL*B3F)*SPU
     *                +LL*(PONULL*B1F+POL*A1F)*SMU)*DQBS)
        ENDIF
        IF (LQ.EQ.-4) THEN
         CALL HSSAB1(X,Y,-LL,POL,HPOL, 1,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU
     *                +LL*(PONULL*B3F+POL*A3F)*SMU)*QBC
     *              -HPOL*((PONULL*A3F+POL*B3F)*SPU
     *                +LL*(PONULL*B1F+POL*A1F)*SMU)*DQBC)
        ENDIF
        IF (LQ.EQ.-5) THEN
         CALL HSSAB1(X,Y,-LL,POL,HPOL, 2,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU
     *                +LL*(PONULL*B3F+POL*A3F)*SMU)*QBB
     *              -HPOL*((PONULL*A3F+POL*B3F)*SPU
     *                +LL*(PONULL*B1F+POL*A1F)*SMU)*DQBB)
        ENDIF
        IF (LQ.EQ.-6) THEN
         CALL HSSAB1(X,Y,-LL,POL,HPOL, 1,A1F,A3F,B1F,B3F)
         HSSGNC=SPN*(((PONULL*A1F+POL*B1F)*SPU
     *                +LL*(PONULL*B3F+POL*A3F)*SMU)*QBT
     *              -HPOL*((PONULL*A3F+POL*B3F)*SPU
     *                +LL*(PONULL*B1F+POL*A1F)*SMU)*DQBT)
        ENDIF
C
        IF (LQ.EQ.0.AND.IPDFOP.EQ.1) THEN
          CALL HSSAB1(X,Y,LL,POL,HPOL,1,A1U,A3U,B1U,B3U)
          CALL HSSAB1(X,Y,LL,POL,HPOL,2,A1D,A3D,B1D,B3D)
          F13U=(PONULL*A1U+POL*B1U)*SPU
     *        - LL*(PONULL*B3U+POL*A3U)*SMU
          F13D=(PONULL*A1D+POL*B1D)*SPU
     *        - LL*(PONULL*B3D+POL*A3D)*SMU
          G13U=(PONULL*A3U+POL*B3U)*SPU
     *        - LL*(PONULL*B1U+POL*A1U)*SMU
          G13D=(PONULL*A3D+POL*B3D)*SPU
     *        - LL*(PONULL*B1D+POL*A1D)*SMU
          CALL HSSAB1(X,Y,-LL,POL,HPOL,1,A1BU,A3BU,B1BU,B3BU)
          CALL HSSAB1(X,Y,-LL,POL,HPOL,2,A1BD,A3BD,B1BD,B3BD)
          F13BU=(PONULL*A1BU+POL*B1BU)*SPU
     *         + LL*(PONULL*B3BU+POL*A3BU)*SMU
          F13BD=(PONULL*A1BD+POL*B1BD)*SPU
     *         + LL*(PONULL*B3BD+POL*A3BD)*SMU
          G13BU=(PONULL*A3BU+POL*B3BU)*SPU
     *         + LL*(PONULL*B1BU+POL*A1BU)*SMU
          G13BD=(PONULL*A3BD+POL*B3BD)*SPU
     *         + LL*(PONULL*B1BD+POL*A1BD)*SMU
          CQP(1) =         F13U*QU  +HPOL*G13U*DQU
          CQP(2) =CQP(1) + F13BU*QBU-HPOL*G13BU*DQBU
          CQP(3) =CQP(2) + F13D*QD  +HPOL*G13D*DQD
          CQP(4) =CQP(3) + F13BD*QBD-HPOL*G13BD*DQBD
          CQP(5) =CQP(4) + F13D*QS  +HPOL*G13D*DQS
          CQP(6) =CQP(5) + F13BD*QBS-HPOL*G13BD*DQBS
          CQP(7) =CQP(6) + F13U*QC  +HPOL*G13U*DQC
          CQP(8) =CQP(7) + F13BU*QBC-HPOL*G13BU*DQBC
          CQP(9) =CQP(8) + F13D*QB  +HPOL*G13D*DQB
          CQP(10)=CQP(9) + F13BD*QBB-HPOL*G13BD*DQBB
          CQP(11)=CQP(10)+ F13U*QT  +HPOL*G13U*DQT
          CQP(12)=CQP(11)+ F13BU*QBT-HPOL*G13BU*DQBT
          DO 21 I=1,12
            CQP(I)=CQP(I)*SXNORM*SP*X
 21       CONTINUE
          HSSGNC=CQP(12)
        ELSEIF (LQ.EQ.0.AND.IPDFOP.EQ.0) THEN
          R1=CHNULL*Y*Y*(1D0+2D0*MEI2/T)
          R2=CHNULL*(1D0-Y+MPRO2*T/SHAT/SHAT)/X
          R3=-LL*Y*(1D0-Y/2D0)
          CQP(12)=0D0
          DO 221 IB1=1,2
          DO 222 IB2=1,2
          CQP(12)=CQP(12)+F1(IB1,IB2)*R1+F2(IB1,IB2)*R2+F3(IB1,IB2)*R3
     *         -HPOL*(G1(IB1,IB2)*R3*2D0+G3(IB1,IB2)*R1+G4(IB1,IB2)*R2)
 222      CONTINUE
 221      CONTINUE
          HSSGNC=8D0*SPN*CQP(12)/T/T
        ELSEIF (LQ.EQ.0.AND.IPDFOP.GE.2) THEN
          R1=CHNULL*Y*Y*(1D0+2D0*MEI2/T)
          R2=CHNULL*(1D0-Y+MPRO2*T/SHAT/SHAT)/X
          R3=-LL*Y*(1D0-Y/2D0)
          SUMME=0D0
          DO 231 IB1=1,2
          DO 232 IB2=1,2
          SUMME=SUMME+F1(IB1,IB2)*R1+F2(IB1,IB2)*R2+F3(IB1,IB2)*R3
     *        -HPOL*(G1(IB1,IB2)*R3*2D0+G3(IB1,IB2)*R1+G4(IB1,IB2)*R2)
 232      CONTINUE
 231      CONTINUE
          HSSGNC=8D0*SPN*SUMME/T/T
          CALL HSSAB1(X,Y,LL,POL,HPOL,1,A1U,A3U,B1U,B3U)
          CALL HSSAB1(X,Y,LL,POL,HPOL,2,A1D,A3D,B1D,B3D)
          F13U=(PONULL*A1U+POL*B1U)*SPU
     *        - LL*(PONULL*B3U+POL*A3U)*SMU
          F13D=(PONULL*A1D+POL*B1D)*SPU
     *        - LL*(PONULL*B3D+POL*A3D)*SMU
          G13U=(PONULL*A3U+POL*B3U)*SPU
     *        - LL*(PONULL*B1U+POL*A1U)*SMU
          G13D=(PONULL*A3D+POL*B3D)*SPU
     *        - LL*(PONULL*B1D+POL*A1D)*SMU
          CALL HSSAB1(X,Y,-LL,POL,HPOL,1,A1BU,A3BU,B1BU,B3BU)
          CALL HSSAB1(X,Y,-LL,POL,HPOL,2,A1BD,A3BD,B1BD,B3BD)
          F13BU=(PONULL*A1BU+POL*B1BU)*SPU
     *         + LL*(PONULL*B3BU+POL*A3BU)*SMU
          F13BD=(PONULL*A1BD+POL*B1BD)*SPU
     *         + LL*(PONULL*B3BD+POL*A3BD)*SMU
          G13BU=(PONULL*A3BU+POL*B3BU)*SPU
     *         + LL*(PONULL*B1BU+POL*A1BU)*SMU
          G13BD=(PONULL*A3BD+POL*B3BD)*SPU
     *         + LL*(PONULL*B1BD+POL*A1BD)*SMU
          CQP(1) =         F13U*QU  +HPOL*G13U*DQU
          CQP(2) =CQP(1) + F13BU*QBU-HPOL*G13BU*DQBU
          CQP(3) =CQP(2) + F13D*QD  +HPOL*G13D*DQD
          CQP(4) =CQP(3) + F13BD*QBD-HPOL*G13BD*DQBD
          CQP(5) =CQP(4) + F13D*QS  +HPOL*G13D*DQS
          CQP(6) =CQP(5) + F13BD*QBS-HPOL*G13BD*DQBS
          CQP(7) =CQP(6) + F13U*QC  +HPOL*G13U*DQC
          CQP(8) =CQP(7) + F13BU*QBC-HPOL*G13BU*DQBC
          CQP(9) =CQP(8) + F13D*QB  +HPOL*G13D*DQB
          CQP(10)=CQP(9) + F13BD*QBB-HPOL*G13BD*DQBB
          CQP(11)=CQP(10)+ F13U*QT  +HPOL*G13U*DQT
          CQP(12)=CQP(11)+ F13BU*QBT-HPOL*G13BU*DQBT
          RNORM=8D0*SUMME/T/T/CQP(12)*SPN
          DO 24 I=1,12
            CQP(I)=CQP(I)*RNORM
 24       CONTINUE
        ENDIF
C
      ENDIF
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C EPNCEMRC.S(NCVS)
C
C                EP DEEP INELASTIC SCATTERING FOR HERA
C
C                    VIRTUAL CORRECTIONS AND
C                     SOFT BREMSSTRAHLUNG
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC 22. 9. 86   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC H. SPIESBERGER CCC
C
      SUBROUTINE HSSAB0(X,Y,POL,HPOL,LQ,A1,A3,B1,B3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     BORN CROSS SECTION FOR DEEP INELASTIC ELECTRON QUARK
C     SCATTERING
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
C
      T=-(SP-MEI2-MPRO2)*X*Y
      DG=1D0/T
      DZ=1D0/(T-MZ2)
C  STRUCTURE FUNCTIONS
C
        F1GG=FLIND(1,LQ+1,1,1)
        F1GZ=FLIND(1,LQ+1,1,2)
        F1ZZ=FLIND(1,LQ+1,2,2)
        F3GG=-FLIND(2,LQ+1,1,1)
        F3GZ=-FLIND(2,LQ+1,1,2)
        F3ZZ=-FLIND(2,LQ+1,2,2)

C...CONSTRUCTING THE CROSS SECTION
C
      A1=F1GG*DG*DG*FLIND(1,1,1,1)
     *   +2D0*F1GZ*DG*DZ*FLIND(1,1,1,2)
     *   +F1ZZ*DZ*DZ*FLIND(1,1,2,2)
C
      B3=-F3GG*DG*DG*FLIND(2,1,1,1)
     *   -2D0*F3GZ*DG*DZ*FLIND(2,1,1,2)
     *   -F3ZZ*DZ*DZ*FLIND(2,1,2,2)
C
      IF (POL.EQ.0D0.AND.HPOL.EQ.0D0) THEN
        B1=0D0
        A3=0D0
        ELSE
        B1=-F1GG*DG*DG*FLIND(2,1,1,1)
     *     -2D0*F1GZ*DG*DZ*FLIND(2,1,1,2)
     *     -F1ZZ*DZ*DZ*FLIND(2,1,2,2)
C
        A3=F3GG*DG*DG*FLIND(1,1,1,1)
     *     +2D0*F3GZ*DG*DZ*FLIND(1,1,1,2)
     *     +F3ZZ*DZ*DZ*FLIND(1,1,2,2)
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE HSSETF(T)
C---FORM FACTORS
C   RECALCULATE EFFECTIVE COUPLING CONSTANTS INCLUDING SELF ENERGIES
C   RUNNING ALPHA, EFFECTIVE SW2, FORM FACTOR KAPPA (FROM PI_Z)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSRGG,HSSRGZ,HSSRZZ,CG,CM,CZ,HSFHFB
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSASY0/ PONULL,CHNULL,ICPASY
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSDELR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSSMC1/ VAFI1(2,3,2)
     &               ,AFIJ1(3,2,2),BFIJ1(3,2,2),FLIND1(2,3,2,2)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSKNCC/ SXNRCC,SX1NCC
      COMMON /HSFRFF/ ALPFFQ,AKAPPA,GMUFFQ,SWEFF2
C
C  SELF ENERGIES
C
      IF (LPAR(7).GE.1) THEN
        CG=HSSRGG(T)
        PIGGG=DREAL(CG)/T
        ELSE
        PIGGG=0D0
      ENDIF
      ALPFFQ=1D0/(1D0+PIGGG)
      EFFQ=SQRT(ALPFFQ)
      IF (LPAR(8).EQ.1) THEN
        CM=HSSRGZ(T)
        SM=DREAL(CM)/T
        ELSE
        SM=0D0
      ENDIF
      AKAPPA=1D0-CW/SW*SM/(1D0+PIGGG)
      SWEFF2=SW2*AKAPPA
      IF (LPAR(9).EQ.1) THEN
        CZ=HSSRZZ(T)
        PIGZZ=DREAL(CZ)/(T-MZ2)
        ELSE
        PIGZZ=0D0
      ENDIF
      GMUFFQ=1D0/(1D0+PIGZZ)
      BFFQ=SQRT(GMUFFQ)
      IF (LPAR(3).GE.1) BFFQ=BFFQ*BTOP4
C
C---DEFINE FERMION GAUGE BOSON COUPLING CONSTANTS
      IF (LPAR(4).EQ.1) THEN
        B0=1D0/4D0/CW/SW
        B=B0
        ELSE
C---NORMALIZED TO G-MU
        B0=MZ/SQRT(AGF0)/4D0
        B=B0
        IF (LPAR(9).GE.1) B=B0*SQRT(1D0-DELTAR)
      ENDIF
      IF (LPAR(2).EQ.1.AND.LPAR(9).GE.1) B=B*BFFQ
      RHO1=B/B0

      VAFI1(2,1,1)=0D0
      VAFI1(2,2,1)=0D0
      VAFI1(2,3,1)=0D0
      VAFI1(2,1,2)=-B
      VAFI1(2,2,2)=B
      VAFI1(2,3,2)=-B
      VAFI1(1,1,1)=EFFQ
      VAFI1(1,2,1)=-2D0/3D0*EFFQ
      VAFI1(1,3,1)=1D0/3D0*EFFQ
      VAFI1(1,1,2)=B*(4D0*SWEFF2-1D0)
      VAFI1(1,2,2)=B*(1D0-8D0*SWEFF2/3D0)
      VAFI1(1,3,2)=B*(4D0*SWEFF2/3D0-1D0)
C
C..VERTEX CORRECTIONS (NO QED PARTS)
      LPRK11=LPAR(11)
      LPAR(11)=0
C..ELECTRON VERTEX
      IF (LPAR(12).EQ.1) THEN
        LF=1
        DO 1 IVA=1,2
        VAFI1(IVA,LF,1)=VAFI1(IVA,LF,1)
     *                   +EFFQ*DREAL(HSFHFB(T,IVA,LF,1,MEI2))
        VAFI1(IVA,LF,2)=VAFI1(IVA,LF,2)
     *                   +RHO1*DREAL(HSFHFB(T,IVA,LF,2,MEI2))
    1   CONTINUE
      ENDIF
      LPAR(11)=LPRK11
C..QUARK VERTEX
      IF (LPAR(13).EQ.1) THEN
        DO 21 LF=2,3
        DO 22 IVA=1,2
        VAFI1(IVA,LF,1)=VAFI1(IVA,LF,1)
     *                   +EFFQ*DREAL(HSFHFB(T,IVA,LF,1,MQI2))
        VAFI1(IVA,LF,2)=VAFI1(IVA,LF,2)
     *                   +RHO1*DREAL(HSFHFB(T,IVA,LF,2,MQI2))
 22     CONTINUE
 21     CONTINUE
      ENDIF
C
      IGAMMA = 1
      IZ = 2
      IEL = 1
      IFU = 2
      IFD = 3
      INDV = 1
      INDA = 2
C
      DO 31 IF=IEL,IFD
        DO 32 IB1=IGAMMA,IZ
          DO 33 IB2=IGAMMA,IZ
          FLIND1(INDV,IF,IB1,IB2)=
     *     2D0*(VAFI1(INDV,IF,IB1)*VAFI1(INDV,IF,IB2)
     *         +VAFI1(INDA,IF,IB1)*VAFI1(INDA,IF,IB2))
 33   CONTINUE
 32   CONTINUE
 31   CONTINUE
      DO 41 IF=IEL,IFD
        DO 42 IB1=IGAMMA,IZ
          DO 43 IB2=IGAMMA,IZ
          FLIND1(INDA,IF,IB1,IB2)=
     *     2D0*(VAFI1(INDV,IF,IB1)*VAFI1(INDA,IF,IB2)
     *         +VAFI1(INDA,IF,IB1)*VAFI1(INDV,IF,IB2))
 43   CONTINUE
 42   CONTINUE
 41   CONTINUE
C
      DO 51 IVB1 = IGAMMA, IZ
       DO 52 IVB2 = IGAMMA, IZ
        DO 53 IFERM = IFU, IFD
        AFIJ1(IFERM,IVB1,IVB2)=FLIND1(INDV,IFERM,IVB1,IVB2)
     *  *(PONULL*FLIND1(INDV,IEL,IVB1,IVB2)
     *    -POLARI*FLIND1(INDA,IEL,IVB1,IVB2))
        BFIJ1(IFERM,IVB1,IVB2)=FLIND1(INDA,IFERM,IVB1,IVB2)
     *  *(PONULL*FLIND1(INDA,IEL,IVB1,IVB2)
     *    -POLARI*FLIND1(INDV,IEL,IVB1,IVB2))
 53   CONTINUE
 52   CONTINUE
 51   CONTINUE
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE HSSAB1(X,Y,LL,POL,HPOL,LQ,RA1,RA3,RB1,RB3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     1-LOOP VIRTUAL CORRECTIONS TO DEEP INELASTIC ELECTRON QUARK
C     SCATTERING
C     LL = +1,-1: POSITRONS, ELECTRONS
C     POL  = +1,-1: RIGHT-, LEFTHANDED LEPTONS
C
C--->
C--->  NEW VERSION 07.02.1991
C--->  FORM FACTOR APPROACH
C--->
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) CMZ2,CMW2
     *       ,HSBRNC,BSQ,HSBCGA,HSBXCV,HSBXCA,HSBXI0,HSBXI5
     *       ,CVEG,HSFHFB
      COMPLEX(kind=8) BSQC,BSQ0
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSCBMS/ CMW2,CMZ2
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSSMC1/ VAFI1(2,3,2)
     &               ,AFIJ1(3,2,2),BFIJ1(3,2,2),FLIND1(2,3,2,2)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSASY0/ PONULL,CHNULL,ICPASY
      COMMON /HSFRFF/ ALPFFQ,AKAPPA,GMUFFQ,SWEFF2
C
      IF (LPAR(2).EQ.0) THEN
        RA1=0D0
        RA3=0D0
        RB1=0D0
        RB3=0D0
        RETURN
      ENDIF
C
      LQF=MOD(LQ,2)
      IF (LQF.EQ.0) THEN
        LQF=3
        CQ=-1D0/3D0
        ELSE
        LQF=2
        CQ=2D0/3D0
      ENDIF
C
      T=-(SP-MEI2-MPRO2)*X*Y
      DG=1D0/T
      DZ=1D0/(T-MZ2)
C
C..COUPLING CONSTANTS INCLUDING RUNNING ALPHA AND WEAK FORM FACTORS
      CALL HSSETF(T)
      VEZ=VAFI1(1,1,2)
      AEZ=VAFI1(2,1,2)
C
C..STRUCTURE FUNCTIONS
      F1GG=FLIND1(1,LQF,1,1)
      F1GZ=FLIND1(1,LQF,1,2)
      F1ZZ=FLIND1(1,LQF,2,2)
      F3GG=-FLIND1(2,LQF,1,1)
      F3GZ=-FLIND1(2,LQF,1,2)
      F3ZZ=-FLIND1(2,LQF,2,2)
      EPEGG=FLIND1(1,1,1,1)
      EPEGZ=FLIND1(1,1,1,2)
      EPEZZ=FLIND1(1,1,2,2)
      EMEGG=-FLIND1(2,1,1,1)
      EMEGZ=-FLIND1(2,1,1,2)
      EMEZZ=-FLIND1(2,1,2,2)

C..SOFT BREMSSTRAHLUNG (INCLUDING QUARK SELF ENERGY)
      IF (LPAR(11).EQ.1) THEN
        BSQ=HSBRNC(X,Y,LL,LQ)
        RBSQ=DREAL(BSQ)*ALP2PI
        if (icpasy.eq.2) then
          BSQC=(HSBRNC(X,Y,LL,LQ)-HSBRNC(X,Y,-LL,LQ))/2d0 
          BSQ0=(HSBRNC(X,Y,LL,LQ)+HSBRNC(X,Y,-LL,LQ))/2d0
          RBSQC=DREAL(BSQC)*ALP2PI
          RBSQ0=DREAL(BSQ0)*ALP2PI
        endif
        CALL HSSAB0(X,Y,POL,HPOL,LQ,RA10,RA30,RB10,RB30)
C..QED VERTEX
        LPRK15=LPAR(15)
        LPAR(15)=0
        CVEG=HSFHFB(T,1,1,1,MEI2)
        LPAR(15)=LPRK15
        DVERTX=DREAL(CVEG)
        RBSQ=RBSQ+2D0*DVERTX
        if (icpasy.eq.2) then
          RBSQC=RBSQC
          RBSQ0=RBSQ0+2D0*DVERTX
        endif
        ELSE
        RBSQ=0D0
        RA10=0D0
        RA30=0D0
        RB10=0D0
        RB30=0D0
      ENDIF

C..BOX DIAGRAMS
      IF (LPAR(14).EQ.1) THEN
        S=SP*X
        U=-(T+S)
        BXGG = -LL*ALP1PI*DREAL(HSBCGA(T,S)-HSBCGA(T,U))
        BXGG5=  ALP1PI*DREAL(HSBCGA(T,S)+HSBCGA(T,U))
        BXGZ = -LL*2D0*ALP1PI*DLOG((MZ2-T)/(-T))*DLOG(-U/S)
     *         -LL*ALP1PI*DREAL(HSBXCV(T,S,CMZ2)-HSBXCV(T,U,CMZ2))
        BXGZ5=  ALP1PI*DREAL( HSBXCA(T,S,CMZ2) + HSBXCA(T,U,CMZ2) )
        if (icpasy.eq.2) then
        BXGGC = -LL*ALP1PI*DREAL(HSBCGA(T,S)-HSBCGA(T,U))
        BXGG5C= 0D0
        BXGZC = -LL*2D0*ALP1PI*DLOG((MZ2-T)/(-T))*DLOG(-U/S)
     *         -LL*ALP1PI*DREAL(HSBXCV(T,S,CMZ2)-HSBXCV(T,U,CMZ2))
        BXGZ5C= 0D0
        BXGG0 = 0D0
        BXGG50= ALP1PI*DREAL(HSBCGA(T,S)+HSBCGA(T,U))
        BXGZ0 = 0D0
        BXGZ50= ALP1PI*DREAL( HSBXCA(T,S,CMZ2) + HSBXCA(T,U,CMZ2) )
        endif
        F1BGG=F1GG*CQ
        F1BGZ=F1GZ*CQ
        F1BZZ=F1ZZ*CQ
        F3BGG=F3GG*CQ
        F3BGZ=F3GZ*CQ
        F3BZZ=F3ZZ*CQ
      ELSE
        BXGG=0D0
        BXGG5=0D0
        BXGZ=0D0
        BXGZ5=0D0
        F1BGG=0D0
        F1BGZ=0D0
        F1BZZ=0D0
        F3BGG=0D0
        F3BGZ=0D0
        F3BZZ=0D0
      ENDIF
      IF (LPAR(16).EQ.1) THEN
        S=SP*X
        U=-(T+S)
        F1BZZ=F1ZZ*CQ
        F3BZZ=F3ZZ*CQ
        BXZZ = -LL*ALP1PI*DREAL(HSBXI0(T,U,CMZ2) - HSBXI0(T,S,CMZ2))
        BXZZ5=     ALP1PI*DREAL(HSBXI5(T,U,CMZ2) + HSBXI5(T,S,CMZ2))
        if (icpasy.eq.2) then
        BXZZC = -LL*ALP1PI*DREAL(HSBXI0(T,U,CMZ2) - HSBXI0(T,S,CMZ2))
        BXZZ5C= 0D0
        BXZZ0 = 0D0
        BXZZ50= ALP1PI*DREAL(HSBXI5(T,U,CMZ2) + HSBXI5(T,S,CMZ2))
        endif
        FGBWW = (VAFI1(1,LQF,1)+VAFI1(2,LQF,1))
     *         *(VAFI1(1,  1,1)+VAFI1(2,  1,1))/4D0/SW2/SW2
        FZBWW = (VAFI1(1,LQF,2)+VAFI1(2,LQF,2))
     *         *(VAFI1(1,  1,2)+VAFI1(2,  1,2))/4D0/SW2/SW2
        FBZZZ= FLIND1(1,LQF,2,2)*VAFI1(1,LQF,2)
     *        +FLIND1(2,LQF,2,2)*VAFI1(2,LQF,2)
        FBZZ5= FLIND1(1,LQF,2,2)*VAFI1(2,LQF,2)
     *        +FLIND1(2,LQF,2,2)*VAFI1(1,LQF,2)
        IF (LL.LT.0) THEN
          IF (MOD(LQ,2).EQ.1) THEN
            BXWW =  ALP1PI*DREAL(-HSBXI0(T,S,CMW2)+HSBXI5(T,S,CMW2))
            ELSE
            BXWW =  ALP1PI*DREAL( HSBXI0(T,U,CMW2)+HSBXI5(T,U,CMW2))
          ENDIF
        ELSE
          IF (MOD(LQ,2).EQ.1) THEN
            BXWW =  ALP1PI*DREAL(-HSBXI0(T,U,CMW2)+HSBXI5(T,U,CMW2))
            ELSE
            BXWW =  ALP1PI*DREAL( HSBXI0(T,S,CMW2)+HSBXI5(T,S,CMW2))
          ENDIF
        ENDIF
        if (icpasy.eq.2) then
          IF (MOD(LQ,2).EQ.1) THEN
            BXWWC=(ALP1PI*DREAL(-HSBXI0(T,U,CMW2)+HSBXI5(T,U,CMW2))
     *           -ALP1PI*DREAL(-HSBXI0(T,S,CMW2)+HSBXI5(T,S,CMW2)))/2D0
            BXWW0=(ALP1PI*DREAL(-HSBXI0(T,U,CMW2)+HSBXI5(T,U,CMW2))
     *           +ALP1PI*DREAL(-HSBXI0(T,S,CMW2)+HSBXI5(T,S,CMW2)))/2D0
            
          ELSE
            BXWWC=(ALP1PI*DREAL(HSBXI0(T,S,CMW2)+HSBXI5(T,S,CMW2))
     *           -ALP1PI*DREAL(HSBXI0(T,U,CMW2)+HSBXI5(T,U,CMW2)))/2D0
            BXWW0=(ALP1PI*DREAL(HSBXI0(T,S,CMW2)+HSBXI5(T,S,CMW2))
     *           +ALP1PI*DREAL(HSBXI0(T,U,CMW2)+HSBXI5(T,U,CMW2)))/2D0
          ENDIF
          BXWWC=BXWWC*LL
        endif
      ELSE
        BXZZ=0D0
        BXZZ5=0D0
        FGBWW=0D0
        FZBWW=0D0
        FBZZZ=0D0
        FBZZ5=0D0
        BXWW=0D0
      ENDIF
C
C..CONSTRUCTING THE CROSS SECTION
      IF (LPAR(3).LT.3) THEN
        A1=F1GG*DG*DG*EPEGG+2D0*F1GZ*DG*DZ*EPEGZ+F1ZZ*DZ*DZ*EPEZZ
     *     +RBSQ*RA10
        ELSE
        A1=(F1GG*DG*DG*EPEGG+2D0*F1GZ*DG*DZ*EPEGZ+F1ZZ*DZ*DZ*EPEZZ)
     *     *(1D0+RBSQ)
      ENDIF
C
      IF (LPAR(14).EQ.1) THEN
        A1 = A1 +
     *     BXGG*EPEGG*DG*DG*F1BGG
     *    +(BXGG+BXGZ)*EPEGZ*DG*DZ*F1BGZ
     *    -(BXGG5+BXGZ5)*EMEGZ*DG*DZ*F3BGZ
     *    +BXGZ*EPEZZ*DZ*DZ*F1BZZ - BXGZ5*EMEZZ*DZ*DZ*F3BZZ
      ENDIF
      IF (LPAR(16).EQ.1) THEN
        A1 = A1 +
     *     BXZZ*FLIND1(1,1,2,2)*DG*DG*F1BZZ
     *    -BXZZ5*DG*DG*FLIND1(2,1,2,2)*F3BZZ
     *    +BXZZ*(FLIND1(1,1,2,2)*VAFI1(1,1,2)
     *           +FLIND1(2,1,2,2)*VAFI1(2,1,2))*DG*DZ*FBZZZ
     *    -BXZZ5*(FLIND1(1,1,2,2)*VAFI1(2,1,2)
     *           +FLIND1(2,1,2,2)*VAFI1(1,1,2))*DG*DZ*FBZZ5
     *    +DG*(DG*FGBWW+DZ*FZBWW)*BXWW
      ENDIF
C
      IF (LPAR(3).LT.3) THEN
        B3=F3GG*DG*DG*EMEGG+2D0*F3GZ*DG*DZ*EMEGZ+F3ZZ*DZ*DZ*EMEZZ
     *     +RBSQ*RB30
        ELSE
        B3=(F3GG*DG*DG*EMEGG+2D0*F3GZ*DG*DZ*EMEGZ+F3ZZ*DZ*DZ*EMEZZ)
     *     *(1D0+RBSQ)
      ENDIF
C
      IF (LPAR(14).EQ.1) THEN
        B3 = B3
     *    -EPEGG*DG*DG*F1BGG*BXGG5
     *    +(BXGG+BXGZ)*EMEGZ*DG*DZ*F3BGZ
     *    -(BXGG5+BXGZ5)*EPEGZ*DG*DZ*F1BGZ
     *    +BXGZ*EMEZZ*DZ*DZ*F3BZZ - BXGZ5*EPEZZ*DZ*DZ*F1BZZ
      ENDIF
      IF (LPAR(16).EQ.1) THEN
        B3 = B3 +
     *     BXZZ*FLIND1(2,1,2,2)*DG*DG*F3BZZ
     *    -BXZZ5*DG*DG*FLIND1(1,1,2,2)*F1BZZ
     *    +BXZZ*(FLIND1(1,1,2,2)*VAFI1(2,1,2)
     *           +FLIND1(2,1,2,2)*VAFI1(1,1,2))*DG*DZ*FBZZ5
     *    -BXZZ5*(FLIND1(1,1,2,2)*VAFI1(1,1,2)
     *           +FLIND1(2,1,2,2)*VAFI1(2,1,2))*DG*DZ*FBZZZ
     *    +DG*(DG*FGBWW+DZ*FZBWW)*BXWW
      ENDIF
      RA1=A1
      RB3=B3
C
      IF (LPAR(3).LT.3) THEN
        B1=F1GG*DG*DG*EMEGG+2D0*F1GZ*DG*DZ*EMEGZ+F1ZZ*DZ*DZ*EMEZZ
     *     +RBSQ*RB10
        ELSE
        B1=(F1GG*DG*DG*EMEGG+2D0*F1GZ*DG*DZ*EMEGZ+F1ZZ*DZ*DZ*EMEZZ)
     *     *(1D0+RBSQ)
      ENDIF
C
      IF (LPAR(14).EQ.1) THEN
        B1 = B1
     *    +(BXGG+BXGZ)*EMEGZ*DG*DZ*F1BGZ
     *    -(BXGG5+BXGZ5)*EPEGZ*DG*DZ*F3BGZ
     *    +BXGZ*EMEZZ*DZ*DZ*F1BZZ - BXGZ5*EPEZZ*DZ*DZ*F3BZZ
      ENDIF
      IF (LPAR(16).EQ.1) THEN
        B1 = B1 +
     *     BXZZ*EMEZZ*DG*DG*F1BZZ - BXZZ5*DG*DG*EPEZZ*F3BZZ
     *    -BXZZ*(EPEZZ*AEZ-EMEZZ*VEZ)*DG*DZ*FBZZZ
     *    +BXZZ5*(EPEZZ*VEZ-EMEZZ*AEZ)*DG*DZ*FBZZ5
     *    -DG*(DG*FGBWW+DZ*FZBWW)*BXWW
      ENDIF
C
      IF (LPAR(3).LT.3) THEN
        A3=F3GG*DG*DG*EPEGG+2D0*F3GZ*DG*DZ*EPEGZ+F3ZZ*DZ*DZ*EPEZZ
     *     +RBSQ*RA30
        ELSE
        A3=(F3GG*DG*DG*EPEGG+2D0*F3GZ*DG*DZ*EPEGZ+F3ZZ*DZ*DZ*EPEZZ)
     *     *(1D0+RBSQ)
      ENDIF
C
      IF (LPAR(14).EQ.1) THEN
        A3 = A3 +
     *     (BXGG+BXGZ)*EPEGZ*DG*DZ*F3BGZ
     *    -(BXGG5+BXGZ5)*EMEGZ*DG*DZ*F1BGZ
     *    +BXGZ*EPEZZ*DZ*DZ*F3BZZ - BXGZ5*EMEZZ*DZ*DZ*F1BZZ
      ENDIF
      IF (LPAR(16).EQ.1) THEN
        A3 = A3 +
     *     BXZZ*EPEZZ*DG*DG*F3BZZ - BXZZ5*DG*DG*EMEZZ*F1BZZ
     *    +BXZZ*(EPEZZ*VEZ-EMEZZ*AEZ)*DG*DZ*FBZZ5
     *    -BXZZ5*(EPEZZ*AEZ-EMEZZ*VEZ)*DG*DZ*FBZZZ
     *    -DG*(DG*FGBWW+DZ*FZBWW)*BXWW
      ENDIF
      RB1=B1
      RA3=A3


C...To construct the charge-odd terms in HSSGNC, need only the following
      if (icpasy.eq.2) then
      IF (LPAR(3).LT.3) THEN
        A1=+RBSQC*RA10
        ELSE
        A1=(F1GG*DG*DG*EPEGG+2D0*F1GZ*DG*DZ*EPEGZ+F1ZZ*DZ*DZ*EPEZZ)
     *     *RBSQC
      ENDIF
C
      IF (LPAR(14).EQ.1) THEN
        A1 = A1 +
     *     BXGGC*EPEGG*DG*DG*F1BGG
     *    +(BXGGC+BXGZC)*EPEGZ*DG*DZ*F1BGZ
     *    -(BXGG5C+BXGZ5C)*EMEGZ*DG*DZ*F3BGZ
     *    +BXGZC*EPEZZ*DZ*DZ*F1BZZ - BXGZ5C*EMEZZ*DZ*DZ*F3BZZ
      ENDIF
      IF (LPAR(16).EQ.1) THEN
        A1 = A1 +
     *     BXZZC*FLIND1(1,1,2,2)*DG*DG*F1BZZ
     *    -BXZZ5C*DG*DG*FLIND1(2,1,2,2)*F3BZZ
     *    +BXZZC*(FLIND1(1,1,2,2)*VAFI1(1,1,2)
     *           +FLIND1(2,1,2,2)*VAFI1(2,1,2))*DG*DZ*FBZZZ
     *    -BXZZ5C*(FLIND1(1,1,2,2)*VAFI1(2,1,2)
     *           +FLIND1(2,1,2,2)*VAFI1(1,1,2))*DG*DZ*FBZZ5
     *    +DG*(DG*FGBWW+DZ*FZBWW)*BXWWC
      ENDIF
C
      IF (LPAR(3).LT.3) THEN
        B3=F3GG*DG*DG*EMEGG+2D0*F3GZ*DG*DZ*EMEGZ+F3ZZ*DZ*DZ*EMEZZ
     *     +RBSQ0*RB30
        ELSE
        B3=(F3GG*DG*DG*EMEGG+2D0*F3GZ*DG*DZ*EMEGZ+F3ZZ*DZ*DZ*EMEZZ)
     *     *(1D0+RBSQ0)
      ENDIF
C
      IF (LPAR(14).EQ.1) THEN
        B3 = B3
     *    -EPEGG*DG*DG*F1BGG*BXGG50
     *    +(BXGG0+BXGZ0)*EMEGZ*DG*DZ*F3BGZ
     *    -(BXGG50+BXGZ50)*EPEGZ*DG*DZ*F1BGZ
     *    +BXGZ0*EMEZZ*DZ*DZ*F3BZZ - BXGZ50*EPEZZ*DZ*DZ*F1BZZ
      ENDIF
      IF (LPAR(16).EQ.1) THEN
        B3 = B3 +
     *     BXZZ0*FLIND1(2,1,2,2)*DG*DG*F3BZZ
     *    -BXZZ50*DG*DG*FLIND1(1,1,2,2)*F1BZZ
     *    +BXZZ0*(FLIND1(1,1,2,2)*VAFI1(2,1,2)
     *           +FLIND1(2,1,2,2)*VAFI1(1,1,2))*DG*DZ*FBZZ5
     *    -BXZZ50*(FLIND1(1,1,2,2)*VAFI1(1,1,2)
     *           +FLIND1(2,1,2,2)*VAFI1(2,1,2))*DG*DZ*FBZZZ
     *    +DG*(DG*FGBWW+DZ*FZBWW)*BXWW0
      ENDIF
      RA1=A1
      RB3=B3
C
      IF (LPAR(3).LT.3) THEN
        B1=+RBSQC*RB10
        ELSE
        B1=(F1GG*DG*DG*EMEGG+2D0*F1GZ*DG*DZ*EMEGZ+F1ZZ*DZ*DZ*EMEZZ)
     *     *RBSQC
      ENDIF
C
      IF (LPAR(14).EQ.1) THEN
        B1 = B1
     *    +(BXGGC+BXGZC)*EMEGZ*DG*DZ*F1BGZ
     *    -(BXGG5C+BXGZ5C)*EPEGZ*DG*DZ*F3BGZ
     *    +BXGZC*EMEZZ*DZ*DZ*F1BZZ - BXGZ5C*EPEZZ*DZ*DZ*F3BZZ
      ENDIF
      IF (LPAR(16).EQ.1) THEN
        B1 = B1 +
     *     BXZZC*EMEZZ*DG*DG*F1BZZ - BXZZ5C*DG*DG*EPEZZ*F3BZZ
     *    -BXZZC*(EPEZZ*AEZ-EMEZZ*VEZ)*DG*DZ*FBZZZ
     *    +BXZZ5C*(EPEZZ*VEZ-EMEZZ*AEZ)*DG*DZ*FBZZ5
     *    -DG*(DG*FGBWW+DZ*FZBWW)*BXWWC
      ENDIF
C
      IF (LPAR(3).LT.3) THEN
        A3=F3GG*DG*DG*EPEGG+2D0*F3GZ*DG*DZ*EPEGZ+F3ZZ*DZ*DZ*EPEZZ
     *     +RBSQ0*RA30
        ELSE
        A3=(F3GG*DG*DG*EPEGG+2D0*F3GZ*DG*DZ*EPEGZ+F3ZZ*DZ*DZ*EPEZZ)
     *     *(1D0+RBSQ0)
      ENDIF
C
      IF (LPAR(14).EQ.1) THEN
        A3 = A3 +
     *     (BXGG0+BXGZ0)*EPEGZ*DG*DZ*F3BGZ
     *    -(BXGG50+BXGZ50)*EMEGZ*DG*DZ*F1BGZ
     *    +BXGZ0*EPEZZ*DZ*DZ*F3BZZ - BXGZ50*EMEZZ*DZ*DZ*F1BZZ
      ENDIF
      IF (LPAR(16).EQ.1) THEN
        A3 = A3 +
     *     BXZZ0*EPEZZ*DG*DG*F3BZZ - BXZZ50*DG*DG*EMEZZ*F1BZZ
     *    +BXZZ0*(EPEZZ*VEZ-EMEZZ*AEZ)*DG*DZ*FBZZ5
     *    -BXZZ50*(EPEZZ*AEZ-EMEZZ*VEZ)*DG*DZ*FBZZZ
     *    -DG*(DG*FGBWW+DZ*FZBWW)*BXWW0
      ENDIF
      RB1=B1
      RA3=A3
      endif
      
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SOFT BREMSSTRAHLUNG FOR ELECTRON-QUARK SCATTERING
C
      FUNCTION HSBRNC(X,Y,LL,LQF)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSBRNC,HSSPEN
     1       ,BEIEF,BQIQF,BEIQI,BEFQF,BEFQI,BEIQF
     2       ,ZEIEF,ZQIQF,ZEIQI,ZEFQF,ZEFQI,ZEIQF
      DIMENSION CQFL(2)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSISGM/ TCUTQ,TCUTQS
      DATA CQFL /0.6666666666666666D0, -0.3333333333333333D0/
C
      CQF=CQFL(LQF)
      GSP=SP-MEI2-MPRO2
      Q2=GSP*X*Y
      S=X*SP
      U=-(S-Q2)
      DELTA2=DELTA*DELTA
      EEI=EELE
      EEI2=EEI*EEI
      EQI=SP/4D0/EEI*X
      EQI2=EQI*EQI
      EEF=EEI*(1D0-Y)+Q2/4D0/EEI
      EEF2=EEF*EEF
      EQF=EEI+EQI-EEF
      EQF2=EQF*EQF
C
      HSBRNC=0D0
      IF (LPAR(12).EQ.1) THEN
        ZEIEF=DCMPLX(1D0-4D0*EEI*EEF/Q2,0D0)
        BEIEF=-HSSPEN(ZEIEF)
        HSBRNC=HSBRNC
     *          -2D0*DLOG(MEI2/Q2)*(1D0+DLOG(4D0*DELTA2/Q2) )
     *          -DLOG(DELTA2/EEI2)-DLOG(DELTA2/EEF2)
     *          +2D0*DREAL(BEIEF)-2D0*PI*PI/3D0
     *          -DLOG(EEI2/EEF2)*DLOG(EEI2/EEF2)/4D0
     *          +DLOG(Q2/4D0/EEI/EEF)
     *             *(DLOG(4D0*EEI2/MEI2)+DLOG(4D0*EEF2/MEI2)
     *               +DLOG(Q2/4D0/EEI/EEF)                )
      ENDIF
      IF (LPAR(13).EQ.1) THEN
        TCQ2=TCUTQ*TCUTQ
        TCQS2=TCUTQS*TCUTQS
        ZQIQF=DCMPLX(1D0-4D0*EQI*EQF/Q2,0D0)
        BQIQF=-HSSPEN(ZQIQF)
        HSBRNC=HSBRNC+CQF*CQF*(
     *           -DLOG(Q2/MQI2)
     *           +2D0*DREAL(BQIQF)
     *           +9D0/2D0-4D0*PI*PI/3D0
     *           -DLOG(EQI2/EQF2)*DLOG(EQI2/EQF2)/4D0
     *           -DLOG(4D0*EQI*EQF/Q2)*DLOG(4D0*EQI*EQF2/Q2)
     *           -(3D0/2D0+DLOG(DELTA2/EQF2))
     *                 *DLOG((TCQS2*EQF2+MQF2)/Q2)
     *           -(3D0/2D0+DLOG(DELTA2/EQI2))
     *                 *DLOG((TCQ2*EQI2+MQI2)/Q2))
      ENDIF
      IF (LPAR(14).EQ.1) THEN
        ZEIQI=DCMPLX(1D0-4D0*EEI*EQI/S,0D0)
        BEIQI=-HSSPEN(ZEIQI)
        ZEFQF=DCMPLX(1D0-4D0*EEF*EQF/S,0D0)
        BEFQF=-HSSPEN(ZEFQF)
        ZEFQI=DCMPLX(1D0-4D0*EEF*EQI/(-U),0D0)
        BEFQI=-HSSPEN(ZEFQI)
        ZEIQF=DCMPLX(1D0-4D0*EEI*EQF/(-U),0D0)
        BEIQF=-HSSPEN(ZEIQF)
        HSBRNC=HSBRNC+2D0*LL*CQF*(
     *           +2D0*DLOG(4D0*DELTA2/Q2)*DLOG((-U)/S)
     *           +DREAL(-BEIQI-BEFQF+BEFQI+BEIQF)      )
      ENDIF
      RETURN
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        FUNCTION HSFONE(Q2,RM,RN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SKALARES EINSCHLEIFENINTEGRAL, DOPPELTLANG, 'REGULAERER ANTEIL'
C       F(Q2,RM,RN)=B0(Q2,RM,RN)-B0(0D0,RM,RN)  'SUBTRAHIERTES F'
C       Q2=QUADRAT DES DIE SCHLEIFE DURCHLAUFENDEN IMPULSES
C       RM,RN: MASSEN DER TEILCHEN AUF BEIDEN ARMEN
C       D O U B L E P R E C I S I O N
C-----------------------------------------------------------------------
C       19.10.83
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX(kind=8) HSFONE
        DOUBLEPRECISION M,N,PI,A,S,T,B,C,D,Q2,RM,RN,U,V,W
        DATA PI/3.1415926535897932384626438D0/
        M=RM
        N=RN
        IF (M .EQ. N ) GOTO 30
        IF (N .EQ. 0.D0) GOTO 310
        IF (M .EQ. 0.D0) GOTO 300
C---->  ALLGEMEINER FALL
        IF (Q2 .NE. 0D0) GOTO 520
        B=0D0
        A=0D0
        GOTO 560
520     U=M*M+N*N
        V=M*M-N*N
        W=M*N
        S=M+N
        T=M-N
        C=DSQRT(DABS(S*S-Q2))
        D=DSQRT(DABS(T*T-Q2))
        B=1D0+(V/Q2-U/V)*DLOG(N/M)
        IF (2D0*W .LE. DABS(Q2-U)) GOTO 550
        B=B-2D0*C*D/Q2*DATAN(D/C)
        A=0D0
        GOTO 560
550     A=C*D/Q2
        B=B-DSIGN(1D0,Q2-U)*A*DLOG((C+D)*(C+D)/(4D0*W))
        A=PI*A
        IF (Q2 .GE. U) GOTO 560
        A=0D0
560     CONTINUE
        HSFONE=DCMPLX(B,A)
        RETURN
C---->  GLEICHE MASSEN
30      IF (Q2 .NE. 0D0) GOTO 40
        B=0D0
        A=0D0
        GOTO 560
40      U=4D0*M*M
        IF (DABS(Q2).LT.U/1D4) THEN
          B=Q2/6D0/M/M
          A=0D0
          GOTO 560
        ENDIF
        V=DSQRT(DABS(1D0-U/Q2))
        IF ((Q2 .GE. 0D0) .AND. (Q2 .LT. U)) GOTO 50
        B=2D0-V*DLOG((V+1D0)*(V+1D0)/U*DABS(Q2))
        A=PI*V
        IF (Q2 .GE. U) GOTO 560
        A=0D0
        GOTO 560
50      B=2D0-2D0*V*DATAN(1D0/V)
        A=0D0
        GOTO 560
C---->  EINE MASSE NULL
300     M=N
310     IF (Q2 .NE. M*M) GOTO 320
        A=0D0
        B=1D0
        GOTO 560
320     B=1D0
        IF(Q2 .EQ. 0D0) B=0D0
        A=B*(1D0-M*M/(Q2+(1D0-B)))
        B=B-A*DLOG(DABS(1D0-Q2/M/M))
        A=PI*A
        IF (Q2 .GT. M*M) GOTO 560
        A=0D0
        GOTO 560
        END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        FUNCTION HSSFZZ(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       HSSFZZ(Q2): BEITRAEGE ZUR Z0-BOSONSELBSTENERGIE
C                   'REGULAERER' ANTEIL
C       D O U B L E P R E C I S I O N
C-----------------------------------------------------------------------
C       09.11.83
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSFZZ,HSFONE,NYANT,ANTZH
     *          ,FQWW,FQZH,FQEE,FQMYMY,FQTATA
     *          ,FQUU,FQCC,FQTT,FQDD,FQSS,FQBB
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
        DMHZ2=MH2-MZ2
        SMHZ2=MH2+MZ2
        FQWW=HSFONE(Q2,MW,MW)
        FQZH=HSFONE(Q2,MZ,MH)
        FQEE=HSFONE(Q2,ME,ME)
        FQMYMY=HSFONE(Q2,MMY,MMY)
        FQTATA=HSFONE(Q2,MTAU,MTAU)
        FQUU=HSFONE(Q2,MU,MU)
        FQCC=HSFONE(Q2,MC,MC)
        FQTT=HSFONE(Q2,MT,MT)
        FQDD=HSFONE(Q2,MD,MD)
        FQSS=HSFONE(Q2,MS,MS)
        FQBB=HSFONE(Q2,MB,MB)
        DLHW=DLOG(MH/MW)
        DLZW=-DLOG(CW)
        DLZH=DLOG(MZ/MH)
        DLEW2=DLOG(ME2/MW2)
        DLYW2=DLOG(MMY2/MW2)
        DLAW2=DLOG(MTAU2/MW2)
        DLUW2=DLOG(MU2/MW2)
        DLCW2=DLOG(MC2/MW2)
        DLTW2=DLOG(MT2/MW2)
        DLDW2=DLOG(MD2/MW2)
        DLSW2=DLOG(MS2/MW2)
        DLBW2=DLOG(MB2/MW2)
        NYANT=0D0
        IF(Q2 .NE. 0D0) GOTO 10
         ANTZH=0.5D0*SMHZ2+MH2*MZ2/DMHZ2*2D0*DLZH
        GOTO 20
10      ANTZH=FQZH/Q2*DMHZ2*DMHZ2
        IF(Q2 .LT. 0D0) NYANT=5D0/3D0-DLOG(-Q2/MW2)
        IF(Q2 .GT. 0D0) NYANT=DCMPLX(5D0/3D0-DLOG(Q2/MW2),PI)
20      S2C224=1D0/(24D0*SW2*CW2)
        S2C208=3D0*S2C224
        S2_L=8D0*SW2*(SW2-0.5D0)+1D0
        S2_2=SW2/0.28125D0*(SW2-0.75D0)+1D0
        S2_1=SW2/1.125D0*(SW2-1.5D0)+1D0
        HSSFZZ=+ALP1PI*(
     G  +DFLOAT(LPAR(15))/(12D0*SW2)*
     G   (-(Q2/1.5D0+(2D1*MW2+1D1*Q2)*FQWW                     )*CW2
     G    +(3D0*MW2*FQWW+Q2/6D0
     G      -MH2*DLHW-MZ2*DLZW
     G      +0.25D0*(1D1*MZ2-2D0*MH2+Q2)*
     G         (1D0+SMHZ2/DMHZ2*DLZH
     G             -DLZW-DLHW      +FQZH   )
     G      +0.25D0*ANTZH                                      )/CW2
     G    +(Q2/6D0+(2D0*MW2+0.25D0*Q2)*FQWW)*(CW2-SW2)*(CW2-SW2)/CW2)
     L  +S2C224*
     L   (S2_L*((2D0*ME2+Q2)*FQEE-Q2/3D0)    -3D0*ME2*FQEE    )
     L  +S2C224*
     L   (S2_L*((2D0*MMY2+Q2)*FQMYMY-Q2/3D0) -3D0*MMY2*FQMYMY )
     L  +S2C224*
     L   (S2_L*((2D0*MTAU2+Q2)*FQTATA-Q2/3D0)-3D0*MTAU2*FQTATA)
     N  +S2C224*Q2*(NYANT+DLEW2)
     N  +S2C224*Q2*(NYANT+DLYW2)
     N  +S2C224*Q2*(NYANT+DLAW2) )
        HSSFZZ=HSSFZZ+ALP1PI*(
     2  +S2C208*
     2   (S2_2*((2D0*MU2+Q2)*FQUU-Q2/3D0)    -3D0*MU2*FQUU    )
     2  +S2C208*
     2   (S2_2*((2D0*MC2+Q2)*FQCC-Q2/3D0)    -3D0*MC2*FQCC    )
     2  +S2C208*
     2   (S2_2*((2D0*MT2+Q2)*FQTT-Q2/3D0)    -3D0*MT2*FQTT    )
     1  +S2C208*
     1   (S2_1*((2D0*MD2+Q2)*FQDD-Q2/3D0)    -3D0*MD2*FQDD    )
     1  +S2C208*
     1   (S2_1*((2D0*MS2+Q2)*FQSS-Q2/3D0)    -3D0*MS2*FQSS    )
     1  +S2C208*
     1   (S2_1*((2D0*MB2+Q2)*FQBB-Q2/3D0)    -3D0*MB2*FQBB    ) )
        RETURN
        END

        FUNCTION HSSFWW(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       HSSFWW(Q2): BEITRAEGE ZUR W+-BOSONSELBSTENERGIE
C                   NICHT RENORMIERT
C       D O U B L E P R E C I S I O N
C-----------------------------------------------------------------------
C       09.11.83  BER 16.08.85HS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSFWW,HSFONE,FQ0W,FQZW,FQHW,FQUD,FQCS,FQTB
     *          ,FQ0E,FQ0MY,FQ0TA
     *          ,ANT1,ANT2,ANT3,ANT4,ANT5,ANT6,ANT7,ANT8,ANT9
     *          ,ANT71,ANT81
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
        FQ0W=HSFONE(Q2,0D0,MW)
        FQZW=HSFONE(Q2,MZ,MW)
        FQHW=HSFONE(Q2,MH,MW)
        ANT71=1D0-DLOG(MZ2/MW2)/SW2+FQZW
        IF (MH2.EQ.MW2) THEN
          ANT81=1D0+MH2/MW2+FQHW
          ELSE
          ANT81=1D0-DLOG(MH2/MW2)*MH2/(MH2-MW2)+FQHW
        ENDIF
        IF (MD2.EQ.MU2) THEN
          DLMUD=0D0
          ELSE
          DLMUD=1D0-(MU2+MD2)/(MU2-MD2)*DLOG(MU/MD)
        ENDIF
        FQ0E=HSFONE(Q2,0D0,ME)
        FQ0MY=HSFONE(Q2,0D0,MMY)
        FQ0TA=HSFONE(Q2,0D0,MTAU)
        FQUD=HSFONE(Q2,MU,MD)
        FQCS=HSFONE(Q2,MC,MS)
        FQTB=HSFONE(Q2,MT,MB)
        IF(Q2 .NE. 0D0) GOTO 10
        ANT9=MW2
        ANT1=ME2
        ANT2=MMY2
        ANT3=MTAU2
        ANT7=0.5D0*(MW2+MZ2)+MW2*MZ2/(MZ2-MW2)*DLOG(MW2/MZ2)
        IF (MH2.EQ.MW2) THEN
          ANT8=0.5D0*(MW2-MH2)
          ELSE
          ANT8=0.5D0*(MW2+MH2)+MW2*MH2/(MH2-MW2)*DLOG(MW2/MH2)
        ENDIF
        IF (MD2.EQ.MU2) THEN
          ANT4=0.5D0*(MU2-MD2)
          ELSE
          ANT4=0.5D0*(MU2+MD2)+MU2*MD2/(MD2-MU2)*DLOG(MU2/MD2)
        ENDIF
        ANT5=0.5D0*(MC2+MS2)+MC2*MS2/(MS2-MC2)*DLOG(MC2/MS2)
        ANT6=0.5D0*(MT2+MB2)+MT2*MB2/(MB2-MT2)*DLOG(MT2/MB2)
        GOTO 100
10      ANT9=2D0*MW2*MW2*FQ0W/Q2
        ANT1=2D0*ME2*ME2*FQ0E/Q2
        ANT2=2D0*MMY2*MMY2*FQ0MY/Q2
        ANT3=2D0*MTAU2*MTAU2*FQ0TA/Q2
        ANT7=(MZ2-MW2)*(MZ2-MW2)*FQZW/Q2
        ANT8=(MH2-MW2)*(MH2-MW2)*FQHW/Q2
        ANT4=(MD2-MU2)*(MD2-MU2)*FQUD/Q2
        ANT5=(MS2-MC2)*(MS2-MC2)*FQCS/Q2
        ANT6=(MB2-MT2)*(MB2-MT2)*FQTB/Q2
100     HSSFWW=+ALP1PI*(
     G +DFLOAT(LPAR(15))/(12D0*SW2)*
     G  (+(-(7D0*(MZ2+MW2)+1D1*Q2)*ANT71
     G     -4D0*MZ2*DLOG(MZ2/MW2)-Q2/1.5D0+2D0*ANT7         )*CW2
     G   -((4D0*MW2+1D1*Q2)*FQ0W+4D0*MW2+Q2/0.09375D0-ANT9  )*SW2
     G   +(1D0-DLOG(MZ2/MW2)/SW2+FQZW       )*3D0*MW2*SW2*SW2/CW2
     G   +( ANT81                                       )*3D0*MW2
     G   -(MZ2*DLOG(MZ2/MW2)+MH2*DLOG(MH2/MW2))*0.5D0+Q2/3D0
     G   -(MW2+MZ2-0.5D0*Q2)*ANT71*0.5D0+ANT7*0.25D0
     G   -(MH2+MW2-0.5D0*Q2)*ANT81*0.5D0+ANT8*0.25D0             )
     L +1D0/(12D0*SW2)*
     L  ((Q2-ME2*0.5D0)*(HSFONE(Q2,0D0,ME)+1D0)
     L        -Q2/3D0-0.25D0*ANT1)
     L +1D0/(12D0*SW2)*
     L  ((Q2-MMY2*0.5D0)*(HSFONE(Q2,0D0,MMY)+1D0)
     L      -Q2/3D0-0.25D0*ANT2)
     L +1D0/(12D0*SW2)*
     L  ((Q2-MTAU2*0.5D0)*(HSFONE(Q2,0D0,MTAU)+1D0)
     L    -Q2/3D0-0.25D0*ANT3) )
       HSSFWW=HSSFWW+ALP1PI*(
     Q +1D0/(4D0*SW2)*
     Q  (+(Q2-(MU2+MD2)*0.5D0)*
     Q     (+FQUD+DLMUD)
     Q   -Q2/3D0         -0.5D0*ANT4)
     Q +1D0/(4D0*SW2)*
     Q  (+(Q2-(MC2+MS2)*0.5D0)*
     Q     (+FQCS+1D0
     Q      -(MC2+MS2)/(MC2-MS2)*DLOG(MC/MS))
     Q   -Q2/3D0          -0.5D0*ANT5)
     Q +1D0/(4D0*SW2)*
     Q  (+(Q2-(MT2+MB2)*0.5D0)*
     Q     (+FQTB+1D0
     Q      -(MT2+MB2)/(MT2-MB2)*DLOG(MT/MB))
     Q   -Q2/3D0          -0.5D0*ANT6))
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION HSSFGZ(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       HSSFGZ(Q2): BEITRAEGE ZUR Z-GAMMA-MISCHUNGSENERGIE
C                   'REGULAERER ANTEIL'
C       D O U B L E P R E C I S I O N
C-----------------------------------------------------------------------
C       09.11.83  BER. 16.08.85 HS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSFGZ,HSFONE,FQEE,FQMYMY,FQTATA
     *          ,FQUU,FQCC,FQTT,FQDD,FQSS,FQBB
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
        FQEE=HSFONE(Q2,ME,ME)
        FQMYMY=HSFONE(Q2,MMY,MMY)
        FQTATA=HSFONE(Q2,MTAU,MTAU)
        FQUU=HSFONE(Q2,MU,MU)
        FQCC=HSFONE(Q2,MC,MC)
        FQTT=HSFONE(Q2,MT,MT)
        FQDD=HSFONE(Q2,MD,MD)
        FQSS=HSFONE(Q2,MS,MS)
        FQBB=HSFONE(Q2,MB,MB)
        HSSFGZ=+ALP1PI*( +DFLOAT(LPAR(15))/4D0/SW2*
     G  (-(+(3D0*SW*CW+SW/(CW*6D0))*Q2
     G     +(4D0*SW*CW+SW/(CW*0.75D0))*MW2)
     G                                *HSFONE(Q2,MW,MW)-Q2*SW/CW/9D0)
     L  +1D0/(3D0*SW*CW)*
     L   (0.25D0-SW2)*  ((2D0*ME2+Q2)*FQEE-Q2/3D0)
     L  +1D0/(3D0*SW*CW)*
     L   (0.25D0-SW2)*  ((2D0*MMY2+Q2)*FQMYMY-Q2/3D0)
     L  +1D0/(3D0*SW*CW)*
     L   (0.25D0-SW2)*  ((2D0*MTAU2+Q2)*FQTATA-Q2/3D0)
     2  +1D0/(2.25D0*SW*CW)*
     2   (0.375D0-SW2)*  ((2D0*MU2+Q2)*FQUU-Q2/3D0)
     2  +1D0/(2.25D0*SW*CW)*
     2   (0.375D0-SW2)*  ((2D0*MC2+Q2)*FQCC-Q2/3D0)
     2  +1D0/(2.25D0*SW*CW)*
     2   (0.375D0-SW2)*  ((2D0*MT2+Q2)*FQTT-Q2/3D0)
     1  +1D0/(9D0*SW*CW)*
     1   (0.75D0-SW2)*  ((2D0*MD2+Q2)*FQDD-Q2/3D0)
     1  +1D0/(9D0*SW*CW)*
     1   (0.75D0-SW2)*  ((2D0*MS2+Q2)*FQSS-Q2/3D0)
     1  +1D0/9D0/SW/CW*(0.75D0-SW2)*((2D0*MB2+Q2)*FQBB-Q2/3D0))
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  RENORMIERTE SELBSTENERGIEN
C  IN DEN UNTERPROGRAMMEN
C     HSSRGG, HSSRGZ, HSSRZZ, HSSRWW
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FUNCTION HSSRGG(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       HSSRGG(Q2): BEITRAEGE ZUR PHOTONSELBSTENERGIE  (RENORMIERT)
C              Q2 = Q**2 VIERERIMPULSQUADRAT
C-----------------------------------------------------------------------
C       09.11.83
C       14.05.91 HS: BURKHARD'S PARAMETRIZATION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSRGG,HSFONE
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      HSSRGG=DCMPLX(0D0,0D0)
      IF (LPAR(7).EQ.1) THEN
C...HADRONIC CONTRIBUTION FROM EFFECTIVE QUARK LOOPS
        HSSRGG=+ALP1PI*(
     G   -DFLOAT(LPAR(15))*((3D0*Q2+4D0*MW2)*HSFONE(Q2,MW,MW)
     G                                       - Q2/1.5D0)/4D0
     L   +(+(Q2+2D0*ME2  )*HSFONE(Q2,ME  ,ME  )
     L     +(Q2+2D0*MMY2 )*HSFONE(Q2,MMY ,MMY )
     L     +(Q2+2D0*MTAU2)*HSFONE(Q2,MTAU,MTAU) - Q2)/3D0
     2   +(+(Q2+2D0*MU2  )*HSFONE(Q2,MU  ,MU  )
     2     +(Q2+2D0*MC2  )*HSFONE(Q2,MC  ,MC  )
     2     +(Q2+2D0*MT2  )*HSFONE(Q2,MT  ,MT  ) - Q2)/2.25D0
     1   +(+(Q2+2D0*MD2  )*HSFONE(Q2,MD  ,MD  )
     1     +(Q2+2D0*MS2  )*HSFONE(Q2,MS  ,MS  )
     1     +(Q2+2D0*MB2  )*HSFONE(Q2,MB  ,MB  ) - Q2)/9D0        )
      ELSEIF(LPAR(7).EQ.2) THEN
C...HADRONIC CONTRIBUTION FROM BURKHARDT'S PARAMETRIZATION
        HSSRGG=+ALP1PI*(
     G   -DFLOAT(LPAR(15))*((3D0*Q2+4D0*MW2)*HSFONE(Q2,MW,MW)
     G                                       - Q2/1.5D0)/4D0
     L   +(+(Q2+2D0*ME2  )*HSFONE(Q2,ME  ,ME  )
     L     +(Q2+2D0*MMY2 )*HSFONE(Q2,MMY ,MMY )
     L     +(Q2+2D0*MTAU2)*HSFONE(Q2,MTAU,MTAU) - Q2)/3D0
     2   +(+(Q2+2D0*MT2  )*HSFONE(Q2,MT  ,MT  ) - Q2/3D0)/2.25D0     )
     H   - DCMPLX(HSHADQ(Q2)*Q2,0D0)
      ELSEIF(LPAR(7).EQ.0) THEN
       HSSRGG=DCMPLX(0D0,0D0)
      ELSE
      ENDIF
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      FUNCTION HSHADQ(S)
C  HADRONIC IRREDUCIBLE QQ SELF-ENERGY: TRANSVERSE
* THIS IS A SLIGHTLY MODIFIED VERSION OF BURKHARDTS CODE. DB. TR.,01/91
C     parametrize the real part of the photon self energy function
C     by  a + b ln(1+C*:S:) , as in my 1981 TASSO note but using
C     updated values, extended using RQCD up to 100 TeV
C     for details see:
C     H.Burkhardt, F.Jegerlehner, G.Penso and C.Verzegnassi
C     in CERN Yellow Report on "Polarization at LEP" 1988
C     H.BURKHARDT, CERN/ALEPH, AUGUST 1988
C     negative values mean t - channel (spacelike)
C     positive values mean s - channel (timelike )
C     in the space like values around 1 GeV are typical for luminosity
C     the values at 92 GeV ( Z mass ) give the light quark contribution
C     to delta r
C     take care of the sign of REPI when using this in different
C     programs
C     Here REPI was chosen to
C     be positive (so that it corresponds directly to delta alpha)
C     often its assumed to be negative.
C
C     the imaginary part is proportional to R (had / mu cross section)
C     and is therefore 0 below threshold ( in all the spacelike region)
C     note also that alpha_s usually has been derived from the measured
C     values of R.
C     Changing alpha_s to values incompatible with current data
C     would imply to be also inconsistent with RE,IM PI
C     defined here
C
C     H.BURKHARDT
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DATA A1,B1,C1/   0.0   ,  0.00835,  1.0  /
      DATA A2,B2,C2/   0.0   ,  0.00238,  3.927 /
      DATA A3,B3,C3/ 0.00165 ,  0.00300,  1.0  /
      DATA A4,B4,C4/ 0.00221 ,  0.00293,  1.0  /
C
      DATA PI/3.141592653589793D0/,ALFAIN/137.0359895D0/,INIT/0/
C
      IF(INIT.EQ.0) THEN
        INIT=1
        ALFA=1./ALFAIN
        ALFAPI=1./PI/ALFAIN
      ENDIF
      T=ABS(S)
      IF(T.LT.0.3**2) THEN
        REPIAA=A1+B1*LOG(1.+C1*T)
      ELSEIF(T.LT.3.**2) THEN
        REPIAA=A2+B2*LOG(1.+C2*T)
      ELSEIF(T.LT.100.**2) THEN
        REPIAA=A3+B3*LOG(1.+C3*T)
      ELSE
        REPIAA=A4+B4*LOG(1.+C4*T)
      ENDIF
C     as imaginary part take -i alfa/3 Rexp
      HSHADQ=REPIAA
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      FUNCTION HSDSGQ(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       DIFFERENCE OF THE RENORMALIZED PHOTON SELF ENERGY DIVIDED BY Q2
C       OF THE PERTURBATIVE RESULT WITH EFFECTIVE QUARK MASSES AND
C       THE PARAMETRIZATION OF BURKHARD
C       14.05.91 HS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) SGEFFQ,HSFONE
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)

      IF (LPAR(7).LT.2) THEN
C...NO CORRECTION
        HSDSGQ=0D0
      ELSE
        SGEFFQ=+ALP1PI*(
     2   +(+(Q2+2D0*MU2  )*HSFONE(Q2,MU  ,MU  )
     2     +(Q2+2D0*MC2  )*HSFONE(Q2,MC  ,MC  ) - 2D0*Q2/3D0)/2.25D0
     1   +(+(Q2+2D0*MD2  )*HSFONE(Q2,MD  ,MD  )
     1     +(Q2+2D0*MS2  )*HSFONE(Q2,MS  ,MS  )
     1     +(Q2+2D0*MB2  )*HSFONE(Q2,MB  ,MB  ) - Q2        )/9D0   )
        HSDSGQ=-HSHADQ(Q2)-DREAL(SGEFFQ)/Q2
      ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      FUNCTION HSSRGZ(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       HSSRGZ(Q2): BEITRAEGE ZUR GAMMA-Z-MISCHUNG (RENORMIERT)
C              Q2 = Q**2 VIEREIMPULSQUADRAT
C-----------------------------------------------------------------------
C       21.10.83
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSRGZ,HSSFGZ,HSSFZZ,HSSFWW,CAFINW,CAFINZ
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
        SQ=ALP4PI*(DLOG(MU2/MD2)+DLOG(MC2/MS2)+DLOG(MT2/MB2))
        SP=ALP4PI*DLOG(MT2/MB2)*(MT2-MB2)/4./SW2/MW2
        CAFINZ=HSSFZZ(MZ2)
        CAFINW=HSSFWW(MW2)
        RDMZ2=DREAL(CAFINZ)
        RDMW2=DREAL(CAFINW)
        MRENK1=CW/SW*(RDMZ2/MZ2 - RDMW2/MW2) +SQ/6./SW/CW + SP*CW/SW
        HSSRGZ=-(HSSFGZ(Q2)+ Q2 *MRENK1)
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      FUNCTION HSSRZZ(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       HSSRZZ(Q2): BEITRAEGE ZUR Z-BOSONSELBSTENERGIE (RENORMIERT)
C              Q2 = Q**2 VIEREIMPULSQUADRAT
C-----------------------------------------------------------------------
C       21.10.83
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) CJFINZ,CAFINZ,CAFINW
     *          ,HSSRZZ,HSSFZZ,HSSFWW
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
        SQ=ALP4PI*(DLOG(MU2/MD2)+DLOG(MC2/MS2)+DLOG(MT2/MB2))
        SP=ALP4PI*DLOG(MT2/MB2)*(MT2-MB2)/4./SW2/MW2
        CAFINZ=HSSFZZ(MZ2)
        CAFINW=HSSFWW(MW2)
        RDMZ2=DREAL(CAFINZ)
        RDMW2=DREAL(CAFINW)
        MRENK0=(CW2/SW2-1D0)*(RDMZ2/MZ2 - RDMW2/MW2) + ALP2PI/3D0
     *         + SQ/3./SW2 + SP*(CW2/SW2-1D0)
        CJFINZ=DCMPLX(RDMZ2,0D0)
        HSSRZZ=(HSSFZZ(Q2)-CJFINZ+(Q2-MZ2)*MRENK0)
        IF (LPAR(7).EQ.2) THEN
          DDALPP=HSDSGQ(MZ2)*(Q2-MZ2)
          HSSRZZ=HSSRZZ+DCMPLX(DDALPP,0D0)
        ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      FUNCTION HSSRWW(Q2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       HSSRWW(Q2): BEITRAEGE ZUR Z-BOSONSELBSTENERGIE (RENORMIERT)
C              Q2 = Q**2 VIEREIMPULSQUADRAT
C-----------------------------------------------------------------------
C       21.10.83
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSRWW,CJFINW,CAFINZ,CAFINW
     *          ,HSSFZZ,HSSFWW
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
        SQ=ALP4PI*(DLOG(MU2/MD2)+DLOG(MC2/MS2)+DLOG(MT2/MB2))
        SP=ALP4PI*DLOG(MT2/MB2)*(MT2-MB2)/4./SW2/MW2
        CAFINZ=HSSFZZ(MZ2)
        CAFINW=HSSFWW(MW2)
        RDMW2=DREAL(CAFINW)
        RDMZ2=DREAL(CAFINZ)
        MRENKW=CW2/SW2*(RDMZ2/MZ2 - RDMW2/MW2) + ALP2PI/3D0
     *         + SQ/3./SW2 + SP*CW2/SW2
        CJFINW=DCMPLX(RDMW2,0D0)
        HSSRWW=(HSSFWW(Q2)-CJFINW+(Q2-MW2)*MRENKW)
        IF (LPAR(7).EQ.2) THEN
          DDALPP=HSDSGQ(MZ2)*(Q2-MW2)
          HSSRWW=HSSRWW+DCMPLX(DDALPP,0D0)
        ENDIF
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION HSCLN(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       KOMPLEXER NATUERLICHER LOGARITHMUS MIT SCHNITT= NEG. REELLE ACHS
C-----------------------------------------------------------------------
C       20.07.83 UE 20.08.85
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX(kind=8) HSCLN,Z
        DOUBLE PRECISION X,Y,R,PHI,HSSIGN
        EXTERNAL HSSIGN
C---->  HSSIGN(A,B)=|A|*SGN(B), NICHT DIE INTRINSIC FUNCTION VERWENDEN!
C       PI=3.1415926535897932384
        X=REAL(Z)
        Y=IMAG(Z)
        R=ABS(Z)
        
        if (x.eq.0d0) then
          PHI=HSSIGN(1.5707963267948966192D0,Y)
          HSCLN=DCMPLX(DLOG(R),PHI)
          return
        elseif (x.lt.0d0) then
          if (y.eq.0d0) then
            PHI=3.1415926535897932384D0
          else
            PHI=DATAN(Y/X)+HSSIGN(3.1415926535897932384D0,Y)
          endif
          HSCLN=DCMPLX(DLOG(R),PHI)
          return
        elseif (x.gt.0d0) then
          PHI=DATAN(Y/X)
          HSCLN=DCMPLX(DLOG(R),PHI)
          return
        else
          write(*,*) ' Error in HSCLN, execution stopped '
          stop
        endif
        
c        IF (X)  10,20,30
c20      PHI=HSSIGN(1.5707963267948966192D0,Y)
c        GOTO 100
c10      IF (Y .EQ. 0D0) GOTO 40
c        PHI=DATAN(Y/X)+HSSIGN(3.1415926535897932384D0,Y)
c        GOTO 100
c40      PHI=3.1415926535897932384D0
c        GOTO 100
c30      PHI=DATAN(Y/X)
c100     HSCLN=DCMPLX(DLOG(R),PHI)
c        RETURN
        END

        FUNCTION HSSIGN(X,Y)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       HSSIGN(X,Y)=|X|*SGN(Y) IM GEGENSATZ ZUR INTRINSIC FUNCTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        DOUBLE PRECISION HSSIGN,X,Y

        if (y.lt.0d0) then
          HSSIGN=-DABS(X)
        elseif (y.eq.0d0) then
          HSSIGN=0D0
        else
          HSSIGN=DABS(X)
        endif
        
c        IF (Y)  10,20,30
c10      HSSIGN=-DABS(X)
c        RETURN
c20      HSSIGN=0D0
c        RETURN
c30      HSSIGN=DABS(X)
        RETURN
        END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        FUNCTION HSSPEN(Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       SPENCE-FUNKTION KOMPLEX, FREI NACH HOLLIK
C-----------------------------------------------------------------------
C       20.07.83    UE 20.08.85
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        COMPLEX(kind=8) HSSPEN,HSCLN,W,SUM,Z,U
        DOUBLE PRECISION RZ,AZ,A1
        DOUBLE PRECISION B(9)
        DATA B /
     1   0.1666666666666666666666666667D0,
     2  -0.0333333333333333333333333333D0,
     3   0.0238095238095238095238095238D0,
     4  -0.0333333333333333333333333333D0,
     5   0.0757575757575757575757575758D0,
     6  -0.2531135531135531135531135531D0,
     7   1.1666666666666666666666666667D0,
     8  -7.09215686274509804D0         ,
     9  54.97117794486215539D0         /
C     BEACHTE:                 B(N)=B2N
C     B(1)=1./6.
C     B(2)=-1./30.
C     B(3)=1./42.
C     B(4)=-1./30.
C     B(5)=5./66.
C     B(6)=-691./2730.
C     B(7)=7./6.
C     B(8)=-3617./510.
C     B(9)=43867./798.
C     B(10)=-174611./330.
C     B(11)=854513./138.
C     PI=3.1415926535897932384
C     PI*PI/6.=1.6449..., PI*PI/3=3.28986...
C
      RZ=DREAL(Z)
      AZ=ABS(Z)
      A1=ABS(1D0-Z)
CH
CH FOR VS FORTRAN:
      IF (AZ.LT.1D-15) THEN
         HSSPEN=DCMPLX(0D0,0D0)
         RETURN
      ENDIF
CH
      IF((RZ .EQ. 1D0) .AND. (DIMAG(Z) .EQ. 0D0)) GOTO 40
      IF(RZ.GT.5D-1) GOTO 20
      IF(AZ.GT.1D0) GOTO 10
      W=-HSCLN(1D0-Z)
      SUM=W-0.25D0*W*W
      U=W
      DO 1 K=1,9
      U=U*W*W/FLOAT(2*K*(2*K+1))
      SUM=SUM+U*B(K)
 1    CONTINUE
      HSSPEN=SUM
      RETURN
10    W=-HSCLN(1D0-1D0/Z)
      SUM=W-0.25*W*W
      U=W
      DO 11 K=1,9
      U=U*W*W/DFLOAT(2*K*(2*K+1))
      SUM=SUM+U*B(K)
11    CONTINUE
      HSSPEN=-SUM-1.64493406684822643D0-.5*HSCLN(-Z)**2
      RETURN
20    IF(A1.GT.1D0) GOTO 30
      W=-HSCLN(Z)
      SUM=W-0.25*W*W
      U=W
      DO 21 K=1,9
      U=U*W*W/FLOAT(2*K*(2*K+1))
      SUM=SUM+U*B(K)
21    CONTINUE
      HSSPEN=-SUM+1.64493406684822643D0-HSCLN(Z)*HSCLN(1.-Z)
      RETURN
30    W=HSCLN(1D0-1./Z)
      SUM=W-0.25*W*W
      U=W
      DO 31 K=1,9
      U=U*W*W/FLOAT(2*K*(2*K+1))
      SUM=SUM+U*B(K)
31    CONTINUE
      HSSPEN=SUM+3.28986813369645287D0+.5*HSCLN(Z-1.)**2
     *      -HSCLN(Z)*HSCLN(1.-Z)
      CONTINUE
      RETURN
40    HSSPEN=1.64493406684822643D0
      RETURN
      END
C=======================================================================
      FUNCTION HSCLM1(Q2,MF2)
C       PHOTONIC VERTEX CORRECTION
C       IR-FINITE PART
C                                 OHNE LOG**2
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSCLM1
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM

      if (q2.lt.0d0) then
        AMB1 = DLOG(-Q2/MF2) + 4d0*(PI*PI/12d0-1D0)
        HSCLM1=DCMPLX(AMB1)
      elseif (q2.eq.0d0) then
        HSCLM1 = (0D0,0D0)
      else
        AMB1 = DLOG(Q2/MF2) + 4d0*(PI*PI/3d0-1D0)
        HSCLM1=DCMPLX(AMB1)
      endif
      return
      
c      IF (Q2) 10,20,30
c10    AMB1 = DLOG(-Q2/MF2)
c     1         +4.*(PI*PI/12. - 1D0)
c      HSCLM1=DCMPLX(AMB1)
c      GOTO 40
c20    HSCLM1 = (0D0,0D0)
c      GOTO 40
c30    AMB1 = DLOG(Q2/MF2)
c     1         +4.*(PI*PI/3. - 1D0)
c      HSCLM1=DCMPLX(AMB1)
c40    RETURN

      END
C=======================================================================
      FUNCTION HSCLM2(Q2,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSCLM2,W,HSCLN,HSSPEN,CM2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      W=CM2/DCMPLX(Q2,0D0)
      IF(REAL(W).GT.0) GOTO 100
      HSCLM2=.5-2.*(2.+W)-(2.*W+3)*HSCLN(-W)+
     1        2.*(1+W)**2*(HSSPEN(1+1/W)-PI*PI/6)
      GOTO 200
100   HSCLM2=.5-2.*(2.+W)-(2.*W+3)*HSCLN(W)+
     1        2.*(1+W)**2*(HSCLN(W)*HSCLN((W+1.)/W)-HSSPEN(-1/W))
     2        -DCMPLX(0D0,1D0)*PI*(3.+2.*W-2.*(1+W)**2*HSCLN((1.+W)/W))
200   CONTINUE
      RETURN
      END
C=======================================================================
      FUNCTION HSCLM3(Q2,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSCLM3,CM2,W,HSCLN,CHI,CLH
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      W=CM2/DCMPLX(Q2,0D0)
      
      IF(REAL(W).GT.0) GOTO 100
      CHI=SQRT(1d0-4d0*W)
      CLH=HSCLN((1+CHI)/(CHI-1))
      HSCLM3=5d0/6d0-2d0/3d0*W+(2d0*W+1d0)/3d0*CHI*CLH
     1        +2d0/3d0*W*(W+2d0)*CLH*CLH
      GOTO 400
      
100   if (q2.lt.4d0*real(cm2)) then
        CHI=SQRT(4d0*W-1d0)
        CLH=ATAN(1d0/DREAL(CHI))
        HSCLM3=5d0/6d0-2d0/3d0*W+2d0/3d0*(2d0*W+1d0)*CHI*CLH
     1        -8d0/3d0*W*(W+2d0)*CLH*CLH
      else
        CHI=SQRT(1d0-4d0*W)
        CLH=HSCLN( (1d0+CHI)/(1d0-CHI) )
        HSCLM3=5d0/6d0-2d0*W/3d0+(2d0*W+1d0)/3d0*CHI*CLH
     1        +2d0/3d0*W*(W+2d0)*(CLH*CLH - PI*PI)
     2        -DCMPLX(0D0,1D0)*PI
     3         *( (2d0*W+1d0)/3d0*CHI+2d0/3d0*W*(W+2d0)*CLH)
      endif
      
c100   IF (Q2 - 4.*REAL(CM2)) 200,300,300
c200   CHI=SQRT(4.*W-1)
c      CLH=ATAN(1./DREAL(CHI))
c      HSCLM3=5./6.-2./3.*W+2./3.*(2.*W+1)*CHI*CLH
c     1        -8./3.*W*(W+2.)*CLH*CLH
c      GOTO 400
c300   CHI=SQRT(1.-4.*W)
c      CLH=HSCLN( (1.+CHI)/(1.-CHI) )
c      HSCLM3=5./6.-2.*W/3.+(2.*W+1.)/3.*CHI*CLH
c     1       +2./3.*W*(W+2.)*(CLH*CLH - PI*PI)
c     2       -DCMPLX(0D0,1D0)*PI*( (2.*W+1.)/3.*CHI+2./3.*W*(W+2.)*CLH)

400   CONTINUE
      RETURN
      END
C=======================================================================
      FUNCTION HSCLM4(Q2,CM1,CM2,MF2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSCLM4,CM1,CM2
     1          ,W1,W2,X1,X2,C12,HSCLN,HSSPEN
      RM1=SQRT(DREAL(CM1))
      RM2=SQRT(DREAL(CM2))
      SM2=(RM1+RM2)*(RM1+RM2)
      DM2=(RM1-RM2)*(RM1-RM2)
      W1=CM1/DCMPLX(Q2,0D0)
      W2=CM2/DCMPLX(Q2,0D0)
      IF (REAL(CM1).EQ.0) GOTO 401
      IF (REAL(CM2).EQ.0) GOTO 402
      IF ( (Q2.GT.DM2).AND.(Q2.LT.SM2) ) GOTO 100
      X1=(1.-W1+W2)/2.+SQRT((1.-W1+W2)*(1.-W1+W2)-4.*W2)/2.
      X2=(1.-W1+W2)/2.-SQRT((1.-W1+W2)*(1.-W1+W2)-4.*W2)/2.
      GOTO 200
100   X1=(1.-W1+W2)/2.+DCMPLX(0D0,1D0)/2.
     1                 *SQRT(4.*W2-(1.-W1+W2)*(1.-W1+W2) )
      X2=(1.-W1+W2)/2.-DCMPLX(0D0,1D0)/2.
     1                 *SQRT(4.*W2-(1.-W1+W2)*(1.-W1+W2) )
200   C12=HSCLN(CM1/CM2)/2.
      HSCLM4=1./6. + (W1+W2)/(W1-W2)*C12 - (W1-W2)/3.*C12
     1       +(W1+W2+1.)/3.*(C12-1.)
     2       +(W1+W2+1.)/3.*(X1*HSCLN(X1/(X1-1.))
     3                       + X2*HSCLN(-X2/(1.-X2)) )
     4       -2./3.*(W1+W2+W1*W2)*HSCLN(X1/(X1-1.))*HSCLN(-X2/(1.-X2))
      GOTO 500
401   W1=W2
      X1=HSCLN(CM2/MF2)
      X2=HSCLN((CM2-Q2)/MF2)
      GOTO 404
402   X1=HSCLN(CM1/MF2)
      X2=HSCLN((CM1-Q2)/MF2)
404   HSCLM4=.5 - W1/3. + (1.-W1*W1)/3.*HSCLN((W1-1.)/W1)
     1       +2./3.*X1 + W1/3.*(X2*X2-X1*X1+2.*HSSPEN(1./(1.-W1)) )
500   RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION HSBXI5(G,H,CM)
      COMPLEX(kind=8) HSBXI5,MY,A2,C,X1,X2,Y1,Y2,HSSPEN,HSCLN,CM
      DOUBLE PRECISION G,H
       MY=CM/G
       A2=DCMPLX(-H/G,0D0)
       C=SQRT(1D0-4D0*MY*(1D0-MY/A2))
       X1=(1D0+C)/2D0
       X2=(1D0-C)/2D0
       C=SQRT(1D0-4D0*MY)
       Y1=(1D0+C)/2D0
       Y2=(1D0-C)/2D0
       C=(A2-.5-2*A2*MY+MY*MY/A2+MY*MY)/(A2-1)/(X2-X1)
     *           *(HSSPEN(X2/(X2-Y2))+HSSPEN(X2/(X2-Y1))
     *            -HSSPEN(X1/(X1-Y2))-HSSPEN(X1/(X1-Y1)))
     *   -(Y2-Y1)*HSCLN(-Y1/Y2)/2-HSCLN(A2/MY)/2
     *   +(A2-MY-.5)/(A2-1)
     *   *(HSSPEN(A2/MY)+HSCLN(A2/MY)*HSCLN(1-A2/MY)+HSCLN(-Y1/Y2)**2)
       HSBXI5=C/(A2-1)
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION HSBXI0(G,H,CM)
      COMPLEX(kind=8) HSBXI0,HSBXI5,MY,A2,C,X1,X2,Y1,Y2,HSSPEN,HSCLN,CM
      DOUBLE PRECISION G,H
       MY=CM/G
       A2=DCMPLX(-H/G,0D0)
       C=SQRT(1D0-4D0*MY*(1D0-MY/A2))
       X1=(1D0+C)/2D0
       X2=(1D0-C)/2D0
       C=SQRT(1D0-4D0*MY)
       Y1=(1D0+C)/2D0
       Y2=(1D0-C)/2D0
       HSBXI0=2/(X2-X1)*(HSSPEN(X2/(X2-Y2))+HSSPEN(X2/(X2-Y1))
     *              -HSSPEN(X1/(X1-Y2))-HSSPEN(X1/(X1-Y1)))
     *   +2*(HSCLN(-Y1/Y2))**2
     *   +HSBXI5(G,H,CM)
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION HSBXCV (RG,RH,CM)
      COMPLEX(kind=8) HSBXCV,HSBXCA,CM,G,H,HSSPEN,HSCLN,IEPS,T1,T2,T3
      DOUBLE PRECISION RG,RH
      H=DCMPLX(RH,0D0)
      G=DCMPLX(RG,0D0)
      IEPS=DCMPLX(0D0,1D-8)
      T1=HSBXCA(RG,RH,CM)
      T2=-2D0*HSCLN(-G/H-IEPS)*HSCLN(CM/(CM-G))
      T3=2D0*HSSPEN((CM+H)/H)
      HSBXCV=T1+T2+T3
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION HSBXCA (RG,RH,CM)
      COMPLEX(kind=8) HSBXCA,CM,G,H,HSSPEN,HSCLN,IEPS
      DOUBLE PRECISION RG,RH
      H=DCMPLX(RH,0D0)
      G=DCMPLX(RG,0D0)
      IEPS=DCMPLX(0D0,1D-8)
      HSBXCA =(G-CM)/(G+H)*(HSCLN(H/(G-CM))-CM/G*HSCLN(CM/(CM-G))+
     *     (G+2D0*H+CM)/(G+H)*
     *(HSSPEN(G/CM)-HSSPEN(-H/CM)+HSCLN(-H/CM)*
     *    HSCLN((CM-G)/(CM+H))))
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     HSBCGA: FORM FACTOR FOR GAMMA-GAMMA-BOX DIAGRAMS
C     IR-FINITE PART ONLY
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION HSBCGA(T,S)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSBCGA,HSCLN,EPS
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      EPS=DCMPLX(0D0,1D-9)
      EPS=HSCLN(S/(T+EPS) )
      HSBCGA=T/2./(T+S)*EPS - T*(T+2.*S)/4./(T+S)/(T+S)*(EPS*EPS+PI*PI)
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     VERTEX CORRECTIONS
C     IR-FINITE PARTS, DOUBLE-LOG CONSTRIBUTIONS SUBTRACTED
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION HSFHFB(T,IVA,LF,LB,MF2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSFHFB,HSCLM1,HSCLM2,HSCLM3,CMW2,CMZ2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSCBMS/ CMW2,CMZ2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      IF ((LPAR(11).EQ.1).AND.
     I    (LPAR(12).EQ.1).OR.(LPAR(13).EQ.1)) THEN
C---PHOTON-EXCHANGE
        HSFHFB=VAFI(IVA,LF,LB)*VAFI(1,LF,1)**2 * ALP4PI*HSCLM1(T,MF2)
        ELSE
        HSFHFB=DCMPLX(0D0,0D0)
      ENDIF
      IF (LPAR(15).EQ.0) RETURN
C---Z-EXCHANGE
      IF (IVA.EQ.1) THEN
        IVAS=2
        ELSE
        IVAS=1
      ENDIF
      FCZ=VAFI(IVA,LF,LB)*(VAFI(1,LF,2)**2+VAFI(2,LF,2)**2)
     *   +VAFI(IVAS,LF,LB)*2D0*VAFI(1,LF,2)*VAFI(2,LF,2)
      HSFHFB=HSFHFB + ALP4PI*FCZ*HSCLM2(T,CMZ2)
C---W-EXCHANGE
C---W-EXCHANGE
      IF (LF.EQ.1.AND.LB.EQ.1) GOTO 10
      LFS=3
      IF (LF.EQ.3) LFS=2
      FCW=(VAFI(1,LFS,LB)+VAFI(2,LFS,LB))/4D0/SW2
      IF (LF.EQ.1.AND.LB.EQ.2) FCW=1D0/8D0/CW/SW/SW2
      HSFHFB=HSFHFB + ALP4PI*FCW*HSCLM2(T,CMW2)
C---TRIPLE GAUGE BOSON VERTEX
   10 CONTINUE
      IF (LB.EQ.1) THEN
        FC3B=3D0/4D0/SW2
        ELSE
        FC3B=-3D0/4D0/SW2*CW/SW
      ENDIF
      IF (MOD(LF,2).EQ.0) FC3B=-FC3B
      HSFHFB=HSFHFB + ALP4PI*FC3B*HSCLM3(T,CMW2)
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND FOR NONRADIATIVE CONTRIBUTION IN CHARGED CURRENT
C   (ARGUMENTS AS REQUIRED BY INTEGRATION ROUTINE IN INIT22;
C    HSCCG1 AND HSCCG2 DIFFER ONLY IN THE LIST OF ARGUMENTS)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCCG1(NDIMEN,ARGUM)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSINTL/ XL,XU
      DIMENSION ARGUM(NDIMEN)
C
      DX=XU-XL
      X=XL+ARGUM(1)*DX
      GS=SP-MEI2-MPRO2
      YMAXX=X*(1D0-4D0*MEI2*MPRO2/GS/GS)/(X*(1D0+X*MPRO2/GS)+MEI2/GS)
      Q2L=Q2MIN
      Q2U=X*GS
      IF(ICUT.EQ.2) THEN
        QQ2MIN=X*(WMIN**2-MPRO2)/(1D0 - X)
        Q2L=MAX(Q2MIN,QQ2MIN)
        Q2U=X*GS
      ELSEIF(ICUT.GE.3) THEN
        QQ2MIN=X*(WMIN**2-MPRO2)/(1D0 - X)
        QQQ2MN=X*YMIN*GS
C---CUT ON ELECTRON SCATTERING ANGLE (MASSES NEGLECTED)
        YMINTH=1D0/(1D0+X*CTHCON)
        QT2MIN=X*YMINTH*GS
C---CUT ON ELECTRON TRANSVERSE MOMENTUM (MASSES NEGLECTED)
        QP2MIN=X*GS/2D0*(1D0-DSQRT(1D0-PTXM0/X))
        QP2MAX=X*GS/2D0*(1D0+DSQRT(1D0-PTXM0/X))
        YP2MAX=QP2MAX/GS/X
        Q2L=MAX(Q2MIN,QQ2MIN,QQQ2MN,QT2MIN,QP2MIN)
        Q2U=X*MIN(YMAX,YMAXX,YP2MAX)*GS
        Q2U=MIN(Q2U,Q2MAX)
      ENDIF
      DQ2=DMAX1(Q2U-Q2L,0D0)
      Q2=Q2L+DQ2*ARGUM(2)
      IF(IPRINT.GT.20) WRITE(LUNTES,'(A,2D15.6)')
     &                      ' HSCCG1: X, Q2',X,Q2
      HSCCG1=HSCC22(X,Q2)*DX*DQ2
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND FOR NONRADIATIVE CONTRIBUTION IN CHARGED CURRENT
C   (ARGUMENTS AS REQUIRED BY ESTMAX FOR ESTIMATION OF LOCAL MAXIMA;
C    BESIDES THE ARGUMENT LIST HSCCG1 = HSCCG2 )
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCCG2(XARG)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      DIMENSION XARG(2)
C
      DX=XMAX-XMIN
      X=XMIN+XARG(1)*DX
      Z=XARG(2)
      GS=SP-MEI2-MPRO2
      YMAXX=X*(1D0-4D0*MEI2*MPRO2/GS/GS)/(X*(1D0+X*MPRO2/GS)+MEI2/GS)
      Q2L=Q2MIN
      Q2U=X*GS
      IF(ICUT.EQ.2) THEN
        QQ2MIN=X*(WMIN**2-MPRO2)/(1D0 - X)
        Q2L=MAX(Q2MIN,QQ2MIN)
        Q2U=X*GS
      ELSEIF(ICUT.GE.3) THEN
        QQ2MIN=X*(WMIN**2-MPRO2)/(1D0 - X)
        QQQ2MN=X*YMIN*GS
C---CUT ON ELECTRON SCATTERING ANGLE (MASSES NEGLECTED)
        YMINTH=1D0/(1D0+X*CTHCON)
        QT2MIN=X*YMINTH*GS
C---CUT ON ELECTRON TRANSVERSE MOMENTUM (MASSES NEGLECTED)
        QP2MIN=X*GS/2D0*(1D0-DSQRT(1D0-PTXM0/X))
        QP2MAX=X*GS/2D0*(1D0+DSQRT(1D0-PTXM0/X))
        YP2MAX=QP2MAX/GS/X
        Q2L=MAX(Q2MIN,QQ2MIN,QQQ2MN,QT2MIN,QP2MIN)
        Q2U=X*MIN(YMAX,YMAXX,YP2MAX)*GS
        Q2U=MIN(Q2U,Q2MAX)
      ENDIF
      DQ2=DMAX1(Q2U-Q2L,0D0)
      Q2=Q2L+DQ2*Z
      IF(IPRINT.GT.20) WRITE(LUNTES,'(A,3D15.6)')
     &                      ' HSCCG2: X, Z, Q2',X,Z,Q2
      HSCCG2=HSCC22(X,Q2)*DQ2*DX
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   D2SIG/DX*DQ2 WITH COMPLETE 1-LOOP SOFT AND VIRTUAL CORRECTIONS
C   FROM H. SPIESBERGER
C
C        NOTE: OUTPUT FROM H.S.  DSIG / DX*DY
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCC22(X,Q2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSKPXY/ XX,Y
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      DATA NEVERR /0/
      LOGICAL LERR
      DATA LERR /.TRUE./
C
      GSP=SP-MEI2-MPRO2
      XX=X
      IF(IPRINT.GT.20)
     &  WRITE(LUNTES,'(A/3(1PD13.5),2(0PF8.3),2I3)')
     *         ' HSCC22: SP, X, Q2, POLARI, HPOLAR, LLEPT, LQUA',
     *                   SP, X, Q2, POLARI, HPOLAR, LLEPT, LQUA
      Y=Q2/X/GSP
      HSCC22=HSSGCC(X,Y,LLEPT,POLARI,HPOLAR,LQUA)/X/SP
C
      IF(HSCC22.LE.0D0) THEN
        HSCC22=0D0
        NEVERR=NEVERR+1
        IF (NEVERR.LT.20) THEN
         WRITE(LUNTES,'(A,/,4(1PD13.5),2I5,2(0PF8.3)/A/2(6(1PD13.5)/))')
     +     ' HSCC22: X, Y, Q2, HSCC22, LLEPT, LQUA, POLARI, HPOLAR',
     +               X, Y, Q2, HSCC22, LLEPT, LQUA, POLARI, HPOLAR, 
     +     ' HSPDFQ: QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT',
     +               QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
         WRITE(LUNTES,'(2A,/,2(6(1PD13.5)/))')
     +     ' HSDPDF: DQU,DQBU,DQD,DQBD,DQS,DQBS,',
     +                                  'DQC,DQBC,DQB,DQBB,DQT,DQBT',
     +               DQU,DQBU,DQD,DQBD,DQS,DQBS,
     +                                   DQC,DQBC,DQB,DQBB,DQT,DQBT
        ELSEIF (LERR) THEN
         LERR=.FALSE.
         WRITE(LUNTES,'(A,I3,A)')
     &     ' ERROR HSCC22 < 0 HAS OCCURED ',NEVERR,
     &     ' TIMES, NO FURTHER WARNINGS ARE PRINTED'
        ELSE
        ENDIF
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND OF KP TERM (INITIAL STATE RADIATION IN CHARGED CURRENT)
C   HS(19/07/94)
C   FROM epcctot.f (28/08/98)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCCKL(X)
C
C  X(1) -->  XX
C  X(2) -->  Q**2  -->  Y
C  X(3) -->  XS
C  X(4) -->  LOG(2*K.P)
C  X(5) -->  TS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSKNCC/ SXNRCC,SX1NCC
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSKPXY/ XX,Y
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSXSLM/ XSMIN,XSCUT
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSPSPC/ IPHSPC
      COMMON /HSWGTC/ IWEIGS
      DIMENSION X(5)
C
C---X-VALUE
      XX=XMIN+(XMAX-XMIN)*X(1)
C---Y-VALUE
      GS=SP-MEI2-MPRO2
      YMAXX=XX*(1D0-4D0*MEI2*MPRO2/GS/GS)/(XX*(1D0+XX*MPRO2/GS)+MEI2/GS)
      IF(ICUT.LT.3) THEN
C---CUT IN EXTERNALLY DEFINED Q**2(MIN)
        Q2L=Q2MIN
        Q2U=XX*GS
        ELSEIF(ICUT.EQ.3) THEN
C---CUT IN Y
        QQ2MIN=XX*YMIN*GS
C---CUT ON ELECTRON SCATTERING ANGLE (MASSES NEGLECTED)
        YMINTH=1D0/(1D0+XX*CTHCON)
        QT2MIN=XX*YMINTH*GS
C---CUT ON ELECTRON TRANSVERSE MOMENTUM (MASSES NEGLECTED)
        QP2MIN=XX*GS/2D0*(1D0-DSQRT(1D0-PTXM0/XX))
        QP2MAX=XX*GS/2D0*(1D0+DSQRT(1D0-PTXM0/XX))
        YP2MAX=QP2MAX/GS/XX
        Q2L=MAX(Q2MIN,QQ2MIN,QT2MIN,QP2MIN)
        Q2U=XX*MIN(YMAX,YMAXX,YP2MAX)*GS
        Q2U=MIN(Q2U,Q2MAX)
        ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSCCKL'
        STOP
      ENDIF
      DQ2=DMAX1(Q2U-Q2L,0D0)
C
C---CUT IN W LATER
C
      Q2=Q2L+X(2)*DQ2
      Y=Q2/XX/GS
C---EXTERNAL WEIGHT
      IACPT=1
      IF (IWEIGS.GT.0) THEN
        CALL HSWGTX(XX,Y,IACPT)
        IF (IACPT.EQ.0) THEN
          HSCCKL=0D0
          RETURN
        ENDIF
      ENDIF
      CALL HSDELO(XX,Y)
      XSMIN=HSCXSM(XX,Y)
      XSMAX=1D0
      XSMINI=XSMIN
      IF(XSMAX.LE.XSMINI) THEN
        HSCCKL=0D0
        RETURN
      ENDIF
C---SUBSTITUTION FOR XS
      SMF=MEI2+MQI2+MQF2
      UMIN=DLOG(XSMINI-XX-SMF/Y/SP)
      UMAX=DLOG(XSMAX-XX-SMF/Y/SP)
      UV=UMIN+(UMAX-UMIN)*X(3)
      XS=XX+SMF/Y/SP+DEXP(UV)
C
      IF(IPRINT.GT.30) THEN
        WRITE(LUNTES,230)SP,XX,Y,XSMIN,XSMAX,XS
230     FORMAT(/,' SP = ',1PD12.3,' X = ',D12.6,'   Y = ',D12.6,/
     F        ,' XSMIN = ',D17.11,'   XSMAX = ',D17.11,'  XS = ',D12.6)
      ENDIF
C---INTEGRATE IN CMSYSTEM (K+QS=0)
      CALL HSCCMS(XX,Y,XS)
      IF(IPHSPC.EQ.1) THEN
        HSCCKL=0D0
        RETURN
      ENDIF
      CALL HSCLAB(XX,Y,XS)
      CALL HSLZK1(ZMIN,ZMAX)
      IF(IPHSPC.EQ.1.OR.ZMIN.GE.ZMAX) THEN
        HSCCKL=0D0
        RETURN
      ENDIF
C---SUBSTITUTION W = LOG(E/P-COSTHETA)
      IF (ZMAX.EQ.1D0) THEN
        WVMIN=DLOG(MEI2/2D0/PEL/PEL)
       ELSE
        WVMIN=DLOG(EEL/PEL-ZMAX)
      ENDIF
      WVMAX=DLOG(EEL/PEL-ZMIN)
      WVV=WVMIN+(WVMAX-WVMIN)*X(4)
      ZP=EEL/PEL-DEXP(WVV)
C---NO SUBSTITUTION FOR PHI
      CALL HSCPHL(ZP,PHPMIN)
      PHPMAX=PI
      PHIP=PHPMIN+(PHPMAX-PHPMIN)*X(5)
C
      SZ2=1D0-ZP*ZP
      IF (SZ2.LT.0D0) THEN
        HSCCKL=0D0
        RETURN
      ENDIF
      SZ=DSQRT(SZ2)
      CPHI=DCOS(PHIP)
      DKP=OMEGA*(EEL-PEL*ZP)
      DKPS=OMEGA*(ES-PS*(SINE*SZ*CPHI+COSE*ZP))
      DKQ=OMEGA*(EQ-PS*(SINE*SZ*CPHI+COSE*ZP)+PEL*ZP)
      DKQS=DKP+DKQ-DKPS
C---CHECK CONDITION ON PHOTON ENERGY
      OMH=(PH*DKQ+PQH*DKP)/(PH*EQH+PQH*EH)
      IF (OMH.LT.DELTA) THEN
        HSCCKL= 0D0
        RETURN
      ENDIF
C
      S=XS*SP
      T=-XX*Y*SP
      U=-XS*(1D0-Y)*SP
      TS=T-2D0*(DKP-DKPS)
      SS=S-2D0*(DKP+DKQ)
      US=U+2D0*(DKPS-DKQ)
C
      CALL HSPVER(XS,-TS)
      IF (HPOLAR.NE.0D0) THEN
        CALL HSDPVR(XS,-TS)
      ENDIF
      CALL HSIROT
      D  = 1D0/(T -MW2)
      DS = 1D0/(TS-MW2)
      CCMSLQ = 2D0*DS*DS*SS*DKQ*DKQS
     *       + 2D0*D*D*S*DKPS*DKP
     *       - D*DS*(  S*SS*US + S*(SS-T)*DKQS - SS*(S-TS)*DKP
     *                  -US*(SS*DKQ-S*DKPS)  )
     *       -2D0*D*DS*DS*(SS*T*DKQ - S*T*DKQS - S*(SS-US)*DKPS - S*SS*T
     *                      - (SS*U - 4D0*SS*DKQ)*DKP ) * DKQS
     *       +2D0*D*D*DS*(SS*TS*DKP - S*TS*DKPS + SS*(S-US)*DKQ- S*SS*TS
     *                      + (S*U + 4D0*S*DKPS)*DKQS ) * DKP
     *       +4D0*D*D*DS*DS * DKP * DKQS
     *         *( - S*SS*T + S*SS*(DKP-DKPS)
     *            + 2D0*S *DKPS*DKQS + 2D0*SS*DKP *DKQ
     *            + 2D0*US*DKPS*DKQ  + 2D0*U *DKP *DKQS
     *            - 2D0*T *DKQ *DKQS - 2D0*TS*DKP *DKPS )
      CCMSLB = - 2D0*DS*DS*U*DKQS*DKQS
     *         - 2D0*D*D*U*DKP*DKP
     *       + D*DS*( -U*US*US + 4D0*U*DKP*DKQS
     *                        - 2D0*U*US*(DKQS-DKP)  )
     *       -2D0*D*DS*DS*(-U*T*DKQS + US*T*DKQ - US*(U-S)*DKPS - U*US*T
     *                      - (U*SS + 4D0*U*DKQS)*DKP ) * DKQS
     *       +2D0*D*D*DS*(U*TS*DKP - US*TS*DKPS + US*(U-SS)*DKQ- U*US*TS
     *                      + (S*U - 4D0*U*DKP)*DKQS ) * DKP
     *       +4D0*D*D*DS*DS * DKP * DKQS
     *         *( - U*US*T + U*US*(DKP-DKPS)
     *            - 2D0*US*DKPS*DKQ  - 2D0*U *DKP *DKQS
     *            - 2D0*S *DKPS*DKQS - 2D0*SS*DKP *DKQ
     *            - 2D0*T *DKQ *DKQS - 2D0*TS*DKP *DKPS )
      CCMSQ  = - MEI2*D*D*SS*SS/DKP
      CCMSB  = - MEI2*D*D*U *U /DKP
      CCMSLQ = CCMSLQ/DKQS
      CCMSLB = CCMSLB/DKQS
C      CCMSHQ = - MQF2*D*D*S*S/DKQS/DKQS*DKP
C      CCMSHB = - MQF2*D*D*U*U/DKQS/DKQS*DKP
      CCMSHQ = 0D0
      CCMSHB = 0D0
C
      IF (LLEPT.EQ.-1) THEN
      CQP(1) = (QU +HPOLAR*DQU)   * (CCMSQ + CCMSLQ + CCMSHQ)
      CQP(2) = 0D0
      CQP(3) = 0D0
      CQP(4) = (QBD-HPOLAR*DQBD) * (CCMSB + CCMSLB + CCMSHB)
      CQP(5) = 0D0
      CQP(6) = (QBS-HPOLAR*DQBS) * (CCMSB + CCMSLB + CCMSHB)
      CQP(7) = (QC +HPOLAR*DQC)  * (CCMSQ + CCMSLQ + CCMSHQ)
      CQP(8) = 0D0
      CQP(9) = 0D0
C      CQP(10)= (QBB-HPOLAR*DQBB) * (CCMSB + CCMSLB + CCMSHB)
Chs: no b-quarks
      CQP(10)= 0D0
      CQP(11)= (QT +HPOLAR*DQT)  * (CCMSQ + CCMSLQ + CCMSHQ)
      CQP(12)= 0D0
      ELSEIF(LLEPT.EQ.1) THEN
      CQP(1) = 0D0
      CQP(2) = (QBU-HPOLAR*DQBU) * (CCMSQ + CCMSLQ + CCMSHQ)
      CQP(3) = (QD +HPOLAR*DQD)  * (CCMSB + CCMSLB + CCMSHB)
      CQP(4) = 0D0
      CQP(5) = (QS +HPOLAR*DQS)  * (CCMSB + CCMSLB + CCMSHB)
      CQP(6) = 0D0
      CQP(7) = 0D0
      CQP(8) = (QBC-HPOLAR*DQBC) * (CCMSQ + CCMSLQ + CCMSHQ)
C      CQP(9) = (QB +HPOLAR*DQB)  * (CCMSB + CCMSLB + CCMSHB)
Chs: no b-quarks
      CQP(9) = 0D0
      CQP(10)= 0D0
      CQP(11)= 0D0
      CQP(12)= (QBT-HPOLAR*DQBT) * (CCMSQ + CCMSLQ + CCMSHQ)
      ENDIF

      DO 2 IFL = 2,12
        CQP(IFL) = CQP(IFL-1) + CQP(IFL)
 2    CONTINUE
      SUMME=CQP(12)
C
      RPOL=(1D0+LLEPT*POLARI)/2D0
      SQGRAM=(EEL+EQ-ES)/OMEGA
      HSCCKL=SUMME*Y*SX1NCC*RPOL/SQGRAM/XS
     *       *(UMAX-UMIN)*(XS-XX-SMF/Y/SP)
     *       *(WVMAX-WVMIN)*(EEL/PEL-ZP)/DKP
     *       *2D0*(PHPMAX-PHPMIN)
     *       *(XMAX-XMIN)*DQ2/(XX*SP)
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   LIMITS FOR AZIMUTH (1/KP-PARAMETRIZATION FOR CC)
C   FROM CPHPLM (epcctot.f, 28/08/98)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSCPHL(Z,PHPMIN)
      IMPLICIT double precision (A-H,M,O-Z)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
C
      DEPS=DELTA/OMEGA*(PQH*EELE+PELE*EQH)
      SIGM=PQH*EEL+PELE*EQ
      TAU=PELE*(PS*COSE-PEL)+PQH*PEL
      IF (((1D0-Z*Z).LT.1D-10).AND.((SIGM-DEPS-TAU*Z).LT.0D0)) THEN
        PHPMIN=PI
        RETURN
      ELSEIF (((1D0-Z*Z).LT.1D-10).AND.((SIGM-DEPS-TAU*Z).GT.0D0)) THEN
        PHPMIN=0D0
        RETURN
      ENDIF
      CMAX=(SIGM-DEPS-TAU*Z)/PELE/PS/SINE/DSQRT(1D0-Z*Z)
      IF (CMAX.LE.(-1D0)) THEN
        PHPMIN=PI
        RETURN
      ENDIF
      IF (CMAX.GE.1D0) THEN
        PHPMIN=0D0
        RETURN
      ELSE
        PHPMIN=DACOS(CMAX)
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND OF KQ TERM (RADIATION FROM THE INITIAL QUARK
C                         IN CHARGED CURRENT SCATTERING)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCCQI(X)
C
C  X(1) -->  XX
C  X(2) -->  Y
C  X(3) -->  XS
C  X(4) -->  LOG(2*K.P)
C  X(5) -->  TS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSKPXY/ XX,Y
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSXSLM/ XSMIN,XSCUT
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSWGTC/ IWEIGS
      DIMENSION X(5)
C
C      DATA ICOUNT /0/

C---EXTERNAL WEIGHT
c     IACPT=1
c     IF (IWEIGS.GT.0) THEN
c       CALL HSWGTX(XX,Y,IACPT)
c       IF (IACPT.EQ.0) THEN
c         HSCCQI=0D0
c         RETURN
c       ENDIF
c     ENDIF

      HSCCQI=0D0
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND OF KQS TERM (IN CHARGED CURRENT SCATTERING)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCCQF(X)
C
C  X(1) -->  XX
C  X(2) -->  Q2   -->  Y
C  X(3) -->  U    -->  XS
C  X(4) -->  A1=(2*K.P)
C  X(5) -->  TS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSKNCC/ SXNRCC,SX1NCC
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSKPXY/ XX,Y
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSPCUT/ PTMIN,PTXM0
      COMMON /HSXSLM/ XSMIN,XSCUT
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSWGTC/ IWEIGS
      DIMENSION X(5)
C
C---X-VALUE
      XX=XMIN+(XMAX-XMIN)*X(1)
      GS=SP-MEI2-MPRO2
      YMAXX=XX*(1D0-4D0*MEI2*MPRO2/GS/GS)/(XX*(1D0+XX*MPRO2/GS)+MEI2/GS)
C---Y-VALUE
      IF(ICUT.LT.3) THEN
C---CUT IN EXTERNALLY DEFINED Q**2(MIN)
        Q2L=Q2MIN
        Q2U=XX*GS
        ELSEIF(ICUT.EQ.3) THEN
C---CUT IN Y
        QQ2MIN=XX*YMIN*GS
C---CUT ON ELECTRON SCATTERING ANGLE (MASSES NEGLECTED)
        YMINTH=1D0/(1D0+XX*CTHCON)
        QT2MIN=XX*YMINTH*GS
C---CUT ON ELECTRON TRANSVERSE MOMENTUM (MASSES NEGLECTED)
        QP2MIN=XX*GS/2D0*(1D0-DSQRT(1D0-PTXM0/XX))
        QP2MAX=XX*GS/2D0*(1D0+DSQRT(1D0-PTXM0/XX))
        YP2MAX=QP2MAX/GS/XX
        Q2L=MAX(Q2MIN,QQ2MIN,QT2MIN,QP2MIN)
        Q2U=XX*MIN(YMAX,YMAXX,YP2MAX)*GS
        Q2U=MIN(Q2U,Q2MAX)
        ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSCCQF'
        STOP
      ENDIF
      DQ2=DMAX1(Q2U-Q2L,0D0)
C---CUT IN W LATER
      Q2=Q2L+X(2)*DQ2
      Y=Q2/(XX*GS)
C
C---EXTERNAL WEIGHT
      IACPT=1
      IF (IWEIGS.GT.0) THEN
        CALL HSWGTX(XX,Y,IACPT)
        IF (IACPT.EQ.0) THEN
          HSCCQF=0D0
          RETURN
        ENDIF
      ENDIF
      CALL HSFIVC(XX,Y)
      CALL HSDELO(XX,Y)
C
      XSMIN=HSXSMN(XX,Y)
      XSMAX=1D0
      XSMINI=XSMIN
      IF(XSMAX.LE.XSMINI) THEN
        HSCCQF=0D0
        RETURN
      ENDIF
C---SUBSTITUTION FOR XS
      M2S =(MEI2+MQI2+MQF2)/Y/SP
      UMIN=Y*SP/MQF2*(DLOG(XSMINI-XX-M2S)
     *               -DLOG(XSMINI-XX-(MEI2+MQI2)/Y/SP) )
      UMAX=Y*SP/MQF2*(DLOG(XSMAX-XX-M2S)
     *               -DLOG(XSMAX-XX-(MEI2+MQI2)/Y/SP) )
      UV=UMIN+(UMAX-UMIN)*X(3)
      EXUV=DEXP(UV*MQF2/Y/SP)
      XS=((XX+(MEI2+MQI2)/Y/SP)*EXUV-XX-M2S)/(EXUV-1D0)
Ch    UMIN=DLOG(XSMINI-XX)
Ch    UMAX=DLOG(XSMAX-XX)
Ch    UV=UMIN+(UMAX-UMIN)*X(3)
Ch    XS=XX+DEXP(UV)
C
      IF(IPRINT.GT.30) THEN
        WRITE(LUNTES,230)SP,XX,Y,XSMIN,XSMAX,XS
230     FORMAT(/,' SP = ',1PD12.3,' X = ',D12.6,'   Y = ',D12.6,/
     F        ,' XSMIN = ',D17.11,'   XSMAX = ',D17.11,'  XS = ',D12.6)
      ENDIF
C
      CALL HSFCMS(XX,Y,XS)
      CALL HSFLAB(XS)
      CALL HSLZK1(ZMIN,ZMAX)
C---NO SUBSTITUTION FOR A1=K.P
      A1MAX=2D0*OMEGA*(EEL-PEL*ZMIN)
      IF (ZMAX.EQ.1D0) THEN
        A1MIN=2D0*OMEGA*MEI2/2D0/EEL
      ELSE
        A1MIN=2D0*OMEGA*(EEL-PEL*ZMAX)
      ENDIF
      A1=A1MIN+(A1MAX-A1MIN)*X(4)
C---NO SUBSTITUTION FOR TS
      CALL HSLTS1(A1,XX,Y,XS,TSMIN,TSMAX,TSM,TSP,CFKP)
      TS=TSMIN+(TSMAX-TSMIN)*X(5)
      SQGRM2=-CFKP*(TS-TSM)*(TS-TSP)
      IF (SQGRM2.LE.0D0) THEN
        HSCCQF=0D0
        RETURN
        ELSE
        SQGRAM=DSQRT(SQGRM2)
      ENDIF
      CALL HSFIV1(XX,Y,XS,A1,TS)
      CALL HSPVER(XS,-TS)
      IF (HPOLAR.NE.0D0) THEN
        CALL HSDPVR(XS,-TS)
      ENDIF
      CALL HSIROT
C---CHECK CONDITION ON PHOTON ENERGY
      OMH=(PH*DKQ+PQH*DKP)/(PH*EQH+PQH*EH)
      IF (OMH.LT.DELTA) THEN
        HSCCQF=0D0
        RETURN
      ENDIF
C
C---MATRIX ELEMENT
      DBOS=1D0/(T-MW2)
      CCMSQ=-MQF2*DBOS*DBOS*S*S/DKQS/DKQS
      CCMSB=-MQF2*DBOS*DBOS*U*U/DKQS/DKQS

      IF (LLEPT.EQ.-1) THEN
        CQP(1)=(QU +HPOLAR*DQU)*CCMSQ
        CQP(2)=0D0
        CQP(3)=0D0
        CQP(4)=(QBD-HPOLAR*DQBD)*CCMSB
        CQP(5)=0D0
        CQP(6)=(QBS-HPOLAR*DQBS)*CCMSB
        CQP(7)=(QC +HPOLAR*DQC)*CCMSQ
        CQP(8)=0D0
        CQP(9)=0D0
C        CQP(10)=(QBB-HPOLAR*DQBB)*CCMSB
Chs: no b-quarks
        CQP(10)= 0D0
        CQP(11)=(QT +HPOLAR*DQT)*CCMSQ
        CQP(12)=0D0
      ELSEIF (LLEPT.EQ.1) THEN
        CQP(1)=0D0
        CQP(2)=(QBU-HPOLAR*DQBU)*CCMSQ
        CQP(3)=(QD +HPOLAR*DQD)*CCMSB
        CQP(4)=0D0
        CQP(5)=(QS +HPOLAR*DQS)*CCMSB
        CQP(6)=0D0
        CQP(7)=0D0
        CQP(8)=(QBC-HPOLAR*DQBC)*CCMSQ
C        CQP(9)=(QB +HPOLAR*DQB)*CCMSB
Chs: no b-quarks
        CQP(9)=0D0
        CQP(10)=0D0
        CQP(11)=0D0
        CQP(12)=(QBT-HPOLAR*DQBT)*CCMSQ
      ENDIF

      DO 2 IFL=2,12
        CQP(IFL)=CQP(IFL-1)+CQP(IFL)
 2    CONTINUE
      SUMME=CQP(12)
C
      RPOL=(1D0+LLEPT*POLARI)/2D0
      HSCCQF=SUMME*Y*SX1NCC*RPOL/SQGRAM/XS
     *      *(UMAX-UMIN)*(XS-XX-M2S)*(XS-XX-(MEI2+MQI2)/Y/SP)
Ch   *      *(UMAX-UMIN)*(XS-XX)
     *      *(A1MAX-A1MIN)*(TSMAX-TSMIN)
     *      *(XMAX-XMIN)*DQ2/(XX*SP)
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     CROSS SECTION AND
C     VIRTUAL AND SOFT BREMSSTRAHLUNG CORRECTIONS
C     FOR DEEP INELASTIC ELECTRON PROTON SCATTERING
C     VIA  CHARGED CURRENT
C     AT HERA
C
C     AUTHOR: H.SPIESBERGER, DEC 1990 (20.03.87)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     CROSS SECTION FOR CHARGED CURRENT ELECTRON PROTON SCATTERING
C
      FUNCTION HSSGCC (X,Y,LL,POL,HPOL,LQ)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSRWW,HSENUW,HSDUWQ,HSDUWA,CMW2,CMZ2,CMW20,CMZ20
     *          ,PIW,CENW,CDUWQ,CDUWA
     *          ,HSWG1L,HSWG2L,CWG2QL,CWG1AL
     *          ,HSWG1R,HSWG2R,CWG1QR,CWG2QR,CWG1AR,CWG2AR
     *          ,HSIWZ1,HSIWZ2
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSKNCC/ SXNRCC,SX1NCC
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSCBMS/ CMW2,CMZ2
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSDELR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSWGTC/ IWEIGS
C
      HSSGCC=0D0

C---EXTERNAL WEIGHT
      IACPT=1
      IF (IWEIGS.GT.0) THEN
        CALL HSWGTX(X,Y,IACPT)
        IF (IACPT.EQ.0) THEN
          HSSGCC=0D0
          RETURN
        ENDIF
      ENDIF
C
C---PREPARE GENERAL FACTORS, CALL PARTON DISTRIBUTIONS
      T=-(SP-MEI2-MPRO2)*X*Y
      RPOL=(1D0+LL*POL)/2D0
      SN0CC=SXNRCC*SP*X/(T-MW2)/(T-MW2)*RPOL
      FQA=(1D0-Y)*(1D0-Y)
      CALL HSDELO(X,Y)
      CALL HSPVER(X,-T)
      IF (HPOL.NE.0D0) THEN
        CALL HSDPVR(X,-T)
      ENDIF
      CALL HSIROT
      DO 10 I=1,12
        CQP(I)=0D0
 10   CONTINUE
 
C
C---BORN CROSS SECTION
      IF (LPAR(2).EQ.0) THEN
        IF (LL.EQ.-1) THEN
          CQP(1) =CQP(1)  + (QU+HPOL*DQU)
          CQP(4) =CQP(4)  + (QBD-HPOL*DQBD)*FQA
          CQP(6) =CQP(6)  + (QBS-HPOL*DQBS)*FQA
          CQP(7) =CQP(7)  + (QC+HPOL*DQC)
C          CQP(10)=CQP(10) + (QBB-HPOL*DQBB)*FQA
Chs: no b-quarks
          CQP(11)=CQP(11) + (QT+HPOL*DQT)
        ELSEIF(LL.EQ.1) THEN
          CQP(2) =CQP(2)  + (QBU-HPOL*DQBU)
          CQP(3) =CQP(3)  + (QD+HPOL*DQD)*FQA
          CQP(5) =CQP(5)  + (QS+HPOL*DQS)*FQA
          CQP(8) =CQP(8)  + (QBC-HPOL*DQBC)
C          CQP(9) =CQP(9)  + (QB+HPOL*DQB)*FQA
Chs: no b-quarks
          CQP(12)=CQP(12) + (QBT-HPOL*DQBT)
        ELSE
          WRITE(LUNOUT,*) ' ERROR IN HSSGCC: WRONG LEPTON CHARGE = ',LL
          STOP
        ENDIF
        IF (LQ.GT.0) THEN
          HSSGCC=CQP(2*LQ-1)*SN0CC
        ELSEIF (LQ.LT.0) THEN
          HSSGCC=CQP(-2*LQ)*SN0CC
        ELSEIF (LQ.EQ.0) THEN
          CQP(1)=CQP(1)*SN0CC
          DO 11 I=2,12
            CQP(I)=CQP(I-1)+CQP(I)*SN0CC
 11       CONTINUE
          HSSGCC=CQP(12)
        ENDIF
        RETURN
      ENDIF
C
C---1-LOOP CORRECTIONS
      IF (LPAR(2).EQ.1) THEN
        DFQ=1D0
        DFA=1D0
        CQU=2D0/3D0
        CQD=-1D0/3D0
        S=X*SP
        U=-X*SP-T
        CMW20=DCMPLX(MW2,-1D-8)
        CMZ20=DCMPLX(MZ2,-1D-8)
C---NORMALIZATION
        SN1CC=SN0CC
C---W SELF ENERGY
        IF (LPAR(10).GT.0) THEN
          PIW=HSSRWW(T)
          RPIW=DREAL(PIW)/(T-MW2)
          DFQ=DFQ-2D0*RPIW
          DFA=DFA-2D0*RPIW
        ENDIF
        IF ((LPAR(4).GT.1).AND.(LPAR(10).GT.0)) THEN
          DFQ=DFQ-2D0*(DELTAR-DRPIW2)
          DFA=DFA-2D0*(DELTAR-DRPIW2)
        ENDIF

C---PHOTONIC CORRECTIONS ACCORDING TO BARDIN'S PRESCRIPTION
        IF (LPAR(11).EQ.1) THEN
          DFQ=DFQ+HSCCBQ(X,Y)
          DFA=DFA+HSCCBA(X,Y)
        ENDIF
C---REMNANT OF PHOTONIC CORRECTIONS
        IF (LPAR(10).EQ.0) GOTO 300
        DFQ=DFQ+HSCCSQ(X,Y)
        DFA=DFA+HSCCSA(X,Y)
C---LEPTONIC VERTEX
        CENW=HSENUW(T)
        DELVNU=2D0*DREAL(CENW)
        DFQ=DFQ+DELVNU
        DFA=DFA+DELVNU
        IF (LPAR(12).EQ.1) THEN
          DFQ=DFQ+ALP2PI*((DLOG(MZ2/MD2)+DLOG(MZ2/MEI2))/2D0+9D0/2D0)
          DFA=DFA+ALP2PI*(-(DLOG(MZ2/MU2)-DLOG(MZ2/MEI2))/2D0)
C---GAMMA W BOXES (LEPTONIC PART, ONLY LOGS)
C         CWG2Q=CFWG2(T,U,CMW20)
C         CWG1A=CFWG1(T,U,CMW20)
          CWG2QL=HSWG2L(T,U,CMW20)
          CWG1AL=HSWG1L(T,U,CMW20)
          DFQ=DFQ+4D0*DREAL(CWG2QL)
          DFA=DFA-4D0*DREAL(CWG1AL)
        ENDIF
C---QUARK VERTEX
        CDUWQ=HSDUWQ(T)
        CDUWA=HSDUWA(T)
        DFQ=DFQ+2D0*DREAL(CDUWQ)
        DFA=DFA+2D0*DREAL(CDUWA)
C---QUARK PROPAGATOR RESIDUUM
        IF (LPAR(13).EQ.1) THEN
          DSQ=ALP2PI*(-DLOG(MZ2/MU2)+DLOG(MZ2/MD2))/2D0
          DFQ=DFQ+CQU*CQU*DSQ
          DFA=DFA+CQD*CQD*DSQ
        ENDIF
C---LEPTON QUARK INTERFERENCE
        IF (LPAR(14).EQ.1) THEN
          DFQ=DFQ-ALP2PI*CQU*(DLOG(MZ2/MD2)+9D0/2D0)
          DFA=DFA-ALP2PI*CQD*(DLOG(MZ2/MU2)+9D0/2D0)
C---GAMMA W BOXES (INTERFERENCE PART)
          CWG1QR=HSWG1L(T,S,CMW20)+HSWG1R(T,S,CMW20)
          CWG2QR=HSWG2L(T,U,CMW20)+HSWG2R(T,U,CMW20)
          CWG1AR=HSWG1L(T,U,CMW20)+HSWG1R(T,U,CMW20)
          CWG2AR=HSWG2L(T,S,CMW20)+HSWG2R(T,S,CMW20)
          DFQ=DFQ-4D0*CQU*(DREAL(CWG1QR)+DREAL(CWG2QR))
          DFA=DFA-4D0*CQD*(DREAL(CWG1AR)+DREAL(CWG2AR))
        ENDIF
C
C---WEAK Z-W BOXES
        IF (LPAR(15).EQ.1) THEN
C---NON-LOG TERMS FROM GAMMA-W BOXES
          CWG2QR=HSWG2R(T,U,CMW20)
          CWG1AR=HSWG1R(T,U,CMW20)
          DFQ=DFQ+4D0*DREAL(CWG2QR)
          DFA=DFA-4D0*DREAL(CWG1AR)
C---W-Z BOXES
          VEZ=VAFI(1,1,2)
          VNZ=VAFI(2,2,2)
          VUZ=VAFI(1,2,2)
          VDZ=VAFI(1,3,2)
          AEZ=-VNZ
          ANZ= VNZ
          ADZ=-VNZ
          AUZ= VNZ
          BWZ1P=ALP2PI*((VEZ+AEZ)*(VUZ+AUZ)+(VNZ+ANZ)*(VDZ+ADZ))
     *         *DREAL(HSIWZ1(T,S,CMW20,CMZ20))
          BWZ2P=ALP2PI*((VEZ+AEZ)*(VDZ+ADZ)+(VNZ+ANZ)*(VUZ+AUZ))
     *         *DREAL(HSIWZ2(T,U,CMW20,CMZ20))
          BWZ1A=ALP2PI*((VEZ+AEZ)*(VDZ+ADZ)+(VNZ+ANZ)*(VUZ+AUZ))
     *         *DREAL(HSIWZ2(T,S,CMW20,CMZ20))
          BWZ2A=ALP2PI*((VEZ+AEZ)*(VUZ+AUZ)+(VNZ+ANZ)*(VDZ+ADZ))
     *         *DREAL(HSIWZ1(T,U,CMW20,CMZ20))
          DFQ=DFQ+4D0*(T-MW2)*(BWZ1P+BWZ2P)
          DFA=DFA+4D0*(T-MW2)*(BWZ1A+BWZ2A)
        ENDIF
C---
  300   CONTINUE
        IF (LL.EQ.-1) THEN
          CQP(1) =CQP(1)  + (QU+HPOL*DQU)*DFQ
          CQP(4) =CQP(4)  + (QBD-HPOL*DQBD)*FQA*DFA
          CQP(6) =CQP(6)  + (QBS-HPOL*DQBS)*FQA*DFA
          CQP(7) =CQP(7)  + (QC+HPOL*DQC)*DFQ
C          CQP(10)=CQP(10) + (QBB-HPOL*DQBB)*FQA*DFA
Chs: no b-quarks
          CQP(11)=CQP(11) + (QT+HPOL*DQT)*DFQ
        ELSEIF(LL.EQ.1) THEN
          CQP(2) =CQP(2)  + (QBU-HPOL*DQBU)*DFQ
          CQP(3) =CQP(3)  + (QD+HPOL*DQD)*FQA*DFA
          CQP(5) =CQP(5)  + (QS+HPOL*DQS)*FQA*DFA
          CQP(8) =CQP(8)  + (QBC-HPOL*DQBC)*DFQ
C          CQP(9) =CQP(9)  + (QB+HPOL*DQB)*FQA*DFA
Chs: no b-quarks
          CQP(12)=CQP(12) + (QBT-HPOL*DQBT)*DFQ
        ELSE
          WRITE(LUNOUT,*) ' ERROR IN HSSGCC: WRONG LEPTON CHARGE = ',LL
          STOP
        ENDIF
C---TOTAL CROSS SECTION
        IF (LQ.GT.0) THEN
          HSSGCC=CQP(2*LQ-1)*SN1CC
        ELSEIF (LQ.LT.0) THEN
          HSSGCC=CQP(-2*LQ)*SN1CC
        ELSEIF (LQ.EQ.0) THEN
          CQP(1)=CQP(1)*SN1CC
          DO 12 I=2,12
            CQP(I)=CQP(I-1)+CQP(I)*SN1CC
 12       CONTINUE
          HSSGCC=CQP(12)
        ENDIF
      ENDIF

      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   SOFT BREMSSTRAHLUNG CORRECTION FOR ELECTRON QUARK SCATTERING
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCCSQ(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C
      HSCCSQ=0D0
      CQU=2D0/3D0
      GSP=SP-MEI2-MPRO2
      S=X*SP
      T=-GSP*X*Y
      U=-X*SP-T
      DLMWME=DLOG(MW2/MEI2)
      DLMWMU=DLOG(MW2/MU2)
      DLMWMD=DLOG(MW2/MD2)
      DLQ2MW=DLOG(-T/MW2)
C
C---LEPTONIC PART
      IF (LPAR(12).EQ.1) THEN
        HSCCSQ=HSCCSQ+ALP2PI
     .  *(-0.5D0*DLMWME*(DLMWME+3D0)-0.5D0*DLMWMD*(DLMWMD+3D0)
     .    -DLOG(-MW2/U)*(DLOG(-MEI2/U)+DLOG(-MD2/U))
     .    +DLQ2MW*(DLQ2MW-2D0*DLOG(T/U)))
      ENDIF
C---QUARKONIC PART
      IF (LPAR(13).EQ.1) THEN
        HSCCSQ=HSCCSQ+ALP2PI*CQU*CQU
     .  *(-0.5D0*DLMWMU*(DLMWMU+3D0)-0.5D0*DLMWMD*(DLMWMD+3D0)
     .    -DLOG(-MW2/T)*(DLOG(-MU2/T)+DLOG(-MD2/T))
     .    -2D0*DLOG(-MD2/T)+DLQ2MW*(DLQ2MW-3D0))
      ENDIF
C---LEPTON QUARK INTERFERENCE
      IF (LPAR(14).EQ.1) THEN
        HSCCSQ=HSCCSQ+ALP2PI*(-CQU)
     .  *(+DLOG(-U/S)*DLOG(ME2/S)+DLOG(-T/S)*DLOG(MU2/S)
     .    -DLOG(MW2*MW2/T/U)*DLOG(-MD2/U)-2D0*DLOG(-MD2/T)
     .    -2D0*DLQ2MW*DLOG(-T/S)-DLOG(-U/S)*DLOG(U/T)
     .    -DLMWMD*(DLMWMD+3D0)+DLQ2MW*(DLQ2MW-3D0+2D0*DLOG(-U/S)))
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   SOFT BREMSSTRAHLUNG CORRECTION FOR ELECTRON ANTI-QUARK SCATTERING
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCCSA(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C
      HSCCSA=0D0
      CQD=-1D0/3D0
      GSP=SP-MEI2-MPRO2
      S=X*SP
      T=-GSP*X*Y
      U=-X*SP-T
      DLMWME=DLOG(MW2/MEI2)
      DLMWMD=DLOG(MW2/MD2)
      DLMWMU=DLOG(MW2/MU2)
      DLQ2MW=DLOG(-T/MW2)
C
C---LEPTONIC PART
      IF (LPAR(12).EQ.1) THEN
        HSCCSA=HSCCSA+ALP2PI
     .  *(-0.5D0*DLMWME*(DLMWME+3D0)-0.5D0*DLMWMU*(DLMWMU+3D0)
     .    -DLOG(-MW2/U)*(DLOG(-MEI2/U)+DLOG(-MU2/U))
     .    +DLQ2MW*(DLQ2MW-2D0*DLOG(T/U)))
      ENDIF
C---QUARKONIC PART
      IF (LPAR(13).EQ.1) THEN
        HSCCSA=HSCCSA+ALP2PI*CQD*CQD
     .  *(-0.5D0*DLMWMU*(DLMWMU+3D0)-0.5D0*DLMWMD*(DLMWMD+3D0)
     .    -DLOG(-MW2/T)*(DLOG(-MU2/T)+DLOG(-MD2/T))
     .    -2D0*DLOG(-MD2/T)+DLQ2MW*(DLQ2MW-3D0))
      ENDIF
C---LEPTON QUARK INTERFERENCE
      IF (LPAR(14).EQ.1) THEN
        HSCCSA=HSCCSA+ALP2PI*(+CQD)
     .  *(+DLOG(-U/S)*DLOG(MEI2/S)+DLOG(-T/S)*DLOG(MD2/S)
     .    -DLOG(MW2*MW2/T/U)*DLOG(-MU2/U)-2D0*DLOG(-MD2/T)
     .    -DLOG(-U/S)*DLOG(U/T)-2D0*DLQ2MW*DLOG(-T/S)
     .    -DLMWMU*(DLMWMU+3D0)+DLQ2MW*(DLQ2MW-3D0+2D0*DLOG(-U/S)))
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCCBQ(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSPEN,ZEQI,ZEQF,ZQIQF,CEQI,CEQF,CQIQF
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C
      HSCCBQ=0D0
      CQU=2D0/3D0
      GSP=SP-MEI2-MPRO2
      S=X*SP
      T=-GSP*X*Y
      U=-X*SP-T
      EEI=EELE
      EQI=X*SP/4D0/EEI
      ENU=EEI*(1D0-Y)-T/4D0/EEI
      EQF=EEI+EQI-ENU
C
      DLMEEE=DLOG(MEI2/4D0/EEI/EEI)
      DLMUEI=DLOG(MU2/4D0/EQI/EQI)
      DLMDEF=DLOG(MD2/4D0/EQF/EQF)
      DLMWME=DLOG(MW2/MEI2)
      DLMWMU=DLOG(MW2/MU2)
      DLMWMD=DLOG(MW2/MD2)
      DLQ2MW=DLOG(-T/MW2)
C
C---LEPTONIC PART
      IF (LPAR(12).EQ.1) THEN
        CEQF=DCMPLX(1D0+4D0*EEI*EQF/U,0D0)
        ZEQF=HSSPEN(CEQF)
        REQF=-DREAL(ZEQF)-DLMEEE*DLMEEE/4D0-DLMDEF*DLMDEF/4D0-PI*PI/3D0
        HSCCBQ=HSCCBQ+ALP2PI
     .  *(2D0*REQF-DLOG(DELTA*DELTA/EQF/EQF)-DLOG(DELTA*DELTA/EEI/EEI)
     .    +0.5D0*DLMWME*(DLMWME+3D0)+0.5D0*DLMWMD*(DLMWMD+3D0)
     .    -DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(-MEI2/U)+DLOG(-MD2/U))
     .    -DLQ2MW*(DLQ2MW-2D0*DLOG(T/U)))
      ENDIF
C---QUARKONIC PART
      IF (LPAR(13).EQ.1) THEN
        CQIQF=DCMPLX(1D0+4D0*EQI*EQF/T,0D0)
        ZQIQF=HSSPEN(CQIQF)
        RQIQF=-DREAL(ZQIQF)-DLMUEI*DLMUEI/4D0-DLMDEF*DLMDEF/4D0
     .        -PI*PI/3D0
        HSCCBQ=HSCCBQ+ALP2PI*CQU*CQU
     .  *(2D0*RQIQF+0.5D0*DLMWMU*(DLMWMU+3D0)+0.5D0*DLMWMD*(DLMWMD+3D0)
     .    -DLOG(DELTA*DELTA/EQI/EQI)-DLOG(DELTA*DELTA/EQF/EQF)
     .    -DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(-MU2/T)+DLOG(-MD2/T))
     .    -DLQ2MW*(DLQ2MW-3D0))
      ENDIF
C---LEPTON QUARK INTERFERENCE
      IF (LPAR(14).EQ.1) THEN
        CEQF=DCMPLX(1D0+4D0*EEI*EQF/U,0D0)
        ZEQF=HSSPEN(CEQF)
        REQF=-DREAL(ZEQF)-DLMEEE*DLMEEE/4D0-DLMDEF*DLMDEF/4D0-PI*PI/3D0
        CQIQF=DCMPLX(1D0+4D0*EQI*EQF/T,0D0)
        ZQIQF=HSSPEN(CQIQF)
        RQIQF=-DREAL(ZQIQF)-DLMUEI*DLMUEI/4D0-DLMDEF*DLMDEF/4D0
     .        -PI*PI/3D0
        CEQI=DCMPLX(1D0-4D0*EEI*EQI/S,0D0)
        ZEQI=HSSPEN(CEQI)
        REQI=-DREAL(ZEQI)-DLMEEE*DLMEEE/4D0-DLMUEI*DLMUEI/4D0-PI*PI/3D0
        HSCCBQ=HSCCBQ+ALP2PI*(-CQU)
     .  *(-2D0*(REQI-RQIQF-REQF)-2D0*DLOG(DELTA*DELTA/EQF/EQF)
     .    +DLMWMD*(DLMWMD+3D0)
     .    +DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(MEI2/S)+DLOG(MU2/S))
     .    -DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(-MU2/T)+DLOG(-MD2/T))
     .    -DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(-MEI2/U)+DLOG(-MD2/U))
     .    -DLQ2MW*(DLQ2MW-3D0+2D0*DLOG(-U/S)))
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   SOFT BREMSSTRAHLUNG CORRECTION FOR ELECTRON ANTI-QUARK SCATTERING
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCCBA(X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSPEN,ZEQI,ZEQF,ZQIQF,CEQI,CEQF,CQIQF
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
C
      HSCCBA=0D0
      CQD=-1D0/3D0
      GSP=SP-MEI2-MPRO2
      S=X*SP
      T=-GSP*X*Y
      U=-X*SP-T
      EEI=EELE
      EQI=X*SP/4D0/EEI
      ENU=EEI*(1D0-Y)-T/4D0/EEI
      EQF=EEI+EQI-ENU
C
      DLMEEE=DLOG(MEI2/4D0/EEI/EEI)
      DLMUEF=DLOG(MU2/4D0/EQF/EQF)
      DLMDEI=DLOG(MD2/4D0/EQI/EQI)
      DLMWME=DLOG(MW2/MEI2)
      DLMWMU=DLOG(MW2/MU2)
      DLMWMD=DLOG(MW2/MD2)
      DLQ2MW=DLOG(-T/MW2)
C
C---LEPTONIC PART
      IF (LPAR(12).EQ.1) THEN
        CEQF=DCMPLX(1D0+4D0*EEI*EQF/U,0D0)
        ZEQF=HSSPEN(CEQF)
        REQF=-DREAL(ZEQF)-DLMEEE*DLMEEE/4D0-DLMUEF*DLMUEF/4D0-PI*PI/3D0
        HSCCBA=HSCCBA+ALP2PI
     .  *(2D0*REQF-DLOG(DELTA*DELTA/EQF/EQF)-DLOG(DELTA*DELTA/EEI/EEI)
     .    +0.5D0*DLMWME*(DLMWME+3D0)+0.5D0*DLMWMU*(DLMWMU+3D0)
     .    -DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(-MEI2/U)+DLOG(-MU2/U))
     .    -DLQ2MW*(DLQ2MW-2D0*DLOG(T/U)))
      ENDIF
C---QUARKONIC PART
      IF (LPAR(13).EQ.1) THEN
        CQIQF=DCMPLX(1D0+4D0*EQI*EQF/T,0D0)
        ZQIQF=HSSPEN(CQIQF)
        RQIQF=-DREAL(ZQIQF)-DLMDEI*DLMDEI/4D0-DLMUEF*DLMUEF/4D0
     .        -PI*PI/3D0
        HSCCBA=HSCCBA+ALP2PI*CQD*CQD
     .  *(2D0*RQIQF-DLOG(DELTA*DELTA/EQI/EQI)-DLOG(DELTA*DELTA/EQF/EQF)
     .    +0.5D0*DLMWMU*(DLMWMU+3D0)+0.5D0*DLMWMD*(DLMWMD+3D0)
     .    -DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(-MU2/T)+DLOG(-MD2/T))
     .    -DLQ2MW*(DLQ2MW-3D0))
      ENDIF
C---LEPTON QUARK INTERFERENCE
      IF (LPAR(14).EQ.1) THEN
        CEQF=DCMPLX(1D0+4D0*EEI*EQF/U,0D0)
        ZEQF=HSSPEN(CEQF)
        REQF=-DREAL(ZEQF)-DLMEEE*DLMEEE/4D0-DLMUEF*DLMUEF/4D0-PI*PI/3D0
        CEQI=DCMPLX(1D0-4D0*EEI*EQI/S,0D0)
        ZEQI=HSSPEN(CEQI)
        REQI=-DREAL(ZEQI)-DLMEEE*DLMEEE/4D0-DLMDEI*DLMDEI/4D0-PI*PI/3D0
        CQIQF=DCMPLX(1D0+4D0*EQI*EQF/T,0D0)
        ZQIQF=HSSPEN(CQIQF)
        RQIQF=-DREAL(ZQIQF)-DLMDEI*DLMDEI/4D0-DLMUEF*DLMUEF/4D0
     .        -PI*PI/3D0
        HSCCBA=HSCCBA+ALP2PI*(+CQD)
     .  *(+2D0*(REQF+RQIQF-REQI)-2D0*DLOG(DELTA*DELTA/EQF/EQF)
     .    -DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(-MEI2/U)+DLOG(-MU2/U))
     .    -DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(-MU2/T)+DLOG(-MD2/T))
     .    +DLOG(4D0*DELTA*DELTA/MW2)*(DLOG(MEI2/S)+DLOG(MD2/S))
     .    +DLMWMU*(DLMWMU+3D0)-DLQ2MW*(DLQ2MW-3D0+2D0*DLOG(-U/S)))
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     1-LOOP FORM FACTORS, WEAK CORRECTIONS
C
C     FOR ELECTRON QUARK SCATTERING:
C      DIRECT BOX:  (T,S)  IN CFWG1
C     CROSSED BOX:  (T,U)  IN CFWG2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSWG1L(GS,GT,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSWG1L,CM2,HSCIR,HSCWLL,HSCWQL
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      HSWG1L=GT/2D0*HSCIR(GT,MU2)
     *       +CM2/2D0*(HSCWLL(GS,CM2)+HSCWQL(GS,CM2,MU2))
      HSWG1L=HSWG1L*ALP2PI
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSWG1R(GS,GT,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSWG1R,CM2,HSD13C,HSCMW,HSCWLR,HSCWQR
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      GU=-(GS+GT)
      HSWG1R=-GU*HSD13C(GT,GS,CM2)+GT/2D0*HSCMW(GT,CM2)
     *       +CM2/2D0*(HSCWLR(GS,CM2)+HSCWQR(GS,CM2))
      HSWG1R=HSWG1R*ALP2PI
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSWG2L(GS,GT,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSWG2L,CM2,HSCWLL,HSCWQL,HSCIR
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      HSWG2L=-CM2/2D0*(HSCWLL(GS,CM2)+HSCWQL(GS,CM2,MD2))
     *       -GT/2D0*HSCIR(GT,MD2)
      HSWG2L=HSWG2L*ALP2PI
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSWG2R(GS,GT,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSWG2R,CM2,HSD13C,HSCMW,HSCWLR,HSCWQR,HSCLN
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      GU=-(GS+GT)
      HSWG2R=((GS-CM2)*(GS-CM2)/2D0/GU-(GT+CM2))*HSD13C(GT,GS,CM2)
     *       -CM2/2D0*(HSCWLR(GS,CM2)+HSCWQR(GS,CM2))
     *       -GT/2D0*HSCMW(GT,CM2)
     *       +(GS-CM2)/2D0/GU*(-HSCLN(-GT/CM2)
     *                         +(GS-CM2)/GS*HSCLN((CM2-GS)/CM2))
      HSWG2R=HSWG2R*ALP2PI
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSD13C(FS,FT,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSD13C,CM2,HSSPEN,HSCLN
      FU=-(FS+FT)
      HSD13C=1D0/FU*(HSSPEN(-FS/CM2)-HSSPEN(FT/CM2)
     *               +HSCLN(-FS/CM2)*HSCLN((CM2+FS)/(CM2-FT)))
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCIR(FS,MQ2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSCIR
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2

      if (fs.lt.0d0) then
        DLNE=DLOG(-MEI2/FS)
        DLNQ=DLOG(-MQ2/FS)
        AIR=-1D0/FS*(DLNE*DLNE/4D0+DLNQ*DLNQ/4D0+PI*PI/6D0)
        HSCIR=DCMPLX(AIR,0D0)
      elseif (fs.eq.0d0) then
        HSCIR=(0D0,0D0)
      else
        DLNE=DLOG(MEI2/FS)
        DLNQ=DLOG(MQ2/FS)
        AIR=-1D0/FS*(DLNE*DLNE/4D0+DLNQ*DLNQ/4D0+2D0*PI*PI/3D0)
        HSCIR=DCMPLX(AIR,0D0)
      endif

c      IF (FS) 1,2,3
c    1 DLNE=DLOG(-MEI2/FS)
c      DLNQ=DLOG(-MQ2/FS)
c      AIR=-1D0/FS*(DLNE*DLNE/4D0+DLNQ*DLNQ/4D0+PI*PI/6D0)
c      HSCIR=DCMPLX(AIR,0D0)
c      RETURN
c    2 HSCIR=(0D0,0D0)
c      RETURN
c    3 DLNE=DLOG(MEI2/FS)
c      DLNQ=DLOG(MQ2/FS)
c      AIR=-1D0/FS*(DLNE*DLNE/4D0+DLNQ*DLNQ/4D0+2D0*PI*PI/3D0)
c      HSCIR=DCMPLX(AIR,0D0)

      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCWLL(FS,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSCWLL,CM2,CSX,HSCLN
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      CSX=DCMPLX(FS,1D-6)
      HSCWLL=-1D0/FS*HSCLN(-CSX/MEI2)*HSCLN(CM2/(CM2-FS))
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCWLR(FS,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSCWLR,CM2,CSX,HSCLN,HSSPEN
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      CSX=DCMPLX(FS,1D-6)
      HSCWLR=-1D0/FS*(HSCLN(-CSX/CM2)*HSCLN(-CSX/CM2)/2D0
     .                +HSSPEN((CSX-CM2)/CSX)-PI*PI/3D0
     .                +DCMPLX(0D0,1D0)*PI*HSCLN((CSX-CM2)/CSX))
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCWQL(FS,CM2,MQ2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSCWQL,HSCLN,CSX,CM2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      CSX=DCMPLX(FS,1D-6)
      HSCWQL=-1D0/FS*HSCLN(-CSX/MQ2)*HSCLN(CM2/(CM2-FS))
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCWQR(FS,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSCWQR,HSSPEN,HSCLN,CSX,CM2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      CSX=DCMPLX(FS,1D-6)
      HSCWQR=-1D0/FS*(HSCLN(-CSX/CM2)*HSCLN(-CSX/CM2)/2D0
     .                +HSSPEN((CSX-CM2)/CSX)-PI*PI/3D0
     .                +DCMPLX(0D0,1D0)*PI*HSCLN((CSX-CM2)/CSX))
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCMW(FS,CM2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSCMW,HSSPEN,CM2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      HSCMW=-1D0/FS*(HSSPEN((FS+CM2)/CM2)-PI*PI/6D0)
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       FUNCTION HSIWZ1(GS,GT,CM1,CM2)
       COMPLEX(kind=8) HSIWZ1,HSD0,HSCMWZ,CM1,CM2
       DOUBLE PRECISION GS,GT
       HSIWZ1=GT/2D0*HSD0(GS,GT,CM1,CM2)-HSCMWZ(GS,CM1,CM2)
       END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
       FUNCTION HSIWZ2(GS,GT,CM1,CM2)
       COMPLEX(kind=8) HSIWZ2,HSD0,HSCMW,HSCMWZ,CM1,CM2,HSCLN,HSFONE
       COMPLEX(kind=8) A,MY1,MY2,C,Y1,Y2,COEFF1,COEFF2,COEFF3
       DOUBLE PRECISION GS,GT,RM1,RM2
       A = DCMPLX(-GT/GS,0D0)
       MY1 = CM1/GS
       MY2 = CM2/GS
       RM1=DSQRT(DREAL(CM1))
       RM2=DSQRT(DREAL(CM2))
       C=SQRT(1D0+(MY1-MY2)*(MY1-MY2)-2D0*(MY1+MY2))
       Y1=(1D0-MY1+MY2+C)/2D0
       Y2=(1D0-MY1+MY2-C)/2D0
       COEFF1=0.25D0/(GS+GT)/(GS+GT)*
     *      ( GS*(GS*GT+2D0*CM1*CM2)
     *       -GT*(GS-CM1)*(GS-CM1) - GT*(GS-CM2)*(GS-CM2)
     *       -2D0*GT*(GS+GT)*(GT+CM1+CM2)                )
       COEFF2=-GT*(GS+2D0*GT+CM1+CM2)/4D0/(GS+GT)/(GS+GT)
       COEFF3=(GS*(GS-CM1-CM2)+2D0*GT*(GS+GT))/4D0/(GS+GT)/(GS+GT)
       HSIWZ2=COEFF1*HSD0(GS,GT,CM1,CM2)
     *     +COEFF2*(HSCMW(GT,CM1)+HSCMW(GT,CM2))
     *     +COEFF3*2D0*HSCMWZ(GS,CM1,CM2)
     *     -0.5D0/(GS+GT)*( 1D0 - HSCLN(-GT/CM1)
     *                      -CM2/(CM1-CM2)*HSCLN(CM2/CM1)
     *                      -HSFONE(GS,RM1,RM2)             )
       END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSIXY(X,Y)
      COMPLEX(kind=8) HSIXY,X,Y,HSSPEN,HSCLN
      HSIXY=HSSPEN((1D0-Y)/(X-Y))-HSSPEN(Y/(Y-X))
     *    + HSCLN((1D0-X)/(Y-X))*(HSCLN(1D0-Y)-HSCLN(X-Y))
     *    - HSCLN(X/(X-Y))*(HSCLN(-Y)-HSCLN(X-Y))
     *    + HSCLN((X-1D0)/X)*HSCLN(X-Y)
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSCMWZ(GS,CM1,CM2)
      DOUBLE PRECISION GS
      COMPLEX(kind=8) HSCMWZ,CM1,CM2,X1,X2,C,HSCLN
      C=CDSQRT((GS-CM1-CM2)*(GS-CM1-CM2)-4D0*CM1*CM2)
      X1=(GS-CM2+CM1+C)/2D0/GS
      X2=(GS-CM2+CM1-C)/2D0/GS
      HSCMWZ=-1D0/GS*HSCLN(X1/(X1-1D0))*HSCLN(X2/(X2-1D0))
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSD0(GS,GT,CM1,CM2)
      COMPLEX(kind=8) HSD0,CM1,CM2,HSSPEN
      COMPLEX(kind=8) A,MY1,MY2,CX,CY,XX,X1X2,Y1Y2,Y1,Y2,X1,X2
      DOUBLE PRECISION GS,GT
      A = DCMPLX(-GT/GS,0D0)
      MY1 = CM1/GS
      MY2 = CM2/GS
      CX=SQRT(1D0+(MY1-MY2)*(MY1-MY2)-2D0*(MY1+MY2)+4D0*MY1*MY2/A)
      CY=SQRT(1D0+(MY1-MY2)*(MY1-MY2)-2D0*(MY1+MY2))
      XX=1D0-MY1+MY2
      X1X2=CM2/GS*(GT+CM1)/GT
      Y1Y2=CM2/GS
      IF (DREAL(XX).GT.0D0) THEN
        X1=(XX+CX)/2D0
        Y1=(XX+CY)/2D0
        X2=X1X2/X1
        Y2=Y1Y2/Y1
        ELSE
        X2=(XX-CX)/2D0
        Y2=(XX-CY)/2D0
        X1=X1X2/X2
        Y1=Y1Y2/Y2
      ENDIF
      HSD0 = 1D0/GS/GT/(X1-X2)*
     *    ( HSSPEN((1D0-X1)/(Y1-X1)) - HSSPEN(-X1/(Y1-X1))
     *     -HSSPEN((1D0-X2)/(Y2-X2)) + HSSPEN(-X2/(Y2-X2))
     *     +HSSPEN((1D0-X1)/(Y2-X1)) - HSSPEN(-X1/(Y2-X1))
     *     -HSSPEN((1D0-X2)/(Y1-X2)) + HSSPEN(-X2/(Y1-X2)) )
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSENUW(T)
C---IR-FINITE PART
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSENUW,HSCLM2,HSCLM4,HSCLN,CMW2,CMZ2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSCBMS/ CMW2,CMZ2
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      HSENUW=DCMPLX(0D0,0D0)
      IF (LPAR(12).GT.0) THEN
      HSENUW=HSENUW+ALP4PI*(HSCLN(MEI2/CMW2)
     1                     +3D0*HSCLM4(T,CMW2,(0D0,0D0),MEI2))
      ENDIF
      IF (LPAR(15).EQ.1) THEN
      HSENUW=HSENUW+ALP4PI*(3D0*(3D0*CW2-1D0)/2D0/SW2
     1           +((2D0*SW2-1D0)/2D0/SW2+3D0*CW2/SW2/SW2)*DLOG(CW2)
     2           +(2D0*SW2-1D0)/4D0/SW2/CW2*HSCLM2(T,CMZ2)
     3           +3D0*CW2/SW2*HSCLM4(T,CMW2,CMZ2,0D0))
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSDUWQ(T)
C---IR-FINITE PART
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSDUWQ,HSCLN,HSCLM1,HSCLM2,HSCLM4,CMW2,CMZ2
      COMPLEX(kind=8) AMBU,AMBD
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSCBMS/ CMW2,CMZ2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      CQU=2D0/3D0
      HSDUWQ=DCMPLX(0D0,0D0)
C
C---PHOTONIC PART
      IF (LPAR(12).EQ.1) THEN
        HSDUWQ=HSDUWQ+ALP4PI*(-HSCLN(CMZ2/MD2)-9D0/2D0
     *                        +3D0*HSCLM4(T,CMW2,(0D0,0D0),MD2))
      ENDIF
      IF (LPAR(14).EQ.1) THEN
        AMBU=HSCLM1(T,MU2)+DLOG(-T/MU2)*DLOG(-T/MU2)
        AMBD=HSCLM1(T,MD2)+DLOG(-T/MD2)*DLOG(-T/MD2)
        HSDUWQ=HSDUWQ+ALP4PI*CQU*(
     *      -(3D0*DLOG(MD/MU)+.5D0*AMBU+.5D0*AMBD)
     *      -SW2/CW2*HSCLM2(T,CMZ2)
     *      +3D0*HSCLM4(T,CMW2,(0D0,0D0),MU2)
     *      -3D0*HSCLM4(T,CMW2,(0D0,0D0),MD2)
     *      +HSCLN(CMZ2/MD2)+9D0/2D0)
      ENDIF
      IF (LPAR(13).EQ.1) THEN
        AMBU=HSCLM1(T,MU2)+DLOG(-T/MU2)*DLOG(-T/MU2)
        AMBD=HSCLM1(T,MD2)+DLOG(-T/MD2)*DLOG(-T/MD2)
        HSDUWQ=HSDUWQ+ALP4PI*CQU*CQU*(
     *             3D0*DLOG(MD/MU)+.5D0*AMBU+.5D0*AMBD
     *             +SW2/CW2*HSCLM2(T,CMZ2))
      ENDIF
C---WEAK PART
      IF (LPAR(15).EQ.1) THEN
        HSDUWQ=HSDUWQ+ALP4PI*(
     *      +(2D0*SW2-1D0)/4D0/SW2/CW2*HSCLM2(T,CMZ2)
     *      +3D0*CW2/SW2*HSCLM4(T,CMW2,CMZ2,0D0)
     *      +3D0/SW2+(1D0/2D0/SW2-3D0*CW2/SW2/SW2)*DLOG(1D0/CW2))
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSDUWA(T)
C---IR-FINITE PART
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSDUWA,HSCLN,HSCLM1,HSCLM2,HSCLM4,CMW2,CMZ2
      COMPLEX(kind=8) AMBU,AMBD
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSCBMS/ CMW2,CMZ2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      HSDUWA=DCMPLX(0D0,0D0)
      CQD=-1D0/3D0
C
C---PHOTONIC PART
      IF (LPAR(12).EQ.1) THEN
        HSDUWA=HSDUWA+ALP4PI*(3D0*HSCLM4(T,CMW2,(0D0,0D0),MU2) )
      ENDIF
      IF (LPAR(14).EQ.1) THEN
        AMBU=HSCLM1(T,MU2)+DLOG(-T/MU2)*DLOG(-T/MU2)
        AMBD=HSCLM1(T,MD2)+DLOG(-T/MD2)*DLOG(-T/MD2)
        HSDUWA=HSDUWA+ALP4PI*CQD*(
     *      +(3D0*DLOG(MD/MU)+.5D0*AMBU+.5D0*AMBD)
     *      +SW2/CW2*HSCLM2(T,CMZ2)
     *      +3D0*HSCLM4(T,CMW2,(0D0,0D0),MU2)
     *      -3D0*HSCLM4(T,CMW2,(0D0,0D0),MD2)
     *      +HSCLN(CMZ2/MD2)+9D0/2D0)
      ENDIF
      IF (LPAR(13).EQ.1) THEN
        AMBU=HSCLM1(T,MU2)+DLOG(-T/MU2)*DLOG(-T/MU2)
        AMBD=HSCLM1(T,MD2)+DLOG(-T/MD2)*DLOG(-T/MD2)
        HSDUWA=HSDUWA+ALP4PI*CQD*CQD*(
     *        3D0*DLOG(MD/MU)+.5D0*AMBU+.5D0*AMBD
     *        +SW2/CW2*HSCLM2(T,CMZ2))
      ENDIF
C---WEAK PART
      IF (LPAR(15).EQ.1) THEN
        HSDUWA=HSDUWA+ALP4PI*(
     *      +(2D0*SW2-1D0)/4D0/SW2/CW2*HSCLM2(T,CMZ2)
     *      +3D0*CW2/SW2*HSCLM4(T,CMW2,CMZ2,0D0)
     *      +3D0/SW2+(1D0/2D0/SW2-3D0*CW2/SW2/SW2)*DLOG(1D0/CW2))
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   INTERFACE FOR CALLS OF PARTON DISTRIBUTIONS FROM HERACLES
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSPVER(X,Q2)
      DOUBLE PRECISION QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      DOUBLE PRECISION X,Q2
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      DIMENSION XPQ(-6:6)
C---
      RX=REAL(X)
      RQ2=REAL(Q2)
C---
C...USE PYSTFU ROUTINES via LYSTFU from LEPTO 6.5
      CALL LYSTFU(2212,RX,RQ2,XPQ)
      QU=XPQ(2)/X
      QD=XPQ(1)/X
      QS=XPQ(3)/X
      QC=XPQ(4)/X
      QB=XPQ(5)/X
      QT=XPQ(6)/X
      QBU=XPQ(-2)/X
      QBD=XPQ(-1)/X
      QBS=XPQ(-3)/X
      QBC=XPQ(-4)/X
      QBB=XPQ(-5)/X
      QBT=XPQ(-6)/X

      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   INTERFACE FOR CALLS OF POLARIZED PARTON DISTRIBUTIONS 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSDPVR(X,Q2)
      DOUBLE PRECISION X,Q2
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      DOUBLE PRECISION DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD,IDPVR
      DOUBLE PRECISION XPPDF,GRAD,UUB,DDB,DGLU
      DIMENSION XPPDF(-3:3),GRAD(-3:3,11)
      LOGICAL FIRST
      DATA FIRST /.TRUE./

      ISET=IDPVR/100
      MODE=MOD(IDPVR,100)
C   USE DSSV
      IF (ISET.EQ.1) THEN
        IF (FIRST) THEN
          CALL DSSVINI
          FIRST=.FALSE.
        ENDIF
        CALL DSSVFIT(X,Q2,DQU,DQD,DQBU,DQBD,DQS,DGLU)
        DQU=(DQU+DQBU)
        DQD=(DQD+DQBD)
C   USE DNS
      ELSEIF (ISET.EQ.2) THEN
        IF (FIRST) THEN
          CALL INIDNS
          FIRST=.FALSE.
        ENDIF
        CALL POLFIT(MODE,X,Q2,DQU,DQD,DQBU,DQBD,DQS,DGLU,G1P,G1N)
        DQU=(DQU+DQBU)
        DQD=(DQD+DQBD)
C   USE DS
      ELSEIF (ISET.EQ.3) THEN
        IF (FIRST) THEN
          CALL INIDS
          FIRST=.FALSE.
        ENDIF
        CALL DS2000(MODE,X,Q2,DQU,DQD,DQBU,DQBD,DQS,DGLU,G1P,G1N,G1D)
        DQU=(DQU+DQBU)
        DQD=(DQD+DQBD)
C   USE GSLO
      ELSEIF (ISET.EQ.4) THEN
        IF (FIRST) THEN
          CALL POLINI
          FIRST=.FALSE.
        ENDIF
        CALL polpar(MODE,X,Q2,DQU,DQD,DGLU,DQBU,DQS)
        DQBD=DQBU
        DQU=(DQU+DQBU)
        DQD=(DQD+DQBD)
C   USE GSNLO
      ELSEIF (ISET.EQ.5) THEN
        IF (FIRST) THEN
          CALL NLOINI
          FIRST=.FALSE.
        ENDIF
        CALL polnlo(MODE,X,Q2,DQU,DQD,DGLU,DQBU,DQBD,DQS)
        DQU=(DQU+DQBU)
        DQD=(DQD+DQBD)
C   USE BB
      ELSEIF (ISET.EQ.6) THEN
        CALL PPDF(MODE,X,Q2,DQU,DUV,DQD,DDV,DGLU,DGL,
     1           DQBU,DQB,G1P,DG1P,G1N,DG1N)
        DQBD=DQBU
        DQS=DQBU
        DQU=(DQU+DQBU)
        DQD=(DQD+DQBD)
C   USE AAC
      ELSEIF (ISET.EQ.7) THEN
        CALL AAC08PDF(Q2,X,MODE,XPPDF,GRAD)
        DQU=XPPDF(1)
        DQD=XPPDF(2)
        DQS=XPPDF(3)
        DQBU=XPPDF(-1)
        DQBD=XPPDF(-2)
        DQBS=XPPDF(-3)
C   USE LSS
      ELSEIF (ISET.EQ.8) THEN
        CALL LSS2010(MODE,X,Q2,UUB,DDB,DQU,DQD,DQBU,DQBD,DQS,DGLU)
C   USE GRSV
      ELSEIF (ISET.EQ.9) THEN
        CALL PARPOL(MODE,X,Q2,DQU,DQD,DQBU,DQBD,DQS,DGLU,G1P,G1N)
      ELSE
        WRITE(LUNOUT,*) ' INVALID POLARIZED PARTON CODE'
        WRITE(LUNOUT,*) ' IDPVR = ',IDPVR
        WRITE(LUNOUT,*) ' ***** EXECUTION STOPPED ***** '
        STOP
      ENDIF

      DQU=DQU/X
      DQD=DQD/X
      DQS=DQS/X
      DQC=0D0
      DQB=0D0
      DQT=0D0
c      DGLU=DGLU/X
      DQBU=DQBU/X
      DQBD=DQBD/X
      DQBS=DQS
      DQBC=0D0
      DQBB=0D0
      DQBT=0D0

      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   ISOSPIN ROTATION FOR NUCLEAR TARGETS (NEEDED FOR IPDFOP=1) 
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSIROT
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSNUCL/ HNA,HNZ,INUMOD

      IF (INT(HNA).EQ.1.AND.INT(HNZ).EQ.1) RETURN
      IF (INT(HNA).EQ.1.AND.INT(HNZ).EQ.0) THEN
C---NEUTRON TARGET
C---unpolarized
        QDUMD=QD
        QDUMU=QU 
        QD=QDUMU
        QU=QDUMD
        QDUMD=QBU
        QDUMU=QBD
        QBU=QDUMU
        QBD=QDUMD
C---polarized
        QDUMD=DQD
        QDUMU=DQU 
        DQD=QDUMU
        DQU=QDUMD
        QDUMD=DQBU
        QDUMU=DQBD
        DQBU=QDUMU
        DQBD=QDUMD
      ELSE
C---OTHER TARGETS
        IF (INUMOD.EQ.1000) THEN
C---SIMPLE SHADOWING MODEL: FIRST CALCULATE DEUTERON MIXING
C   PER NUCLEON
C---unpolarized
          HNRAT=HSNRAT(X)
          QDUMD=QD
          QDUMU=QU
          QD=HNRAT*(QDUMU+QDUMD)/2D0
          QU=HNRAT*QD
          QDUMD=QBU
          QDUMU=QBD
          QBD=HNRAT*(QDUMU+QDUMD)/2D0
          QBU=HNRAT*QBD
C---polarized
          QDUMD=DQD
          QDUMU=DQU
          DQD=HNRAT*(QDUMU+QDUMD)/2D0
          DQU=HNRAT*DQD
          QDUMD=DQBU
          QDUMU=DQBD
          DQBD=HNRAT*(QDUMU+QDUMD)/2D0
          DQBU=HNRAT*DQBD
        ELSE 
C---WEIGHTED AVERAGE OF UP AND DOWN, CORRECTION INCLUDED IN PDFs
C---unpolarized
          QDUMD=QD
          QDUMU=QU
          QU=(HNZ*QDUMU+(HNA-HNZ)*QDUMD)/HNA
          QD=((HNA-HNZ)*QDUMU+HNZ*QDUMD)/HNA
          QDUMD=QBU
          QDUMU=QBD
          QBU=(HNZ*QDUMD+(HNA-HNZ)*QDUMU)/HNA
          QBD=((HNA-HNZ)*QDUMD+HNZ*QDUMU)/HNA
C---polarized
          QDUMD=DQD
          QDUMU=DQU
          DQU=(HNZ*QDUMU+(HNA-HNZ)*QDUMD)/HNA
          DQD=((HNA-HNZ)*QDUMU+HNZ*QDUMD)/HNA
          QDUMD=DQBU
          QDUMU=DQBD
          DQBU=(HNZ*QDUMD+(HNA-HNZ)*QDUMU)/HNA
          DQBD=((HNA-HNZ)*QDUMD+HNZ*QDUMU)/HNA
        ENDIF
      ENDIF

      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
************************************************************************
C
C   SUBROUTINES FOR THE INCLUSION OF THE LONGITUDINAL STRUCTURE
C   FUNCTION IN HERACLES, since VERSION 4.5
C
C   TAKEN / MODIFIED FROM LEPTO 6.1 BY G.INGELMAN
C
************************************************************************

      SUBROUTINE HSLUFL(XA,Q2A,F2EM,FL)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C...COMMON BLOCKS FROM HERACLES (NOTE: NAMES HAVE PARTLY BEEN CHANGED)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSCUTS/ XMINH,XMAXH,Q2MINH,Q2MAXH,YMINH,YMAXH,WMINH,GMIN
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMAX,NPYMIN
      REAL            PYSTOP,PYSLAM
C...COMMON BLOCKS FROM LEPTO
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      REAL            CUT            ,PARL    ,X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      REAL            PARI    ,EWQC       ,QC   ,ZL     ,ZQ     ,PQ
      COMMON /FLGRID/ NFX,NFQ,XR(2),QR(2),FLQT(41,16),FLGT(41,16),
     &FLMT(41,16)
      REAL                    XR   ,QR   ,FLQT       ,FLGT       ,
     &FLMT
      REAL FLQ,FLG,FLM
      LOGICAL LFIRST,LFIRS1
      DATA LFIRST,LFIRS1 /2*.TRUE./

      IF (LFIRS1) THEN
      LFIRS1=.FALSE.
      QC(1)=-.33333
      QC(2)=.66667
      QC(3)=-.33333
      QC(4)=.66667
      QC(5)=-.33333
      QC(6)=.66667
      QC(7)=-.33333
      QC(8)=.66667
C...ONLY WARNINGS PRINTED, EXECUTION NOT STOPPED IF ERROR IN LEPTO
ckc..
      LST3=LST(3)
      LST(3)=1
C...F_L INCLUDED
      LST(11)=IFLOPT
C...To be sure that LNSTRF is initialized
      CALL HSPVER(0.1D0,100D0)
      ENDIF

      FLQ=0D0
      FLG=0D0
      FLM=0D0
      FLT=0D0
      X=SNGL(XA)
      Q2=SNGL(Q2A)
      IF (LQCD.EQ.1.OR.LTM.EQ.1) THEN
C...F_L FROM INTERPOLATION
        IF (LFIRST) THEN
C...INITIALIZATION OF GRID
        LFIRST=.FALSE.
        PARL(21)=SNGL(SP-MEI2-MPRO2)
        NFX=41
        NFQ=16
        DO 11 I=1,NFX
        DO 12 J=1,NFQ
        FLQT(I,J)=0.
        FLGT(I,J)=0.
 12     CONTINUE
 11     CONTINUE
        CALL FLTABL
        ENDIF
        CALL FLIPOL(FLQ,FLG,FLM)
      ELSEIF (LQCD.EQ.2.OR.LTM.EQ.2) THEN
C...F_L FROM INTEGRATION EVENT-BY-EVENT
        IF (LFIRST) THEN
C...INITIALIZATION
        LFIRST=.FALSE.
C...PROTON TARGET
        PARL(1)=1.
        PARL(2)=1.
        PARI(11)=(PARL(1)-PARL(2))/PARL(1)
C...KINEMATIC LIMITS
        XR(1)=SNGL(XMINH)
        XR(2)=SNGL(XMAXH)
        QR(1)=SNGL(Q2MINH)
        ENDIF
        CALL FLINTG(FLQ,FLG,FLM)
      ELSE
        RETURN
      ENDIF

      IF (LHT.GE.1) FLT=8D0*PARL(19)/Q2A*F2EM
      FL=DBLE(FLQ+FLG+FLM)+FLT
      LST(3)=LST3
      RETURN
      END

C*********************************************************************
C...translate input from ALFAS input

      SUBROUTINE DIALFS
      COMMON /HSALFS/ PAR111,PAR112,PARL11,PARL19,MST111,MST115
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U

      MSTU(111)=MST111
      MSTU(115)=MST115
      PARU(111)=PAR111
      PARU(112)=PAR112
      RETURN
      END

C*********************************************************************
C...Translate input from FLONG input flag

      SUBROUTINE DIFLOP
      COMMON /HSALFS/ PAR111,PAR112,PARL11,PARL19,MST111,MST115
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U

      PARL(11)=PARL11
      PARL(19)=PARL19
      RETURN
      END

C
C***********************************************************************
C
C **********************************************************************
C

      SUBROUTINE FLTABL

C...Integrates the longitudinal structure function, store on grid
C...in x, Q**2.
C   (changed: 12.01.94, HS)
C   (changed: 29.03.95, HS)
chs..use the true value to determine limits of grid also in the case of
c....initial state radiation with reduced effective S stored on PARL(21)

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
C.HS FOR USE IN HERACLES+DJANGO since VERSION 4.5:
C... TRANSFER KINEMATIC LIMITS FROM /HSCUTS/
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      DOUBLE PRECISION MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      DOUBLE PRECISION SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMINH,XMAXH,Q2MINH,Q2MAXH,YMINH,YMAXH,WMINH,GMIN
      DOUBLE PRECISION XMINH,XMAXH,Q2MINH,Q2MAXH,YMINH,YMAXH,WMINH,GMIN
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
C.HS  COMMON /LINTRL/ PSAVE(3,4,5),KSAVE(4),XMIN,XMAX,YMIN,YMAX,
C.HS &Q2MIN,Q2MAX,W2MIN,W2MAX,ILEP,INU,IG,IZ
      COMMON /LINTEG/ NTOT,NPASS
      COMMON /FLGRID/ NFX,NFQ,XR(2),QR(2),FLQT(41,16),FLGT(41,16),
     &FLMT(41,16)
      EXTERNAL FLQINT,FLGINT,FLTINT
C.HS(12/1/94)
C-HS(10.01.94)
C...Use always fixed center-of-mass energy as upper limit
      PARL21=PARL(21)
      PARL(21)=SNGL(SP-MEI2-MPRO2)
      DO 1 IO=1,15
        IF (INT3(IO).NE.0.OR.ISAM3(IO).NE.0) GOTO 2
    1 CONTINUE
      XMAX=SNGL(XMAXH)
      Q2MIN=SNGL(Q2MINH)
      GOTO 3
    2 CONTINUE
      XMAX=0.999
      Q2MIN=1.0
    3 CONTINUE
      XMIN=SNGL(XMINH)

      LQCD=MOD(LST(11),10)
      LTM=MOD(LST(11)/10,10)
      LHT=LST(11)/100
      IF(LST(3).GE.3) WRITE(LUNOUT,1000) LST(11),LQCD,LTM,LHT
      IF(LQCD.LT.1.AND.LTM.LT.1) GOTO 900
      CALL LTIMEX(T1)
      DO 11 IX=1,NFX
      DO 12 IQ=1,NFQ
      FLQT(IX,IQ)=0.
      FLGT(IX,IQ)=0.
      FLMT(IX,IQ)=0.
 12   CONTINUE
 11   CONTINUE
      QR(1)=Q2MIN
      XR(1)=XMIN
      XR(2)=XMAX
      DO 500 IX=1,NFX
      X=10**(ALOG10(XR(1))+(ALOG10(XR(2))-ALOG10(XR(1)))*(IX-1)/(NFX-1))
      QR(2)=X*PARL(21)
      IF(QR(1).GT.QR(2)) GOTO 500
      LQ=0
      DO 400 IQ=1,NFQ
      Q2=10**(ALOG10(QR(1))+(ALOG10(QR(2))-ALOG10(QR(1)))*
     &(IQ-1)/(NFQ-1))
CTEST IF(LQ.GT.0) GOTO 500
      IF(Q2.GT.PARL(21)) LQ=LQ+1
      Y=Q2/(PARL(21)*X)
      IF(Y.LT.0.0.OR.Y.GT.1.0) LQ=LQ+1
      PARL(25)=ULALPS(Q2)
      IF(LQCD.EQ.1) THEN
C...Quark part.
        ACCUR=PARL(11)
        IT=0
  100   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLQINT,EPS,FLQ)
        IF(FLQ.LT.1) THEN
          ACCUR=FLQ*PARL(11)
          IF(IT.LT.2) GOTO 100
        ENDIF
        FLQT(IX,IQ)=FLQ
C...Gluon part.
        ACCUR=PARL(11)
        IT=0
  200   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLGINT,EPS,FLG)
        IF(FLG.LT.1.) THEN
          ACCUR=FLG*PARL(11)
          IF(IT.LT.2) GOTO 200
        ENDIF
        FLGT(IX,IQ)=FLG
      ENDIF
      IF(LTM.EQ.1) THEN
C...Target mass  part.
        ACCUR=PARL(11)
        IT=0
  300   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLTINT,EPS,FLM)
        IF(FLM.LT.1) THEN
          ACCUR=FLM*PARL(11)
          IF(IT.LT.2) GOTO 300
        ENDIF
        FLMT(IX,IQ)=FLM
      ENDIF
  400 CONTINUE
  500 CONTINUE
      CONTINUE
      CALL LTIMEX(T2)
      IF(LST(3).GE.3) WRITE(LUNOUT,1100) T2-T1
  900 CONTINUE
      PARL(21)=PARL21
      RETURN

 1000 FORMAT(' Initialisation for FL; QCD, target mass, higher twist: ',
     &/,' LST(11) =',I5,' --> LQCD, LTM, LHT =',3I3)
 1100 FORMAT(' FL integrations performed if LQCD=1 and/or LTM=1, ',
     &'results on grid.'/,' Time for FL integrations is ',F7.1,' sec.')
      END

C **********************************************************************

      SUBROUTINE FLIPOL(FLQ,FLG,FLM)

C...QCD and target mass contributions to longitudinal structure function
C...from interpolation on x,Q2 grid.
chs..transfer true center-of-mass energy from HERACLES, use the true
c....value to determine limits of grid also in the case of initial
c....state radiation with reduced effective S stored on PARL(21)

      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      DOUBLE PRECISION MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      DOUBLE PRECISION SP,EELE,PELE,EPRO,PPRO
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /FLGRID/ NFX,NFQ,XR(2),QR(2),FLQT(41,16),FLGT(41,16),
     &FLMT(41,16)
      DATA NOUT/0/,NWARN/10/

      LQCD=MOD(LST(11),10)
      LTM=MOD(LST(11)/10,10)
      LHT=LST(11)/100
      XP=X
      Q2P=Q2

C...Use always fixed center-of-mass energy as upper limit
      PARL21=PARL(21)
      PARL(21)=SNGL(SP-MEI2-MPRO2)
C...NOTE: tiny mismatch between present x-value and those on grid.
      QR(2)=X*PARL(21)
      IF(QR(1).GT.QR(2)) GOTO 900
      IF(X.LT.XR(1).OR.X.GT.XR(2).OR.
     &Q2.LT.QR(1).OR.Q2.GT.QR(2)) THEN
C...x and/or Q2 outside grid limits, write warning for first NWARN cases
        IF(LST(2).GE.0) THEN
          NOUT=NOUT+1
          IF(LST(3).GE.1.AND.NOUT.LE.NWARN)
     &       WRITE(LUNOUT,1000) X,Q2,NWARN
        ENDIF
        IF(X.LT.XR(1)) XP=XR(1)
        IF(X.GT.XR(2)) XP=XR(2)
        IF(Q2.LT.QR(1)) Q2P=QR(1)
        IF(Q2.GT.QR(2)) Q2P=QR(2)
      ENDIF

      IX=INT((ALOG10(XP)-ALOG10(XR(1)))/
     &       (ALOG10(XR(2))-ALOG10(XR(1)))*(NFX-1)+1)
      IQ=INT((ALOG10(Q2P)-ALOG10(QR(1)))/
     &       (ALOG10(QR(2))-ALOG10(QR(1)))*(NFQ-1)+1)
      IX=MIN(IX,NFX-1)
      IQ=MIN(IQ,NFQ-1)
      Q2L=10**(ALOG10(QR(1))+(ALOG10(QR(2))-ALOG10(QR(1)))*
     &(IQ-1)/(NFQ-1))
      Q2H=10**(ALOG10(QR(1))+(ALOG10(QR(2))-ALOG10(QR(1)))*
     &(IQ  )/(NFQ-1))
      XL=10**(ALOG10(XR(1))+(ALOG10(XR(2))-ALOG10(XR(1)))*
     &(IX-1)/(NFX-1))
      XH=10**(ALOG10(XR(1))+(ALOG10(XR(2))-ALOG10(XR(1)))*
     &(IX  )/(NFX-1))
      QD=1D0
      IF (Q2H.NE.Q2L) QD=(Q2P-Q2L)/(Q2H-Q2L)
      XD=(XP-XL)/(XH-XL)

      IF(LQCD.EQ.1) THEN
        X1P=(FLQT(IX+1,IQ)-FLQT(IX,IQ))*XD+FLQT(IX,IQ)
        X2P=(FLQT(IX+1,IQ+1)-FLQT(IX,IQ+1))*XD+FLQT(IX,IQ+1)
        FLQ=(X2P-X1P)*QD+X1P
        X1P=(FLGT(IX+1,IQ)-FLGT(IX,IQ))*XD+FLGT(IX,IQ)
        X2P=(FLGT(IX+1,IQ+1)-FLGT(IX,IQ+1))*XD+FLGT(IX,IQ+1)
        FLG=(X2P-X1P)*QD+X1P
      ENDIF
      IF(LTM.EQ.1) THEN
        X1P=(FLMT(IX+1,IQ)-FLMT(IX,IQ))*XD+FLMT(IX,IQ)
        X2P=(FLMT(IX+1,IQ+1)-FLMT(IX,IQ+1))*XD+FLMT(IX,IQ+1)
        FLM=(X2P-X1P)*QD+X1P
      ENDIF
  900 PARL(21)=PARL21
      RETURN
ckc..format for 'x=' and 'Q2=' changed
 1000 FORMAT(' Warning: x=',E9.3,' or Q2=',E9.3,' outside grid,',
     &' for FL interpolation',/,10X,'value on grid limit used.',
     &' Only first',I5,' warnings printed.',/)
      END

C **********************************************************************

      SUBROUTINE FLINTG(CFLQ,CFLG,CFLM)

C...Event-by-event calculation of contribution to longitudinal
C...structure function from QCD and target mass effects.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTEG/ NTOT,NPASS
      EXTERNAL FLQINT,FLGINT,FLTINT

      LQCD=MOD(LST(11),10)
      LTM=MOD(LST(11)/10,10)
      LHT=LST(11)/100
      PARL(25)=ULALPS(Q2)
      IF(LQCD.EQ.2) THEN
C...FL from QCD, quark and gluon contributions.
        ACCUR=PARL(11)
        IT=0
  100   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLQINT,EPS,CFLQ)
        IF(CFLQ.LT.1) THEN
          ACCUR=CFLQ*PARL(11)
          IF(IT.LT.2) GOTO 100
        ENDIF
        ACCUR=PARL(11)
        IT=0
  200   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLGINT,EPS,CFLG)
        IF(CFLG.LT.1.) THEN
          ACCUR=CFLG*PARL(11)
          IF(IT.LT.2) GOTO 200
        ENDIF
      ENDIF
      IF(LTM.EQ.2) THEN
        ACCUR=PARL(11)
        IT=0
  300   IT=IT+1
        NTOT=0
        NPASS=0
        EPS=ACCUR
        CALL GADAP(X,1.,FLTINT,EPS,CFLM)
        IF(CFLM.LT.1.) THEN
          ACCUR=CFLM*PARL(11)
          IF(IT.LT.2) GOTO 300
        ENDIF
      ENDIF

      RETURN
      END

C **********************************************************************

      FUNCTION FLQINT(Z)

C...Quark contribution integrand to QCD longitudinal structure function.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTEG/ NTOT,NPASS
      DIMENSION XPQ(-6:6)
      DATA PI/3.14159/
      NTOT=NTOT+1
      CALL LNSTRF(Z,Q2,XPQ)
      FLQINT=0.
      DO 100 I=-LST(12),LST(12)
      IF(I.EQ.0) GOTO 100
      FLQINT=FLQINT+QC(IABS(I))**2*XPQ(I)
  100 CONTINUE
      FLQINT=4./3.*PARL(25)/PI*(X/Z)**2*FLQINT/Z
      NPASS=NPASS+1

      RETURN
      END

C **********************************************************************

      FUNCTION FLGINT(Z)
 
C...Gluon contribution integrand to QCD longitudinal structure function.

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTEG/ NTOT,NPASS
      DIMENSION XPQ(-6:6)
      DATA PI/3.14159/
      NTOT=NTOT+1
      CALL LNSTRF(Z,Q2,XPQ)
      FLGINT=20./9.*PARL(25)/PI*(X/Z)**2*(1.-X/Z)/Z*XPQ(0)
      NPASS=NPASS+1

      RETURN
      END

C **********************************************************************

      FUNCTION FLTINT(Z)

C...Integrand for target mass correction contribution to
C...quark longitudinal structure function

      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
      COMMON /LINTER/ PARI(40),EWQC(2,2,8),QC(8),ZL(2,4),ZQ(2,8),PQ(17)
      COMMON /LINTEG/ NTOT,NPASS
      DIMENSION XPQ(-6:6)
      DATA PM2/0.8804/
      NTOT=NTOT+1
      CALL LNSTRF(Z,Q2,XPQ)
      FLTINT=0.
      DO 100 I=-LST(12),LST(12)
      IF(I.EQ.0) GOTO 100
      FLTINT=FLTINT+QC(IABS(I))**2*XPQ(I)
  100 CONTINUE
      FLTINT=4.*PM2/Q2*(X/Z)**2*X*FLTINT
      NPASS=NPASS+1

      RETURN
      END

C#######################################################################
C
C   One- and two-dimensional adaptive Gaussian integration routines.
C
C **********************************************************************

      SUBROUTINE GADAP(A0,B0,F,EPS,SUM)
C
C   PURPOSE           - INTEGRATE A FUNCTION F(X)
C   METHOD            - ADAPTIVE GAUSSIAN
C   USAGE             - CALL GADAP(A0,B0,F,EPS,SUM)
C   PARAMETERS  A0    - LOWER LIMIT (INPUT,REAL)
C               B0    - UPPER LIMIT (INPUT,REAL)
C               F     - FUNCTION F(X) TO BE INTEGRATED. MUST BE
C                       SUPPLIED BY THE USER. (INPUT,REAL FUNCTION)
C               EPS   - DESIRED RELATIVE ACCURACY. IF SUM IS SMALL EPS
C                       WILL BE ABSOLUTE ACCURACY INSTEAD. (INPUT,REAL)
C               SUM   - CALCULATED VALUE FOR THE INTEGRAL (OUTPUT,REAL)
C   PRECISION         - SINGLE
C   REQ'D PROG'S      - F
C   AUTHOR            - T. JOHANSSON, LUND UNIV. COMPUTER CENTER, 1973
C   REFERENCE(S)      - THE AUSTRALIAN COMPUTER JOURNAL,3 P.126 AUG. -71
C
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON/GADAP1/ NUM,IFU
      EXTERNAL F
      DIMENSION A(300),B(300),F1(300),F2(300),F3(300),S(300),N(300)
    1 FORMAT(16H GADAP:I TOO BIG)
C     DSUM(F1F,F2F,F3F,AA,BB)=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)
C     is now external: GADSUM
      IF(EPS.LT.1.0E-8) EPS=1.0E-8
      RED=1.3
      L=1
      I=1
      SUM=0.
      C=SQRT(15.)/5.
      A(1)=A0
      B(1)=B0
      F1(1)=F(0.5*(1+C)*A0+0.5*(1-C)*B0)
      F2(1)=F(0.5*(A0+B0))
      F3(1)=F(0.5*(1-C)*A0+0.5*(1+C)*B0)
      IFU=3
      S(1)=GADSUM(F1(1),F2(1),F3(1),A0,B0)
  100 CONTINUE
      L=L+1
      N(L)=3
      EPS=EPS*RED
      A(I+1)=A(I)+C*(B(I)-A(I))
      B(I+1)=B(I)
      A(I+2)=A(I)+B(I)-A(I+1)
      B(I+2)=A(I+1)
      A(I+3)=A(I)
      B(I+3)=A(I+2)
      W1=A(I)+(B(I)-A(I))/5.
      U2=2.*W1-(A(I)+A(I+2))/2.
      F1(I+1)=F(A(I)+B(I)-W1)
      F2(I+1)=F3(I)
      F3(I+1)=F(B(I)-A(I+2)+W1)
      F1(I+2)=F(U2)
      F2(I+2)=F2(I)
      F3(I+2)=F(B(I+2)+A(I+2)-U2)
      F1(I+3)=F(A(I)+A(I+2)-W1)
      F2(I+3)=F1(I)
      F3(I+3)=F(W1)
      IFU=IFU+6
      IF(IFU.GT.5000) GOTO 130
      S(I+1)=GADSUM(F1(I+1),F2(I+1),F3(I+1),A(I+1),B(I+1))
      S(I+2)=GADSUM(F1(I+2),F2(I+2),F3(I+2),A(I+2),B(I+2))
      S(I+3)=GADSUM(F1(I+3),F2(I+3),F3(I+3),A(I+3),B(I+3))
      SS=S(I+1)+S(I+2)+S(I+3)
      I=I+3
      IF(I.GT.300)GOTO 120
      SOLD=S(I-3)
      IF(ABS(SOLD-SS).GT.EPS*(1.+ABS(SS))/2.) GOTO 100
      SUM=SUM+SS
      I=I-4
      N(L)=0
      L=L-1
  110 CONTINUE
      IF(L.EQ.1) GOTO 130
      N(L)=N(L)-1
      EPS=EPS/RED
      IF(N(L).NE.0) GOTO 100
      I=I-1
      L=L-1
      GOTO 110
  120 WRITE(LUNOUT,1)
  130 RETURN
      END

      FUNCTION GADSUM(F1F,F2F,F3F,AA,BB)
      GADSUM=5./18.*(BB-AA)*(F1F+1.6*F2F+F3F)
      RETURN
      END


C######################################################################
C
C   Various routines to give structure function parametrizations.
C
C ********************************************************************

c...hs taken from LEPTO 6.5 (modified)
      SUBROUTINE LYSTFU(KF,X,Q2,XPQ)

C...Interface to PYSTFU in PYTHIA 5.7 to get parton density distributions,
C...i.e. momentum weighted probability distributions xq(x,Q2), xg(x,Q2).
C...Also gives intrinsic charm and beauty distributions.
      COMMON /LEPTOU/ CUT(14),LST(40),PARL(30),XLP,YLP,W2LP,Q2LP,ULP
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
c...hs: upgrade to Pythia6
      double precision parp,pari
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON /ARSTRF/ KFSAVE(2),XSAVE(2),XQ2SAV(2),XPQSAV(2,-6:6)
c...hs
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMAX,NPYMIN
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD,IDPVR
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      DOUBLE PRECISION POLARI
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      DOUBLE PRECISION SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      DOUBLE PRECISION MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      DIMENSION XPQ(-6:6)
      double precision dxpqy(-25:25)
      DOUBLE PRECISION HSLOQS,DQ2,DX
C...hs
      LOGICAL LFIRST
      DATA LFIRST/.TRUE./

      IF (LFIRST) THEN
        LFIRST=.FALSE.
C...Transfer HERACLES switches
        MSTU(11)=LUNOUT
        LST(15)=ICODE
        LST(16)=ILIB
        LST(12)=NPYMAX
        PMAS(6,1)=PYSTOP
C...Initialize PYTHIA for parton densities.
        IF(LST(15).GT.0) THEN
C...Set switches and parameters for parton densities in PYSTFU.
          MSTP(51)=LST(15)
          MSTP(52)=LST(16)
          MSTP(58)=LST(12)
        ENDIF
        ROOTS=REAL(SQRT(SP-MEI2-MPRO2))
        IF (LLEPT.LE.0) THEN
          CALL PYINIT('NONE','e-','p',ROOTS)
         ELSE
          CALL PYINIT('NONE','e+','p',ROOTS)
        ENDIF
C...Reset parameters from input
        CALL DIALFS
        CALL DIFLOP
        MSTU(112)=NPYMAX
        PARL(26)=REAL(PARP(1))
      ENDIF

C...Reset arrays etc.
      DO 100 KFL=-6,6
      XPQ(KFL)=0.0
      XPQSAV(1,KFL)=0.
 100  CONTINUE
      XSAVE(1)=X
      XQ2SAV(1)=Q2
      KFSAVE(1)=KF
C...Check x and particle species.
      IF(X.LE.0..OR.X.GE.1.) THEN
        WRITE(MSTU(11),5000) X
        RETURN
      ENDIF

      IF(LST(15).EQ.-4.OR.LST(15).EQ.-5) THEN
C...Intrinsic charm/bottom quark distribution in the proton...
        IF(Q2.LT.1.) RETURN
C...from Phys. Lett 93B (1980) 451,
C...Amount of intrinsic charm PARL(12)=BETA^2
        XPQ(4)=X**3*1800.*PARL(12)*
     &         ((1.-X)/3.*(1.+10.*X+X**2)+2.*X*(1.+X)*LOG(X))
C...plus first order QCD-correction parametrized with polynomia
        IF(X.LT.0.9) THEN
          XCORR=0.22024E-1*X-0.77833E-1*X**2-0.47292*X**3+
     &          2.104*X**4-2.1698*X**5-0.84891*X**6+1.8882*X**7+
     &          0.8989*X**8-2.1072*X**9+0.76351*X**10
        ELSE
          XCORR=-1.
        ENDIF
C...and a Q2 dependence on that
CJR        XCORR=1.125*XCORR*0.190424*EXP(1.15*LOG(LOG(Q2)))
        IF(Q2.GT.1) THEN
           XCORR=1.125*XCORR*0.190424*EXP(1.15*LOG(LOG(Q2)))
        ELSE
           XCORR=1.125*XCORR*0.190424
        ENDIF
C...smooth cut-off of the structure function !
        XPQ(4)=MAX(XPQ(4)+XCORR,XPQ(4)/Q2)
        XPQ(-4)=XPQ(4)
        IF(LST(15).EQ.-5) THEN
C...Intrinsic bottom assumed to have the same shape as zeroth
C...approximation but suppressed by (mc/mb)**2=0.1 approximately
          XPQ(5)=XPQ(4)*0.1
          XPQ(-5)=XPQ(5)
          XPQ(4)=0.
          XPQ(-4)=0.
        ENDIF
      ELSE
C...Parton densities from PYSTFU in PYTHIA 6
        DQ2=DBLE(Q2)
        DX=DBLE(X)
        CALL PYPDFU(KF,DX,DQ2,DXPQY)
        DO 110 KFL=-6,6
          XPQ(KFL)=REAL(DXPQY(KFL))
 110    CONTINUE
      ENDIF

C...H.S.
C...APPLY LOW Q2 SUPPRESSION
      IF (ILQMOD.EQ.1) THEN
        DQ2=DBLE(Q2)
        DX=DBLE(X)
        DCORRP=REAL(HSLOQS(DQ2))
        DO 200 IF=1,6
          XPQ(IF)=XPQ(IF)*DCORRP
          XPQ(-IF)=XPQ(-IF)*DCORRP
 200    CONTINUE
      ENDIF

      DO 120 KFL=-6,6
        XPQSAV(1,KFL)=XPQ(KFL)
 120  CONTINUE

C...Formats for error printouts.
 5000 FORMAT(' Error in LYSTFU: x =',1P,E12.4,' outside physical range')

      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSLOQS(Q2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      HSLOQS=1D0-DEXP(-3.37D0*Q2)
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   H. SPIESBERGER 22.03.91
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   SUBROUTINE FOR STRUCTURE FUNCTIONS
C   CALCULATED FROM PARTON DISTRIBUTION FUNCTIONS
C   UPDATED: 28.08.92 FOR LOW Q2
C   UPDATED: 13.08.93 FOR LOW FL
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   FOR CONVENTIONS OF THE USE OF PARTON DISTRIBUTION FUNCTIONS
C   AND STRUCTURE FUNCTIONS AT LOW Q2 SEE THE MANUAL
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSWPDF
C---CHECK CONSISTENCY OF PARTON DISTRIBUTION OR STRUCTURE FUNCTION
C   INPUT AND WRITE CHOSEN OPTIONS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD,IDPVR
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMAX,NPYMIN
      REAL            PYSTOP,PYSLAM
      COMMON /HSALFS/ PAR111,PAR112,PARL11,PARL19,MST111,MST115
      REAL            PAR111,PAR112,PARL11,PARL19
      INTEGER         MST111,MST115
      DOUBLE PRECISION VALUE(20)
c      CHARACTER*20 PARM(20)
      CHARACTER PARM*20
      DIMENSION PARM(20)
      DATA VALUE/20*0D0/,PARM/20*' '/
c
      CHARACTER CHPSET*2
      CHARACTER LHAPRM*20
      CHARACTER LHAPTH*232
      CHARACTER LHAPATH*232
      COMMON /LHAPDFC/ LHAPATH
      SAVE /LHAPDFC/
C
C---OPTION FOR STRUCTURE FUNCTION INPUT
C---CHECK CONSISTENCY WITH OTHER INPUT DEFINITIONS
      IOPCHQ=0
      IF (IPDFOP.EQ.0.OR.IFLOPT.GE.1) THEN
        DO 3 IL=13,16
        IF (LPAR(IL).NE.0) THEN
          LPAR(IL)=0
          IOPCHQ=1
        ENDIF
    3   CONTINUE
        IF (IPDFOP.EQ.0)
     *  WRITE(LUNOUT,'(10X,A,I4)')
     *        ' STRUCTURE FUNCTION INPUT: IPDFOP = ',IPDFOP
        IF (IOPCHQ.NE.0) THEN
          WRITE(LUNOUT,'(4(10X,A,/))')
     *    ' WARNING: ONLY LEPTONIC CORRECTIONS AND SELF ENRGIES',
     *    ' CAN BE APPLIED FOR STRUCTURE FUNCTION INPUT',
     *    ' OR IF F_L IS INCLUDED ',
     *    ' LPAR(13...16) ARE SET TO 0'
        ENDIF
        IOPCHQ=0
        IF(INT3(4).GT.0 .OR. ISAM3(4).GT.0) THEN
          INT3(4)=0
          ISAM3(4)=0
          IOPCHQ=1
        ENDIF
        IF (IOPCHQ.NE.0) THEN
          WRITE(LUNOUT,'(3(10X,A/))')
     *    ' WARNING: QUARKONIC BREMSSTRAHLUNG CANNOT BE SEPARATED',
     *    ' FOR STRUCTURE FUNCTION INPUT OR IF F_L IS INCLUDED',
     *    ' INT3(4) AND ISAM3(4) ARE SET TO 0'
        ENDIF
      ENDIF

C---Initialization of LHAPDF 
      IF (ILIB.EQ.2) THEN
        LHAPTH=LHAPATH
c        CALL SetPDFPath(LHAPTH)
c        WRITE(LUNOUT,'(2(10X,A/))')
c     *  ' USING PATH TO PDFSETS: ',
c     *  LHAPATH
        ROOTS=SQRT(SP-MEI2-MPRO2)
C...Initialize parton distributions: PDFLIB.
        PARM(1)='NPTYPE'
        VALUE(1)=1
        PARM(2)='NGROUP'
        VALUE(2)=ICODE/1000
        PARM(3)='NSET'
        VALUE(3)=MOD(ICODE,1000)
        PARM(4)='TMAS'
        VALUE(4)=PYSTOP
c   CALL PDFSET via PYINIT 
        ROOTS=SQRT(SP-MEI2-MPRO2)
        IF (LLEPT.LE.0) THEN
          CALL PYINIT('NONE','e-','p',ROOTS)
         ELSE
          CALL PYINIT('NONE','e+','p',ROOTS)
        ENDIF
      ENDIF

C---WRITE OPTIONS
      IF (ILIB.EQ.1) THEN
        WRITE(LUNOUT,'(/A/A/)')
     *   '           PARTON DISTRIBUTIONS TAKEN FROM PYSTFU   *****'
     *  ,'           VIA LYSTFU FROM LEPTO 6.5                *****'
        IF (ICODE.NE.0) THEN
        WRITE(LUNOUT,'(A,I5)')
     *   '           WITH IDENTIFICATION CODE = ',ICODE
        ELSE
        WRITE(LUNOUT,'(A)')
     *   '           SEE VALUE OF ILQMOD AND THE MANUAL'
        ENDIF
      ELSEIF (ILIB.EQ.2) THEN
        WRITE(LUNOUT,'(/A)')
     *   '           PARTON DISTRIBUTIONS TAKEN FROM PDFLIB   *****'
        WRITE(LUNOUT,'(/A)')
     *   '           OR FROM LHAPDF VIA CALLS THROUGH PYTHIA  *****'
        WRITE(LUNOUT,'(A,I5)')
     *   '           WITH IDENTIFICATION CODE IVAL = ',ICODE
      ELSEIF (ILIB.EQ.3) THEN
        WRITE(LUNOUT,'(/A/A/)')
     *   ' *****  WARNING: WRONG CODE FOR ILIB: ILIB=3 NOT    *****'
     *  ,'                 ALLOWED (IS OBSOLETE)              *****'
      ENDIF
      WRITE(LUNOUT,'(/A/)')
     *      ' *****  LOW Q2 MODEL FOR STRUCTURE FUNCTIONS:    *****'
      WRITE(LUNOUT,'(A,I2)')
     *   '           ILQMOD = ',ILQMOD
      WRITE(LUNOUT,'(/8(A/))')
     * '           ILQMOD =  0: UNMODIFIED PARTON DISTRIBUTIONS'
     *,'                  =  1: LOW Q2 SUPPRESSED PDF    '
     *,'                  =  2: BRASSE AND STEIN WITH PDF           '
     *,'                  =  3: ALLM(1997) WITH PDF  '
     *,'                  =  4: BADELEK AND KWIECINSKI WITH PDF '
     *,'                  =  5: DONNACHIE AND LANDSHOFF WITH PDF '
     *,'                       (FOR DETAILS, SEE THE MANUAL )'
     *,'                  = 10: STRUCTURE FUNCTIONS FROM USER ROUTINE' 
      IF (ILQMOD.GE.4.AND.ILQMOD.NE.10) THEN
        WRITE(LUNOUT,'(//3(A/))')
     *' ***** WARNING: MAKE SURE THAT PARTON DISTRIBUTION FUNCTIONS'
     *,'                ARE CHOSEN CONSISTENTLY WITH THE LOW Q2 '
     *,'                BEHAVIOUR OF THE F_2 PARAMETRIZATION '
      ENDIF
      IF (ILQMOD.EQ.4.AND.XMIN.LT.1D-5) THEN
        WRITE(LUNOUT,'(//3(A/))')
     *' ***** WRONG COMBINATION OF STRUCTURE FUNCTION CODE WITH '
     *,'       KINEMATIC LIMITS: BADELEK-KWIECINSKI NOT VALID FOR '
     *,'       X-VALUES BELOW 1E-5. EXECUTION STOPPED '
        STOP
      ENDIF

C---FOR CHARGED CURRENT: NO FL, REMOVE B QUARKS
      IF (INT2(2).GT.0
     *   .OR.INT3(7).GT.0.OR.INT3(8).GT.0.OR.INT3(9).GT.0
     *   .OR.ISAM2(2).GT.0
     *   .OR.ISAM3(7).GT.0.OR.ISAM3(8).GT.0.OR.ISAM3(9).GT.0
     *   ) THEN
        IF (IFLOPT.NE.0) THEN
          IFLOPT=0
          WRITE(LUNOUT,'(//A/)')
     *' ***** LONGITUDINAL STRUCTURE FUNCTION NOT INCLUDED FOR CC *****'
        ENDIF
        IF (NPYMAX.GE.5) THEN
          NPYMAX=4
          WRITE(LUNOUT,'(//A/)')
     *' ***** NUMBER OF FLAVORS: NPYMAX RESET TO 4 (NO TOP IN CC) *****'
        ENDIF
      ENDIF

C---OPTIONS FOR LONGITUDINAL STRUCTURE FUNCTION
C---TO BE USED TOGETHER WITH PARTON DISTRIBUTIONS
      IF (IFLOPT.EQ.0) THEN
        WRITE(LUNOUT,'(//A/)')
     *    ' *****  LONGITUDINAL STRUCTURE FUNCTION NOT INCLUDED *****'
      ELSEIF (IFLOPT.GE.1) THEN
        IF (ILQMOD.GE.2) THEN
          IFLOPT=0
          LQCD=0
          LTM=0
          LHT=0
          WRITE(LUNOUT,'(//4(A/))')
     *      ' *****  WARNING: LONGITUDINAL STRUCTURE FUNCTION *****'
     *     ,'           NOT INCLUDED. '
     *     ,'           INCONSISTENT INPUT: IFLOPT > 0 AND '
     *     ,'           ILQMOD > 1, IFLOPT SET TO 0.'
        ELSE
          LQCD=MOD(IFLOPT,10)
          LTM=MOD(IFLOPT/10,10)
          LHT=IFLOPT/100
          IF (ILQMOD.LE.1) IPDFOP=2
          WRITE(LUNOUT,'(//A/,3(A,I3,/))')
     *      ' *****  LONGITUDINAL STRUCTURE FUNCTION INCLUDED *****'
     *     ,'           QCD CONTRIBUTION TO F_L: LQCD = ',LQCD
     *     ,'           TARGET MASS EFFECTS:      LTM = ', LTM
     *     ,'           HIGHER TWIST:             LHT = ', LHT
          WRITE(LUNOUT,'(//A/,2(A,I5,/),2(A,F10.4,/))')
     *      ' *****  DETERMINATION OF ALPHA_S NEEDED FOR F_L: *****'
     *     ,'           ORDER OF ALPHA_S IN ULALPS: MST111 = ',MST111
     *     ,'           TREATMENT OF SINGULARITY:   MST115 = ',MST115
     *     ,'           FIX ALPHA_S VALUE:          PAR111 = ',PAR111
     *     ,'           LAMBDA IN RUNNING ALPHA_S:  PAR112 = ',PAR112
     *     ,'           ACCUARCY IN FL-INTEGRATION: PARL11 = ',PARL11
     *     ,'           PARAMETER FOR HIGHER TWIST: PARL19 = ',PARL19
          WRITE(LUNOUT,'(//A/A/A)')
     *     ' *****  NOTE: LOW Q2 BEHAVIOUR OF F_L IS DETERMINED ****'
     *    ,'              BY THE VALUE OF ILQMOD, '
     *    ,'              SEE MANUAL FOR DETAILS '
          IF (IPDFOP.LT.2) WRITE(LUNOUT,'(//A/A/A)')
     *      ' *****  WARNING: WITH THIS OPTION NO SEPARATION ****'
     *     ,'        OF FLAVORS FOR THE TOTAL CROSS SECTION '
     *     ,'        ---> DJANGO CAN NOT RUN '
        ENDIF
      ENDIF

C---NUCLEAR TARGET
      IF (INT(HNA).NE.1.OR.INT(HNZ).NE.1) THEN
Cv4.6.19
c        IF (ILQMOD.LE.1) IPDFOP=2
        IF (ILQMOD.LE.1) IPDFOP=1
        WRITE(LUNOUT,'(//A/)')
     *    ' *****  NUCLEAR TARGET  *****'
        WRITE(LUNOUT,'(/2(A,F5.0,/))')
     *    '           A-NUCLEUS = ',HNA
     *   ,'           Z-NUCLEUS = ',HNZ
        WRITE(LUNOUT,'(//A/)')
     *    ' *****  MODEL FOR NUCLEAR CORRECTION / SHADOWING  *****'
        IF ((ICODE.GE.100050).AND.(ICODE.LE.100169)) THEN
          INUMOD=0
          WRITE(LUNOUT,'(/2(A/))')
     *    ' *****  NO EXTRA CORRECTION FOR HKN PARTONS       *****'
     *   ,'        INUMOD IS SET TO 0                             '
        ENDIF
        WRITE(LUNOUT,'(A,I4)')
     *   '           INUMOD = ',INUMOD
        WRITE(LUNOUT,'(/6(A/))')
     * '           INUMOD =  0:    NO NUCLEAR CORRECTION '
     *,'                  =  1000: SIMPLE Q2-INDEP. NUCL. SHADOWING '
     *,'                  =  2000: EKS98 '
     *,'                  =  3000: EPS08 '
     *,'                  >  4000: EPS09 '
     *,'                     (4000 + 100*iorder + ipset) '
        IF (INUMOD.EQ.0) THEN
        ELSEIF (INUMOD.EQ.1000) THEN
        ELSEIF (INUMOD.EQ.2000) THEN
          Call setlhaparm('EKS98')
        ELSEIF (INUMOD.EQ.3000) THEN
          Call setlhaparm('EPS08')
        ELSEIF (INUMOD.GT.4000) THEN
          IF ((INUMOD-4000)/100.EQ.1) THEN
            ipset=inumod-4100
            write(chpset,'(I2)') ipset
            lhaprm='EPS09LO,'//chpset
            Call setlhaparm(lhaprm)
          ELSEIF ((INUMOD-4000)/100.EQ.2) THEN
            ipset=inumod-4200
            write(chpset,'(I2)') ipset
            lhaprm='EPS09NLO,'//chpset
            Call setlhaparm(lhaprm)
          ENDIF
        ENDIF
      ENDIF
      
C---Polarized nucleon PDFs
      IF (HPOLAR.NE.0D0) THEN
        WRITE(LUNOUT,'(//A/)')
     *    ' *****  POLARIZED PARTON DISTRIBUTIONS  *****'
        WRITE(LUNOUT,'(A,I4)')
     *   '           IDPVR = ',IDPVR
        IPDSET=IDPVR/100
        IPDMOD=MOD(IDPVR,100)
c
        IF (IPDSET.EQ.1) THEN
        WRITE(LUNOUT,'(/2(A/))')
     * '           DSSV '
     *,'           de Florian, Sassot, Stratmann, Vogelsang (2008) '
c
        ELSEIF (IPDSET.EQ.2) THEN
        WRITE(LUNOUT,'(2(A,I4))')
     *   '           IPDSET = ',IPDSET,'   IPDMOD = ',IPDMOD
        IF (IPDMOD.GT.4.or.IPDMOD.LT.1) THEN
          WRITE(LUNOUT,*) ' INVALID POLARIZED PARTON CODE'
          WRITE(LUNOUT,*) ' ***** EXECUTION STOPPED ***** '
        STOP
        ELSE
          WRITE(LUNOUT,'(/2(A/))')
     *   '           DNS '
     *  ,'           de Florian, Navarro, Sassot (2005) '
        ENDIF
c
        ELSEIF (IPDSET.EQ.3) THEN
        WRITE(LUNOUT,'(2(A,I4))')
     *   '           IPDSET = ',IPDSET,'   IPDMOD = ',IPDMOD
        IF (IPDMOD.GT.12.or.IPDMOD.LT.1) THEN
          WRITE(LUNOUT,*) ' INVALID POLARIZED PARTON CODE'
          WRITE(LUNOUT,*) ' ***** EXECUTION STOPPED ***** '
        ELSE
          WRITE(LUNOUT,'(/2(A/))')
     *   '           DS2000 '
     *  ,'           de Florian, Sassot (2000) '
        ENDIF
c
        ELSEIF (IPDSET.EQ.4) THEN
        WRITE(LUNOUT,'(2(A,I4))')
     *   '           IPDSET = ',IPDSET,'   IPDMOD = ',IPDMOD
        IF (IPDMOD.GT.2) THEN
          WRITE(LUNOUT,*) ' INVALID POLARIZED PARTON CODE'
          WRITE(LUNOUT,*) ' ***** EXECUTION STOPPED ***** '
        ELSE
          WRITE(LUNOUT,'(/2(A/))')
     *   '           GSLO '
     *  ,'           Gehrmann, Stirling (1996) '
        ENDIF
c
        ELSEIF (IPDSET.EQ.5) THEN
        WRITE(LUNOUT,'(2(A,I4))')
     *   '           IPDSET = ',IPDSET,'   IPDMOD = ',IPDMOD
        IF (IPDMOD.GT.2) THEN
          WRITE(LUNOUT,*) ' INVALID POLARIZED PARTON CODE'
          WRITE(LUNOUT,*) ' ***** EXECUTION STOPPED ***** '
        ELSE
          WRITE(LUNOUT,'(/2(A/))')
     *   '           GSNLO '
     *  ,'           Gehrmann, Stirling (1996) '
        ENDIF
c
        ELSEIF (IPDSET.EQ.6) THEN
        WRITE(LUNOUT,'(2(A,I4))')
     *   '           IPDSET = ',IPDSET,'   IPDMOD = ',IPDMOD
        IF (IPDMOD.GT.4.or.IPDMOD.LT.1) THEN
          WRITE(LUNOUT,*) ' INVALID POLARIZED PARTON CODE'
          WRITE(LUNOUT,*) ' ***** EXECUTION STOPPED ***** '
        ELSE
          WRITE(LUNOUT,'(/2(A/))')
     *   '           BB '
     *  ,'           Bluemlein, Boettcher (2002) '
        ENDIF
c
        ELSEIF (IPDSET.EQ.7) THEN
        WRITE(LUNOUT,'(2(A,I4))')
     *   '           IPDSET = ',IPDSET,'   IPDMOD = ',IPDMOD
        IF (IPDMOD.GT.2.or.IPDMOD.LT.1) THEN
          WRITE(LUNOUT,*) ' INVALID POLARIZED PARTON CODE'
          WRITE(LUNOUT,*) ' ***** EXECUTION STOPPED ***** '
        ELSE
          WRITE(LUNOUT,'(/2(A/))')
     *   '           AAC '
     *  ,'           Hirai, Kumano (2008) '
        ENDIF
c
        ELSEIF (IPDSET.EQ.8) THEN
        WRITE(LUNOUT,'(2(A,I4))')
     *   '           IPDSET = ',IPDSET,'   IPDMOD = ',IPDMOD
        IF (IPDMOD.GT.2.or.IPDMOD.LT.1) THEN
          WRITE(LUNOUT,*) ' INVALID POLARIZED PARTON CODE'
          WRITE(LUNOUT,*) ' ***** EXECUTION STOPPED ***** '
        ELSE
          WRITE(LUNOUT,'(/2(A/))')
     *   '           LSS '
     *  ,'           Leader, Sidorov, Stamenov (2010) '
        ENDIF
c
        ELSEIF (IPDSET.EQ.9) THEN
        WRITE(LUNOUT,'(2(A,I4))')
     *   '           IPDSET = ',IPDSET,'   IPDMOD = ',IPDMOD
        IF (IPDMOD.GT.4.or.IPDMOD.LT.1) THEN
          WRITE(LUNOUT,*) ' INVALID POLARIZED PARTON CODE'
          WRITE(LUNOUT,*) ' ***** EXECUTION STOPPED ***** '
        ELSE
          WRITE(LUNOUT,'(/2(A/))')
     *   '           GRSV '
     *  ,'           Gluck, Reya, Stratmann, Vogelsang (2001) '
        ENDIF
c
        ELSE
        WRITE(LUNOUT,'(//3(A/))')
     * ' ***** UNKNOWN CODE FOR POLARIZED PARTON DISTRIBUTIONS '
     *,'       EXECUTION STOPPED '
        STOP        
      ENDIF
      ENDIF
      
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSCPDF(IMOD)
C---Write final statistics for calls to HSDVPR with x, Q2 out of range
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD,IDPVR
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
c
      integer ixoutmn,iqoutmn,ixoutmx,iqoutmx
      DOUBLE PRECISION XMIND, XMAXD, Q2MIND, Q2MAXD
      COMMON /HSDPEX/ ixoutmn, ixoutmx, iqoutmn, iqoutmx,
     *                XMIND, XMAXD, Q2MIND, Q2MAXD
c
      IF (IMOD.EQ.1) THEN
c...write after integration and reset counters for event sampling
       IF (ixoutmn.GT.0) THEN
        WRITE(LUNTES,'(2(/A),1PE13.4,A,I8,A,/A,1PE13.4)')
     *     ' *** WARNING: X OUT OF RANGE FOR POLARIZED PDF: ',
     *     '     X<',XMIND, '   OCCURED ', ixoutmn, 
     *     '     TIMES DURING INTEGRATION ',
     *     '     X WAS SET TO THE LOWER LIMIT XMIND =',XMIND
        ixoutmn=0
       ENDIF
       IF (ixoutmx.GT.0) THEN
        WRITE(LUNTES,'(2(/A),1PE13.4,A,I8,A,/A,1PE13.4)')
     *     ' *** WARNING: X OUT OF RANGE FOR POLARIZED PDF: ',
     *     '     X>',XMAXD, '   OCCURED ', ixoutmx, 
     *     '     TIMES DURING INTEGRATION ',
     *     '     X WAS SET TO THE UPPER LIMIT XMAXD =',XMAXD
        ixoutmx=0
       ENDIF
       IF (iqoutmn.GT.0) THEN
        WRITE(LUNTES,'(2(/A),1PE13.4,A,I8,A,/A,1PE13.4)')
     *     ' *** WARNING: Q2 OUT OF RANGE FOR POLARIZED PDF: ',
     *     '     Q2<',Q2MIND, '   OCCURED ', iqoutmn, 
     *     '     TIMES DURING INTEGRATION ',
     *     '     Q2 WAS SET TO THE LOWER LIMIT Q2MIND =',Q2MIND
        iqoutmn=0
       ENDIF
       IF (iqoutmx.GT.0) THEN
        WRITE(LUNTES,'(2(/A),1PE13.4,A,I8,A,/A,1PE13.4)')
     *     ' *** WARNING: Q2 OUT OF RANGE FOR POLARIZED PDF: ',
     *     '     Q2>',Q2MAXD, '   OCCURED ', iqoutmx, 
     *     '     TIMES DURING INTEGRATION ',
     *     '     Q2 WAS SET TO THE UPPER LIMIT Q2MAXD =',Q2MAXD
        iqoutmx=0
       ENDIF
      ELSEIF (IMOD.EQ.2) THEN
c...write after event sampling
       IF (ixoutmn.GT.0) THEN
        WRITE(LUNTES,'(2(/A),1PE13.4,A,I8,A,/A,1PE13.4)')
     *     ' *** WARNING: X OUT OF RANGE FOR POLARIZED PDF: ',
     *     '     X<',XMIND, '   OCCURED ', ixoutmn, 
     *     '     TIMES DURING EVENT SAMPLING ',
     *     '     X WAS SET TO THE LOWER LIMIT XMIND =',XMIND
       ENDIF
       IF (ixoutmx.GT.0) THEN
        WRITE(LUNTES,'(2(/A),1PE13.4,A,I8,A,/A,1PE13.4)')
     *     ' *** WARNING: X OUT OF RANGE FOR POLARIZED PDF: ',
     *     '     X>',XMAXD, '   OCCURED ', ixoutmx, 
     *     '     TIMES DURING EVENT SAMPLING ',
     *     '     X WAS SET TO THE UPPER LIMIT XMAXD =',XMAXD
       ENDIF
       IF (iqoutmn.GT.0) THEN
        WRITE(LUNTES,'(2(/A),1PE13.4,A,I8,A,/A,1PE13.4)')
     *     ' *** WARNING: Q2 OUT OF RANGE FOR POLARIZED PDF: ',
     *     '     Q2<',Q2MIND, '   OCCURED ', iqoutmn, 
     *     '     TIMES DURING EVENT SAMPLING ',
     *     '     Q2 WAS SET TO THE LOWER LIMIT Q2MIND =',Q2MIND
       ENDIF
       IF (iqoutmx.GT.0) THEN
        WRITE(LUNTES,'(2(/A),1PE13.4,A,I8,A,/A,1PE13.4)')
     *     ' *** WARNING: Q2 OUT OF RANGE FOR POLARIZED PDF: ',
     *     '     Q2>',Q2MAXD, '   OCCURED ', iqoutmx, 
     *     '     TIMES DURING EVENT SAMPLING ',
     *     '     Q2 WAS SET TO THE UPPER LIMIT Q2MAXD =',Q2MAXD
       ENDIF
      ENDIF

      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE FIUSER(X,Q2,ZF1,ZF2,IPDFR)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      
      ZF1=0D0
      ZF2=0D0
      IPDFR=0

C--->
C    insert here your conditions on x,Q2 for range of applicability 
C    of your parametrization. If outside, set IPDFR = 1 and return. 
C    The calling routine will then use parton distributions as specified
C    by ILIB and ICODE to calculate structure functions
C 
c      IF (X.LT.XFMIN.OR.X.GT.XFMAX.
c     &.OR.Q2.LT.QFMIN.OR.Q2.GT.QFMAX) THEN
c        IPDFR=1
c        RETURN
c      ENDIF 
      IF (X.LT.0D0.OR.Q2.LT.0D0) THEN
         IPDFR=1
         RETURN
      ENDIF
C--->
      
C--->
C    insert here your parametrization for F1 and F2 (electromagnetic
C    part, no contribution from Z exchange).
C    Note: F1 = (F2-FL)/(2*X)
C
c      ZF1=
c      ZF2=
C--->
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSSTRF(X,Q2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSFIJK/ F1(2,2),F2(2,2),F3(2,2)
      COMMON /HSGIJK/ G1(2,2),G3(2,2),G4(2,2)
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMAX,NPYMIN
      REAL            PYSTOP,PYSLAM
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD,IDPVR
C...HSRADF USED IN SUBROUTINE STRFBS (CALL TO STEIN OR BRASSE)
      COMMON /HSRADC/ AMP , AMP2,RPI ,RPI2 , ALFA, AML , AML2
      COMMON /HSRADF/ AP,AP2,AP4,AL2,AL4,ALPS,CALPI,W2PIT,W2TR,TTR
      DIMENSION DSBOS(2)
      LOGICAL LFIRST
      DATA LFIRST/.TRUE./

      IF (LFIRST) THEN
        LFIRST=.FALSE.
        W2TR=4D0
        TTR=6D0
        AMP2=MPRO2
        AP=2D0*MPRO
        W2PIT=1.15184D0
        W2MIN=1.2321D0
        S=SP
        IF (ILQMOD.EQ.4.AND.ICODE.EQ.0) THEN
          ILIB=2
          ICODE=3034
        ENDIF
        IF (ILQMOD.EQ.5.AND.ICODE.EQ.0) THEN
          ILIB=2
          ICODE=3025
        ENDIF
      ENDIF
      DSBOS(1)=1D0
      DSBOS(2)=Q2/(Q2+MZ2)
      ZF1=0D0
      ZF2=0D0

      IF (ILQMOD.EQ.0.OR.ILQMOD.EQ.1) THEN
C---CALCULATE STRUCTURE FUNCTIONS FROM PARTON DISTRIBUTIONS
        GOTO 1000
      ELSEIF (ILQMOD.EQ.2) THEN
C---FOR LOW Q2 ONLY GAMMA EXCHANGE, PARAMETRIZATION BY BRASSE AND
C   STEIN, COMBINED WITH PARTON DISTRIBUTION FUNCTIONS FOR LARGE Q2
        IF (Q2.LT.TTR) THEN
          DO 11 IB1=1,2
          DO 12 IB2=1,2
            F1(IB1,IB2)=0D0
            F2(IB1,IB2)=0D0
            F3(IB1,IB2)=0D0
 12       CONTINUE
 11       CONTINUE
          CALL STRFBS(X,Q2,ZF1,ZF2)
          F1(1,1)=ZF1
          F2(1,1)=ZF2
          RETURN
         ELSE
          GOTO 1000
        ENDIF
      ELSEIF (ILQMOD.EQ.3.AND.ICODE.NE.0) THEN
C---FOR LOW Q2 ONLY GAMMA EXCHANGE, PARAMETRIZATION BY ALLM (1997) 
C   COMBINED WITH PARTON DISTRIBUTION FUNCTIONS FOR LARGE Q2
C   AND BRASSE PARAMETRIZATION FOR LOW W 
        DO 31 IB1=1,2
        DO 32 IB2=1,2
          F1(IB1,IB2)=0D0
          F2(IB1,IB2)=0D0
          F3(IB1,IB2)=0D0
 32     CONTINUE
 31     CONTINUE
        W2=(1D0-X)/X*Q2+MPRO2
        IF(Q2.LT.TTR) THEN
          IF (W2.LT.W2TR) THEN
            CALL STRFBS(X,Q2,ZF1,ZF2)
            ELSE
            CALL HSSTAL(X,Q2,ZF1,ZF2)
          ENDIF
         ELSEIF (Q2.GE.5000D0.OR.X.GE.0.85D0) THEN
          GOTO 1000
         ELSE
          CALL HSSTAL(X,Q2,ZF1,ZF2)
        ENDIF
        F1(1,1)=ZF1
        F2(1,1)=ZF2
        RETURN
      ELSEIF (ILQMOD.EQ.3.AND.ICODE.EQ.0) THEN
        WRITE(LUNOUT,'(//3(A/))')
     *' ***** WRONG STRUCTURE FUNCTION CODE:'
     *,'       ILQMOD = 3 AND ICODE = 0'
     *,'       EXECUTION STOPPED '
        STOP
      ELSEIF (ILQMOD.EQ.4) THEN
C---BADELEK KWIECINSKI
        DO 41 IB1=1,2
        DO 42 IB2=1,2
          F1(IB1,IB2)=0D0
          F2(IB1,IB2)=0D0
          F3(IB1,IB2)=0D0
 42     CONTINUE
 41     CONTINUE
        IF ((X.GT.0.1D0.AND.Q2.GT.6D0).OR.Q2.GT.998.8D0) THEN
          GOTO 1000
          ELSE
          ANU=Q2/X/2D0
          IF (ANU.LT.10D0) THEN
            CALL STRFBS(X,Q2,ZF1,ZF2)
            ELSE
            CALL HSSTBK(X,Q2,ZF1,ZF2)
          ENDIF
        ENDIF
        F1(1,1)=ZF1
        F2(1,1)=ZF2
        GOTO 2000
      ELSEIF (ILQMOD.EQ.5) THEN
C---DONNACHIE LANDSHOFF with pdf's
        DO 51 IB1=1,2
        DO 52 IB2=1,2
          F1(IB1,IB2)=0D0
          F2(IB1,IB2)=0D0
          F3(IB1,IB2)=0D0
 52     CONTINUE
 51     CONTINUE
        IF (Q2.LT.10D0) THEN
          CALL HSSTDL(X,Q2,ZF1,ZF2)
          F1(1,1)=ZF1
          F2(1,1)=ZF2
          GOTO 2000
         ELSE
          GOTO 1000
        ENDIF
      ELSEIF (ILQMOD.EQ.10) THEN
C---STRUCTURE FUNCTIONS FROM USER ROUTINE
        DO 61 IB1=1,2
        DO 62 IB2=1,2
          F1(IB1,IB2)=0D0
          F2(IB1,IB2)=0D0
          F3(IB1,IB2)=0D0
 62     CONTINUE
 61     CONTINUE
        CALL FIUSER(X,Q2,ZF1,ZF2,IPDFR)
        IF (IPDFR.EQ.1) GOTO 1000
        F1(1,1)=ZF1
        F2(1,1)=ZF2
        GOTO 2000
      ELSE
        WRITE(LUNOUT,200) ILQMOD
        STOP
  200   FORMAT(/,'          WRONG VALUE FOR ILQMOD: ',I8,/
     F          ,' ******** EXECUTION STOPPED IN HSSTRF ',/)
      ENDIF

1000  CONTINUE
      CALL HSPVER(X,Q2)
      IF (HPOLAR.NE.0D0) THEN
        CALL HSDPVR(X,Q2)
      ENDIF
      IF (INT(HNA).EQ.1.AND.INT(HNZ).EQ.1) THEN
C---PROTON TARGET
C---unpolarized
        QUT=QU+QC+QT
        QDT=QD+QS+QB
        QUBT=QBU+QBC+QBT
        QDBT=QBD+QBS+QBB
C---polarized
        DQUT=DQU+DQC+DQT
        DQDT=DQD+DQS+DQB
        DQUBT=DQBU+DQBC+DQBT
        DQDBT=DQBD+DQBS+DQBB
      ELSEIF (INT(HNA).EQ.1.AND.INT(HNZ).EQ.0) THEN
C---NEUTRON TARGET
C---unpolarized
        QDT=QU+QC+QT
        QUT=QD+QS+QB
        QDBT=QBU+QBC+QBT
        QUBT=QBD+QBS+QBB
C---polarized
        DQDT=DQU+DQC+DQT
        DQUT=DQD+DQS+DQB
        DQDBT=DQBU+DQBC+DQBT
        DQUBT=DQBD+DQBS+DQBB
      ELSE
C---OTHER TARGETS
        IF (INUMOD.EQ.1000) THEN
C---SIMPLE SHADOWING MODEL: FIRST CALCULATE DEUTERON STRUCTURE FUNCTIONS
C   PER NUCLEON
C---unpolarized
          HNRAT=HSNRAT(X)
          QDT=HNRAT*(QU+QC+QT+QD+QS+QB)/2D0
          QUT=HNRAT*QDT
          QDBT=HNRAT*(QBU+QBC+QBT+QBD+QBS+QBB)/2D0
          QUBT=HNRAT*QDBT
C---polarized
          DQDT=HNRAT*(DQU+DQC+DQT+DQD+DQS+DQB)/2D0
          DQUT=HNRAT*DQDT
          DQDBT=HNRAT*(DQBU+DQBC+DQBT+DQBD+DQBS+DQBB)/2D0
          DQUBT=HNRAT*DQDBT
        ELSE 
C---WEIGHTED AVERAGE OF UP AND DOWN, CORRECTION INCLUDED IN PDFs
C---unpolarized
          QUT=(HNZ*(QU+QC+QT)+(HNA-HNZ)*(QD+QS+QB))/HNA
          QDT=((HNA-HNZ)*(QU+QC+QT)+HNZ*(QD+QS+QB))/HNA
          QUBT=(HNZ*(QBU+QBC+QBT)+(HNA-HNZ)*(QBD+QBS+QBB))/HNA
          QDBT=((HNA-HNZ)*(QBU+QBC+QBT)+HNZ*(QBD+QBS+QBB))/HNA
C---polarized
          DQUT=(HNZ*(DQU+DQC+DQT)+(HNA-HNZ)*(DQD+DQS+DQB))/HNA
          DQDT=((HNA-HNZ)*(DQU+DQC+DQT)+HNZ*(DQD+DQS+DQB))/HNA
          DQUBT=(HNZ*(DQBU+DQBC+DQBT)+(HNA-HNZ)*(DQBD+DQBS+DQBB))/HNA
          DQDBT=((HNA-HNZ)*(DQBU+DQBC+DQBT)+HNZ*(DQBD+DQBS+DQBB))/HNA
        ENDIF
      ENDIF

      DO 101 IB1=1,2
       DO 102 IB2=1,2
C---unpolarized
         F1(IB1,IB2) = (AFIJ(2,IB1,IB2)*(QUT+QUBT)
     &                 +AFIJ(3,IB1,IB2)*(QDT+QDBT))
     &                 /4D0*DSBOS(IB1)*DSBOS(IB2)  /2D0
         F2(IB1,IB2) = 2D0*X*F1(IB1,IB2)
         F3(IB1,IB2) = (BFIJ(2,IB1,IB2)*(QUT-QUBT)
     &                  +BFIJ(3,IB1,IB2)*(QDT-QDBT))
     &                   /4D0*DSBOS(IB1)*DSBOS(IB2)
C---polarized
         G1(IB1,IB2) = (BGIJ(2,IB1,IB2)*(DQUT+DQUBT)
     &                 +BGIJ(3,IB1,IB2)*(DQDT+DQDBT))
     &                 /4D0*DSBOS(IB1)*DSBOS(IB2)  /2D0
         G3(IB1,IB2) = (AGIJ(2,IB1,IB2)*(DQUT-DQUBT)
     &                  +AGIJ(3,IB1,IB2)*(DQDT-DQDBT))
     &                   /4D0*DSBOS(IB1)*DSBOS(IB2) /2D0
         G4(IB1,IB2) = 2D0*X*G3(IB1,IB2)
 102  CONTINUE
 101  CONTINUE

C---INCLUDE LONGITUDINAL STRUCTURE FUNCTION
 2000 CONTINUE
      IF (IFLOPT.GT.0) THEN
      FL=0D0
      F2EM=F2(1,1)
      Q2L=Q2
      IF (Q2L.LT.TTR) THEN
        Q2L=TTR
        XP21=SNGL(X)*SNGL(SP-MPRO2-MEI2)
        IF (Q2L.GE.XP21) Q2L=XP21
      ENDIF
      CALL HSLUFL(X,Q2L,F2EM,FL)
      F1(1,1)=(F2(1,1)-FL)/2D0/X
      IF (F1(1,1).LT.0D0) F1(1,1)=0D0
      ENDIF

C---APPLY LOW Q2 SUPPRESSION
C   MOVED TO LYSTFU
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   SUBROUTINE FOR STRUCTURE FUNCTIONS INCLUDING SELF ENERGIES
C   AND LEPTONIC QED CORRECTIONS
C   CALCULATED FROM PARTON DISTRIBUTION FUNCTIONS
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSSTR1(X,Q2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSRGG,HSSRGZ,HSSRZZ,CG,CM,CZ
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSFIJK/ F1(2,2),F2(2,2),F3(2,2)
      COMMON /HSGIJK/ G1(2,2),G3(2,2),G4(2,2)
      COMMON /HSSMCP/ VAFI(2,3,2),AFIJ(3,2,2),BFIJ(3,2,2),FLIND(2,3,2,2)
      COMMON /HSSMCG/ AGIJ(3,2,2),BGIJ(3,2,2)
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSASY0/ PONULL,CHNULL,ICPASY
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGSW/  SW,CW,SW2,CW2
     *              ,MW,MZ,MH,ME,MMY,MTAU,MU,MD,MS,MC,MB,MT
     *              ,MW2,MZ2,MH2,ME2,MMY2,MTAU2,MU2,MD2,MS2,MC2,MB2,MT2
      COMMON /HSDELR/ DELTAR,AGF0,DRHOT,DALPMZ,XGMT,ALPQCD,BTOP4,DRPIW2
      COMMON /HSPDFO/ IPDFOP,IFLOPT,LQCD,LTM,LHT
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HYSTFU/ PYSTOP,PYSLAM,NPYMAX,NPYMIN
      REAL            PYSTOP,PYSLAM
      COMMON /HSSTRP/ ICODE,ILIB,ILQMOD,IDPVR
C...HSRADF USED IN SUBROUTINE STRFBS (CALL TO STEIN OR BRASSE)
      COMMON /HSRADC/ AMP , AMP2,RPI ,RPI2 , ALFA, AML , AML2
      COMMON /HSRADF/ AP,AP2,AP4,AL2,AL4,ALPS,CALPI,W2PIT,W2TR,TTR
      DIMENSION VAFI1(2,3,2),AFIJ1(3,2,2),BFIJ1(3,2,2),FLIND1(2,3,2,2),
     *          AGIJ1(3,2,2),BGIJ1(3,2,2)
      DIMENSION DSBOS(2),DSBS1(2)
      LOGICAL LFIRST
      DATA LFIRST/.TRUE./

      IF (LFIRST) THEN
        LFIRST=.FALSE.
        W2TR=4D0
        TTR=6D0
        AMP2=MPRO2
        AP=2D0*MPRO
        W2PIT=1.15184D0
        W2MIN=1.2321D0
        S=SP
        IF (ILQMOD.EQ.4.AND.ICODE.EQ.0) THEN
          ILIB=2
          ICODE=3034
        ENDIF
        IF (ILQMOD.EQ.5.AND.ICODE.EQ.0) THEN
          ILIB=2
          ICODE=3025
        ENDIF
      ENDIF
      ZF1=0D0
      ZF2=0D0
      T=-Q2
      SHAT=SP-MEI2-MPRO2
      Y=Q2/SHAT/X
      DSBOS(1)=1D0
      DSBOS(2)=Q2/(Q2+MZ2)

      IF (ILQMOD.EQ.0.OR.ILQMOD.EQ.1) THEN
C---CALCULATE STRUCTURE FUNCTIONS FROM PARTON DISTRIBUTIONS
        GOTO 1000
      ELSEIF (ILQMOD.EQ.2) THEN
C---FOR LOW Q2 ONLY GAMMA EXCHANGE, PARAMETRIZATION BY BRASSE AND
C   STEIN, COMBINED WITH PARTON DISTRIBUTION FUNCTIONS FOR LARGE Q2
        IF (Q2.LT.TTR) THEN
          DO 41 IB1=1,2
          DO 42 IB2=1,2
            F1(IB1,IB2)=0D0
            F2(IB1,IB2)=0D0
            F3(IB1,IB2)=0D0
 42       CONTINUE
 41       CONTINUE
          CALL STRFBS(X,Q2,AF1,AF2)
          F1(1,1)=AF1
          F2(1,1)=AF2
C..PHOTON SELF ENERGY
          IF (LPAR(7).GE.1) THEN
            CG=HSSRGG(T)
            PIGGG=DREAL(CG)/T
            ELSE
            PIGGG=0D0
          ENDIF
          DVACGG=1D0/(1D0+PIGGG)
C..LEPTONIC QED VERTEX CORRECTIONS INCLUDING SOFT BREMSSTRAHLUNG
          IF (LPAR(12).GE.1) THEN
            DVRTXV=HSDQDV(X,Q2)*ALP2PI
            DVRTXS=HSDQDS(Q2)*ALP2PI
            DVRTX1=DVRTXV-(Q2+4D0*MEI2)/(Q2-2D0*MEI2)*DVRTXS
            DVRTX2=DVRTXV+
     *         (2D0*(1D0-Y)+Y*Y/2D0)/(1D0-Y-MPRO2*Q2/SHAT/SHAT)*DVRTXS
            ELSE
            DVRTX1=0D0
            DVRTX2=0D0
          ENDIF
C
          IF (LPAR(3).LT.3) THEN
          F1(1,1)=F1(1,1)*(DVRTX1+DVACGG*DVACGG)
          F2(1,1)=F2(1,1)*(DVRTX2+DVACGG*DVACGG)
          ELSE
          F1(1,1)=F1(1,1)*(1D0+DVRTX1)*DVACGG*DVACGG
          F2(1,1)=F2(1,1)*(1D0+DVRTX2)*DVACGG*DVACGG
          ENDIF
          RETURN
         ELSE
          GOTO 1000
        ENDIF
      ELSEIF (ILQMOD.EQ.3.AND.ICODE.NE.0) THEN
C---FOR LOW Q2 ONLY GAMMA EXCHANGE, PARAMETRIZATION BY ABRAMOWICZ,
C   LEVY, LEVIN, MAOR, COMBINED WITH PARTON DISTRIBUTION FUNCTIONS FOR
C   LARGE Q2 FROM MORFIN AND TUNG
        DO 61 IB1=1,2
        DO 62 IB2=1,2
          F1(IB1,IB2)=0D0
          F2(IB1,IB2)=0D0
          F3(IB1,IB2)=0D0
 62     CONTINUE
 61     CONTINUE
        W2=(1D0-X)/X*Q2+MPRO2
        IF(Q2.LT.TTR) THEN
          IF (W2.LT.W2TR) THEN
            CALL STRFBS(X,Q2,ZF1,ZF2)
            ELSE
            CALL HSSTAL(X,Q2,ZF1,ZF2)
          ENDIF
         ELSEIF (Q2.GE.5000D0.OR.X.GE.0.85D0) THEN
          GOTO 1000
         ELSE
          CALL HSSTAL(X,Q2,ZF1,ZF2)
        ENDIF
        F1(1,1)=ZF1
        F2(1,1)=ZF2
C..PHOTON SELF ENERGY
        IF (LPAR(7).GE.1) THEN
          CG=HSSRGG(T)
          PIGGG=DREAL(CG)/T
          ELSE
          PIGGG=0D0
        ENDIF
        DVACGG=1D0/(1D0+PIGGG)
C..LEPTONIC QED VERTEX CORRECTIONS INCLUDING SOFT BREMSSTRAHLUNG
        IF (LPAR(12).GE.1) THEN
          DVRTXV=HSDQDV(X,Q2)*ALP2PI
          DVRTXS=HSDQDS(Q2)*ALP2PI
          DVRTX1=DVRTXV-(Q2+4D0*MEI2)/(Q2-2D0*MEI2)*DVRTXS
          DVRTX2=DVRTXV+
     *          (2D0*(1D0-Y)+Y*Y/2D0)/(1D0-Y-MPRO2*Q2/SHAT/SHAT)*DVRTXS
          ELSE
          DVRTX1=0D0
          DVRTX2=0D0
        ENDIF
C
        IF (LPAR(3).LT.3) THEN
        F1(1,1)=F1(1,1)*(DVRTX1+DVACGG*DVACGG)
        F2(1,1)=F2(1,1)*(DVRTX2+DVACGG*DVACGG)
        ELSE
        F1(1,1)=F1(1,1)*(1D0+DVRTX1)*DVACGG*DVACGG
        F2(1,1)=F2(1,1)*(1D0+DVRTX2)*DVACGG*DVACGG
        ENDIF
        RETURN
      ELSEIF (ILQMOD.EQ.3.AND.ICODE.EQ.0) THEN
        WRITE(LUNOUT,'(//3(A/))')
     *' ***** WRONG STRUCTURE FUNCTION CODE:'
     *,'       ILQMOD = 3 AND ICODE = 0'
     *,'       EXECUTION STOPPED '
        STOP
      ELSEIF (ILQMOD.EQ.4) THEN
C---FOR LOW Q2 ONLY GAMMA EXCHANGE,
C   PARAMETRIZATION BY BADELEK AND KWIECINSKI
C   COMBINED WITH BRASSE STEIN FOR SMALL HADRONIC MASS
        DO 81 IB1=1,2
        DO 82 IB2=1,2
          F1(IB1,IB2)=0D0
          F2(IB1,IB2)=0D0
          F3(IB1,IB2)=0D0
 82     CONTINUE
 81     CONTINUE
        IF ((X.GT.0.1D0.AND.Q2.GT.6D0).OR.Q2.GT.998.8D0) THEN
          GOTO 1000
          ELSE
          ANU=Q2/X/2D0
          IF (ANU.LT.10D0) THEN
            CALL STRFBS(X,Q2,ZF1,ZF2)
            ELSE
            CALL HSSTBK(X,Q2,ZF1,ZF2)
          ENDIF
        ENDIF
        F1(1,1)=ZF1
        F2(1,1)=ZF2
C..PHOTON SELF ENERGY
        IF (LPAR(7).GE.1) THEN
          CG=HSSRGG(T)
          PIGGG=DREAL(CG)/T
          ELSE
          PIGGG=0D0
        ENDIF
        DVACGG=1D0/(1D0+PIGGG)
C..LEPTONIC QED VERTEX CORRECTIONS INCLUDING SOFT BREMSSTRAHLUNG
        IF (LPAR(12).GE.1) THEN
          DVRTXV=HSDQDV(X,Q2)*ALP2PI
          DVRTXS=HSDQDS(Q2)*ALP2PI
          DVRTX1=DVRTXV-(Q2+4D0*MEI2)/(Q2-2D0*MEI2)*DVRTXS
          DVRTX2=DVRTXV+
     *          (2D0*(1D0-Y)+Y*Y/2D0)/(1D0-Y-MPRO2*Q2/SHAT/SHAT)*DVRTXS
          ELSE
          DVRTX1=0D0
          DVRTX2=0D0
        ENDIF
C
        IF (LPAR(3).LT.3) THEN
        F1(1,1)=F1(1,1)*(DVRTX1+DVACGG*DVACGG)
        F2(1,1)=F2(1,1)*(DVRTX2+DVACGG*DVACGG)
        ELSE
        F1(1,1)=F1(1,1)*(1D0+DVRTX1)*DVACGG*DVACGG
        F2(1,1)=F2(1,1)*(1D0+DVRTX2)*DVACGG*DVACGG
        ENDIF
        F2EM=F2(1,1)
        GOTO 2000
      ELSEIF (ILQMOD.EQ.5) THEN
C---FOR LOW Q2 ONLY GAMMA EXCHANGE,
C   PARAMETRIZATION BY DONNACHIE AND LANDSHOFF
C   COMBINED WITH PARTON DISTRIBUTION FUNCTIONS FOR LARGE Q2
        DO 91 IB1=1,2
        DO 92 IB2=1,2
          F1(IB1,IB2)=0D0
          F2(IB1,IB2)=0D0
          F3(IB1,IB2)=0D0
 92     CONTINUE
 91     CONTINUE
        IF (Q2.GT.10D0) THEN
          GOTO 1000
          ELSE
          CALL HSSTDL(X,Q2,ZF1,ZF2)
        ENDIF
        F1(1,1)=ZF1
        F2(1,1)=ZF2
C..PHOTON SELF ENERGY
        IF (LPAR(7).GE.1) THEN
          CG=HSSRGG(T)
          PIGGG=DREAL(CG)/T
          ELSE
          PIGGG=0D0
        ENDIF
        DVACGG=1D0/(1D0+PIGGG)
C..LEPTONIC QED VERTEX CORRECTIONS INCLUDING SOFT BREMSSTRAHLUNG
        IF (LPAR(12).GE.1) THEN
          DVRTXV=HSDQDV(X,Q2)*ALP2PI
          DVRTXS=HSDQDS(Q2)*ALP2PI
          DVRTX1=DVRTXV-(Q2+4D0*MEI2)/(Q2-2D0*MEI2)*DVRTXS
          DVRTX2=DVRTXV+
     *          (2D0*(1D0-Y)+Y*Y/2D0)/(1D0-Y-MPRO2*Q2/SHAT/SHAT)*DVRTXS
          ELSE
          DVRTX1=0D0
          DVRTX2=0D0
        ENDIF
C
        IF (LPAR(3).LT.3) THEN
        F1(1,1)=F1(1,1)*(DVRTX1+DVACGG*DVACGG)
        F2(1,1)=F2(1,1)*(DVRTX2+DVACGG*DVACGG)
        ELSE
        F1(1,1)=F1(1,1)*(1D0+DVRTX1)*DVACGG*DVACGG
        F2(1,1)=F2(1,1)*(1D0+DVRTX2)*DVACGG*DVACGG
        ENDIF
        F2EM=F2(1,1)
        GOTO 2000
      ELSEIF (ILQMOD.EQ.10) THEN
C---FOR LOW Q2 ONLY GAMMA EXCHANGE,
C   STRUCTURE FUNCTIONS FROM USER ROUTINE
C   PARAMETRIZATION BY DONNACHIE AND LANDSHOFF
        DO 101 IB1=1,2
        DO 102 IB2=1,2
          F1(IB1,IB2)=0D0
          F2(IB1,IB2)=0D0
          F3(IB1,IB2)=0D0
 102    CONTINUE
 101    CONTINUE
        CALL FIUSER(X,Q2,ZF1,ZF2,IPDFR)
        IF (IPDFR.EQ.1) GOTO 1000
        F1(1,1)=ZF1
        F2(1,1)=ZF2
C..PHOTON SELF ENERGY
        IF (LPAR(7).GE.1) THEN
          CG=HSSRGG(T)
          PIGGG=DREAL(CG)/T
          ELSE
          PIGGG=0D0
        ENDIF
        DVACGG=1D0/(1D0+PIGGG)
C..LEPTONIC QED VERTEX CORRECTIONS INCLUDING SOFT BREMSSTRAHLUNG
        IF (LPAR(12).GE.1) THEN
          DVRTXV=HSDQDV(X,Q2)*ALP2PI
          DVRTXS=HSDQDS(Q2)*ALP2PI
          DVRTX1=DVRTXV-(Q2+4D0*MEI2)/(Q2-2D0*MEI2)*DVRTXS
          DVRTX2=DVRTXV+
     *          (2D0*(1D0-Y)+Y*Y/2D0)/(1D0-Y-MPRO2*Q2/SHAT/SHAT)*DVRTXS
          ELSE
          DVRTX1=0D0
          DVRTX2=0D0
        ENDIF
C
        IF (LPAR(3).LT.3) THEN
        F1(1,1)=F1(1,1)*(DVRTX1+DVACGG*DVACGG)
        F2(1,1)=F2(1,1)*(DVRTX2+DVACGG*DVACGG)
        ELSE
        F1(1,1)=F1(1,1)*(1D0+DVRTX1)*DVACGG*DVACGG
        F2(1,1)=F2(1,1)*(1D0+DVRTX2)*DVACGG*DVACGG
        ENDIF
        F2EM=F2(1,1)
        GOTO 2000
      ELSE
        WRITE(LUNOUT,200) ILQMOD
        STOP
  200   FORMAT(/,'          WRONG VALUE FOR ILQMOD: ',I8,/
     F          ,' ******** EXECUTION STOPPED IN HSTRF1 ',/)
      ENDIF

1000  CONTINUE
      CALL HSPVER(X,Q2)
      IF (HPOLAR.NE.0D0) THEN
        CALL HSDPVR(X,Q2)
      ENDIF
      IF (INT(HNA).EQ.1.AND.INT(HNZ).EQ.1) THEN
C---PROTON TARGET
C---unpolarized
        QUT=QU+QC+QT
        QDT=QD+QS+QB
        QUBT=QBU+QBC+QBT
        QDBT=QBD+QBS+QBB
C---polarized
        DQUT=DQU+DQC+DQT
        DQDT=DQD+DQS+DQB
        DQUBT=DQBU+DQBC+DQBT
        DQDBT=DQBD+DQBS+DQBB
      ELSEIF (INT(HNA).EQ.1.AND.INT(HNZ).EQ.0) THEN
C---NEUTRON TARGET
C---unpolarized
        QDT=QU+QC+QT
        QUT=QD+QS+QB
        QDBT=QBU+QBC+QBT
        QUBT=QBD+QBS+QBB
C---polarized
        DQDT=DQU+DQC+DQT
        DQUT=DQD+DQS+DQB
        DQDBT=DQBU+DQBC+DQBT
        DQUBT=DQBD+DQBS+DQBB
      ELSE
C---OTHER TARGETS
        IF (INUMOD.EQ.1000) THEN
C---SIMPLE SHADOWING MODEL: FIRST CALCULATE DEUTERON STRUCTURE FUNCTIONS
C   PER NUCLEON
C---unpolarized
          HNRAT=HSNRAT(X)
          QDT=HNRAT*(QU+QC+QT+QD+QS+QB)/2D0
          QUT=HNRAT*QDT
          QDBT=HNRAT*(QBU+QBC+QBT+QBD+QBS+QBB)/2D0
          QUBT=HNRAT*QDBT
C---polarized
          DQDT=HNRAT*(DQU+DQC+DQT+DQD+DQS+DQB)/2D0
          DQUT=HNRAT*DQDT
          DQDBT=HNRAT*(DQBU+DQBC+DQBT+DQBD+DQBS+DQBB)/2D0
          DQUBT=HNRAT*DQDBT
        ELSE 
C---WEIGHTED AVERAGE OF UP AND DOWN, CORRECTION INCLUDED IN PDFs
C---unpolarized
          QUT=(HNZ*(QU+QC+QT)+(HNA-HNZ)*(QD+QS+QB))/HNA
          QDT=((HNA-HNZ)*(QU+QC+QT)+HNZ*(QD+QS+QB))/HNA
          QUBT=(HNZ*(QBU+QBC+QBT)+(HNA-HNZ)*(QBD+QBS+QBB))/HNA
          QDBT=((HNA-HNZ)*(QBU+QBC+QBT)+HNZ*(QBD+QBS+QBB))/HNA
C---polarized
          DQUT=(HNZ*(DQU+DQC+DQT)+(HNA-HNZ)*(DQD+DQS+DQB))/HNA
          DQDT=((HNA-HNZ)*(DQU+DQC+DQT)+HNZ*(DQD+DQS+DQB))/HNA
          DQUBT=(HNZ*(DQBU+DQBC+DQBT)+(HNA-HNZ)*(DQBD+DQBS+DQBB))/HNA
          DQDBT=((HNA-HNZ)*(DQBU+DQBC+DQBT)+HNZ*(DQBD+DQBS+DQBB))/HNA
        ENDIF
      ENDIF

C
C..SELF ENERGIES
      IF (LPAR(7).GE.1) THEN
        CG=HSSRGG(T)
        PIGGG=DREAL(CG)/T
        ELSE
        PIGGG=0D0
      ENDIF
      IF (LPAR(8).EQ.1.AND.LPAR(17).EQ.0) THEN
        CM=HSSRGZ(T)
        SM=DREAL(CM)/T
        AKAPPA=1D0-CW/SW*SM/(1D0+PIGGG)
        SWEFF2=SW2*AKAPPA
        ELSE
        SM=0D0
        AKAPPA=1D0
        SWEFF2=SW2
      ENDIF
      IF (LPAR(9).EQ.1.AND.LPAR(17).EQ.0) THEN
        CZ=HSSRZZ(T)
        PIGZZ=DREAL(CZ)/(T-MZ2)
        ELSE
        PIGZZ=0D0
      ENDIF
      GMUFFQ=1D0/(1D0+PIGZZ)
      BFFQ=SQRT(GMUFFQ)*BTOP4
C
C..REDEFINE FERMION GAUGE BOSON COUPLING CONSTANTS TO INCLUDE SELF
C..ENERGIES
      IF (LPAR(4).EQ.1) THEN
        B0=1D0/4D0/CW/SW
        B=B0
        ELSE
C..NORMALIZED TO G-MU
        B0=MZ/SQRT(AGF0)/4D0
        B=B0
        IF (LPAR(9).GE.1) B=B0*SQRT(1D0-DELTAR)
      ENDIF
      IF (LPAR(2).EQ.1.AND.LPAR(9).GE.1) B=B*BFFQ
      RHO1=B/B0
      DSBS1(1)=DSBOS(1)/(1D0+PIGGG)
      DSBS1(2)=DSBOS(2)*GMUFFQ*RHO1

      VAFI1(2,1,1)=0D0
      VAFI1(2,2,1)=0D0
      VAFI1(2,3,1)=0D0
      VAFI1(2,1,2)=-B
      VAFI1(2,2,2)=B
      VAFI1(2,3,2)=-B
      VAFI1(1,1,1)=1D0
      VAFI1(1,2,1)=-2D0/3D0
      VAFI1(1,3,1)=1D0/3D0
      VAFI1(1,1,2)=B*(4D0*SWEFF2-1D0)
      VAFI1(1,2,2)=B*(1D0-8D0*SWEFF2/3D0)
      VAFI1(1,3,2)=B*(4D0*SWEFF2/3D0-1D0)
 
C..LEPTONIC QED VERTEX CORRECTIONS INCLUDING SOFT BREMSSTRAHLUNG
      IF (LPAR(12).GE.1) THEN
        DVRTXV=HSDQDV(X,Q2)*ALP2PI
        DVRTXS=HSDQDS(Q2)*ALP2PI
        DVRTX1=DVRTXV-(Q2+4D0*MEI2)/(Q2-2D0*MEI2)*DVRTXS
        DVRTX2=DVRTXV+
     *         (2D0*(1D0-Y)+Y*Y/2D0)/(1D0-Y-MPRO2*Q2/SHAT/SHAT)*DVRTXS
        ELSE
        DVRTX1=0D0
        DVRTX2=0D0
      ENDIF

      INDV = 1
      INDA = 2
      IEL=1
      DO 11 IF=1,3
        DO 12 IB1=1,2
          DO 13 IB2=1,2
          FLIND1(INDV,IF,IB1,IB2)=
     *     2D0*(VAFI1(INDV,IF,IB1)*VAFI1(INDV,IF,IB2)
     *         +VAFI1(INDA,IF,IB1)*VAFI1(INDA,IF,IB2))
          FLIND1(INDA,IF,IB1,IB2)=
     *     2D0*(VAFI1(INDV,IF,IB1)*VAFI1(INDA,IF,IB2)
     *         +VAFI1(INDA,IF,IB1)*VAFI1(INDV,IF,IB2))
 13       CONTINUE
 12       CONTINUE
 11       CONTINUE
      DO 21 IVB1=1,2
       DO 22 IVB2=1,2
        DO 23 IFERM=2,3
        AFIJ1(IFERM,IVB1,IVB2)=FLIND1(INDV,IFERM,IVB1,IVB2)
     *  *(PONULL*FLIND1(INDV,IEL,IVB1,IVB2)
     *    -POLARI*FLIND1(INDA,IEL,IVB1,IVB2))
        BFIJ1(IFERM,IVB1,IVB2)=FLIND1(INDA,IFERM,IVB1,IVB2)
     *  *(PONULL*FLIND1(INDA,IEL,IVB1,IVB2)
     *    -POLARI*FLIND1(INDV,IEL,IVB1,IVB2))
        AGIJ1(IFERM,IVB1,IVB2)=FLIND1(INDA,IFERM,IVB1,IVB2)
     *  *(PONULL*FLIND1(INDV,IEL,IVB1,IVB2)
     *    -POLARI*FLIND1(INDA,IEL,IVB1,IVB2))
        BGIJ1(IFERM,IVB1,IVB2)=FLIND1(INDV,IFERM,IVB1,IVB2)
     *  *(PONULL*FLIND1(INDA,IEL,IVB1,IVB2)
     *    -POLARI*FLIND1(INDV,IEL,IVB1,IVB2))
 23   CONTINUE
 22   CONTINUE
 21   CONTINUE
C
      DO 31 IB1=1,2
       DO 32 IB2=1,2
       IF (LPAR(3).LT.3) THEN
C---unpolarized
       F1(IB1,IB2)=
     &   (AFIJ(2,IB1,IB2)*(QUT+QUBT)+AFIJ(3,IB1,IB2)*(QDT+QDBT))
     &     /8D0*DSBOS(IB1)*DSBOS(IB2)*DVRTX1
     &  +(AFIJ1(2,IB1,IB2)*(QUT+QUBT)+AFIJ1(3,IB1,IB2)*(QDT+QDBT))
     &     /8D0*DSBS1(IB1)*DSBS1(IB2)
       F2(IB1,IB2)=2D0*X*(
     &   (AFIJ(2,IB1,IB2)*(QUT+QUBT)+AFIJ(3,IB1,IB2)*(QDT+QDBT))
     &     /8D0*DSBOS(IB1)*DSBOS(IB2)*DVRTX2
     &  +(AFIJ1(2,IB1,IB2)*(QUT+QUBT)+AFIJ1(3,IB1,IB2)*(QDT+QDBT))
     &     /8D0*DSBS1(IB1)*DSBS1(IB2))
       F3(IB1,IB2)=
     &   (BFIJ(2,IB1,IB2)*(QUT-QUBT)+BFIJ(3,IB1,IB2)*(QDT-QDBT))
     &     /4D0*DSBOS(IB1)*DSBOS(IB2)*DVRTX1
     &  +(BFIJ1(2,IB1,IB2)*(QUT-QUBT)+BFIJ1(3,IB1,IB2)*(QDT-QDBT))
     &     /4D0*DSBS1(IB1)*DSBS1(IB2)
C---polarized
       G1(IB1,IB2)=
     &   (BGIJ(2,IB1,IB2)*(DQUT+DQUBT)+BGIJ(3,IB1,IB2)*(DQDT+DQDBT))
     &     /8D0*DSBOS(IB1)*DSBOS(IB2)*DVRTX1
     &  +(BGIJ1(2,IB1,IB2)*(DQUT+DQUBT)+BGIJ1(3,IB1,IB2)*(DQDT+DQDBT))
     &     /8D0*DSBS1(IB1)*DSBS1(IB2)
       G4(IB1,IB2)=2D0*X*(
     &   (AGIJ(2,IB1,IB2)*(DQUT-DQUBT)+AGIJ(3,IB1,IB2)*(DQDT-DQDBT))
     &     /8D0*DSBOS(IB1)*DSBOS(IB2)*DVRTX2
     &  +(AGIJ1(2,IB1,IB2)*(DQUT-DQUBT)+AGIJ1(3,IB1,IB2)*(DQDT-DQDBT))
     &     /8D0*DSBS1(IB1)*DSBS1(IB2))
       G3(IB1,IB2)=
     &   (AGIJ(2,IB1,IB2)*(DQUT-DQUBT)+AGIJ(3,IB1,IB2)*(DQDT-DQDBT))
     &     /4D0*DSBOS(IB1)*DSBOS(IB2)*DVRTX1
     &  +(AGIJ1(2,IB1,IB2)*(DQUT-DQUBT)+AGIJ1(3,IB1,IB2)*(DQDT-DQDBT))
     &     /4D0*DSBS1(IB1)*DSBS1(IB2)
      ELSE
C---unpolarized
       F1(IB1,IB2)=
     &   (AFIJ1(2,IB1,IB2)*(QUT+QUBT)+AFIJ1(3,IB1,IB2)*(QDT+QDBT))
     &     /8D0*DSBS1(IB1)*DSBS1(IB2)*(1D0+DVRTX1)
       F2(IB1,IB2)=2D0*X*
     &   (AFIJ1(2,IB1,IB2)*(QUT+QUBT)+AFIJ1(3,IB1,IB2)*(QDT+QDBT))
     &     /8D0*DSBS1(IB1)*DSBS1(IB2)*(1D0+DVRTX2)
       F3(IB1,IB2)=
     &   (BFIJ1(2,IB1,IB2)*(QUT-QUBT)+BFIJ1(3,IB1,IB2)*(QDT-QDBT))
     &     /4D0*DSBS1(IB1)*DSBS1(IB2)*(1D0+DVRTX1)
C---polarized
       G1(IB1,IB2)=
     &   (BGIJ1(2,IB1,IB2)*(DQUT+DQUBT)+BGIJ1(3,IB1,IB2)*(DQDT+DQDBT))
     &     /8D0*DSBS1(IB1)*DSBS1(IB2)*(1D0+DVRTX1)
       G4(IB1,IB2)=2D0*X*
     &   (AGIJ1(2,IB1,IB2)*(DQUT-DQUBT)+AGIJ1(3,IB1,IB2)*(DQDT-DQDBT))
     &     /8D0*DSBS1(IB1)*DSBS1(IB2)*(1D0+DVRTX2)
       G3(IB1,IB2)=
     &   (AGIJ1(2,IB1,IB2)*(DQUT-DQUBT)+AGIJ1(3,IB1,IB2)*(DQDT-DQDBT))
     &     /4D0*DSBS1(IB1)*DSBS1(IB2)*(1D0+DVRTX1)
      ENDIF
 32   CONTINUE
 31   CONTINUE
      F2EM=X*(AFIJ(2,1,1)*(QUT+QUBT)+AFIJ(3,1,1)*(QDT+QDBT))/4D0

C---INCLUDE LONGITUDINAL STRUCTURE FUNCTION
 2000 CONTINUE
      IF (IFLOPT.GT.0) THEN
        FL=0D0
        Q2L=Q2
        IF (Q2L.LT.TTR) THEN
          Q2L=TTR
          XP21=SNGL(X)*SNGL(SP-MPRO2-MEI2)
          IF (Q2L.GE.XP21) Q2L=XP21
        ENDIF
        CALL HSLUFL(X,Q2L,F2EM,FL)
        IF (LPAR(3).LT.3) THEN
          F1(1,1)=(F2EM-FL)/2D0/X*(DVRTX1+DSBS1(1)*DSBS1(1))
          ELSE
          F1(1,1)=(F2EM-FL)/2D0/X*(1D0+DVRTX1)*DSBS1(1)*DSBS1(1)
        ENDIF
        IF (F1(1,1).LT.0D0) F1(1,1)=0D0
      ENDIF

C---APPLY LOW Q2 SUPPRESSION
C   MOVED TO LYSTFU
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SOFT BREMSSTRAHLUNG FOR ELECTRON-QUARK SCATTERING
C     NO RESTRICTION ON Q2, NO QUARKONIC AND NO INTERFERENCE PARTS
C
      FUNCTION HSDQDV(X,Q2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSPEN,HSFONE,Z1,Z2,Z3,Z4
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
C
      SHAT=SP-MPRO2-MEI2
      DELTA2=DELTA*DELTA
      EEI=EELE
      PEI=PELE
      EEI2=EEI*EEI
      AMFEL=SHAT/2D0-Q2/2D0/X
      BMFEL=Q2/2D0+MEI2
      EEF=(BMFEL*PPRO+AMFEL*PELE)/(PELE*EPRO+PPRO*EELE)
      PEF=SQRT((EEF-MEF)*(EEF+MEF))
      EEF2=EEF*EEF
      Y=Q2/SHAT/X
C
      HSDQDV=0D0
      IF (LPAR(12).EQ.1) THEN
        TAU=DSQRT(Q2*(Q2+4*MEI2))
        BETA=(TAU+Q2)*(TAU+Q2)/4D0/MEI2/Q2
        VAU=(TAU+(BETA+1D0)*Q2)/2D0/(BETA*EEI-EEF)
        Z1=DCMPLX(TAU*(TAU+Q2)/2D0/Q2/MEI2,1D-6)
        Z2=DCMPLX(2D0*TAU/(TAU+Q2),-1D-6)
        DC0=2D0*(Q2+2D0*MEI2)/TAU*DREAL(HSSPEN(Z1)-HSSPEN(Z2))
     *   +2D0-3*DREAL(HSFONE(-Q2,MEI,MEI))
        Z1=DCMPLX((EEI+PEI)*VAU/BETA/MEI2,1D-6)
        Z2=DCMPLX(VAU/(EEI+PEI)/BETA,-1D-6)
        Z3=DCMPLX((EEF+PEF)*VAU/MEI2,1D-6)
        Z4=DCMPLX(VAU/(EEF+PEF),-1D-6)
        DR1=DREAL(HSSPEN(Z1))
        DR2=DREAL(HSSPEN(Z2))
        DR3=DREAL(HSSPEN(Z3))
        DR4=DREAL(HSSPEN(Z4))
        A1S=DABS(VAU*(EEI+PEI)-BETA*MEI2)
        A2S=DABS(VAU*(EEF+PEF)-MEF2)
        A1=A1S*A1S
        A2=A2S*A2S
        DL13=DLOG(BETA*EEI/(EEI+PEI))*
     *   (3D0*DLOG(EEI)+DLOG(EEI+PEI))/2D0-DLOG(EEI/(EEF+PEF))*
     *   (3D0*DLOG(EEI)+DLOG((EEF+PEF)/A2))/2D0
        DL2=DLOG(BETA/EEI*(EEI+PEI))*DLOG(BETA*EEI*(EEI+PEI)
     *   /(VAU-BETA*(EEI+PEI))/(VAU-BETA*(EEI+PEI)))/2D0
        DREG=DR1+DR2-DR3-DR4+DL2+DL13
        IF (A1S.GT.1D-5) THEN
          DL4=DLOG((EEF+PEF)/EEI)*DLOG(EEI*(EEF+PEF)
     *     /(VAU-(EEF+PEF))/(VAU-(EEF+PEF)))/2D0
          DREG=DREG-DL4
          DL13A=4D0*DLOG(MEI/EEI)*DLOG(BETA*A2S/A1S*(EEF+PEF)
     *     /(EEI+PEI))/2D0+DLOG(BETA*EEI/(EEI+PEI))*DLOG(BETA/A1)/2D0
          DREG=DREG+DL13A
     *     +DLOG(EEI/VAU)*DLOG((VAU*(EEF+PEF)-MEF2)
     *     *(VAU-(EEF+PEF))*(EEI+PEI)/(VAU*(EEI+PEI)-BETA*MEI2)
     *     /(VAU-BETA*(EEI+PEI))/(EEF+PEF))
        ELSEIF (A1S.NE.0D0) THEN
          DL4=DLOG((EEF+PEF)/EEI)*DLOG(EEI*(EEF+PEF)
     *     /(VAU-(EEF+PEF))/(VAU-(EEF+PEF)))/2D0
          DREG=DREG-DL4
          DL13A=2D0*DLOG(MEI/EEI)*DLOG(BETA*A2S*(EEF+PEF)
     *     /(EEI+PEI))+DLOG(BETA*EEI/(EEI+PEI))*DLOG(BETA)/2D0
          DREG=DREG+DL13A
     *     +DLOG(EEI/VAU)*DLOG(DABS((VAU*(EEF+PEF)-MEF2)
     *     *(VAU-(EEF+PEF))*(EEI+PEI)/(VAU-BETA*(EEI+PEI))/(EEF+PEF)))
     *     +A1S/BETA/MEI2*DLOG(A1S)
        ELSE
          DL4=DLOG((EEF+PEF)/EEI)*DLOG(EEI*(EEF+PEF))/2D0
          DREG=DREG-DL4
          DL13A=2D0*DLOG(MEI/EEI)*DLOG(BETA*A2S*(EEF+PEF)/(EEI+PEI))
     *     +DLOG(BETA*EEI/(EEI+PEI))*DLOG(BETA)/2D0
          DREG=DREG+DL13A
     *     +DLOG(EEI/VAU)*DLOG(DABS((VAU*(EEF+PEF)-MEF2)
     *     *(EEI+PEI)/(VAU-BETA*(EEI+PEI))/(EEF+PEF)))
          DARG2=1D0-(EEF+PEF)/VAU
          IF (DABS(DARG2).GT.1D-5) THEN
            DREG=DREG+DLOG((VAU-(EEF+PEF))*(VAU-(EEF+PEF)))
     *       *DLOG((EEF+PEF)/VAU)/2D0
            ELSEIF (DARG2.NE.0D0) THEN
            DREG=DREG-DARG2*DLOG(DABS(DARG2)*VAU)
            ELSE
          ENDIF
        ENDIF
        DR00=2D0*(Q2+2D0*MEI2)/TAU*(DREG
     *           +DABS(DLOG(EEI/(EEI+PEI)))**2D0
     *           -DABS(DLOG(EEI/(EEF+PEF)))**2D0)
        HSDQDV = HSDQDV
     *         + 2D0*DLOG(MEI2/4D0/DELTA2)
     *           *(1D0-(Q2+2D0*MEI2)/TAU*DLOG(BETA))
     *           +2D0*EEI/PEI*DLOG((EEI+PEI)/MEI)
     *           +2D0*EEF/PEF*DLOG((EEF+PEF)/MEF)
     *           +DC0+DR00
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     SOFT BREMSSTRAHLUNG FOR ELECTRON-QUARK SCATTERING
C     NO RESTRICTION ON Q2, NO QUARKONIC AND NO INTERFERENCE PARTS
C     (SCALAR PART)
C
      FUNCTION HSDQDS(Q2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSFONE
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
C
      HSDQDS=0D0
      IF (LPAR(12).EQ.1) THEN
        HSDQDS = HSDQDS + 2D0*MEI2/(Q2+4D0*MEI2)
     *          *(DREAL(HSFONE(-Q2,MEI,MEI))-2D0)
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C    THE REST OF THE CODE CONTAINS PARAMETRIZATIONS FOR
C    STRUCTURE FUNCTIONS
C    (COURTESY TO D.Y.BARDIN, G.LEVMAN)
C
      SUBROUTINE STRFBS(X,Q2,AF1,AF2)
C
C   CALLED FOR IPART > 2000
C   CALCULATION OF STRUCTURE FUNCTIONS F1 AND F2
C   FOR LARGE Q2 THE PARTON MODEL RELATIONS ARE USED WITH
C   PARAMETRIZATIONS DETERMINED BY (IPART-2000)
C   FOR SMALL Q2 THE PARAMETRIZATIONS BY STEIN ET AL. AND BRASSE
C   ARE USED (SUBROUTINES FROM A.AKHUNDOV, D. BARDIN ET AL.: TERAD91)
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C...DECLARATIONS FOR HERACLES
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
C...DECLARATIONS FOR STRUCTURE FUNCTIONS PARAMETRIZATIONS
C...SUBROUTINES TAKEN FROM A.AKHUNDOV ET AL. (TERAD91)
      COMMON /HSRADC/ AMP , AMP2,RPI ,RPI2 , ALFA, AML , AML2
      COMMON /HSRADF/ AP,AP2,AP4,AL2,AL4,ALPS,CALPI,W2PIT,W2TR,TTR
      COMMON /HSSINC/ SS,SAP,SAMP2,SW2PIT,SW2MIN,SW2TR
      REAL            SS,SAP,SAMP2,SW2PIT,SW2MIN,SW2TR
C...LOCAL DECLARATIONS
      LOGICAL LFIRST
      DATA W2MIN/1.2321D0/
      DATA LFIRST/.TRUE./
C...
C...INITIALIZATION
      IF (LFIRST) THEN
        LFIRST=.FALSE.
CHS                                                   *** CHECK THIS ***
        IMODEL=0
CHS
*
* SETTING OF CONSTANTS
        RPI =PI
        RPI2=PI*PI
        AML=MEI
        AMP=MPRO
        AML2=MEI2
        AMP2=MPRO2
        ALFA =ALPHA
        ALFAI=ALPHA
        CALPI=1D0/ALFA/PI
        AP  = 2D0*AMP
        AP2 = 2D0*AMP2
        AP4 = 4D0*AMP2
        AL2 = 2D0*AML2
        AL4 = 4D0*AML2
        ALPS= 4D0*AML2*AMP2
*
        S=SP
        SS=SNGL(S-AMP2-AML2)
        SAP = SNGL(AP)
        SAMP2=SNGL(AMP2)
*
C       W2TR =  4D0
C       TTR   = 6D0
C       W2TR =  0D0
C       TTR   = 0D0
* THE TWO FOREGOING PARAMETERS MAY BE USED TO SUBDIVIDE THE 2-DIMENSIO-
* NAL INTEGRATION REGION INTO UP TO 3 PARTS. FOR W2TR=TTR=0D0, THERE
* IS ONLY ONE REGION. CHOOSING THE PARAMETER TTR (IN GEV**2) > 0, ONE
* SPLITS THE REGION INTO TWO PARTS, ONE OF THEM WITH SMALL T=Q'**2
* WHERE THE STRUCTURE FUNCTIONS SHOULD HAVE A SPECIAL BEHAVIOUR.
* IN ADDITION, THIS REGION OF SMALL T CAN BE SPLIT INTO 2 SUBREGIONS
* BY A CHOICE OF THE PARAMETER W2TR (IN GEV**2) > M.PROTON**2.
* THIS CAN BE USED FOR A PROPER TREATMENT OF THE RESONANCE REGION.
*
*     PION THRESHOLD
        IF(IMODEL.EQ.0) THEN
          W2PIT = 1.15184D0
            ELSE
          W2PIT = AMP2*(1D0+1D-10)
        ENDIF
        SW2PIT=SNGL(W2PIT)
        SW2MIN=SNGL(W2MIN)
        SW2TR =SNGL(W2TR )
      ENDIF

C...CALCULATION OF STRUCTURE FUNCTIONS

      W2 = (1D0-X)/X*Q2 + MPRO2
      IF (W2.LT.W2TR) THEN
        CALL BRASSE(Q2,W2,ZF1,ZF2)
        ELSE
        CALL STEIN(Q2,W2,ZF1,ZF2)
      ENDIF

      AF1=ZF1
      AF2=ZF2
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE STEIN(Q2,W2,F1,F2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSRADC/ AMP , AMP2, PI , PI2 , ALFA, AML , AML2
      COMMON /HSRADF/ AP,AP2,AP4,AL2,AL4,ALPS,CALPI,W2PIT,W2TR,TTR
      DIMENSION AN(5)
      DATA AN /1.0621D0, -2.2594D0, 10.54D0, -15.8277D0, 6.7931D0/
C
      X=Q2/(Q2+W2-AMP2)
      X1 = 1D0-X
C
C     STEIN'S FIT OF W2 IN THE LOW X REGION
C
      OS=1D0+W2/Q2
      GE2=1D0/((1D0+.61D0*Q2)*(1D0+2.31D0*Q2)*(1D0+.04D0*Q2))**2
      GM2 = 7.7841d0*GE2
      TAU=Q2/AP4
      W2EL=(GE2+TAU*GM2)/(1D0+TAU)
      ANU=(W2+Q2-AMP2)/AP
      SR=0D0
      DO 2 I=1,5
        SR=SR+AN(I)*(1D0-1D0/OS)**(I+2)
 2    CONTINUE
      F2=(1D0-W2EL)*SR
C
C     F1 = 2.*AMP*W1
C
C   SLAC VALUES OF R
C
      R=.18D0
C
      ANU2=ANU**2
      F1=AP*(1D0+ANU2/Q2)/ANU/(1D0+R)*F2
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE BRASSE(ZT,ZAMF2,ZF1,ZF2)
      IMPLICIT DOUBLE PRECISION (Z)
      COMMON /HSSINC/ S,AP,AMP2,W2PIT,W2MIN,W2TR
C
      FMF2=SNGL(ZAMF2)
      IF(FMF2.LE.W2MIN)   CALL WABC(W2MIN)
      IF(FMF2.GT.W2MIN.AND.FMF2.LE.W2TR) CALL WABC(FMF2)
      IF(FMF2.GT.W2TR )   CALL WABC(W2TR )
      T=SNGL(ZT)
      CALL RF12(T,FMF2,RF1,RF2)
      ZF1=DBLE(RF1)
      ZF2=DBLE(RF2)
      END

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      SUBROUTINE RF12(TB,HM,RF1,RF2)
      COMMON /HSSINC/ S,AP,AMP2,W2PIT,W2MIN,W2TR
C
      CW=3471.16/389383.4
C
C      IF(HM-1.6)1,3,3
      IF ((HM-1.6).LT.0.) THEN
         GOTO 1
      ELSE
         GOTO 3
      ENDIF
 1    CONTINUE
C 1    IF(TB-2.0)2,3,3
      IF ((TB-2.0).LT.0.) THEN
         GOTO 2
      ELSE
         GOTO 3
      ENDIF
2     SRF=0D0
      GO TO 4
3     SRF=.18
4     R=SRF
      XB=S-TB+AMP2-HM
      SN=TB*AMP2/S/XB
      CN=1.-SN
      TN=SN/CN
      ANB=(S-XB)/AP
      ANB2=ANB**2
      VEPS=1./(1.+2.*TN*(1.+ANB2/TB))
C      IF(HM-W2MIN)5,5,6
      IF ((HM-W2MIN).LE.0.) THEN
         GOTO 5
      ELSE
         GOTO 6
      ENDIF
5     STOT=((HM-W2PIT)/(W2MIN-W2PIT))*
     &     SVTOT(TB,W2MIN,VEPS)
      GO TO 7
   6  STOT=SVTOT(TB,HM,VEPS)
7     RF1 =CW*(HM-AMP2)/(1.+VEPS*R)*STOT
      RF2=ANB/AP*TB/(TB+ANB2)*(1.+R)*RF1
      END

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      FUNCTION SVTOT(TB,HM,VEPS)
      COMMON /HSWABC/ A1,B1,C1,A2,B2,C2,A3,B3,C3
      COMMON /HSSINC/ S,AP,AMP2,W2PIT,W2MIN,W2TR
C
C      IF(VEPS-.9)2,1,1
      IF ((VEPS-0.9).LT.0.) THEN
         GOTO 2
      ELSE
         GOTO 1
      ENDIF
1     A=A1
      B=B1
      C=C1
      GO TO 5
C 2    IF(VEPS-.6)4,4,3
 2    CONTINUE
      IF ((VEPS-0.6).LE.0.) THEN
         GOTO 4
      ELSE
         GOTO 3
      ENDIF
3     A=A2
      B=B2
      C=C2
      GO TO 5
4     A=A3
      B=B3
      C=C3
5     GD2=1./(1.+TB/.71)**4
      ANB=(HM+TB-AMP2)/AP
      XMQ=SQRT(ANB**2+TB)
      XMQ0=(HM-AMP2)/AP
      AL=LOG(XMQ/XMQ0)
      SVTOT=GD2*EXP(A+B*AL+C*ABS(AL)**3)
      END

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      SUBROUTINE WABC(W2)
C
C     RESONANCE REGION PARAMETRIZATION
C
      DIMENSION ARW(56),A1(56),B1(56),C1(56),A2(56),B2(56),C2(56),
     *A3(56),B3(56),C3(56)
      COMMON /HSWABC/ AA1,BB1,CC1,AA2,BB2,CC2,AA3,BB3,CC3
C
      DATA ARW /
     +1.110,1.125,1.140,1.155,1.170,1.185,1.200,1.215,1.230,1.245,1.260,
     +1.275,1.290,1.305,1.320,1.335,1.350,1.365,1.380,1.395,1.410,1.425,
     +1.440,1.455,1.470,1.485,1.500,1.515,1.530,1.545,1.560,1.575,1.590,
     +1.605,1.620,1.635,1.650,1.665,1.680,1.695,1.710,1.725,1.740,1.755,
     +1.77,1.79,1.81,1.83,1.85,1.87,1.89,1.91,1.93,1.95,1.97,1.99/
      DATA A1 /
     +5.045,5.126,5.390,5.621,5.913,5.955,6.139,6.178,6.125,5.999,5.769,
     +5.622,5.431,5.288,5.175,5.131,5.003,5.065,5.045,5.078,5.145,5.156,
     +5.234,5.298,5.371,5.457,5.543,5.519,5.465,5.384,5.341,5.328,5.275,
     +5.296,5.330,5.375,5.428,5.478,5.443,5.390,5.333,5.296,5.223,5.159,
     +5.146,5.143,5.125,5.158,5.159,5.178,5.182,5.195,5.160,5.195,5.163,
     +5.172/
      DATA B1 /
     +0.798,1.052,1.213,1.334,1.397,1.727,1.750,1.878,1.887,1.927,2.041,
     +2.089,2.148,2.205,2.344,2.324,2.535,2.464,2.564,2.610,2.609,2.678,
     +2.771,2.890,2.982,3.157,3.188,3.315,3.375,3.450,3.477,3.471,3.554,
     +3.633,3.695,3.804,3.900,4.047,4.290,4.519,4.709,4.757,4.840,5.017,
     +5.015,5.129,5.285,5.322,5.546,5.623,5.775,5.894,6.138,6.151,6.301,
     +6.542 /
      DATA C1 /
     +0.043,0.024,0.000,-.013,-.023,-.069,-.060,-.080,-.065,-.056,-.065,
     +-.056,-.043,-.034,-.054,-.018,-.046,-.015,-.029,-.048,-.032,-.046,
     +-.084,-.115,-.105,-.159,-.164,-.181,-.203,-.220,-.245,-.264,-.239,
     +-.302,-.299,-.318,-.388,-.393,-.466,-.588,-.622,-.568,-.574,-.727,
     +-.665,-.704,-.856,-.798,-1.048,-.980,-1.021,-1.092,-1.313,-1.341,
     +-1.266,-1.473 /
      DATA A2 /
     +-.050,1.082,4.119,3.898,5.990,6.033,6.160,6.219,6.117,5.959,5.451,
     +5.675,5.417,5.238,5.084,4.913,4.458,5.012,5.128,5.140,5.261,5.370,
     +5.416,5.466,5.508,5.578,5.629,5.623,5.555,5.503,5.471,5.419,5.390,
     +5.423,5.396,5.451,5.400,5.446,5.448,5.421,5.308,5.248,5.193,5.099,
     +5.076,5.054,5.064,5.028,5.080,5.076,5.173,5.011,5.136,5.055,5.035,
     +5.059 /
      DATA B2 /
     +4.849,4.120,2.244,2.824,1.257,1.548,1.730,1.805,1.866,1.914,2.272,
     +1.840,2.102,2.221,2.368,2.622,3.146,2.537,2.473,2.533,2.610,2.461,
     +2.570,2.700,2.799,2.929,3.016,3.128,3.215,3.271,3.288,3.328,3.390,
     +3.376,3.567,3.600,3.945,4.046,4.124,4.270,4.667,4.818,4.856,5.069,
     +5.112,5.236,5.331,5.609,5.628,5.819,5.461,6.391,5.995,6.455,6.565,
     +7.063 /
      DATA C2 /
     +-.285,-.212,-.100,-.155,-.006,-.025,-.066,-.086,-.071,-.054,-.084,
     +0.012,-.037,-.046,-.058,-.085,-.131,-.031,-.031,-.033,-.093,0.000,
     +-.052,-.083,-.093,-.104,-.096,-.163,-.173,-.192,-.191,-.202,-.239,
     +-.168,-.286,-.239,-.408,-.402,-.355,-.400,-.603,-.642,-.637,-.773,
     +-.768,-.879,-.941,-1.075,-1.118,-1.223,-.812,-1.581,-1.219,-1.605,
     +-1.055,-4.060 /
      DATA A3 /        -
     +3.024,3.796,6.003,6.339,6.071,6.834,6.166,6.239,6.149,5.988,5.155,
     +5.727,5.454,5.262,5.230,5.153,4.365,5.168,5.162,5.174,5.335,5.378,
     +5.396,5.436,5.564,5.615,5.574,5.560,5.470,5.394,5.317,5.346,5.311,
     +5.314,5.323,5.271,5.315,5.373,5.433,5.235,5.249,5.180,5.071,5.072,
     +5.071,5.036,4.985,4.976,5.021,5.000,5.007,4.980,5.025,5.031,5.018,
     +5.108/
      DATA B3 /
     +6.613,1.935,0.774,0.653,1.017,0.709,1.699,1.768,1.929,1.943,2.703,
     +1.892,1.932,2.244,2.040,1.980,3.206,2.109,2.587,2.479,2.493,2.437,
     +2.581,2.324,2.699,2.971,3.207,3.325,3.375,3.466,3.563,3.366,3.477,
     +3.279,3.691,3.965,4.148,4.106,4.035,4.687,4.738,4.782,5.218,5.050,
     +5.062,4.979,5.532,5.659,5.677,6.081,6.110,6.366,6.433,6.079,6.650,
     +6.644 /
      DATA C3 /
     +-.350,-.005,0.036,0.054,0.053,0.083,-.054,-.059,-.087,-.060,-.158,
     +-.017,0.064,-.050,0.071,0.169,-.155,0.159,-.106,-.046,-.071,0.005,
     +-.089,0.171,-.098,-.182,-.216,-.199,-.192,-.195,-.256,-.152,-.274,
     +-.017,-.373,-.509,-.537,-.425,-.286,-.652,-.631,-.555,-.941,-.719,
     +-.645,-.474,-.960,-1.124,-1.088,-1.359,-1.350,-1.574,-1.921,-.926,
     +-1.845,-2.098 /
C
      W=SQRT(W2)
C
      CALL PARINV(W,ARW,A1,56,AA1)
      CALL PARINV(W,ARW,B1,56,BB1)
      CALL PARINV(W,ARW,C1,56,CC1)
      CALL PARINV(W,ARW,A2,56,AA2)
      CALL PARINV(W,ARW,B2,56,BB2)
      CALL PARINV(W,ARW,C2,56,CC2)
      CALL PARINV(W,ARW,A3,56,AA3)
      CALL PARINV(W,ARW,B3,56,BB3)
      CALL PARINV(W,ARW,C3,56,CC3)
      END

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      SUBROUTINE PARINV(W,ARW,A,N,AA)
      DIMENSION ARW(56),A(56)
      AA=DIVDIF(A,ARW,N,W,1)
      END

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   Interface to the parametrization of ALLM for F2 (1997 update)
C   H.Abramowicz, A.Levy, DESY 97-251 (Dec.1997) (hep-ph/9712415)
C
      SUBROUTINE HSSTAL(X,Q2,F1R,F2R)
      DOUBLE PRECISION X,Q2,F1R,F2R
      double precision xmp2,w2,gdfit10
      character target*1
      SX=SNGL(X)
      SQ2=SNGL(Q2)
c-hs replaced old parametrization by new fit from HERMES
c-hs (21.01.2011)
c      F2E=F2ALLM(SX,SQ2)
      target='p'
      XMP2=0.8802d0
      W2=Q2*(1d0/X-1d0)+XMP2
      F2R=gdfit10(X,Q2,W2,target)
c-hs
      F1R=F2R/2D0/X
      RETURN
      END

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      FUNCTION f2allm(x,q2)             

      REAL M02,M12,LAM2,M22
      COMMON/ALLM/SP,AP,BP,SR,AR,BR,S,XP,XR,F2P,F2R
C  POMERON
      PARAMETER (
     , S11   =   0.28067, S12   =   0.22291, S13   =   2.1979,
     , A11   =  -0.0808 , A12   =  -0.44812, A13   =   1.1709,
     , B11   =   0.60243**2, B12   =   1.3754**2, B13   =   1.8439,
     , M12   =  49.457 )
 
C  REGGEON
      PARAMETER (
     , S21   =   0.80107, S22   =   0.97307, S23   =   3.4942,
     , A21   =   0.58400, A22   =   0.37888, A23   =   2.6063,
     , B21   =   0.10711**2, B22   =   1.9386**2, B23   =   0.49338,
     , M22   =   0.15052 )
C
      PARAMETER ( M02=0.31985, LAM2=0.065270, Q02=0.46017 +LAM2 )
      PARAMETER ( ALFA=112.2, XMP2=0.8802)
C
      W2=q2*(1./x -1.)+xmp2
      W=sqrt(w2)
C
      IF(Q2.EQ.0.) THEN
       S=0.
       Z=1.
C
C   POMERON
C
       XP=1./(1.+(W2-XMP2)/(Q2+M12))
       AP=A11
       BP=B11
       SP=S11
       F2P=SP*XP**AP
C
C   REGGEON
C
       XR=1./(1.+(W2-XMP2)/(Q2+M22))
       AR=A21
       BR=B21
       SR=S21
       F2R=SR*XR**AR

      ELSE

       S=LOG(LOG((Q2+Q02)/LAM2)/LOG(Q02/LAM2))
       Z=1.-X
C
C   POMERON
C
       XP=1./(1.+(W2-XMP2)/(Q2+M12))
       AP=A11+(A11-A12)*(1./(1.+S**A13)-1.)
       BP=B11+B12*S**B13
       SP=S11+(S11-S12)*(1./(1.+S**S13)-1.)
       F2P=SP*XP**AP*Z**BP
C
C   REGGEON
C
       XR=1./(1.+(W2-XMP2)/(Q2+M22))
       AR=A21+A22*S**A23
       BR=B21+B22*S**B23
       SR=S21+S22*S**S23
       F2R=SR*XR**AR*Z**BR
C
      ENDIF

c      CIN=ALFA/(Q2+M02)*(1.+4.*XMP2*Q2/(Q2+W2-XMP2)**2)/Z
c      SIGal=CIN*(F2P+F2R)
c      f2allm=sigal/alfa*(q2**2*(1.-x))/(q2+4.*xmp2*x**2)
      f2allm = q2/(q2+m02)*(F2P+F2R)
 
      RETURN
      END

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

*----------------------------------------------------
*     created by Lara De Nardo, October 28, 2010
*     output is F2
*----------------------------------------------------

c-hs  renamed: name may not contain '-'
c     HS (21.01.201)
      DOUBLE PRECISION FUNCTION gdfit10(dx,dQ2,dW2,target)
      
      IMPLICIT NONE

      CHARACTER TARGET*1
      DOUBLE PRECISION PARAM(23)
      DOUBLE PRECISION F2P,F2R
      DOUBLE PRECISION XP,XR,CP,CR,AP,AR,BP,BR,t,dQ2,dx,dW2

c-hs  this was undefined in the original version
c     added by HS (21.01.2011)
      double precision aprmass
      aprmass=938.28D-3
c-hs
      
      if (target.eq.'p') then

         param(  1)=       0.504875
         param(  2)=      34.046472
         param(  3)=       0.066144
         param(  4)=       1.528063
         param(  5)=       0.065270
         param(  6)=      -0.119775
         param(  7)=      -0.490448
         param(  8)=       1.287539
         param(  9)=       0.984561
         param( 10)=       2.604510
         param( 11)=       1.404179
         param( 12)=       0.358903
         param( 13)=       0.114424
         param( 14)=       1.120841
         param( 15)=       0.336422 
         param( 16)=       1.113998
         param( 17)=       0.598313
         param( 18)=      -7.908678
         param( 19)=      12.523172
         param( 20)=       0.112027
         param( 21)=       1.051316
         param( 22)=       2.739961
         param( 23)=       1.818436

      elseif(target.eq.'d') then

         param(  1 )=               0.455872      
         param(  2 )=               0.531614      
         param(  3 )=               2.719053      
         param(  4 )=              17.494116      
         param(  5 )=               0.065270      
         param(  6 )=              -0.135367      
         param(  7 )=              -0.909263      
         param(  8 )=               0.677523      
         param(  9 )=              -6.934579      
         param(  10)=               0.400460      
         param(  11)=               0.067910      
         param(  12)=               0.175499      
         param(  13)=               0.706134      
         param(  14)=               2.062847      
         param(  15)=               0.654724      
         param(  16)=               6.036018      
         param(  17)=               1.173878      
         param(  18)=              -1.196548      
         param(  19)=               9.159108      
         param(  20)=               0.222962      
         param(  21)=               0.575077      
         param(  22)=              31.152266      
         param(  23)=               1.488237      
      else
         write(*,*) 'target unknown- STOP'
         stop
      endif
      
      t=dlog((dQ2+param(4))/param(5)) / log(param(4)/param(5))
      t=dlog(t)

      if (dQ2.ne.0.d0) then         
         XP=(dQ2+param(2))/(dQ2/dx+param(2))
         XR=(dQ2+param(3))/(dQ2/dx+param(3))
      else
         XP=1.d0+(dW2-aprmass*aprmass)/(dQ2+param(2))
         XR=1.d0+(dW2-aprmass*aprmass)/(dQ2+param(3))
         XP=1.d0/XP
         XR=1.d0/XR
      endif

      AP=param( 6)+(param( 6)-param( 7))*(1.d0/(1.d0+t**param( 8))-1.d0)
      CP=param(12)+(param(12)-param(13))*(1.d0/(1.d0+t**param(14))-1.d0)

      BP=param( 9)+param(10)*t**param(11)
      AR=param(15)+param(16)*t**param(17)
      BR=param(18)+param(19)*t**param(20)
      CR=param(21)+param(22)*t**param(23)

      F2P=CP*(XP**AP)*(1.d0-dX)**BP
      F2R=CR*(XR**AR)*(1.d0-dX)**BR
      
      gdfit10=(F2P+F2R)*dQ2/(dQ2+param(1))
      
      return
      end

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      FUNCTION RSLAC(X,Q2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      TETA=1.+12.*(Q2/(Q2+1.))*(.015625/(.015625+X*X))
C     Q2TH=5.*(1.-X)**5
      R1=.0672/LOG(Q2/.04)*TETA+.4671/(Q2**4+12.97458)**.25
C     R2=.0635/LOG(Q2/.04)*TETA+.5747/Q2-.3534/(Q2*Q2+.09)
C     R3=.0599/LOG(Q2/.04)*TETA+.5088/SQRT((Q2-Q2TH)**2+4.4441)
      RSLAC=R1
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   Interface to the parametrization of BK for F2, F1
C   B.Badelek and J.Kwiecinski, Nucl.Phys. B295 (1992) 263
C   See comments in HSF2BK below
C
      SUBROUTINE HSSTBK(DX,DQ2,ZF1,ZF2)
      DOUBLE PRECISION DX,DQ2,ZF1,ZF2
      LOGICAL LFIRST
      DATA LFIRST /.TRUE./
c     Example of usage of JKBB F2 routine 'f2pd';
c     for the conditions of usage read the instruction in f2pd.
c     When implementing f2pd in your program use "DATA" statement
c     instead of "READ" to store matrix h(2,20,20)
C     F2PD RENAMED TO HSF2BK (H.S.)
c
      common//h(2,20,20)
c
c     choose "mode, x, Q2" and read in the data from the "mode" file
c
C     PRINT*,'MODE, X, Q2='
C     READ*,MODE,X,Q2
      IF (LFIRST) THEN
      LFIRST=.FALSE.
      lode=21
      mode=3
      CALL HSBKIN
C     DO 200 I=1,2
C     DO 200 IQ=1,20
C     READ(MODE,*)H(I,IQ, 1),H(I,IQ, 2),H(I,IQ, 3),H(I,IQ, 4),H(I,IQ, 5)
C     READ(MODE,*)H(I,IQ, 6),H(I,IQ, 7),H(I,IQ, 8),H(I,IQ, 9),H(I,IQ,10)
C     READ(MODE,*)H(I,IQ,11),H(I,IQ,12),H(I,IQ,13),H(I,IQ,14),H(I,IQ,15)
C     READ(MODE,*)H(I,IQ,16),H(I,IQ,17),H(I,IQ,18),H(I,IQ,19),H(I,IQ,20)
C200   CONTINUE
      ENDIF
c
c     and calculate proton & deuteron F2 (f2p, f2d respectively)
c     for a defined x,Q2 by calling the subroutine f2pd
c
      X=REAL(DX)
      Q2=REAL(DQ2)
      CALL HSF2BK(MODE,X,Q2,F2P,F2D)
      ZF2=DBLE(F2P)
      ZF1=ZF2/2D0/DX
      RETURN
      end

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C

      SUBROUTINE HSF2BK(MODE,X,Q2,F2P,F2D)
C     SUBROUTINE F2PD(MODE,X,Q2,F2P,F2D)
c
c --- A code to calculate the nucleon structure function F2 in the
c --- low q2 and low x region within the GVMD inspired model described in
c --- B.Badelek and J.Kwiecinski, B295 (1992) 263.
c --- The structure functions F2 corresponding to the QCD improved
c --- parton model are obtained from the interpolation formula based on
c --- Tchebyshev polynomials. The vector meson contributions are calculated
c --- directly. Parton model contributions are based on the MRS
c --- parton distributions, A.D.Martin, W.J.Stirling and R.G.Roberts,
c --- PR D47(1993) 145, Phys.Lett B306(1993) 145.
c
c >>>>>>  Model is valid for
c                             nu > 10 GeV,
c                             10**(-5) < x < 0.1,
c                             Q2< 1000 GeV2                  <<<<<<
c
c --- The following data files exist:
c --- mode=1. f2mod1.dat, corresponds to the unshadowed D- parametrisation;
c --- mode=2. f2mod2.dat, to the shadowed D- parametrisation and R=2GeV**(-1);
c --- mode=3. f2mod3.dat, to the shadowed D- parametrisation and R=5GeV**(-1);
c --- mode=4. f2mod4.dat, to the unshadowed D0 set of partons;
c --- mode=5. f2mod5.dat, to the D-';
c --- modes 1-4 have calculations done in the LO and mode=5 in the NLO of QCD.
c --- OBS 1:  mode=1,2,3,4 have a narrower validity range: 10**(-4) < x < 0.1
c --- OBS 2:  x validity range applies in fact to the {\bar x} variable,
c             cf.the publication for the definition; validity in the Bjorken x
c             extends to even lower x values (e.g. to x=0 for photoproduction).
c --- User has to attach a data file containing the interpolation
c --- coefficients for a chosen parametrisation.
c
c --- Subroutine parameters:
c --- mode: an integer defining the input file (e.g.mode=2 for file f2mod2.dat),
c --- q2,x: the usual kinematic variables,
c --- f2p=F2 for the proton and f2d=F2 for the deuteron (returned values).
c --- Example of the usage of the subroutine is given above.
c
c --- History:
c     Oct. 25th, 1993; introduced mode=5 which corresponds to the MRS D-'
c                      parametrisation; NLO evolution (LO before);
c                      extension of model validity down to x=1.E-5
c                      (old limit: x=1.E-4)
c
      common//h(2,20,20)
      dimension alama(5),ftwo(2)
      data nmax, pmass / 20, 0.938272/
      data w02, ymax, q20, q2f / 1.2, 9.2103, 1., 1000./
      data alama/3*0.2304, 0.1732, 0.2304/
c
      qbar=q2+w02
      cbar=q2/qbar
      viu=q2/x
      viug=viu/(2.*pmass)
      viubar=viu+w02
      xbar=qbar/viubar
c
      alam=alama(mode)*alama(mode)
      ql0=alog(q20/alam)
      qlf=alog(q2f/alam)
      qlp=ql0+qlf
      qlm=qlf-ql0
      aq=acos((2.*alog(qbar/alam)-qlp)/qlm)
      if (mode.gt.4) ymax=11.5129
      ax=acos((2.*alog(1./xbar)-ymax)/ymax)
c
c
      do 500 j=1,2
      f2=0.
      do 400 nx=1,nmax
      anx=float(nx)-1.
      tx=cos(anx*ax)
      do 300 nq=1,nmax
      anq=float(nq)-1.
      tq=cos(anq*aq)
      f2=f2+tx*tq*h(j,nq,nx)
300   continue
400   continue
      ftwo(j)=4./float(nmax)/float(nmax)*f2
500   continue
c
c --- vector meson contribution
c
      call sigvmes(viug,sigro,sigfi)
      call vmesnuc(q2,sigro,sigfi,fv)
c
      f2p=cbar*ftwo(1)+fv
      f2d=cbar*ftwo(2)+fv
c
      return
      end
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      subroutine sigvmes(viu,sigro,sigfi)
c
c-----calculates sigro, sigfi at given energy viu using the energy dependent
c-----piN and KN total cross-sections
c
      data      amk,     ampi,pmbnag,sigpi0, sigk0,   czad,sigpi1,sigpi2
     *    /0.493646, 0.139567,  2.56, 13.52, 14.17, 0.0102, 22.77, 2.35/
      data sigk1, alfpi1, alfpi2, alfk1
     *    /14.75,  0.369,   0.37, 0.515/
c
      gampi=viu/ampi
      gamk=viu/amk
      gapa=alog(gampi)
      gaka=alog(gamk)
      spip=sigpi1*viu**(-alfpi1)
      spip=spip+sigpi0*(1.+czad*gapa*gapa)
      spim=sigpi2*viu**(-alfpi2)
      skp=sigk1*viu**(-alfk1)
      skp=skp+sigk0*(1.+czad*gaka*gaka)
      sigro=pmbnag*spip
      sigfi=pmbnag*(2.*skp-spip-spim)
      return
      end
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      subroutine vmesnuc(q2,sigro,sigfi,fvmes2)
c
c-----calculates the vector meson contribution to the nucleon (!) f2.
c
      dimension vmass(3),vcoupl(3),sigin(3)
      data pi/3.1415926/, vmass/0.7683, 0.78195, 1.019412/
      data vcoupl/1.98, 21.07, 13.83/
c
      sigin(1)=sigro
      sigin(2)=sigro
      sigin(3)=sigfi
      sumv=0.
      do 100 l=1,3
      sv=vmass(l)*vmass(l)
      denv=pi*vcoupl(l)*(q2+sv)*(q2+sv)
      sumv=sumv+sv*sv*sigin(l)/denv
100   continue
      fvmes2=q2/4./pi*sumv
      return
      end
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C...INPUT FOR BADELEK/KWIECINSKI PARAMETRIZATION
C   NEW VERSION VALID FOR EXTENDED KINEMATIC RANGE (DEC. 93)

      SUBROUTINE HSBKIN
      COMMON //H1(2,20,20)
      DIMENSION H(2,20,20)
      DATA (H(1,1,N),N=1,20)/
     .316.3296644740772,   485.8975656910705,  259.6136606319242,
     .117.3307786807664,    32.16818321591835,
     . 10.64768876117214,    2.848065566825837, -0.7215608659036665,
     .  1.187429709239199,  -0.7932390539536663,
     .  0.4714339692858772, -0.1831410742125502, 1.0170687772727414E-03,
     .  8.4179584992745884E-02,-0.1117787295555190,
     .  9.1866157500536292E-02,-6.6935139980532398E-02,
     .  4.6752933908006442E-02,-2.0748727933986192E-02,
     .  1.1678824404468792E-02/
      DATA (H(1,2,N),N=1,20)/
     .254.3089293501125,   418.7477285501654,  230.8904850564271,
     . 93.64447846743268,   32.48257976526893,
     .  6.879896644486573,   2.936391479892251, -4.9063167715756220E-02,
     . -2.4767278893035123E-02, 0.3955014763248987,
     . -0.4533908636720629,  0.3829556835451522,-0.2505685158029276,
     .  0.1182806533021176, -2.6630772372001804E-02,
     . -2.9048052301551008E-02,4.5916476409516690E-02,
     . -4.2200780634649624E-02,3.1875262968723642E-02,
     . -1.6984452043858281E-02/
      DATA (H(1,3,N),N=1,20)/
     . -5.205517203742071,  -7.613923888380566, -2.139391763920144,
     . -0.8778185462976319, -0.6064977299189174,
     .  0.1964926533236207, -0.3638238753020571, 0.2493255700900238,
     . -0.1134190355714483, -1.6085746861335640E-02,
     .  9.3563287149740028E-02,-0.1181647883494726,
     .  0.1050190308667975, -7.1955618163729553E-02,
     .  3.9320613351923985E-02,
     . -1.2009815438625413E-02,-2.7604838214640864E-03,
     .  9.1599608346141909E-03,-9.3558304809775159E-03,
     .  6.7374042415057438E-03/
      DATA (H(1,4,N),N=1,20)/
     .  1.790899045505892,   3.053839222466831,   1.716760795147843,
     .  0.9389845862757266,  0.3552022319675952,
     .  3.0907987227430154E-02, 0.1089572387552265,
     . -8.7460054261984849E-02, 6.0878745078748803E-02,
     . -1.9558925664357729E-02,
     . -1.3279064772305173E-02, 3.0702202030882964E-02,
     . -3.4410989446204764E-02, 2.8108807989698516E-02,
     . -1.8858381016765550E-02,
     .  9.2738410819645021E-03,-2.9301436541826585E-03,
     . -1.0451584670259943E-03, 2.0371195301416922E-03,
     . -2.1740196816710108E-03/
      DATA (H(1,5,N),N=1,20)/
     . -1.125432508063247,  -1.868590239295427, -1.068872737618474,
     . -0.5180069312658861, -0.1529311946681504,
     . -3.1686389286933733E-02,-3.8870789237635146E-02,
     .  3.3116940661247520E-02,-2.6234759141101501E-02,
     .  1.1450631650217661E-02,
     .  1.7056679752367680E-03,-8.8379478159300044E-03,
     .  1.1203033071420430E-02,-9.9215808824230287E-03,
     .  7.0031753409176388E-03,
     . -4.1486685438180383E-03, 1.8415654248478509E-03,
     . -6.4946774876729796E-06,-2.4179162736877685E-04,
     .  6.9752055132026684E-04/
      DATA (H(1,6,N),N=1,20)/
     .  0.7859860842825152,  1.268601328031006,  0.7198376963883954,
     .  0.3337274972420689,  8.3403187833354106E-02,
     .  2.6098164746585593E-02, 2.0937160406365934E-02,
     . -1.6449241435357990E-02, 1.2500976388980837E-02,
     . -4.9348555980043545E-03,
     . -1.2847024081706998E-03, 3.7857560154279305E-03,
     . -4.2505327971354598E-03, 3.5652230238075717E-03,
     . -2.2208002904650896E-03,
     .  1.4900325694985228E-03,-6.9601213564879949E-04,
     . -8.4776803337891749E-05,-7.4556062576086525E-05,
     . -2.5171560143955646E-04/
      DATA (H(1,7,N),N=1,20)/
     . -0.3458275355070987, -0.5560796396633200, -0.3178436561855143,
     . -0.1432259256486573, -3.1191476091947454E-02,
     . -1.3781369721746755E-02,-8.7187118647218875E-03,
     .  7.6580822379931938E-03,-5.8729154797341954E-03,
     .  2.2832713852619971E-03,
     .  5.9895272062973058E-04,-1.6405192637641120E-03,
     .  1.7203622377909946E-03,-1.3468184125611596E-03,
     .  7.4291499573687678E-04,
     . -4.9990725304443526E-04, 2.2481399476398239E-04,
     .  6.2691047846582758E-05, 4.9961656278537477E-05,
     .  9.0458495053711477E-05/
      DATA (H(1,8,N),N=1,20)/
     . -0.1085455474187870, -0.1655087767157225,-8.9218106237787184E-02,
     . -4.3775058549888573E-02,-1.2630632154353540E-02,
     . -1.2351333903456096E-03,-2.5217085547781231E-03,
     . -4.0266867779742427E-04, 1.2853944564771603E-03,
     . -1.2973691451185036E-03,
     .  8.2431922683277483E-04,-1.3689526356241208E-04,
     . -3.2085154615627093E-04, 4.3498770490028007E-04,
     . -4.6939545426076912E-04,
     .  2.4900431886548259E-04,-1.5150459753039519E-04,
     .  1.0209815802393671E-04, 1.4812796521087268E-05,
     . -7.7146939394380492E-06/
      DATA (H(1,9,N),N=1,20)/
     .  0.3438552509916081,  0.5389308757927751,  0.3013583141188885,
     .  0.1389264048545785,  3.2556417603769212E-02,
     .  1.0326741528094614E-02, 8.2785611717296069E-03,
     . -3.8201384005378077E-03, 1.4216865195395041E-03,
     .  7.3202674486942658E-04,
     . -1.6676915387988519E-03, 1.1644228111795703E-03,
     . -4.2808616876797423E-04,-1.5871321459509406E-05,
     .  4.2473000724665851E-04,
     . -2.0456601380010781E-04, 1.7306727773199913E-04,
     . -2.1646329348452131E-04,-4.2131609511711239E-05,
     . -3.1138105795789064E-05/
      DATA (H(1,10,N),N=1,20)/
     . -0.2666284983331815, -0.4193036547186229, -0.2361279192119010,
     . -0.1072196731354234, -2.3509715113688884E-02,
     . -9.1084669244087511E-03,-6.4858604017878814E-03,
     .  3.6723705418607977E-03,-1.7590346279415437E-03,
     . -2.4900690073546540E-04,
     .  1.2769022045740275E-03,-1.0431994895220903E-03,
     .  4.9766983254701150E-04,-1.0532742411831749E-04,
     . -2.6526758757535447E-04,
     .  1.4025948676539830E-04,-1.3392153823798330E-04,
     .  1.7337473250470669E-04, 2.5950665371892346E-05,
     .  2.9269620241400226E-05/
      DATA (H(1,11,N),N=1,20)/
     . -2.1641655953025585E-03,-2.7280006177934300E-03,
     . -7.6672594403688454E-05,-1.2905786933555744E-03,
     . -1.8107522578205785E-03,
     .  7.6144454771083392E-04, 7.9082559269037676E-05,
     . -4.9312363482609612E-04, 3.6864351208980488E-04,
     . -1.0675557950985583E-04,
     . -7.1443297748211056E-05, 1.2751039991787558E-04,
     . -9.5767267168048481E-05, 3.7682934529567924E-05,
     .  6.3722345818114810E-06,
     . -2.7053378158035031E-05, 2.6923218234557537E-05,
     . -1.6266385457451837E-05, 9.6887873322577205E-06,
     . -2.0352265252877557E-06/
      DATA (H(1,12,N),N=1,20)/
     .  0.2318584783547141,     0.3640874804517640,
     .  0.2033123115140910,     9.3531266922188376E-02,
     .  2.2266520869560401E-02,
     .  7.0983081770853293E-03, 5.4250929901558780E-03,
     . -2.6499440057001076E-03, 1.1986084342017536E-03,
     .  2.4374299263204615E-04,
     . -9.6913977487532241E-04, 7.4712172998354589E-04,
     . -3.4173275196392854E-04, 7.7988600742127338E-05,
     .  1.9221963234209065E-04,
     . -7.1059036124745883E-05, 7.3105232343063069E-05,
     . -1.2575779554997444E-04,-3.4568904139992196E-05,
     . -2.3956753201861310E-05/
      DATA (H(1,13,N),N=1,20)/
     . -0.2551126027054218,    -0.4011298971323523,
     . -0.2249831008736070,    -0.1025452462911845,
     . -2.3367491092718078E-02,
     . -8.4353929707860059E-03,-6.0355086331645017E-03,
     .  3.2918676298121958E-03,-1.6157986916136473E-03,
     . -1.6732176165457908E-04,
     .  1.1080362119742959E-03,-9.1425250635234992E-04,
     .  4.5147914913725408E-04,-1.2111266682289053E-04,
     . -2.0931542116845973E-04,
     .  9.2481020975093035E-05,-9.6705156647598560E-05,
     .  1.4932068209479170E-04, 3.1787023323479374E-05,
     .  2.8251344672428605E-05/
      DATA (H(1,14,N),N=1,20)/
     .  7.8075093211252500E-02, 0.1230915801936174,
     .  6.9635656708936525E-02, 3.1196152781560376E-02,
     .  6.4933981650981401E-03,
     .  2.9483403022493872E-03, 1.8859674419684562E-03,
     . -1.2311906604126394E-03, 6.7465998701432196E-04,
     . -1.3952585110695949E-05,
     . -3.5970934183974964E-04, 3.3258202427922218E-04,
     . -1.8329315501507838E-04, 5.9591646797594392E-05,
     .  6.0591765883761565E-05,
     . -3.4972881170466804E-05, 3.7860821889947154E-05,
     . -5.1447163347695427E-05,-6.3823604161686913E-06,
     . -9.7702386617232954E-06/
      DATA (H(1,15,N),N=1,20)/
     .  0.1438319789368712,     0.2258070543013799,
     .  0.1259919888371545,     5.7927903372604848E-02,
     .  1.3837074640170796E-02,
     .  4.4006172533111952E-03, 3.3508418267598121E-03,
     . -1.6264681715208084E-03, 7.2990997529672634E-04,
     .  1.5856831464849014E-04,
     . -6.0311392583147860E-04, 4.6228050942256978E-04,
     . -2.0937636824327253E-04, 4.5772587830595487E-05,
     .  1.2175076251223997E-04,
     . -4.5887893430150516E-05, 4.6800637373385023E-05,
     . -7.8846348562105379E-05,-2.0974026719137824E-05,
     . -1.4884327106003285E-05/
      DATA (H(1,16,N),N=1,20)/
     . -0.2352880318790727,    -0.3698028202050644,
     . -0.2070853502842321,    -9.4516184750886790E-02,
     . -2.1804421877123153E-02,
     . -7.6657171316025243E-03,-5.5279236751719818E-03,
     .  2.9434083313969099E-03,-1.4231371166309967E-03,
     . -1.7553506755579768E-04,
     .  1.0117564418676352E-03,-8.2273032895695508E-04,
     .  3.9978820619553932E-04,-1.0373709599029986E-04,
     . -1.9443139412460696E-04,
     .  8.3282390876259020E-05,-8.6901792439193437E-05,
     .  1.3622111793675436E-04, 3.0072538256514326E-05,
     .  2.5791484490536205E-05/
      DATA (H(1,17,N),N=1,20)/
     .  0.1340924543603120,     0.2109023755965740,
     .  0.1183732053308312,     5.3772887774393171E-02,
     .  1.2123308000737075E-02,
     .  4.5401943808229436E-03, 3.1665582152355960E-03,
     . -1.7804878814693143E-03, 8.9478244297346624E-04,
     .  6.9306597750372190E-05,
     . -5.8572627300951668E-04, 4.9315828092124907E-04,
     . -2.4876897171992880E-04, 6.9666256248087523E-05,
     .  1.0909715545727975E-04,
     . -5.0493693698681596E-05, 5.3342317261339052E-05,
     . -8.0306058696077234E-05,-1.5569325013133603E-05,
     . -1.5231393378166980E-05/
      DATA (H(1,18,N),N=1,20)/
     .  6.9332981920138960E-02, 0.1088471578287316,
     .  6.0712918472132376E-02, 2.7907682015571262E-02,
     .  6.6789757841758525E-03,
     .  2.1203903268691982E-03, 1.6112807904989852E-03,
     . -7.8054412276863090E-04, 3.5061304367886554E-04,
     .  7.5895844327584269E-05,
     . -2.8959938417618394E-04, 2.2195241775818285E-04,
     . -1.0058448101794640E-04, 2.2138700937261226E-05,
     .  5.8425362973301124E-05,
     . -2.1867229809486621E-05, 2.2394644758000066E-05,
     . -3.7938907355849618E-05,-1.0128219344834088E-05,
     . -7.1864603050829368E-06/
      DATA (H(1,19,N),N=1,20)/
     . -0.2094335438293607,    -0.3291457167597506,
     . -0.1842463995640689,    -8.4099819707091857E-02,
     . -1.9455355646261942E-02,
     . -6.8077143644552556E-03,-4.9091727524013329E-03,
     .  2.6028771767728045E-03,-1.2566865356079515E-03,
     . -1.5769782415207982E-04,
     .  8.9727152319202583E-04,-7.2825324762207682E-04,
     .  3.5325983584691317E-04,-9.1551284431029190E-05,
     . -1.7270491484759411E-04,
     .  7.3384593677747381E-05,-7.6736948361805104E-05,
     .  1.2092337784544062E-04, 2.6905347839721798E-05,
     .  2.2944945547767758E-05/
      DATA (H(1,20,N),N=1,20)/
     .  0.1764241083215092,     0.2773161943401925,
     .  0.1553194347090352,     7.0814195640269395E-02,
     .  1.6292501141539833E-02,
     .  5.7896162067337869E-03, 4.1404075170619767E-03,
     . -2.2254584723765772E-03, 1.0853581943524741E-03,
     .  1.2300838389880883E-04,
     . -7.5875689530633292E-04, 6.2123569222668283E-04,
     . -3.0428300054637019E-04, 8.0505847821641867E-05,
     .  1.4493110464378429E-04,
     . -6.2785749290831476E-05, 6.5866331486169732E-05,
     . -1.0272471275959145E-04,-2.2158585089963858E-05,
     . -1.9501578701913792E-05/
      DATA (H(2,1,N),N=1,20)/
     . 315.3555650808810,      487.1523194474932,   259.6465107900499,
     . 116.4940461036096,      33.21540403363760,
     .  9.831397315937661,      3.218098014278122,
     . -0.6817154010258749,     0.9200946379397350,
     . -0.4866104630931506,
     .  0.2383907776459113,    -5.5988745449716497E-02,
     . -3.9178743903936125E-02, 7.2407546263734249E-02,
     . -7.9561173025909792E-02,
     .  5.8214987013573908E-02,-4.0877887493026279E-02,
     .  2.9363377988324933E-02,-1.1066386701897477E-02,
     .  6.9549005586796583E-03/
      DATA (H(2,2,N),N=1,20)/
     .254.3818289740736,      418.5275186398955,   231.1498494922002,
     . 93.49937336431225,      32.41140777734330,
     .  7.152067200268769,      2.587646605884788,
     .  0.2299337074089342,    -0.1499072195943851,
     .  0.3716051170379300,
     . -0.3395322792371598,     0.2481939288576829,
     . -0.1407742546302435,     5.1804775191776264E-02,
     .  9.2148510769340189E-04,
     . -2.9919145279846772E-02, 3.4634000756206363E-02,
     . -2.7546874461335447E-02, 2.0486183161881020E-02,
     . -9.8246558595267586E-03/
      DATA (H(2,3,N),N=1,20)/
     . -5.222753288349486,     -7.569075206092052,
     . -2.195918358368663,     -0.8299124712324074,
     . -0.6210124599979283,
     .  0.1616486397900743,    -0.2883521656087849,
     .  0.1622022683048591,    -4.7157807034083544E-02,
     . -4.4675383118636211E-02,
     .  8.5796938584630368E-02,-8.8589980315254132E-02,
     .  6.9274423787098518E-02,-4.2103020391994220E-02,
     .  1.9617402469487285E-02,
     . -2.7872258929548117E-03,-5.1104801020305181E-03,
     .  7.2200152134024414E-03,-6.7037539882323981E-03,
     .  4.1539590637957897E-03/
      DATA (H(2,4,N),N=1,20)/
     .  1.795655239403132,      3.042305197786875,
     .  1.731510650345955,      0.9245925080343382,
     .  0.3635697359701548,
     .  3.4770334331769240E-02, 9.1904811816968064E-02,
     . -6.2496362250700219E-02, 3.7414104058495009E-02,
     . -4.6536547847940985E-03,
     . -1.6938004389858978E-02, 2.5626445644071382E-02,
     . -2.4666614104332649E-02, 1.8140589733304925E-02,
     . -1.0864519578660768E-02,
     .  4.4550266994742695E-03,-5.5627123898424732E-04,
     . -1.4020650074635405E-03, 1.7417533968980771E-03,
     . -1.4283687571557594E-03/
      DATA (H(2,5,N),N=1,20)/
     . -1.126853083480120,     -1.865303194677924,
     . -1.073102970239374,     -0.5135741672826667,
     . -0.1562121942761721,
     . -3.1632271398583714E-02,-3.4828615756586772E-02,
     .  2.5849725963632295E-02,-1.8438367321697530E-02,
     .  5.5489037460329903E-03,
     .  4.2475786663977816E-03,-8.3539999305615269E-03,
     .  8.6737311523375181E-03,-6.8302270862958138E-03,
     .  4.1674577161149381E-03,
     . -2.2165421159097054E-03, 6.7951339675775276E-04,
     .  3.7790330032298208E-04,-3.2587438691072173E-04,
     .  5.0242900463093502E-04/
      DATA (H(2,6,N),N=1,20)/
     .  0.7864312243380061,      1.267607980471678,
     .  0.7211224905345020,      0.3323237488854016,
     .  8.4593623825662339E-02,
     .  2.5804678659169554E-02, 1.9951615107583772E-02,
     . -1.4277573204395454E-02, 9.9344292404631553E-03,
     . -2.7679879559438641E-03,
     . -2.4419729832566975E-03, 3.9342797751438191E-03,
     . -3.6125020449757010E-03, 2.6267220260744386E-03,
     . -1.2588652671254420E-03,
     .  7.8114114261018274E-04,-2.1923412548235279E-04,
     . -2.8171486406146958E-04, 5.5576500129959715E-06,
     . -2.0435510342071948E-04/
      DATA (H(2,7,N),N=1,20)/
     . -0.3459714042534403,    -0.5557681331266985,
     . -0.3182485416066988,    -0.1427718698835976,
     . -3.1613247839279500E-02,
     . -1.3618499107698309E-02,-8.4773043724060808E-03,
     .  6.9945056072709231E-03,-5.0253027746090357E-03,
     .  1.5091899675752218E-03,
     .  1.0657727843081655E-03,-1.7743351883985245E-03,
     .  1.5658858287012425E-03,-1.0636890094487570E-03,
     .  4.2176733617930925E-04,
     . -2.4920864200303785E-04, 4.2659968471616453E-05,
     .  1.4690131523422416E-04, 8.2876816011666943E-06,
     .  8.0107476865693879E-05/
      DATA (H(2,8,N),N=1,20)/
     . -0.1084980810772992,    -0.1656088807985460,
     . -8.9087376386210287E-02,-4.3924074265384052E-02,
     . -1.2482577849853671E-02,
     . -1.3065344936064699E-03,-2.5793676278996006E-03,
     . -1.9668365825810749E-04, 1.0041118427874141E-03,
     . -1.0240666464075246E-03,
     .  6.4550834159022521E-04,-6.7816159392562099E-05,
     . -2.8648780311074406E-04, 3.4992809259848387E-04,
     . -3.6282081763720117E-04,
     .  1.6185832607705037E-04,-8.4279742396267165E-05,
     .  6.8736558185330570E-05, 3.3157017208563402E-05,
     . -5.9345829973611333E-06/
      DATA (H(2,9,N),N=1,20)/
     .  0.3438393579778490,     0.5389635985445374,
     .  0.3013154171255550,     0.1389757547543309,
     .  3.2504722326641548E-02,
     .  1.0355376472440389E-02, 8.2913668159813562E-03,
     . -3.8847017032601163E-03, 1.5153449687520053E-03,
     .  6.3614404340656736E-04,
     . -1.6010333670451369E-03, 1.1339006548746402E-03,
     . -4.3426387082593887E-04, 9.4917241579444531E-06,
     .  3.8947712947007999E-04,
     . -1.7455910817432596E-04, 1.4873062177360617E-04,
     . -2.0375990653644737E-04,-4.9620779049301657E-05,
     . -3.1206664618377173E-05/
      DATA (H(2,10,N),N=1,20)/
     . -0.2666231040709893,    -0.4193145044508198,
     . -0.2361136868986530,    -0.1072360910235140,
     . -2.3491759328072356E-02,
     . -9.1194205360529169E-03,-6.4882354386895188E-03,
     .  3.6927084796933584E-03,-1.7902600447898725E-03,
     . -2.1552806369962549E-04,
     .  1.2524701436853841E-03,-1.0306463554775114E-03,
     .  4.9803908648807727E-04,-1.1278777536861013E-04,
     . -2.5364485080527975E-04,
     .  1.2999604384653926E-04,-1.2521969763617665E-04,
     .  1.6865130501965398E-04, 2.8883403189738144E-05,
     .  2.9123344756608840E-05/
      DATA (H(2,11,N),N=1,20)/
     . -2.1660351807644072E-03,-2.7243364001492458E-03,
     . -8.1435272540315466E-05,-1.2851179429450414E-03,
     . -1.8169381076812112E-03,
     .  7.6549364276813217E-04, 7.9308966343408952E-05,
     . -4.9954504267920127E-04, 3.7904763654213230E-04,
     . -1.1838963598946141E-04,
     . -6.2601422280735362E-05, 1.2255928196267087E-04,
     . -9.5336761775792050E-05, 3.9830327528115821E-05,
     .  2.5565310038369642E-06,
     . -2.3563899255824005E-05, 2.3840774617620322E-05,
     . -1.4538768861509441E-05, 8.5705688071448483E-06,
     . -1.9310428292372785E-06/
      DATA (H(2,12,N),N=1,20)/
     .  0.2318591581225479,     0.3640861941488209,
     .  0.2033139209296905,     9.3529466274897100E-02,
     .  2.2268614225486147E-02,
     .  7.0968688618638989E-03, 5.4251737911641374E-03,
     . -2.6479100630039140E-03, 1.1951494982517930E-03,
     .  2.4776160118547664E-04,
     . -9.7230098779385421E-04, 7.4901618532146762E-04,
     . -3.4205697579785480E-04, 7.7389287135517535E-05,
     .  1.9346583836977588E-04,
     . -7.2237973141837467E-05, 7.4188032169166135E-05,
     . -1.2638124806418650E-04,-3.4150899972068985E-05,
     . -2.4009607069959709E-05/
      DATA (H(2,13,N),N=1,20)/
     . -0.2551128802029598,    -0.4011294003881104,
     . -0.2249836551117676,    -0.1025446709226060,
     . -2.3368165547998303E-02,
     . -8.4349197469735719E-03,-6.0355616009768629E-03,
     .  3.2912148917943074E-03,-1.6146531788159599E-03,
     . -1.6869615025303481E-04,
     .  1.1091480240063034E-03,-9.1495592957547284E-04,
     .  4.5164275442913885E-04,-1.2095138684732984E-04,
     . -2.0972035015917538E-04,
     .  9.2876388617450081E-05,-9.7082065154004003E-05,
     .  1.4954255380035639E-04, 3.1633576812386239E-05,
     .  2.8274623990525821E-05/
      DATA (H(2,14,N),N=1,20)/
     .  7.8075235376370969E-02, 0.1230913449954655,
     .  6.9635856975744840E-02, 3.1195988326370767E-02,
     .  6.4935814396565525E-03,
     .  2.9482200390526297E-03, 1.8859683953848872E-03,
     . -1.2309708030527523E-03, 6.7428314752218884E-04,
     . -1.3492726150342946E-05,
     . -3.6008765732794689E-04, 3.3283120260798764E-04,
     . -1.8336057445542679E-04, 5.9548638500715326E-05,
     .  6.0723096244574487E-05,
     . -3.5104120561733622E-05, 3.7990281705310375E-05,
     . -5.1524502850817283E-05,-6.3273339575380029E-06,
     . -9.7794431641205986E-06/
      DATA (H(2,15,N),N=1,20)/
     .  0.1438318812767813,     0.2258072042788636,
     .  0.1259919072418009,     5.7927929457090424E-02,
     .  1.3837061884805410E-02,
     .  4.4006090730575940E-03, 3.3508729547767256E-03,
     . -1.6265535015216533E-03, 7.3003221412595660E-04,
     .  1.5842369175691944E-04,
     . -6.0299655557464298E-04, 4.6220200068652175E-04,
     . -2.0935573839129444E-04, 4.5786152006797103E-05,
     .  1.2170767562643602E-04,
     . -4.5845125254479122E-05, 4.6757378864284945E-05,
     . -7.8820587763527001E-05,-2.0992870871451951E-05,
     . -1.4881256790120993E-05/
      DATA (H(2,16,N),N=1,20)/
     . -0.2352879476493046,    -0.3698029441522917,
     . -0.2070853081732636,    -9.4516163586214274E-02,
     . -2.1804469324534519E-02,
     . -7.6656614736123753E-03,-5.5279712518447083E-03,
     .  2.9434525948627788E-03,-1.4231753036239052E-03,
     . -1.7549895555224489E-04,
     .  1.0117316711568970E-03,-8.2271541058852500E-04,
     .  3.9978841328645314E-04,-1.0374450655089577E-04,
     . -1.9441647123489869E-04,
     .  8.3269069670810831E-05,-8.6888408161558546E-05,
     .  1.3621365092123308E-04, 3.0078193672955222E-05,
     .  2.5790964380796467E-05/
      DATA (H(2,17,N),N=1,20)/
     .  0.1340923727789033,     0.2109024936138492,
     .  0.1183731759744742,     5.3772849717184081E-02,
     .  1.2123377974822865E-02,
     .  4.5401199766745579E-03, 3.1666142102052796E-03,
     . -1.7805202206368379E-03, 8.9479294892310841E-04,
     .  6.9308027891171892E-05,
     . -5.8573458849584875E-04, 4.9316714692979459E-04,
     . -2.4877815808025492E-04, 6.9673037398030878E-05,
     .  1.0909108718445778E-04,
     . -5.0490121639002391E-05, 5.3339261680800218E-05,
     . -8.0305072546186575E-05,-1.5570187753166633E-05,
     . -1.5231913605119485E-05/
      DATA (H(2,18,N),N=1,20)/
     .  6.9333064756997220E-02, 0.1088470386760217,
     .  6.0712944150708274E-02, 2.7907726960858526E-02,
     .  6.6788959412898632E-03,
     .  2.1204735436808044E-03, 1.6112198446343334E-03,
     . -7.8051463146365807E-04, 3.5061166757438840E-04,
     .  7.5881089234121166E-05,
     . -2.8957883517329916E-04, 2.2193434253996885E-04,
     . -1.0057132561036559E-04, 2.2131486682815538E-05,
     .  5.8428725478745573E-05,
     . -2.1867576938304706E-05, 2.2394106721910335E-05,
     . -3.7937567302082411E-05,-1.0129117078619442E-05,
     . -7.1855069091663770E-06/
      DATA (H(2,19,N),N=1,20)/
     . -0.2094336295158042,    -0.3291455937230821,
     . -0.1842464246997682,    -8.4099868361146289E-02,
     . -1.9455269978078640E-02,
     . -6.8078031852446442E-03,-4.9091081178749376E-03,
     .  2.6028476632793045E-03,-1.2566882207296329E-03,
     . -1.5767791759710169E-04,
     .  8.9724593667888015E-04,-7.2823118838405930E-04,
     .  3.5324469912063817E-04,-9.1543537856009860E-05,
     . -1.7270751195724047E-04,
     .  7.3383856366975696E-05,-7.6735115913733591E-05,
     .  1.2092115497396109E-04, 2.6906925838931395E-05,
     .  2.2943793016116745E-05/
      DATA (H(2,20,N),N=1,20)/
     .  0.1764241977012391,     0.2773160660648814,
     .  0.1553194604372024,     7.0814247132839963E-02,
     .  1.6292410696607291E-02,
     .  5.7897098459946868E-03, 4.1403394690611027E-03,
     . -2.2254279622666902E-03, 1.0853609872971168E-03,
     .  1.2298599076304767E-04,
     . -7.5872864510608314E-04, 6.2121142480548227E-04,
     . -3.0426659804017821E-04, 8.0497599697570645E-05,
     .  1.4493354353069614E-04,
     . -6.2784615612362718E-05, 6.5863971103739356E-05,
     . -1.0272209329061569E-04,-2.2160470633835152E-05,
     . -1.9500311080302588E-05/
      DO 203 I=1,2
      DO 202 IQ=1,20
      DO 201 N=1,20
        H1(I,IQ,N)=H(I,IQ,N)
 201  CONTINUE
 202  CONTINUE
 203  CONTINUE
      END

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   Interface to the parametrization of DL for F2
C   A.Donnachie and P.Landshoff, DAMTP 93-23 and M/C-TH 93/11
C   Published in Z.Phys.C61 (1994) (hep-ph/9305319) 

      SUBROUTINE HSSTDL(X,Q2,ZF1,ZF2)
      DOUBLE PRECISION X,Q2,ZF1,ZF2,HSF2DL
      ZF2=HSF2DL(Q2,X)
      ZF1=ZF2/2D0/X
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c this function calculates the proton structure function at all x for
c q**2 < 10 GeV**2 using the parametrisation of Donnachie and Landshoff
c in DAMTP 93-23 and M/C-TH 93/11
C
C   RENAMED TO HSF2DL BY H.S.
c
c-----------------------------------------------------------------------
      FUNCTION HSF2DL(Q2,X)
      implicit double precision (a-h,o-z)
      c=0.219744
      b=0.278516
      x0=0.071348
      d=15.8769
      fm2=0.302408
      xi1=x*(1.+16./q2)
      xi2=x*(1.+1.7/q2)
      if (xi1.ge.1.) f1=0.
      if (xi1.lt.1.) f1=0.027/(xi1**0.0808)*(1.-xi1)**7*q2/(q2+6.25)
      if (xi2.ge.1.) f2=0.
      if (xi2.lt.1.) f2=2.0*c/(9.*xi2**0.0808)*(1.-xi2)**7*q2/(q2+1.)
      aa=1./(1./(2.038*c)-0.173)
      bb=0.489*b
      phi=q2/(q2+bb)
      phis=q2/(q2+aa)
      xi=x*(1.+0.28/q2)
      if (xi.ge.1.0) then
      HSF2DL=0.
      return
      else
      ht=d*x**2*(1.-xi)**2/(1.+q2/fm2)
      if(xi.lt.x0) then
      f3=10.*c/(9.*xi**0.0808)*phis
      s=f3+f2+f1+ht
      HSF2DL=S+B*XI**0.4525*PHI
      else
      ru=3.0*x0/(1.-x0)+0.4525
      rd=4.0*x0/(1.-x0)+0.4525
      buu=2.0/(x0**ru*((1.-x0)**3/0.4525-1./ru+3.0*x0/(1.+ru)
     &-3.0*x0**2/(2.+ru)+x0**3/(3.+ru))+1./ru-3.0/(1.+ru)
     &+3.0/(2.+ru)-1.0/(3.+ru))
      bu=buu*x0**(ru-0.4525)*(1.-x0)**3
      bdd=1.0/(x0**rd*((1.-x0)**4/0.4525-1./rd+4.0*x0/(1.+rd)
     &-6.0*x0**2/(2.+rd)+4.0*x0**3/(3.+rd)-x0**4/(4.+rd))+1/rd
     &-4.0/(1.+rd)+6.0/(2.+rd)-4.0/(3.+rd)+1.0/(4.+rd))
      bd=bdd*x0**(rd-0.4525)*(1.-x0)**4
      u=buu*xi**ru*(1.-xi)**3*phi
      d=bdd*xi**rd*(1.-xi)**4*phi
      rs=9.0*x0/(1.-x0)+0.4525
      rss=7.*x0/(1.-x0)-0.0808
      bss=(b-4.*bu/9.-bd/9.)*x0**(0.4525-rs)/(1.-x0)**9
      css=c/(x0**(0.0808+rss)*(1.-x0)**7)
      f3=(bss*xi**rs*(1.-xi)**2*phi+(10./9.)*css*xi**rss*phis)
     &*(1.-xi)**7
      s=f3+f2+f1+ht
      HSF2DL=4.*U/9.+D/9.+S
      return
      endif
      endif
      return
      end
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C
C***********************************************************************
C
C     NEW SUBROUTINES FOR THE INCLUSION OF
C     ELASTIC AND QUASI-ELASTIC EP SCATTERING
C     AND NUCLEAR SHADOWING
C
C
C     NEW ROUTINES:      IN HSOURCEE: HSELG1 (F)
C                                     HSELG2 (F)
C                                     HSEL22 (F)
C                                     HSELK1 (F)
C                                     HSELK2 (F)
C                                     HSELCO (F)
C                                     HSSGEL (F)
C                                     HSXMAX (F)
C                                     HSNRAT (F)
C                                     HSINIL (S)
C                                     HSFIE0 (S)
C                                     HSFIEL (S)
C                                     HSDELX (S)
C                                     D01AJF (S)
C
C***********************************************************************
C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSINIL(FUN,EPSO,NBIN2,NDO2,SIG2,SIG2E,XX2)
C
C   INITIALIZATION FOR ELASTIC EP EVENT GENERATION
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C      INTEGER MINPTS,MAXPTS
      EXTERNAL FUN
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSINTL/ XL,XU
      DIMENSION XX2(50,1)
      DATA BLOW/0D0/, BUP/1D0/
C
C   PARAMETERS
C
      ACCURA=EPSO
      NDO2=NBIN2
      SIG2=0D0
      SIG2E=0D0
      IFAIL=1
      NDIMEN=1
      CALL D01AJF(FUN,BLOW,BUP,ACCURA,RESULT,ACCFIN,IFAIL)
      IF(IFAIL.NE.0.OR.IPRINT.GT.1)
     &  WRITE(LUNTES,'(A)')
     &  ' D01AJF DID NOT MEET REQUIRED ACCURACY IN ELASTIC EP '
      SIG2=RESULT
      SIG2E=ACCFIN
C
      IF(IPRINT.GT.1) THEN
        WRITE(LUNOUT,'(///A,5X,1PE12.4,A,1PE12.4,A)')
     *        ' CROSS SECTION VALUE SIG2L (WITH ERROR ESTIMATE):',
     *        SIG2, ' +/- ', SIG2E, '  NB'
        WRITE(LUNTES,'(A,1PD10.1)') ' RELATIVE ACCURACY REQUIRED:',
     &        ACCURA
      ENDIF
      G=0D0
      DG=1D0/DFLOAT(NDO2)
      DO 5 I=1,NDO2
        XX2(I,1)=DFLOAT(I)*DG
  5   CONTINUE
      IF(IPRINT.GT.1) THEN
        WRITE(LUNTES,'(A,/,4(5(1PD15.5)/))')
     &                     ' XX2(I)', (XX2(IG,1),IG=1,NDO2)
      ENDIF
      IF(IPRINT.GT.1) WRITE(LUNTES,'(A)') ' HSINIL FINISHED'
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSELG1(ARGUM)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSINTL/ XL,XU
C
      X=1D0
      GSP=SP-MEI2-MPRO2
      Q2MNY=YMIN*GSP
      GL=-1D0/MAX(Q2MIN,Q2MNY)
      YMAXX=(1D0-4D0*MEI2*MPRO2/GSP/GSP)/(1D0+X*MPRO2/GSP+MEI2/GSP)
c      ymaxx=1d0
      ymax1=q2max/gsp
      GU=-1D0/(MIN(YMAX,ymax1,YMAXX)*GSP)
      DG=GU-GL
      G=GL+ARGUM*DG
      Q2=-1D0/G
      IF(IPRINT.GT.20) WRITE(LUNTES,'(A,4D15.6)')
     &                      ' HSELG1: G, Q2',G,Q2
      HSELG1=Q2**2*HSEL22(Q2)*DG
c      write(*,*) q2,hselg1
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSELG2(XARG)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      DIMENSION XARG(1)
C
      X=1D0
      Z=XARG(1)
      GSP=SP-MEI2-MPRO2
      Q2MNY=YMIN*GSP
      GL=-1D0/MAX(Q2MIN,Q2MNY)
      YMAXX=(1D0-4D0*MEI2*MPRO2/GSP/GSP)/(1D0+X*MPRO2/GSP+MEI2/GSP)
      ymax1=q2max/gsp
      GU=-1D0/(MIN(YMAX,ymax1,YMAXX)*GSP)
      DG=GU-GL
      G=GL+Z*DG
      Q2=-1D0/G
      IF(IPRINT.GE.20) WRITE(LUNTES,'(A,3D15.6)')
     &                      ' HSELG2: Z, G, Q2',Z,G,Q2
      HSELG2=Q2**2*HSEL22(Q2)*DG
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSEL22(Q2)
C
C   D2SIG*DQ2 FOR ELASTIC EP INCLUDING SOFT AND VIRTUAL CORRECTIONS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSKPXY/ XX,Y
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
C      DATA NEVERR /0/
C      LOGICAL LERR
C      DATA LERR /.TRUE./
C
      GSP=SP-MEI2-MPRO2
      XX=1D0
      Y=Q2/GSP
      IF(IPRINT.GT.20)
     &  WRITE(LUNTES,'(A/2(1PD13.5),F8.3,2I3)')
     *         ' HSEL22: SP, Q2, POLARI,LLEPT,LQUA',
     *         SP,Q2,POLARI,LLEPT,LQUA
      HSEL22=HSSGEL(Q2,LLEPT,POLARI)/SP
C
C     IF(HSEL22.LT.0D0) THEN
C       NEVERR=NEVERR+1
C       IF (NEVERR.LT.20) THEN
C        WRITE(LUNTES,'(A,/,2(1PD13.5),2I3,F8.3)')
C    +     ' HSEL22: Q2, HSEL22, LLEPT, LQUA, POLARI',
C    +      Q2, HSEL22, LLEPT, LQUA, POLARI
C       ELSEIF (LERR) THEN
C        LERR=.FALSE.
C        WRITE(LUNTES,'(A,I3,A)')
C    &     ' ERROR HSEL22 < 0 HAS OCCURED ',NEVERR,
C    &     ' TIMES, NO FURTHER WARNINGS ARE PRINTED'
C       ELSE
C       ENDIF
C       HSEL22=0D0
C     ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND OF KP TERM (ELASTIC RADIATIVE TAIL)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSELK1(X)
C
C  X(1) -->  XX
C  X(2) -->  G=-1/Q**2   -->  Q**2=-1/G-->  Y
C  X(3) -->  LOG(2*K.P)
C  X(4) -->  TS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSIRCX/ XIRDEL
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSKPXY/ XX,Y
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSXSLM/ XSMIN,XSCUT
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSPSPC/ IPHSPC
      DIMENSION X(4)
      COMPLEX(kind=8) HSSRGG,CG
C
C---QUASI-ELASTIC SCATTERING
      XS=1D0
C---CHOOSE VALUE OF Y
      GS=SP-MEI2-MPRO2
      YMAXX=(1D0-4D0*MEI2*MPRO2/GS/GS)/(1D0+2D0*MEI*MPRO/GS)
      IF (ICUT.LT.3) THEN
        GMIN=-1D0/(Q2MIN/XMAX/GS)
        GMAX=-1D0/DMIN1(1D0,YMAXX)
        ELSEIF(ICUT.EQ.3) THEN
C---CUT IN Y
        GMIN=-1D0/DMAX1(Q2MIN/XMAX/GS,YMIN)
        GMAX=-1D0/DMIN1(1D0,YMAXX,YMAX)
c...cut on lepton scattering angle (fixed target, me=0)
        elseif (icut.eq.4) then
        q2min=4d0*eele*eele*dsin(themin/180d0*pi/2d0)**2
     *       /(1d0+2d0*eele/mpro*dsin(themin/180d0*pi/2d0)**2)
        q2max=4d0*eele*eele*dsin(themax/180d0*pi/2d0)**2
     *       /(1d0+2d0*eele/mpro*dsin(themax/180d0*pi/2d0)**2)
        gmin=-1d0/q2min
        gmax=-1d0/q2max
        ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSELK1'
        STOP
      ENDIF
      GACT=GMIN+X(1)*(GMAX-GMIN)
      Y=-1D0/GACT
      XA=1D0
      CALL HSDELX(XA,Y)
C---X-VALUE
      XXMAX1=HSXMAX(Y)
      XXHH1=1D0-Y-4D0*MEI2*MPRO2/GS/GS
      XXMNY1=(XXHH1+DSQRT(XXHH1*XXHH1-4D0*Y*Y*MEI2*MPRO2/GS/GS))
     &      /2D0/Y/MPRO2*GS
      XXMINY=MEI2/MPRO2/XXMNY1
      XXMIN=DMAX1(XMIN,Q2MIN/Y/GS,XXMINY)
      XXMAX=DMIN1(XMAX,XXMAX1)
      XX=XXMIN+(XXMAX-XXMIN)*X(2)
      Q2=XX*Y*GS
      CALL HSFIVC(XX,Y)
C
      IF(IPRINT.GT.30) THEN
        WRITE(LUNTES,230)SP,XX,Y,Q2,XSMIN,XSCUT,DELEPS,DELTA
230     FORMAT(//,
     F         ' ***************************************************',/
     F        ,' SP = ',1PD12.3,/
     F        ,' X = ',D12.6,'   Y = ',D12.6,/
     F        ,' Q2 = ',D12.6,/
     F        ,' XSMIN = ',D17.11,'   XSCUT = ',D17.11,/
     F        ,' DELEPS = ',D12.6,'   DELTA = ',D14.8)
      ENDIF
C
      CALL HSFCMS(XX,Y,XS)
      IF(IPHSPC.EQ.1) THEN
        HSELK1=0D0
        RETURN
      ENDIF
      CALL HSFLAB(XS)
      CALL HSLZK1(ZMIN,ZMAX)
      IF(IPHSPC.EQ.1) THEN
        HSELK1=0D0
        RETURN
      ENDIF
C---SUBSTITUTION V = LN(A1)
      A1MAX=2D0*OMEGA*(EEL-PEL*ZMIN)
      IF ((ZMAX.GE.0.9999D0).AND.(EEL/MEI.GT.1D3)) THEN
        A1MIN=2D0*OMEGA*MEI2/2D0/EEL
      ELSE
        A1MIN=2D0*OMEGA*(EEL-PEL*ZMAX)
      ENDIF
      if (iphspc.eq.1) then
        hselk1=0d0
        return
      endif
      if ((a1max.le.a1min).or.(a1min.le.0d0).or.(a1max.le.0d0)) then
        hselk1=0d0
        return
      endif
      VMAX=DLOG(A1MAX)
      VMIN=DLOG(A1MIN)
      V=VMIN+(VMAX-VMIN)*X(3)
      A1=DEXP(V)
C---NO SUBSTITUTION FOR TS
      CALL HSLTS1(A1,XX,Y,XS,TSMIN,TSMAX,TSM,TSP,CFKP)
      IF(IPHSPC.EQ.1)THEN
       HSELK1=0D0
       RETURN
      ENDIF
      TS=TSMIN+(TSMAX-TSMIN)*X(4)
      SQGRM2=-CFKP*(TS-TSM)*(TS-TSP)
      IF (SQGRM2.LE.0D0) THEN
        HSELK1=0D0
        RETURN
      ENDIF
      SQGRAM=DSQRT(DABS(SQGRM2))
      CALL HSFIV1(XX,Y,XS,A1,TS)
C
      A2=2D0*DKPS
      
      RUNALP=1D0
      IF (LPAR(3).GE.3) THEN
        CG=HSSRGG(TS)
        RUNALP=1D0/(1D0+DREAL(CG)/TS)
      ENDIF
      R1=4D0*GX*(
     &     -(T+6D0*MEI2)/(A1-TS)/(A1-TS)/A1
     &     +1D0/(A1-TS)/A1
     &     + 2D0/(A1+A2)/A1
     &     -8D0*MEI2*MEI2/(A1*A2+TS*TS)/(A1+A2)/A1
     &     +4D0*MEI2*MEI2/(A1*A1+TS*TS)/A1/A1
     &     -2D0*MEI2/(A1-TS)/A1/A1 )
C
      FAC1 = -(GS*GS + GU*GU - GX*(GS+GU) -4D0*MEI2*MPRO2)
      FAC2 = -GS*GX + MPRO2*T
      FAC3 = (-2D0*GS*GU + GX*(GS+GU-GX))*2D0*MEI2
      FAC4 = GU*(GX-GU)*2D0*MEI2
      R2=4D0*(
     &     -FAC1/(A1+A2-TS)*(1D0/(A1+A2)+1D0/(A1-TS))/A1
     &     +FAC2/(A1-TS)/(A1-TS)/A1
     &     +FAC3/(A1*A2+TS*TS)/(A1+A2)/A1
     &     -2D0*MPRO2/(A1+A2)/A1
     &     +2D0*MEI2*MPRO2/(A1-TS)/(A1-TS)/A1
     &     +2D0*MEI2*MPRO2/(A1-TS)/A1/A1
     &     -MPRO2/(A1-TS)/A1
     &     +FAC4/(A1*A1+TS*TS)/A1/A1  )
      DO 20 IFL=1,12
        CQP(IFL)=0D0
 20   CONTINUE
      CALL HSFIE0(-TS,F1EL,F2EL)
      CQP(12)=(F1EL*R1+F2EL*R2)*RUNALP*RUNALP
      SUMME=CQP(12)
      HSELK1=SUMME*Y*2D0*SX1NRM/SQGRAM
     *      *(VMAX-VMIN)*A1*(TSMAX-TSMIN)
     *      *(XXMAX-XXMIN)*(GMAX-GMIN)/GACT**2
ctest
c      if (iprint.eq.-7) then
c         write(*,*) ' runalp = ',runalp
c         write(*,*) ' f1el,f2el = ',f1el,f2el
c         write(*,*) ' fac1-4 = ',fac1,fac2,fac3,fac4
c         write(*,*) ' r1,r2 = ',r1,r2
c         write(*,*) ' hselk1 = ',hselk1
c      endif
ctest
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND OF KPS TERM (ELASTIC RADIATIVE TAIL)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSELK2(X)
C
C  X(1) -->  XX
C  X(2) -->  G=-1/Q**2   -->  Q**2=-1/G-->  Y
C  X(3) -->  LOG(2*K.PS)
C  X(4) -->  TS
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSKPXY/ XX,Y
      COMMON /HSXSLM/ XSMIN,XSCUT
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSIRCX/ XIRDEL
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSPSPC/ IPHSPC
      DIMENSION X(4)
      COMPLEX(kind=8) HSSRGG,CG

C---QUASI-ELASTIC SCATTERING
      XS=1D0
C---CHOOSE VALUE OF Y
      GS=SP-MEI2-MPRO2
      YMAXX=(1D0-4D0*MEI2*MPRO2/GS/GS)/(1D0+2D0*MEI*MPRO/GS)
      IF (ICUT.LT.3) THEN
        GMINL=DLOG(Q2MIN/XMAX/GS)
        GMAXL=DLOG(DMIN1(1D0,YMAXX))
        ELSEIF(ICUT.EQ.3) THEN
C---CUT IN Y
        GMINL=DLOG(DMAX1(Q2MIN/XMAX/GS,YMIN))
        GMAXL=DLOG(DMIN1(1D0,YMAXX,YMAX))
c...cut on lepton scattering angle (fixed target, me=0)
        elseif (icut.eq.4) then
        q2min=4d0*eele*eele*dsin(themin/180d0*pi/2d0)**2
     *       /(1d0+2d0*eele/mpro*dsin(themin/180d0*pi/2d0)**2)
        q2max=4d0*eele*eele*dsin(themax/180d0*pi/2d0)**2
     *       /(1d0+2d0*eele/mpro*dsin(themax/180d0*pi/2d0)**2)
        gmin=-1d0/q2min
        gmax=-1d0/q2max
        ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSELK2'
        STOP
      ENDIF
      GACT=GMINL+X(1)*(GMAXL-GMINL)
      Y=DEXP(GACT)
      XA=1D0
      CALL HSDELX(XA,Y)
C---X-VALUE
      XXMAX1=HSXMAX(Y)
      XXHH1=1D0-Y-4D0*MEI2*MPRO2/GS/GS
      XXMNY1=(XXHH1+DSQRT(XXHH1*XXHH1-4D0*Y*Y*MEI2*MPRO2/GS/GS))
     &      /2D0/Y/MPRO2*GS
      XXMINY=MEI2/MPRO2/XXMNY1
      XXMIN=DMAX1(XMIN,Q2MIN/Y/GS,XXMINY)
      XXMAX=DMIN1(XMAX,XXMAX1)
      GXMIN=-1D0/XXMIN**2
      GXMAX=-1D0/XXMAX**2
      GX=GXMIN+(GXMAX-GXMIN)*X(2)
      XX=SQRT(-1D0/GX)
      Q2=XX*Y*GS
      CALL HSFIVC(XX,Y)
C
      IF(IPRINT.GT.30) THEN
        WRITE(LUNTES,230)SP,XX,Y,XSMIN,XSCUT,DELEPS,DELTA
230     FORMAT(' ***************************************************',/
     F        ,' SP = ',1PD12.3,/
     F        ,' X = ',D12.6,'   Y = ',D12.6,/
     F        ,' XSMIN = ',D17.11,'   XSCUT = ',D17.11,/
     F        ,' DELEPS = ',D12.6,'   DELTA = ',D14.8,/
     F       ,' ***************************************************',//)
      ENDIF
C
      CALL HSFCMS(XX,Y,XS)
      IF(IPHSPC.EQ.1) THEN
        HSELK2=0D0
        RETURN
      ENDIF
      CALL HSFLAB(XS)
      CALL HSLZK2(ZMIN,ZMAX)
      IF(IPHSPC.EQ.1) THEN
        HSELK2=0D0
        RETURN
      ENDIF
C---SUBSTITUTION V = LN(A2)
      A2MAX=2D0*OMEGA*(ES-PS*ZMIN)
      IF ((ZMAX.GE.0.9999D0).AND.(ES/MEF.GT.1D3)) THEN
        A2MIN=2D0*OMEGA*MEF2/2D0/ES
      ELSE
        A2MIN=2D0*OMEGA*(ES-PS*ZMAX)
      ENDIF
      VMAX=DLOG(A2MAX)
      VMIN=DLOG(A2MIN)
      V=VMIN+(VMAX-VMIN)*X(3)
      A2=DEXP(V)
C---NO SUBSTITUTION FOR = TS
      CALL HSLTS2(A2,XX,Y,XS,TSMIN,TSMAX,TSM,TSP,CFKPS)
      IF(IPHSPC.EQ.1)THEN
       HSELK2=0D0
       RETURN
      ENDIF
      GS00=2D1
      IF ((ABS(TSMIN).GT.GS00).OR.(ABS(TSMAX).GT.GS00)) THEN
        HSELK2=0D0
        RETURN
      ENDIF
c      write(*,*) ' '
c      write(*,*) ' tsmin = ',tsmin
c      write(*,*) ' tsmax = ',tsmax
      TSMINA=DLOG((-GS00-TSMIN)/TSMIN)/GS00
      TSMAXA=DLOG((-GS00-TSMAX)/TSMAX)/GS00
      TSA=TSMINA+(TSMAXA-TSMINA)*X(4)
      TS=-GS00/(1D0+DEXP(GS00*TSA))
c      write(*,*) ' tsmina = ',tsmina
c      write(*,*) ' tsmaxa = ',tsmaxa
c      write(*,*) ' tsa = ',tsa
c      write(*,*) ' ts = ',ts
      SQGRM2=-CFKPS*(TS-TSM)*(TS-TSP)
      IF (SQGRM2.LE.0D0) THEN
        HSELK2=0D0
        RETURN
      ENDIF
      SQGRAM=DSQRT(DABS(SQGRM2))
      CALL HSFIV2(XX,Y,XS,A2,TS)
c      write(*,*) ' xx = ',xx
c      write(*,*) ' y = ',y
c      write(*,*) ' xs = ',xs
c      write(*,*) ' a2 = ',a2
c      write(*,*) ' ts = ',ts
C
      A1=2D0*DKP
      RUNALP=1D0
      IF (LPAR(3).GE.3) THEN
        CG=HSSRGG(TS)
        RUNALP=1D0/(1D0+DREAL(CG)/TS)
      ENDIF
      R1=4D0*GX*(
     &      (T+2D0*MEI2)/(A2-TS)/(A2-TS)/A2
     &     -1D0/(A2-TS)/A2
     &     + 2D0/(A1+A2)/A2
     &     -8D0*MEI2*MEI2/(A1*A2+TS*TS)/(A1+A2)/A2
     &     +4D0*MEI2*MEI2/(A2*A2+TS*TS)/A2/A2
     &     -2D0*MEI2/(A2-TS)/A2/A2 )
      FAC1 =-(GS*GS + GU*GU - GX*(GS+GU) -4D0*MEI2*MPRO2)
      FAC2 = GU*GX - MPRO2*T
      FAC3 = (-2D0*GS*GU + GX*(GS+GU-GX))*2D0*MEI2
      FAC4 = GS*(GX-GS)*2D0*MEI2
      R2=4D0*(
     &     -FAC1/(A1+A2-TS)*(1D0/(A1+A2)+1D0/(A2-TS))/A2
     &     +FAC2/(A2-TS)/(A2-TS)/A2
     &     +FAC3/(A1*A2+TS*TS)/(A1+A2)/A2
     &     -2D0*MPRO2/(A1+A2)/A2
     &     +2D0*MEI2*MPRO2/(A2-TS)/(A2-TS)/A2
     &     +2D0*MEI2*MPRO2/(A2-TS)/A2/A2
     &     +MPRO2/(A2-TS)/A2
     &     +FAC4/(A2*A2+TS*TS)/A2/A2 )
      DO 20 IFL=1,12
        CQP(IFL)=0D0
 20   CONTINUE
      CALL HSFIE0(-TS,F1EL,F2EL)
c      write(*,*) ' f1el = ',f1el
c      write(*,*) ' f2el = ',f2el
c      write(*,*) ' r1 = ',r1
c      write(*,*) ' r2 = ',r2
      CQP(12)=(F1EL*R1+F2EL*R2)*RUNALP*RUNALP
      SUMME=CQP(12)
c      write(*,*) ' summe = ',summe
c      write(*,*) ' sqgram = ',sqgram
      HSELK2=SUMME*Y*2D0*SX1NRM/SQGRAM
     *      *(VMAX-VMIN)*A2*(TSMAXA-TSMINA)*TS*(-GS00-TS)
     *      *(GXMAX-GXMIN)*xx*xx*xx/2d0*(GMAXL-GMINL)*Y
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   INTEGRAND OF TS TERM (ELASTIC RADIATIVE TAIL)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSELCO(X)
C
C  X(1) -->  XX
C  X(2) -->  G=-1/Q**2   -->  Q**2=-1/G-->  Y
C  X(3) -->  LOG(-TS)
C  X(4) -->  2*K.P
C
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSTCUT/ THEMIN,THEMAX,CTHMIN,CTHCON
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSKPXY/ XX,Y
      COMMON /HSXSLM/ XSMIN,XSCUT
      COMMON /HSCMSP/ EQ,PQ,EEL,PEL,ES,PS,COSE,SINE,OMEGA
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSIRCX/ XIRDEL
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSPSPC/ IPHSPC
      DIMENSION X(4)
      COMPLEX(kind=8) HSSRGG,CG
C
C---QUASI-ELASTIC SCATTERING
      XS=1D0
C---CHOOSE VALUE OF Y
      GS=SP-MEI2-MPRO2
      YMAXX=(1D0-4D0*MEI2*MPRO2/GS/GS)/(1D0+2D0*MEI*MPRO/GS)
      IF (ICUT.LT.3) THEN
        GMIN=-1D0/(Q2MIN/XMAX/GS)
        GMAX=-1D0/DMIN1(1D0,YMAXX)
        ELSEIF(ICUT.EQ.3) THEN
C---CUT IN Y
        GMIN=-1D0/DMAX1(Q2MIN/XMAX/GS,YMIN)
        GMAX=-1D0/DMIN1(1D0,YMAXX,YMAX)
c...cut on lepton scattering angle (fixed target, me=0)
        elseif (icut.eq.4) then
        q2min=4d0*eele*eele*dsin(themin/180d0*pi/2d0)**2
     *       /(1d0+2d0*eele/mpro*dsin(themin/180d0*pi/2d0)**2)
        q2max=4d0*eele*eele*dsin(themax/180d0*pi/2d0)**2
     *       /(1d0+2d0*eele/mpro*dsin(themax/180d0*pi/2d0)**2)
        gmin=-1d0/q2min
        gmax=-1d0/q2max
        ELSE
        WRITE(LUNOUT,'(/A,I5/A)') ' WRONG VALUE OF ICUT:',ICUT,
     *                            ' STOP IN HSELCO'
        STOP
      ENDIF
      GACT=GMIN+X(1)*(GMAX-GMIN)
      Y=-1D0/GACT
      XA=1D0
      CALL HSDELX(XA,Y)
C---X-VALUE
      XXMAX1=HSXMAX(Y)
      XXHH1=1D0-Y-4D0*MEI2*MPRO2/GS/GS
      XXMNY1=(XXHH1+DSQRT(XXHH1*XXHH1-4D0*Y*Y*MEI2*MPRO2/GS/GS))
     &      /2D0/Y/MPRO2*GS
      XXMINY=MEI2/MPRO2/XXMNY1
      XXMIN=DMAX1(XMIN,Q2MIN/Y/GS,XXMINY)
      XXMAX=DMIN1(XMAX,XXMAX1)
      XX=XXMIN+(XXMAX-XXMIN)*X(2)
      Q2=XX*Y*GS
      CALL HSFIVC(XX,Y)
C
      IF(IPRINT.GT.30) THEN
        WRITE(LUNTES,230)SP,XX,Y,XSMIN,XSCUT,DELEPS,DELTA
230     FORMAT(' ***************************************************',/
     F        ,' SP = ',1PD12.3,/
     F        ,' X = ',D12.6,'   Y = ',D12.6,/
     F        ,' XSMIN = ',D17.11,'   XSCUT = ',D17.11,/
     F        ,' DELEPS = ',D12.6,'   DELTA = ',D14.8,/
     F       ,' ***************************************************',//)
      ENDIF
C
      XS=1D0
      CALL HSFCMS(XX,Y,XS)
      IF(IPHSPC.EQ.1) THEN
        HSELCO=0D0
        RETURN
      ENDIF
      CALL HSFLAB(XS)
      CALL HSLZTS(ZMIN,ZMAX)
      IF(IPHSPC.EQ.1) THEN
        HSELCO=0D0
        RETURN
      ENDIF
      TSMIN=TP+2D0*OMEGA*(ES-EEL-PQ*ZMAX)
      IF(ZMIN.LT.-0.9999D0)THEN
        TSMAX=-XX*XX*XS*MPRO2/(XS-XX+XS*XS*MPRO2/Y/GS)
      ELSE
        TSMAX=TP+2D0*OMEGA*(ES-EEL-PQ*ZMIN)
      ENDIF
C---SUBSTITUTION R = LN(-TS)
      RMAX=DLOG(-TSMIN)
      RMIN=DLOG(-TSMAX)
      if (rmax.le.rmin) then
        hselco=0d0
        return
      endif
      R=RMIN+(RMAX-RMIN)*X(3)
      TS=-DEXP(R)
C---NO SUBSTITUTION FOR A1
      CALL HSL1TS(TS,XX,Y,XS,A1MIN,A1MAX,A1M,A1P,CFKP)
      IF(IPHSPC.EQ.1)THEN
       HSELCO=0D0
       RETURN
      ENDIF
      if ((a1max.le.a1min).or.(a1min.le.0d0).or.(a1max.le.0d0)) then
        hselco=0d0
        return
      endif
      A1=A1MIN+(A1MAX-A1MIN)*X(4)
      SQGRM2=-CFKP*(A1-A1M)*(A1-A1P)
      SQGRAM=DSQRT(DABS(SQGRM2))
      CALL HSFIV1(XX,Y,XS,A1,TS)
C
      A2 = 2D0*DKPS
      RUNALP=1D0
      IF (LPAR(3).GE.3) THEN
        CG=HSSRGG(TS)
        RUNALP=1D0/(1D0+DREAL(CG)/TS)
      ENDIF
      R1=4D0*GX*(
     &     +(T+6D0*MEI2)/(A1-TS)/(A1-TS)/TS
     &     -(T+4D0*MEI2)/(A1-TS)/TS/TS
     &     -(T+2D0*MEI2)/(A2-TS)/(A2-TS)/TS
     &     -(T+4D0*MEI2)/(A2-TS)/TS/TS
     &     +2D0/TS/TS
     &     -1D0/(A1-TS)/TS
     &     +1D0/(A2-TS)/TS
     &     -8D0*MEI2*MEI2/(A1*A2+TS*TS)/TS/TS
     &     +4D0*MEI2*MEI2/(A1*A1+TS*TS)/TS/TS
     &     +4D0*MEI2*MEI2/(A2*A2+TS*TS)/TS/TS )
      FAC1 =-(GS*GS + GU*GU - GX*(GS+GU) -4D0*MEI2*MPRO2)
      FAC2 = -GS*GX + MPRO2*T
      FAC3 = GU*GX - MPRO2*T
      FAC4 = (-2D0*GS*GU + GX*(GS+GU-GX))*2D0*MEI2
      FAC5 = GU*(GX-GU)*2D0*MEI2
      FAC6 = GS*(GX-GS)*2D0*MEI2
      R2=4D0*(
     &      FAC1/(A1+A2-TS)*(1D0/(A1-TS)+1D0/(A2-TS))/TS
     &     -FAC2/(A1-TS)/(A1-TS)/TS
     &     +FAC2/(A1-TS)/TS/TS
     &     -FAC3/(A2-TS)/(A2-TS)/TS
     &     +FAC3/(A2-TS)/TS/TS
     &     +FAC4/(A1*A2+TS*TS)/TS/TS
     &     -2D0*MEI2*MPRO2/(A1-TS)/(A1-TS)/TS
     &     -2D0*MEI2*MPRO2/(A2-TS)/(A2-TS)/TS
     &     -2D0*MPRO2/TS/TS
     &     +MPRO2/(A1-TS)/TS
     &     -MPRO2/(A2-TS)/TS
     &     +FAC5/(A1*A1+TS*TS)/TS/TS
     &     +FAC6/(A2*A2+TS*TS)/TS/TS)
      DO 20 IFL=1,12
        CQP(IFL)=0D0
 20   CONTINUE
      CALL HSFIE0(-TS,F1EL,F2EL)
      CQP(12)=(F1EL*R1+F2EL*R2)*RUNALP*RUNALP
      SUMME=CQP(12)
      HSELCO=SUMME*Y*2D0*SX1NRM/SQGRAM
     *      *(A1MAX-A1MIN)*(RMAX-RMIN)*(-TS)
     *      *(XXMAX-XXMIN)*(GMAX-GMIN)/GACT**2
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   ELASTIC FORM FACTORS FOR THE PROTON
C   TAKEN FROM STEIN ET AL. (PRD 12 (1975) 1884)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFIE0(Q2,F1EL,F2EL)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSELEP/ IDIPOL
      DIMENSION HDIP(0:5)
      DATA AK /2.7927D0/
      DATA Q20 /0.71D0/
      DATA HDIP /1.0007D0, 1.01807D0, 1.05584D0,
     *           0.836380D0, 0.6864584D0, 0.672830D0/
C...for Larin's parametrization
      dimension ffae2(8),ffam2(8)
      data ffae2 /-0.48818d0, 5.20168d0, -32.22855d0, 112.32519d0, 
     *           -231.64144d0, 283.78979d0, -193.25343d0, 56.33867d0/
      data ffam2 /0.27814d0, -5.6776d0, 36.5010d0, -112.99d0, 
     *           196.640d0, -197.640d0, 108.090d0, -24.988d0/
C...for the Friedrich-Walcher parametrization
      dimension ffae3(7),ffam3(7)
      data ffae3 /1.041d0, 0.765d0, -0.041d0, 6.2d0, -0.23d0, 
     *            0.07d0, 0.27d0/
      data ffam3 /1.002d0, 0.749d0, -0.002d0, 6.0d0, -0.13d0, 
     *            0.35d0, 0.21d0/

ctest
c test for dvcs with cuts on t_1
c      logical lfirst
c      data lfirst /.true./
c      if (lfirst) then 
c        write(*,*) ' +++++++++++++++++++++++++++++++++++++++++++ '
c        write(*,*) ' +++++++++++++++++++++++++++++++++++++++++++ '
c        write(*,*) ' +++++++++++++++++++++++++++++++++++++++++++ '
c        write(*,*) ' +++++++++++++++++++++++++++++++++++++++++++ '
c        write(*,*) ' '
c        write(*,*) ' TEST FOR DVCS: CUTS on T1 '
c        write(*,*) ' '
c        write(*,*) ' +++++++++++++++++++++++++++++++++++++++++++ '
c        write(*,*) ' +++++++++++++++++++++++++++++++++++++++++++ '
c        write(*,*) ' +++++++++++++++++++++++++++++++++++++++++++ '
c        write(*,*) ' +++++++++++++++++++++++++++++++++++++++++++ '
c        lfirst=.false.
c      endif
c      if (Q2.le.0.01d0.or.Q2.ge.1.0d0) then
c        f1el=0d0
c        f2el=0d0
c        return
c      endif
ctest
      TAU=Q2/4D0/MPRO2
      GEP=1D0/(1D0+Q2/Q20)**2
      GEM=GEP*AK
      PDEVE=1D0
      PDEVM=1D0
      IF (IDIPOL.EQ.1) THEN
C...DEVIATION FROM DIPOLE FORM FACTOR
        Q=DSQRT(Q2)
        PDEV=0D0
        DO 1 I=0,5
        PDI=1D0
        DO 2 J=0,5
        IF (J.EQ.I) GOTO 2
        PDI=PDI*(Q-DFLOAT(J))/DFLOAT(I-J)
    2   CONTINUE
        PDEV=PDEV+PDI*HDIP(I)
    1   CONTINUE
        PDEVE=PDEV
        PDEVM=PDEV
        GEP=GEP*PDEVE
        GEM=GEP*PDEVM*AK
C...Larin's parametrization
       ELSEIF (IDIPOL.EQ.2) THEN
        pdeve=1d0
        pdevm=1d0
        do 3 i=1,8
        pdeve=pdeve+ffae2(i)*q2**dfloat(i)
        pdevm=pdevm+ffam2(i)*q2**dfloat(i)
    3   continue
        GEP=GEP*PDEVE
        GEM=GEP*PDEVM*AK
       ELSEIF (IDIPOL.EQ.3) THEN
C...Friedrich-Walcher parametrization
        gep=ffae3(1)/(1d0+q2/ffae3(2))**2
     .     +ffae3(3)/(1d0+q2/ffae3(4))**2
     .     +ffae3(5)*q2*(exp(-0.5d0*((dsqrt(q2)-ffae3(6))/ffae3(7))**2)
     .                  +exp(-0.5d0*((dsqrt(q2)+ffae3(6))/ffae3(7))**2))
        gem=ffam3(1)/(1d0+q2/ffam3(2))**2
     .     +ffam3(3)/(1d0+q2/ffam3(4))**2
     .     +ffam3(5)*q2*(exp(-0.5d0*((dsqrt(q2)-ffam3(6))/ffam3(7))**2)
     .                  +exp(-0.5d0*((dsqrt(q2)+ffam3(6))/ffam3(7))**2))
        gem=gem*ak
      ENDIF
C...STRUCTURE FUNCTIONS W1EL, W2EL -> F1,F2
      W1EL=TAU*GEM*GEM
      W2EL=(GEP*GEP+TAU*GEM*GEM)/(1D0+TAU)
      F1EL=W1EL/(2D0*TAU)
      F2EL=W2EL
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   ELASTIC FORM FACTORS FOR THE PROTON
C   TAKEN FROM STEIN ET AL. (PRD 12 (1975) 1884)
C   INCLUDING VIRTUAL&SOFT CORRECTIONS AND VACUUM POLARIZATION
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSFIEL(Q2,F1EL,F2EL)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMPLEX(kind=8) HSSRGG,CG
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSELEP/ IDIPOL
      DIMENSION HDIP(0:5)
      DATA AK /2.7927D0/
      DATA Q20 /0.71D0/
      DATA HDIP /1.0007D0, 1.01807D0, 1.05584D0,
     *           0.836380D0, 0.6864584D0, 0.672830D0/
C...for Larin's parametrization
      dimension ffae2(8),ffam2(8)
      data ffae2 /-0.48818d0, 5.20168d0, -32.22855d0, 112.32519d0, 
     *           -231.64144d0, 283.78979d0, -193.25343d0, 56.33867d0/
      data ffam2 /0.27814d0, -5.6776d0, 36.5010d0, -112.99d0, 
     *           196.640d0, -197.640d0, 108.090d0, -24.988d0/
C...for the Friedrich-Walcher parametrization
      dimension ffae3(7),ffam3(7)
      data ffae3 /1.041d0, 0.765d0, -0.041d0, 6.2d0, -0.23d0, 
     *            0.07d0, 0.27d0/
      data ffam3 /1.002d0, 0.749d0, -0.002d0, 6.0d0, -0.13d0, 
     *            0.35d0, 0.21d0/

      TAU=Q2/4D0/MPRO2
      GEP=1D0/(1D0+Q2/Q20)**2
      GEM=GEP*AK
      PDEVE=1D0
      PDEVM=1D0
      IF (IDIPOL.EQ.1) THEN
C...DEVIATION FROM DIPOLE FORM FACTOR
        Q=DSQRT(Q2)
        PDEV=0D0
        DO 1 I=0,5
        PDI=1D0
        DO 2 J=0,5
        IF (J.EQ.I) GOTO 2
        PDI=PDI*(Q-DFLOAT(J))/DFLOAT(I-J)
    2   CONTINUE
        PDEV=PDEV+PDI*HDIP(I)
    1   CONTINUE
        PDEVE=PDEV
        PDEVM=PDEV
        GEP=GEP*PDEVE
        GEM=GEP*PDEVM*AK
C...Larin's parametrization
       ELSEIF (IDIPOL.EQ.2) THEN
        pdeve=1d0
        pdevm=1d0
        do 3 i=1,8
        pdeve=pdeve+ffae2(i)*q2**dfloat(i)
        pdevm=pdevm+ffam2(i)*q2**dfloat(i)
    3   continue
        GEP=GEP*PDEVE
        GEM=GEP*PDEVM*AK
C...Friedrich-Walcher parametrization
       ELSEIF (IDIPOL.EQ.3) THEN
        gep=ffae3(1)/(1d0+q2/ffae3(2))**2
     .     +ffae3(3)/(1d0+q2/ffae3(4))**2
     .     +ffae3(5)*q2*(exp(-0.5d0*((dsqrt(q2)-ffae3(6))/ffae3(7))**2)
     .                  +exp(-0.5d0*((dsqrt(q2)+ffae3(6))/ffae3(7))**2))
        gem=ffam3(1)/(1d0+q2/ffam3(2))**2
     .     +ffam3(3)/(1d0+q2/ffam3(4))**2
     .     +ffam3(5)*q2*(exp(-0.5d0*((dsqrt(q2)-ffam3(6))/ffam3(7))**2)
     .                  +exp(-0.5d0*((dsqrt(q2)+ffam3(6))/ffam3(7))**2))
        gem=gem*ak
C...Test: F1=0, F2=1
       ELSEIF (IDIPOL.EQ.4) THEN
        gem=0d0
        gep=1d0
      ENDIF
C...STRUCTURE FUNCTIONS W1EL, W2EL -> F1,F2
      W1EL=TAU*GEM*GEM
      W2EL=(GEP*GEP+TAU*GEM*GEM)/(1D0+TAU)
      F1EL=W1EL/(2D0*TAU)
      F2EL=W2EL
C...Test
      if (idipol.eq.4) then
        f1el=0d0
        f2el=1d0
      endif
      
C...INCLUDE VIRTUAL AND SOFT PHOTON CORRECTIONS
      IF (LPAR(2).GT.0) THEN
      SHAT=SP-MEI2-MPRO2
      X=1D0
      Y=Q2/SHAT
C..PHOTON SELF ENERGY
      IF (LPAR(7).GE.1) THEN
        T=-Q2
        CG=HSSRGG(T)
        PIGGG=DREAL(CG)/T
        ELSE
        PIGGG=0D0
      ENDIF
      DVACGG=1D0/(1D0+PIGGG)
C..LEPTONIC QED VERTEX CORRECTIONS INCLUDING SOFT BREMSSTRAHLUNG
      IF (LPAR(12).GE.1) THEN
        DVRTXV=HSDQDV(X,Q2)*ALP2PI
        DVRTXS=HSDQDS(Q2)*ALP2PI
        DVRTX1=DVRTXV-(Q2+4D0*MEI2)/(Q2-2D0*MEI2)*DVRTXS
        DVRTX2=DVRTXV+
     *     (2D0*(1D0-Y)+Y*Y/2D0)/(1D0-Y-MPRO2*Q2/SHAT/SHAT)*DVRTXS
        ELSE
        DVRTX1=0D0
        DVRTX2=0D0
      ENDIF
      IF (LPAR(3).LT.3) THEN
        F1EL=F1EL*(DVRTX1+DVACGG*DVACGG)
        F2EL=F2EL*(DVRTX2+DVACGG*DVACGG)
        ELSE
        F1EL=F1EL*(1D0+DVRTX1)*DVACGG*DVACGG
        F2EL=F2EL*(1D0+DVRTX2)*DVACGG*DVACGG
      ENDIF
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   CROSS SECTION FOR ELASTIC LEPTON PROTON SCATTERING:
CCCCCCC 01.12.95    CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC H. SPIESBERGER CCC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSSGEL(Q2,LL,POL)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSPARL/ LPAR(20),LPARIN(12)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSCUMS/ CQP(12)
      COMMON /HSKNST/ PI,ALPHA,ALP1PI,ALP2PI,ALP4PI,E,GF,SXNORM,SX1NRM
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      LOGICAL LFIRST
      DATA LFIRST /.TRUE./

      IF (LFIRST) THEN
         IF (POL.NE.0D0) THEN
            WRITE(LUNOUT,'(5X,A)') ' ****** Warning: '
            WRITE(LUNOUT,'(5X,A)') ' ****** POL non-zero in HSSGEL '
            WRITE(LUNOUT,'(5X,A)') ' ****** is not implemented '
            WRITE(LUNOUT,'(5X,A)') ' ****** Execution continued '
         ENDIF
         IF (LL.NE.-1) THEN
            WRITE(LUNOUT,'(5X,A)') ' ****** Warning: '
            WRITE(LUNOUT,'(5X,A)') ' ****** LL not -1 in HSSGEL '
            WRITE(LUNOUT,'(5X,A)') ' ****** is not implemented '
            WRITE(LUNOUT,'(5X,A)') ' ****** Execution continued '
         ENDIF
         LFIRST=.FALSE.
      ENDIF

      SHAT=SP-MEI2-MPRO2
      X=1D0
      Y=Q2/SHAT
      T=-Q2
      SPN=SXNORM*SP
      CALL HSFIVC(X,Y)
      CALL HSDELX(X,Y)
      CALL HSFIEL(-T,F1EL,F2EL)
      DO 10 I=1,12
        CQP(I)=0D0
 10   CONTINUE
      HSSGEL=0D0
      R1=Y*Y*(1D0+2D0*MEI2/T)
      R2=(1D0-Y+MPRO2*T/SHAT/SHAT)/X
      CQP(12)=F1EL*R1+F2EL*R2
      HSSGEL=8D0*SPN*CQP(12)/T/T
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE HSDELX(XA,Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSIRCX/ XIRDEL
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      IF (IOPEGM.GT.0) THEN
        DELTA=EGMIN
        RETURN
      ELSE
        X=XA-XIRDEL
        SIGMA=EPRO/EELE
        XMY=Y-(1D0-X*Y)*SIGMA
        XPY=Y+(1D0-X*Y)*SIGMA
        OMEGA=2D0*EELE*SIGMA*Y*(1D0-X)
     *        /(XPY+DSQRT(4D0*X*Y*SIGMA*(1D0-Y)+XMY*XMY))
        EQUA=X*SIGMA*EELE
        EES=EELE*(1D0-Y+X*Y*SIGMA)
C       DELTA=DMIN1(EELE,EES,OMEGA)*DELEPS
C       DELTA=DMIN1(EELE,EES,OMEGA)
        DELTA=OMEGA
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   XSMIN FROM OMEGA-MAX(XS)=DELTA
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSXMAX(Y)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSGIKP/ GS,GU,GX,TP,UP
      COMMON /HSIRCT/ DELEPS,DELTA,EGMIN,IOPEGM
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
C
      XS=1D0
      X1=XMIN
      X2=1D0
      DO 90 N=1,70,1
        X3=(X1+X2)/2D0
        CALL HSFIVC(X3,Y)
        OM3=HSOMAX(X3,Y,XS)
        IF (OM3.LT.DELTA) THEN
          X2=X3
        ELSE
          X1=X3
        ENDIF
90    CONTINUE
      HSXMAX=X3

C...CHECK ON ES > ME
C     IF (X3.LT.1D-6) THEN
C     BXS=GU*(TP+MEI2)/(GU*GU-4D0*MEI2*MPRO2)
C     DXS=1D0-(GU*GU-4D0*MEI2*MPRO2)/GU/GU
C    *        *(TP*(TP-2D0*MEI2)-7D0*MEI2*MEI2)/(TP+MEI2)/(TP+MEI2)
C     IF (DXS.LT.0D0) RETURN
C     XEMIN=-BXS*(1D0+DSQRT(DXS))
C     XEMAX=-BXS*(1D0-DSQRT(DXS))
C     HSXMAX=DMIN1(HSXMAX,XEMIN,XEMAX)
C     ENDIF

      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   NUCLEAR SHADOWING
C   ASSUME Q2-INDEPENDENT SHADOWING (ANTI-SHADOWING) FOR HEAVY NUCLEI
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION HSNRAT(X)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      LOGICAL LFIRST
      DATA LFIRST/.TRUE./
      DATA AN1/0.130D0/
      DATA AN2/0.456D0/
      DATA AN3/0.773D0/

      HSNRAT=1D0
      IF (INUMOD.NE.1000) RETURN
      
      IF (LFIRST) THEN
        LFIRST=.FALSE.
        HNA3=HNA**(1D0/3D0)
        HMI=1D0-1D0/HNA3-1.145D0/HNA3/HNA3+0.93D0/HNA+0.88D0/HNA/HNA3
     *         -0.59D0/HNA/HNA3/HNA3
        HM1=HMI*AN1
        H1M2=1D0+HMI*AN2
        HM3=HMI*AN3
      ENDIF
      IF (INT(HNA).EQ.2.AND.INT(HNZ).EQ.1) THEN
        HSNRAT=1D0
       ELSE
        HSNRAT=X**HM1*(1D0-HM3*X)*H1M2
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     INTEGRATION SUBROUTINE
C     SUBSTITUTE FOR THE NAGLIB ROUTINE D01FCF
C     Note: The function to be integrated is identified by 
C     INTEGER IFUN = 1: Neutral Current
C                  = 2: Charged Current
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE DX1FCF(BLOW,BUP,IFUN,ACCREQ,ACCFIN,RESULT,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BLOW(2),BUP(2),BLOWC(2),BUPC(2)
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSD01L/ BLOWC,BUPC,ARG1,ACC,INCCC
      EXTERNAL DFNCII
      DO 1 I=1,2
        BLOWC(I)=BLOW(I)
        BUPC(I)=BUP(I)
 1    CONTINUE
      ACC=ACCREQ
      INCCC=IFUN
      RESULT=GAUSK1(DFNCII,BLOWC(2),BUPC(2),ACCREQ)
      ACCFIN=ACCREQ
      IFAIL=0
      RETURN
      END
C
      FUNCTION DFNCII(B1)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      DIMENSION BLOWC(2),BUPC(2)
      EXTERNAL DFNC00
      COMMON /HSD01L/ BLOWC,BUPC,ARG1,ACCREQ,INCCC
      ARG1=B1
      DFNCII=GAUSK2(DFNC00,BLOWC(1),BUPC(1),ACCREQ)
      END
C
      FUNCTION DFNC00(B2)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      DIMENSION ARG(2),BLOWC(2),BUPC(2)
      COMMON /HSD01L/ BLOWC,BUPC,ARG1,ACCREQ,INCCC
      DFNC00=0D0
      NDIM=2
      ARG(1)=ARG1
      ARG(2)=B2
      IF (INCCC.EQ.1) THEN
        DFNC00=HSNCG1(NDIM,ARG)
      ELSEIF (INCCC.EQ.2) THEN
        DFNC00=HSCCG1(NDIM,ARG)
      ENDIF
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     INTEGRATION SUBROUTINE
C     SUBSTITUTE FOR THE NAGLIB ROUTINE D01AJF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE D01AJF(F,BLOW,BUP,ACCREL,RESULT,ACCFIN,IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      EXTERNAL HSELG1
      EXTERNAL F
      BLOWC=BLOW
      BUPC=BUP
      ACCREQ=ACCREL
C      RESULT=GAUSK1(HSELG1,BLOWC,BUPC,ACCREQ)
      RESULT=GAUSK1(F,BLOWC,BUPC,ACCREQ)
      ACCFIN=ACCREL*RESULT
      IFAIL=0
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.
C     GAUSK1(DGAUSS) IS SET EQUAL TO THE APPROXIMATE VALUE OF THE
C     INTEGRAL OF THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY
C     PARAMETER EPS.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION GAUSK1(F,A,B,EPS)
      DOUBLE PRECISION GAUSK1,F,A,B,EPS
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST
      EXTERNAL F
      DATA W / 0.10122 85362 90376 259D0,
     1         0.22238 10344 53374 471D0,
     2         0.31370 66458 77887 287D0,
     3         0.36268 37833 78361 983D0,
     4         0.27152 45941 17540 949D-1,
     5         0.62253 52393 86478 929D-1,
     6         0.95158 51168 24927 848D-1,
     7         0.12462 89712 55533 872D0,
     8         0.14959 59888 16576 732D0,
     9         0.16915 65193 95002 538D0,
     A         0.18260 34150 44923 589D0,
     B         0.18945 06104 55068 496D0/
      DATA X / 0.96028 98564 97536 232D0,
     1         0.79666 64774 13626 740D0,
     2         0.52553 24099 16328 986D0,
     3         0.18343 46424 95649 805D0,
     4         0.98940 09349 91649 933D0,
     5         0.94457 50230 73232 576D0,
     6         0.86563 12023 87831 744D0,
     7         0.75540 44083 55003 034D0,
     8         0.61787 62444 02643 748D0,
     9         0.45801 67776 57227 386D0,
     A         0.28160 35507 79258 913D0,
     B         0.95012 50983 76374 402D-1/
C******************************************************************
C..START.
      GAUSK1=0.0D0
      IF(B.EQ.A) RETURN
      CONST=0.005D0/(B-A)
      BB=A
C..COMPUTATIONAL LOOP.
    1 AA=BB
      BB=B
    2    C1=0.5D0*(BB+AA)
         C2=0.5D0*(BB-AA)
         S8=0.0D0
         DO 3 I=1,4
            U=C2*X(I)
            S8=S8+W(I)*(F(C1+U)+F(C1-U))
    3    CONTINUE
         S8=C2*S8
         S16=0.0D0
         DO 4 I=5,12
            U=C2*X(I)
            S16=S16+W(I)*(F(C1+U)+F(C1-U))
    4    CONTINUE
         S16=C2*S16
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5
         BB=C1
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2
      GAUSK1=0.0D0
      RETURN
    5 GAUSK1=GAUSK1+S16
      IF(BB.NE.B) GO TO 1
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ADAPTIVE DOUBLE PRECISION GAUSSIAN QUADRATURE.
C     GAUSK2(DGAUSS) IS SET EQUAL TO THE APPROXIMATE VALUE OF THE
C     INTEGRAL OF THE FUNCTION F OVER THE INTERVAL (A,B), WITH ACCURACY
C     PARAMETER EPS.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      FUNCTION GAUSK2(F,A,B,EPS)
      DOUBLE PRECISION GAUSK2,F,A,B,EPS
      DOUBLE PRECISION W(12),X(12),AA,BB,C1,C2,U,S8,S16,CONST
      EXTERNAL F
      DATA W / 0.10122 85362 90376 259D0,
     1         0.22238 10344 53374 471D0,
     2         0.31370 66458 77887 287D0,
     3         0.36268 37833 78361 983D0,
     4         0.27152 45941 17540 949D-1,
     5         0.62253 52393 86478 929D-1,
     6         0.95158 51168 24927 848D-1,
     7         0.12462 89712 55533 872D0,
     8         0.14959 59888 16576 732D0,
     9         0.16915 65193 95002 538D0,
     A         0.18260 34150 44923 589D0,
     B         0.18945 06104 55068 496D0/
      DATA X / 0.96028 98564 97536 232D0,
     1         0.79666 64774 13626 740D0,
     2         0.52553 24099 16328 986D0,
     3         0.18343 46424 95649 805D0,
     4         0.98940 09349 91649 933D0,
     5         0.94457 50230 73232 576D0,
     6         0.86563 12023 87831 744D0,
     7         0.75540 44083 55003 034D0,
     8         0.61787 62444 02643 748D0,
     9         0.45801 67776 57227 386D0,
     A         0.28160 35507 79258 913D0,
     B         0.95012 50983 76374 402D-1/
C******************************************************************
C..START.
      GAUSK2=0.0D0
      IF(B.EQ.A) RETURN
      CONST=0.005D0/(B-A)
      BB=A
C..COMPUTATIONAL LOOP.
    1 AA=BB
      BB=B
    2    C1=0.5D0*(BB+AA)
         C2=0.5D0*(BB-AA)
         S8=0.0D0
         DO 3 I=1,4
            U=C2*X(I)
            S8=S8+W(I)*(F(C1+U)+F(C1-U))
    3    CONTINUE
         S8=C2*S8
         S16=0.0D0
         DO 4 I=5,12
            U=C2*X(I)
            S16=S16+W(I)*(F(C1+U)+F(C1-U))
    4    CONTINUE
         S16=C2*S16
         IF( ABS(S16-S8) .LE. EPS*(1.+ABS(S16)) ) GO TO 5
         BB=C1
         IF( 1.D0+ABS(CONST*C2) .NE. 1.D0) GO TO 2
      GAUSK2=0.0D0
      RETURN
    5 GAUSK2=GAUSK2+S16
      IF(BB.NE.B) GO TO 1
      RETURN
      END
