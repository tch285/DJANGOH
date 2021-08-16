CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C
C   This is a collection of program codes                                    C
C   for the evaluation of polarized parton distribution functions            C
C   to be used in the package DJANGOH.4.6.10                                 C
C   (by T. Martini, H. Spiesberger, 7.8.2013)                                C
C                                                                            C
C   Origin: http://ribf.riken.jp/~marco/DSSV/                                C
C           http://hepdata.cedar.ac.uk/pdfs                                  C
C   'DSSV.f'                                                                 C
C   'ds2000.f'                                                               C
C   'aac08.f'                                                                C
C   'grsv2000pdf_g1.f'                                                       C
C   'LSS2010_PDFs.f'                                                         C
C   'polnlo.f'                                                               C
C   'polpar.f'                                                               C
C   'ppdf.f'                                                                 C
C   'readgrid.f'                                                             C
C                                                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C=============================================================================
C=============================================================================


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C
C       ------ DSSV POLARIZED PARTON DISTRIBUTIONS ------                    C
C                                                                            C
C       REFERENCE:                                                           C
C           "GLOBAL ANALYSIS OF HELICITY PARTON DENSITIES                    C
C            AND THEIR UNCERTAINTIES",                                       C
C            arXiv:0804:0422 [hep-ph]                                        C
C                                                                            C
C       INPUT:                                                               C
C           X     :  BJORKEN-X BETWEEN 10**(-5)< X < 1.0                     C
C           Q2    :  SCALE**2 IN GeV**2 BETWEEN  1 < Q**2 < 10**5 GeV**2     C
C                                                                            C
C       OUTPUT:                                                              C
C           DUV   :  X * U VALENCE DISTRIBUTION                              C
C           DDV   :  X * D VALENCE DISTRIBUTION                              C
C           DUBAR :  X * UBAR DISTRIBUTION                                   C
C           DDBAR :  X * DBAR DISTRIBUTION                                   C
C           DSTR  :  X * STRANGE DISTRIBUTION                                C
C           DGLU  :  X * GLUON DISTRIBUTION                                  C
C                                                                            C
C       IMPORTANT:                                                           C
C           ALWAYS X*DISTRIBUTION IS RETURNED !!!                            C
C           ALL PDFs ARE IN THE MSbar SCHEME                                 C
C                                                                            C
C           BEFORE CALLING THE SUBROUTINE 'DSSVFIT' FOR THE FIRST TIME,      C
C           THE SUBROUTINE 'DSSVINI' MUST BE CALLED (ONLY ONCE) TO           C
C           INITIALIZE THE GRIDS !!                                          C
C                                                                            C
C       IN CASE OF PROBLEMS, DOUBTS, ETC, PLEASE E-MAIL US:                  C
C           D. de Florian  deflo@df.uba.ar                                   C
C           R. Sassot      sassot@df.uba.ar                                  C
C           M. Stratmann   marco@ribf.riken.jp                               C
C           W. Vogelsang   vogelsan@quark.phy.bnl.gov                        C
C                                                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C---------------------------------------------------------------------
      SUBROUTINE DSSVFIT(Xin,Q2in,DUV,DDV,DUBAR,DDBAR,DSTR,DGLU)
C---------------------------------------------------------------------
      IMPLICIT NONE
C...
      INTEGER NPART, NX, NQ, NARG
      PARAMETER (NPART=6, NX=47, NQ=30, NARG=2)
C...
      INTEGER NA(NARG)
      DOUBLE PRECISION XUF(NX,NQ), XDF(NX,NQ), XUBF(NX,NQ), XDBF(NX,NQ), 
     1          XSF(NX,NQ), XGF(NX,NQ),
     2          XT(NARG), 
     3          ARRF(100) 
chs     3          ARRF(NX+NQ) 
      DOUBLE PRECISION X, Q2
      DOUBLE PRECISION DUV, DDV, DUBAR, DDBAR, DSTR, DGLU
      DOUBLE PRECISION DFINT
C...      
      COMMON/ DSSVGRID / XUF, XDF, XUBF, XDBF, XSF, XGF, NA, ARRF
C

C...Added for interface to DJANGOH (TM&HS, 10.01.2013)
      DOUBLE PRECISION Xin, Q2in
      integer ixoutmn,iqoutmn,ixoutmx,iqoutmx
      DOUBLE PRECISION XMIND, XMAXD, Q2MIND, Q2MAXD
      COMMON /HSDPEX/ ixoutmn, ixoutmx, iqoutmn, iqoutmx,
     *                XMIND, XMAXD, Q2MIND, Q2MAXD
      logical FIRST
      data FIRST /.true./

C...  CHECK IF X AND Q2 VALUES ARE WITHIN RANGE OF THE GRID: 
C
      IF (FIRST) THEN
        ixoutmn=0
        ixoutmx=0
        iqoutmn=0
        iqoutmx=0
        FIRST=.FALSE.
      ENDIF
      X=Xin
      Q2=Q2in
      XMIND=1.0D-5
      XMAXD=1.0D0
      Q2MIND=0.8D0
      Q2MAXD=1.D6
      if (Xin.lt.XMIND) then 
        X=XMIND
        ixoutmn=ixoutmn+1
      endif
      if (Xin.gt.XMAXD) then
        X=XMAXD
        ixoutmx=ixoutmx+1
       endif
       if (Q2in.lt.Q2MIND) then 
         Q2=Q2MIND
         iqoutmn=iqoutmn+1
       endif
       if (Q2in.gt.Q2MAXD) then
         Q2=Q2MAXD
         iqoutmx=iqoutmx+1
       endif
C
C...  INTERPOLATION AND OUTPUT:
C
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      DUV   = DFINT(NARG,XT,NA,ARRF,XUF)  * (1.D0-X)**3 * X
      DDV   = DFINT(NARG,XT,NA,ARRF,XDF)  * (1.D0-X)**4 * X
      DUBAR = DFINT(NARG,XT,NA,ARRF,XUBF) * (1.D0-X)**8 * X**0.5
      DDBAR = DFINT(NARG,XT,NA,ARRF,XDBF) * (1.D0-X)**8 * X**0.5
      DSTR  = DFINT(NARG,XT,NA,ARRF,XSF)  * (1.D0-X)**8 * X**0.5
      DGLU  = DFINT(NARG,XT,NA,ARRF,XGF)  * (1.D0-X)**5 * X**0.5
C...
chs 60   RETURN
      return
      END
C
C---------------------------
      SUBROUTINE DSSVINI
C---------------------------
      IMPLICIT NONE
C...
      INTEGER NPART, NX, NQ, NARG
      PARAMETER (NPART=6, NX=47, NQ=30, NARG=2)
C...
      INTEGER NA(NARG)
      DOUBLE PRECISION PARTON (NPART,NQ,NX-1)
      DOUBLE PRECISION QS(NQ), XB(NX)
      DOUBLE PRECISION XUF(NX,NQ), XDF(NX,NQ), XUBF(NX,NQ), XDBF(NX,NQ), 
     1          XSF(NX,NQ), XGF(NX,NQ), ARRF(100) 
chs     1          XSF(NX,NQ), XGF(NX,NQ), ARRF(NX+NQ) 
      DOUBLE PRECISION XB0, XB1
      INTEGER M, N
      INTEGER IQ, IX
      integer lunpl1,NextUn
C...
      COMMON/ DSSVGRID / XUF, XDF, XUBF, XDBF, XSF, XGF, NA, ARRF
C
C...  BJORKEN-X AND Q**2 VALUES OF THE GRID :
C
      DATA QS / 0.8D0, 1.0D0, 1.25d0, 1.5D0, 2.d0, 2.5D0, 
     1     4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2     1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3     3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4, 
     4     1.0D5, 1.8D5, 3.2D5, 5.8D5, 1.0D6  /
      DATA XB / 
     1           1.0D-5, 4.D-5, 6.7D-5, 1.0D-4, 1.4D-4, 2.0D-4,
     2           3.0D-4, 4.5D-4, 6.7D-4, 1.0D-3, 1.4D-3, 2.0D-3,
     3           3.0D-3, 4.5D-3, 6.7D-3, 1.0D-2, 1.4D-2, 2.0D-2,
     4           3.0D-2, 4.5D-2, 0.06, 0.08, 0.1, 0.125, 0.15,
     5           0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325,
     6           0.35, 0.375, 0.4,  0.45, 0.5, 0.55, 0.6,
     7           0.65,  0.7,  0.75,  0.8,  0.85, 0.9, 0.95, 1.0/ 

C...
      lunpl1=NextUn()
      OPEN(lunpl1,FILE='polpdf-gridfiles/DSSV-GRID.NLO',STATUS='OLD')
C...
      DO 151 M = 1, NX-1
         DO 152 N = 1, NQ
            READ(lunpl1,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1           PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M)
 90         FORMAT (6(1PE12.4))
 152     CONTINUE
 151  CONTINUE
C...
      CLOSE(lunpl1)
C
C.... ARRAYS FOR THE INTERPOLATION SUBROUTINE :
C
      DO 10 IQ = 1, NQ
         DO 20 IX = 1, NX-1
            XB0 = XB(IX) 
            XB1 = 1.D0-XB(IX)
            XUF(IX,IQ)  = PARTON(1,IQ,IX) / (XB1**3 * XB0)
            XDF(IX,IQ)  = PARTON(2,IQ,IX) / (XB1**4 * XB0)
            XUBF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**8 * XB0**0.5) 
            XDBF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**8 * XB0**0.5) 
            XSF(IX,IQ)  = PARTON(5,IQ,IX) / (XB1**8 * XB0**0.5) 
            XGF(IX,IQ)  = PARTON(6,IQ,IX) / (XB1**5 * XB0**0.5)
 20      CONTINUE
         XUF(NX,IQ)  = 0.D0
         XDF(NX,IQ)  = 0.D0
         XUBF(NX,IQ) = 0.D0
         XDBF(NX,IQ) = 0.D0
         XSF(NX,IQ)  = 0.D0
         XGF(NX,IQ)  = 0.D0
 10   CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
         ARRF(IX) = DLOG(XB(IX))
 30   CONTINUE
      DO 40 IQ = 1, NQ
         ARRF(NX+IQ) = DLOG(QS(IQ))
 40   CONTINUE
chs 16   CONTINUE
C...
      RETURN
      END
*
*...CERN LIBRARY ROUTINE E104 (INTERPOLATION) :
*
      FUNCTION DFINT(NARG,ARG,NENT,ENT,TABLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARG(2),NENT(2),ENT(100),TABLE(1200)
chs      DIMENSION ARG(5),NENT(5),ENT(73),TABLE(1200)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
      DO 5 I=1,NARG
         NCOMB(I)=1
         JB=JA-1+NENT(I)
         DO 2 J=JA,JB
            IF (ARG(I).LE.ENT(J)) GO TO 3
 2       CONTINUE
         J=JB
 3       IF (J.NE.JA) GO TO 4
         J=J+1
 4       JR=J-1
         D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
         IENT(I)=J-JA
         KD=KD+IENT(I)*M
         M=M*NENT(I)
         JA=JB+1
 5    CONTINUE
      DFINT=0.D0
 10   FAC=1.D0
      IADR=KD
      IFADR=1
      DO 15 I=1,NARG
         IF (NCOMB(I).EQ.0) GO TO 12
         FAC=FAC*(1.D0-D(I))
         GO TO 15
 12      FAC=FAC*D(I)
         IADR=IADR-IFADR
         IFADR=IFADR*NENT(I)
 15   CONTINUE
      DFINT=DFINT+FAC*TABLE(IADR)
      IL=NARG
 40   IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
      DO 50  K=IL,NARG
         NCOMB(K)=1
 50   CONTINUE
      GO TO 10
 80   IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END


C=============================================================================
C=============================================================================


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C
C       POLARIZED PARTON DISTRIBUTIONS (NLO and LO)                          C
C       FROM HEP-PH/0007068 - PRD62(2000)094025                              C
C       D. de Florian and R. Sassot                                          C
C                                                                            C
C    The distributions are obtained by calling the subroutine                C
C                                                                            C
C     DS2000(SET,X,Q2,DUV,DDV,DUBAR,DDBAR,DSTR,DGLU,G1P,G1N,G1D)             C
C                                                                            C
C Where                                                                      C
C        SET: 1 SET i+   NLO (MSbar)                                         C
C             2 SET i-   NLO (MSbar)                                         C
C             3 SET ii+  NLO (MSbar)                                         C
C             4 SET ii-  NLO (MSbar)                                         C
C             5 SET iii+ NLO (MSbar)                                         C
C             6 SET iii- NLO (MSbar)                                         C
C             7 SET i+   LO                                                  C
C             8 SET i-   LO                                                  C
C             9 SET ii+  LO                                                  C
C            10 SET ii-  LO                                                  C
C            11 SET iii+ LO                                                  C
C            12 SET iii- LO                                                  C
C                                                                            C
C            X= x_bjorken                                                    C
C            Q2=Q^2                                                          C
C                    DUV :    X * U VALENCE DISTRIBUTION                     C
C                    DDV :    X * D VALENCE DISTRIBUTION                     C
C                    DUBAR :  X * UBAR DISTRIBUTION                          C
C                    DDBAR :  X * DBAR DISTRIBUTION                          C
C                    DSTR :   X * STR= X * STRBAR DISTRIBUTION               C
C                    DGLU :   X * GLUON DISTRIBUTION                         C
C                    G1P :    X * POLARIZED STRUCTURE FUNCTION (PROTON)	     C
C                    G1N :    X * POLARIZED STRUCTURE FUNCTION (NEUTRON)     C
C                    G1D :    X * POLARIZED STRUCTURE FUNCTION(DEUTERON)     C
C                                                                            C
C       REMEMBER: ALWAYS X*DISTRIBUTION !!!                                  C
C                                                                            C
C       BEFORE CALLING THE SUBRUTINE `DS2000` FOR THE FIRST TIME, THE        C
C       SUBROUTINE `INIDS` MUST BE CALLED (ONLY ONCE!) TO READ THE GRIDS.    C
C       (See the test program above as an example)                           C
C                                                                            C
C       RANGE OF VALIDITY OF THE INTERPOLATION:                              C
C       10**(-4)< X < 1.0                                                    C
C       1 < Q**2 < 5*10**4                                                   C
C                                                                            C
C       IN CASE OF PROBLEMS, DOUBTS, ETC, PLEASE REPORT TO                   C
C        daniel.deflorian@cern.ch                                            C
C        sassot@df.uba.ar                                                    C
C                                                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        SUBROUTINE DS2000(MODE,Xin,Q2in,DUV,DDV,DUBAR,DDBAR,
     #                             DSTR,DGLU,G1P,G1N,G1D)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION GU(76,26,12), GD(76,26,12),GUB(76,26,12),
     #            GDB(76,26,12),GS(76,26,12), GG(76,26,12), 
     #   XQ(2), GP(76,26,12), GN(76,26,12)
        COMMON/ GRIDDS / GU,GD,GUB,GDB,GS,GG,GP,GN  

C...Added for interface to DJANGOH (TM&HS, 10.01.2013)
      DOUBLE PRECISION Xin, Q2in
      integer ixoutmn,iqoutmn,ixoutmx,iqoutmx
      DOUBLE PRECISION XMIND, XMAXD, Q2MIND, Q2MAXD
      COMMON /HSDPEX/ ixoutmn, ixoutmx, iqoutmn, iqoutmx,
     *                XMIND, XMAXD, Q2MIND, Q2MAXD
      logical FIRST
      data FIRST /.true./

      IF (FIRST) THEN
        ixoutmn=0
        ixoutmx=0
        iqoutmn=0
        iqoutmx=0
        FIRST=.FALSE.
      ENDIF
      X=Xin
      Q2=Q2in
      XMIND=1.0D-4
      XMAXD=1.0D0
      Q2MIND=1.0D0
      Q2MAXD=5.D4
      if (Xin.lt.XMIND) then 
        X=XMIND
        ixoutmn=ixoutmn+1
      endif
      if (Xin.gt.XMAXD) then
        X=XMAXD
        ixoutmx=ixoutmx+1
      endif
      if (Q2in.lt.Q2MIND) then 
        Q2=Q2MIND
        iqoutmn=iqoutmn+1
      endif
      if (Q2in.gt.Q2MAXD) then
        Q2=Q2MAXD
        iqoutmx=iqoutmx+1
      endif

        XQ(1) = DLOG(X)
        XQ(2) = DLOG(Q2)
        X3=(1.D0-X)**3.D0
        X4=(1.D0-X)**4.D0 
        X5=X**0.5D0
        X6=X**0.6D0
        X7=X**0.7D0
        DUV = PERINOLA(XQ,GU(1,1,MODE)) * X3* X6
        DDV = PERINOLA(XQ,GD(1,1,MODE)) * X4 * X7 
        DUBAR = PERINOLA(XQ,GUB(1,1,MODE)) * X3 * X5
        DDBAR = PERINOLA(XQ,GDB(1,1,MODE)) * X3 * X5
        DSTR = PERINOLA(XQ,GS(1,1,MODE))  * X3 * X5
        DGLU = PERINOLA(XQ,GG(1,1,MODE))  * X3 * X5
        G1P = PERINOLA(XQ,GP(1,1,MODE))  * X3 * X5
        G1N = PERINOLA(XQ,GN(1,1,MODE))  * X4 * X5
        G1D=(G1P+G1N)*0.5D0*(1.D0-1.5D0*0.058D0)
        RETURN
        END


      FUNCTION PERINOLA(ARG,TABLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARG(2),NENT(2),ENT(102),TABLE(1976)
      DIMENSION D(8),NCOMB(8),IENT(8)
      COMMON/ XARRAY / ENT
      NARG=2
      NENT(1)=76
      NENT(2)=26
      KD=1
      M=1
      JA=1
      DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
      DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
      JA=JB+1
 5    CONTINUE
      PERINOLA=0.D0
   10 FAC=1.D0
      IADR=KD
      IFADR=1
      DO 151 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.D0-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
 151  CONTINUE
      PERINOLA=PERINOLA+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
      DO 50  K=IL,NARG
        NCOMB(K)=1
 50   CONTINUE
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END


      SUBROUTINE INIDS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XARRAY(102), GU(76,26,12), GD(76,26,12), 
     #  GUB(76,26,12),GDB(76,26,12),GS(76,26,12), GG(76,26,12),
     #  GP(76,26,12), GN(76,26,12) 
       COMMON/ XARRAY / XARRAY
       COMMON/ GRIDDS / GU,GD,GUB,GDB,GS,GG,GP,GN
       integer lunplx
       dimension lunplx(12)
       
c       OPEN(UNIT=80,FILE='polpdf-gridfiles/i+NLO.GRID',STATUS='OLD')
c       OPEN(UNIT=81,FILE='polpdf-gridfiles/i-NLO.GRID',STATUS='OLD')
c       OPEN(UNIT=82,FILE='polpdf-gridfiles/ii+NLO.GRID',STATUS='OLD')
c       OPEN(UNIT=83,FILE='polpdf-gridfiles/ii-NLO.GRID',STATUS='OLD')
c       OPEN(UNIT=84,FILE='polpdf-gridfiles/iii+NLO.GRID',STATUS='OLD')
c       OPEN(UNIT=85,FILE='polpdf-gridfiles/iii-NLO.GRID',STATUS='OLD')
c       OPEN(UNIT=86,FILE='polpdf-gridfiles/i+LO.GRID',STATUS='OLD')
c       OPEN(UNIT=87,FILE='polpdf-gridfiles/i-LO.GRID',STATUS='OLD')
c       OPEN(UNIT=88,FILE='polpdf-gridfiles/ii+LO.GRID',STATUS='OLD')
c       OPEN(UNIT=89,FILE='polpdf-gridfiles/ii-LO.GRID',STATUS='OLD')
c       OPEN(UNIT=90,FILE='polpdf-gridfiles/iii+LO.GRID',STATUS='OLD')
c       OPEN(UNIT=91,FILE='polpdf-gridfiles/iii-LO.GRID',STATUS='OLD')

       lunplx(1)=NextUn()
       OPEN(lunplx(1),FILE='polpdf-gridfiles/i+NLO.GRID',STATUS='OLD')
       lunplx(2)=NextUn()
       OPEN(lunplx(2),FILE='polpdf-gridfiles/i-NLO.GRID',STATUS='OLD')
       lunplx(3)=NextUn()
       OPEN(lunplx(3),FILE='polpdf-gridfiles/ii+NLO.GRID',STATUS='OLD')
       lunplx(4)=NextUn()
       OPEN(lunplx(4),FILE='polpdf-gridfiles/ii-NLO.GRID',STATUS='OLD')
       lunplx(5)=NextUn()
       OPEN(lunplx(5),FILE='polpdf-gridfiles/iii+NLO.GRID',STATUS='OLD')
       lunplx(6)=NextUn()
       OPEN(lunplx(6),FILE='polpdf-gridfiles/iii-NLO.GRID',STATUS='OLD')
       lunplx(7)=NextUn()
       OPEN(lunplx(7),FILE='polpdf-gridfiles/i+LO.GRID',STATUS='OLD')
       lunplx(8)=NextUn()
       OPEN(lunplx(8),FILE='polpdf-gridfiles/i-LO.GRID',STATUS='OLD')
       lunplx(9)=NextUn()
       OPEN(lunplx(9),FILE='polpdf-gridfiles/ii+LO.GRID',STATUS='OLD')
       lunplx(10)=NextUn()
       OPEN(lunplx(10),FILE='polpdf-gridfiles/ii-LO.GRID',STATUS='OLD')
       lunplx(11)=NextUn()
       OPEN(lunplx(11),FILE='polpdf-gridfiles/iii+LO.GRID',STATUS='OLD')
       lunplx(12)=NextUn()
       OPEN(lunplx(12),FILE='polpdf-gridfiles/iii-LO.GRID',STATUS='OLD')
       
c       DO IFILE=80,91
        DO if=1,12
        IFILE=lunplx(if)
        DO  K = 1, 76 
        DO  J = 1, 26
c        I=IFILE-79
        I=if
        READ(IFILE,40) GU(K,J,I), GD(K,J,I), GUB(K,J,I), GDB(K,J,I),
     #                 GS(K,J,I), GG(K,J,I), GP(K,J,I), GN(K,J,I)
        END DO
        END DO
       CLOSE(IFILE)
       END DO
  40   FORMAT (8(1PE11.4))
       DO LX1=0,25,1
        XX=10.D0**(-4.D0+LX1/25.D0*3.D0)
        XARRAY(LX1+1) = DLOG(XX)
       END DO
       DO LX1=1,50,1
        XX=0.1D0+LX1/50.D0*0.9D0
        XARRAY(LX1+26) = DLOG(XX)
       END DO
       DO LQ=0,25,1
        QQ=0.6D0*10**(5.D0*LQ/25.D0)
        XARRAY(LQ+77) = DLOG(QQ)
      END DO

      close(80)
      close(81)
      close(82)
      close(83)
      close(84)
      close(85)
      close(86)
      close(87)
      close(88)
      close(89)
      close(90)
      close(91)
      
      RETURN
      END 


C=============================================================================
C=============================================================================


C *********************************************************************
C  aac08.f  Version 1.2                                2008/AUG/05
C
C  [Package for the AAC polarized PDFs]
C  Reference: 
C         Determination of gluon polarization 
C             from deep inelastic scattering and collider data
C         M. Hirai, and S. Kumano
C         arXiv:0808.0413.
C
C *********************************************************************
C ---------------------------------------------------------------------
C  SUBROUTINE AAC08PDF(Q2,X,ISET,XPPDF,GRAD):
C
C   Subroutine AAC08PDF returns the values of polarized PDFs
C   and their gradient terms at specified Q^2 and x point 
C   by interpolating the grid data.
C   [ Log(Q^2): LINEAR INTERPOLATION, x: CUBIC SPLINE INTERPOLATION ]
C
C   INPUT:
C     Q2, X ... Q^2 and x values at which the functions are calculated. 
C               Available range: 10^-9 <= X <= 1.0,
C                           1.0 GeV^2 <= Q^2 <= 10^8 GeV^2.
C     ISET=1: Positive type PDFs and their gradient terms of Set-B
C          2: Node type PDFs and their gradient terms
C
C   OUTPUT: Arrays XPPDF(-3:3) & GRAD(I,J)
C
C     XPPDF(I) --> AAC08 polarized PDFs.
C      I = -3 ... s-bar quark
C          -2 ... d-bar quark 
C          -1 ... u-bar quark 
C           0 ... gluon D_g(x)
C           1 ... u quark 
C           2 ... d quark 
C           3 ... s quark 
C
C     GRAD(I,J) --> Gradient terms of AAC08 PDFs
C      I is the same index as the one in XPPDF(I).
C      J indicates the parameter index for a gradient term dXPPDF(I)/da_J
C      (a_J = parameter).
C
C      J= 1..4: g   (delta, nu, kappa, mu)
C         5, 6: d_v (delta, mu)
C         7, 8: u_v (delta, mu)
C        9..11: qb  (delta, nu, kappa)
C
C   For example, the above J=1 indicates d XPPDF_g/d delta_g.
C
C   NOTE: The returned values are PPDF multiplied by x.
C
C      *  Error matrix can be used by declaring a common block:
C         COMMON/ERRM/EM(11,11). This matrix is defined as
C         the inverse matrix of Hessian multiplied by Delta chi^2:
C         EM(i,j)=Delta chi^2*H_ij^-1.
C         The values of Delta chi^2 is 12.647
C *********************************************************************
      SUBROUTINE AAC08PDF(Q2in,Xin,ISET,XPPDF,GRAD)
C ---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NQ=33, NX=117, ND=74, NPDF=6, NSET=3, IFILE=20)
      DIMENSION IREAD(2), QG(NQ), XG(NX),PDFJ1(ND), PDFJ2(ND)
     +         ,XPPDF(-3:3), GRADFF(6,11), GRAD(-3:3,11)
     +         ,BXG(NX,NQ,ND), CXG(NX,NQ,ND), DXG(NX,NQ,ND)
     +         ,PDFG(NX,NQ,ND), EMI(11,11,2)

      COMMON/ERRM/EM(11,11)
      SAVE IREAD, NPAR, BXG, CXG, DXG 
      DATA IREAD /1, 1/ 
      DATA EPS/1.D-12/

C...Added for interface to DJANGOH (TM&HS, 10.01.2013)
      DOUBLE PRECISION Xin, Q2in
      integer ixoutmn,iqoutmn,ixoutmx,iqoutmx
      DOUBLE PRECISION XMIND, XMAXD, Q2MIND, Q2MAXD
      COMMON /HSDPEX/ ixoutmn, ixoutmx, iqoutmn, iqoutmx,
     *                XMIND, XMAXD, Q2MIND, Q2MAXD
      logical FIRST
      data FIRST /.true./
      integer lunplx
      dimension lunplx(2)

C      INCLUDE 'EM.inc'
      DATA ((EMI(I,J,1),J=1,11),I=1,11)/
     + 1.808521E+02,-1.416464E+01, 9.510544E+00, 6.525852E-01,
     + 2.820281E-01,-6.032619E-01,
     + 2.175284E-01, 6.551146E-01,
     +-2.453518E+00, 8.296432E-01,-2.124696E-01,
     +-1.416464E+01, 1.454405E+01, 9.523191E-01,-8.878194E+00,
     +-2.655870E-02, 2.959398E-01,
     +-3.781453E-02,-1.593522E-01,
     +-2.617929E-01,-2.213225E-01, 1.669404E-02,
     + 9.510544E+00, 9.523191E-01, 2.314401E+00, 1.568228E-01,
     + 1.353229E-02, 2.086755E-01,
     +-3.250279E-02,-1.215377E-01,
     +-4.034393E-01,-1.467052E-01, 9.814072E-03,
     + 6.525852E-01,-8.878194E+00, 1.568228E-01, 7.069673E+00,
     +-4.527626E-03, 3.275573E-02,
     +-5.273799E-03,-1.745286E-02,
     +-3.035280E-02,-1.682051E-02, 1.631463E-03,
     + 2.820281E-01,-2.655870E-02, 1.353229E-02,-4.527626E-03,
     + 6.171736E-01,-1.229288E+00,
     +-5.817620E-02,-3.528513E-02,
     + 2.023520E-01, 2.567341E-02, 3.237632E-03,
     +-6.032619E-01, 2.959398E-01, 2.086755E-01, 3.275573E-02,
     +-1.229288E+00, 3.857335E+00,
     +-3.187044E-02,-2.276460E-01,
     +-2.225872E+00, 5.425563E-02,-7.929669E-02,
     + 2.175284E-01,-3.781453E-02,-3.250279E-02,-5.273799E-03,
     +-5.817620E-02,-3.187044E-02,
     + 1.231818E-01, 2.757046E-01,
     + 4.704684E-01,-9.244957E-03, 1.403817E-02,
     + 6.551146E-01,-1.593522E-01,-1.215377E-01,-1.745286E-02,
     +-3.528513E-02,-2.276460E-01,
     + 2.757046E-01, 7.714670E-01,
     + 1.656757E+00,-7.031732E-02, 5.956737E-02,
     +-2.453518E+00,-2.617929E-01,-4.034393E-01,-3.035280E-02,
     + 2.023520E-01,-2.225872E+00,
     + 4.704684E-01, 1.656757E+00,
     + 9.105840E+00,-1.155936E+00, 4.350568E-01,
     + 8.296432E-01,-2.213225E-01,-1.467052E-01,-1.682051E-02,
     + 2.567341E-02, 5.425563E-02,
     +-9.244957E-03,-7.031732E-02,
     +-1.155936E+00, 2.782340E-01,-7.019085E-02,
     +-2.124696E-01, 1.669404E-02, 9.814072E-03, 1.631463E-03,
     + 3.237632E-03,-7.929669E-02,
     + 1.403817E-02, 5.956737E-02,
     + 4.350568E-01,-7.019085E-02, 2.352342E-02/

      DATA ((EMI(I,J,2),J=1,11),I=1,11)/
     + 4.426450E+01,-5.463504E+00,-9.674955E-01,-2.036167E-01,
     + 2.377636E-01,-2.099402E+00,
     + 1.631463E-01, 6.842027E-01,
     +-2.946751E-01, 2.681164E+00,-5.096741E-01,
     +-5.463504E+00, 7.373201E+00,-3.224985E+00,-1.441758E+01,
     +-1.059819E-02, 2.579988E-01,
     +-1.129377E-02,-7.196143E-02,
     +-4.527626E-01,-1.859109E-01, 7.107614E-03,
     +-9.674955E-01,-3.224985E+00, 1.707345E+00, 7.208790E+00,
     +-1.884403E-02, 5.627915E-02,
     +-6.361441E-03,-1.871756E-02,
     + 2.643223E-01,-1.441758E-01, 4.059687E-02,
     +-2.036167E-01,-1.441758E+01, 7.208790E+00, 3.237632E+01,
     +-2.769693E-02, 7.917022E-02,
     +-5.109388E-03,-1.871756E-02,
     + 3.161750E-02,-6.121148E-02, 1.105348E-02,
     + 2.377636E-01,-1.059819E-02,-1.884403E-02,-2.769693E-02,
     + 8.650548E-01,-1.859109E+00,
     +-5.729091E-02,-1.795874E-02,
     + 5.058800E-01, 1.183759E-03, 1.618816E-02,
     +-2.099402E+00, 2.579988E-01, 5.627915E-02, 7.917022E-02,
     +-1.859109E+00, 5.956737E+00,
     +-4.932330E-02,-3.262926E-01,
     +-3.591748E+00, 1.909697E-01,-1.593522E-01,
     + 1.631463E-01,-1.129377E-02,-6.361441E-03,-5.109388E-03,
     +-5.729091E-02,-4.932330E-02,
     + 9.953189E-02, 2.149990E-01,
     + 4.198804E-01,-4.755272E-03, 1.115465E-02,
     + 6.842027E-01,-7.196143E-02,-1.871756E-02,-1.871756E-02,
     +-1.795874E-02,-3.262926E-01,
     + 2.149990E-01, 6.032619E-01,
     + 1.618816E+00,-7.739964E-02, 6.386735E-02,
     +-2.946751E-01,-4.527626E-01, 2.643223E-01, 3.161750E-02,
     + 5.058800E-01,-3.591748E+00,
     + 4.198804E-01, 1.618816E+00,
     + 1.196406E+01,-1.972932E+00, 8.005551E-01,
     + 2.681164E+00,-1.859109E-01,-1.441758E-01,-6.121148E-02,
     + 1.183759E-03, 1.909697E-01,
     +-4.755272E-03,-7.739964E-02,
     +-1.972932E+00, 6.386735E-01,-1.934991E-01,
     +-5.096741E-01, 7.107614E-03, 4.059687E-02, 1.105348E-02,
     + 1.618816E-02,-1.593522E-01,
     + 1.115465E-02, 6.386735E-02,
     + 8.005551E-01,-1.934991E-01, 6.702910E-02/
     
C Q2 AND X GRID.
      DATA QG /
     +  1.000000D+00, 1.467799D+00, 2.154435D+00,
     +  3.162278D+00, 4.641589D+00, 6.812921D+00,
     +  1.000000D+01, 1.467799D+01, 2.154435D+01,
     +  3.162278D+01, 4.641589D+01, 6.812921D+01,
     +  1.000000D+02, 1.778279D+02, 3.162278D+02, 5.623413D+02,
     +  1.000000D+03, 1.778279D+03, 3.162278D+03, 5.623413D+03,
     +  1.000000D+04, 1.778279D+04, 3.162278D+04, 5.623413D+04,
     +  1.000000D+05, 1.778279D+05, 3.162278D+05, 5.623413D+05,
     +  1.000000D+06, 4.641589D+06, 
     +  1.000000D+07, 4.641589D+07,  
     +  1.000000D+08  /

      DATA XG / 
     +  1.000000D-09, 1.333521D-09, 1.778279D-09, 2.371374D-09,
     +  3.162278D-09, 4.216965D-09, 5.623413D-09, 7.498942D-09,
     +  1.000000D-08, 1.333521D-08, 1.778279D-08, 2.371374D-08,
     +  3.162278D-08, 4.216965D-08, 5.623413D-08, 7.498942D-08,
     +  1.000000D-07, 1.333521D-07, 1.778279D-07, 2.371374D-07,
     +  3.162278D-07, 4.216965D-07, 5.623413D-07, 7.498942D-07,
     +  1.000000D-06, 1.333521D-06, 1.778279D-06, 2.371374D-06,
     +  3.162278D-06, 4.216965D-06, 5.623413D-06, 7.498942D-06,
     +  1.000000D-05, 1.333521D-05, 1.778279D-05, 2.371374D-05,
     +  3.162278D-05, 4.216965D-05, 5.623413D-05, 7.498942D-05,
     +  1.000000D-04, 1.333521D-04, 1.778279D-04, 2.371374D-04,
     +  3.162278D-04, 4.216965D-04, 5.623413D-04, 7.498942D-04,
     +  1.000000D-03, 1.154782D-03, 1.333521D-03, 1.539927D-03,
     +  1.778279D-03, 2.053525D-03, 2.371374D-03, 2.738420D-03,
     +  3.162278D-03, 3.651741D-03, 4.216965D-03, 4.869675D-03,
     +  5.623413D-03, 6.493816D-03, 7.498942D-03, 8.659643D-03,
     +  1.000000D-02, 1.154782D-02, 1.333521D-02, 1.539927D-02,
     +  1.778279D-02, 2.053525D-02, 2.371374D-02, 2.738420D-02,
     +  3.162278D-02, 3.651741D-02, 4.216965D-02, 4.869675D-02,
     +  5.623413D-02, 6.493816D-02, 7.498942D-02, 8.659643D-02,
     +  1.000000D-1, 1.250000D-1, 1.500000D-1, 1.750000D-1,
     +  2.000000D-1, 2.250000D-1, 2.500000D-1, 2.750000D-1,
     +  3.000000D-1, 3.250000D-1, 3.500000D-1, 3.750000D-1,
     +  4.000000D-1, 4.250000D-1, 4.500000D-1, 4.750000D-1, 
     +  5.000000D-1, 5.250000D-1, 5.500000D-1, 5.750000D-1,
     +  6.000000D-1, 6.250000D-1, 6.500000D-1, 6.750000D-1,
     +  7.000000D-1, 7.250000D-1, 7.500000D-1, 7.750000D-1,
     +  8.000000D-1, 8.250000D-1, 8.500000D-1, 8.750000D-1,
     +  9.000000D-1, 9.250000D-1, 9.500000D-1, 9.750000D-1,
     +  1.000000D+0 /

C CALCULATE SPLINE COEFFICIENTS.
      IF((IREAD(1).NE.1).AND.(IREAD(2).EQ.ISET)) GO TO 20

C READ GRID DATA AND CALCULATE SPLINE COEFFICIENTS.
      IF(ISET.EQ.1) THEN
c        OPEN(UNIT=IFILE,FILE='polpdf-gridfiles/aacnlo08_posi.grd',
        lunplx(1)=NextUn()
        OPEN(lunplx(1),FILE='polpdf-gridfiles/aacnlo08_posi.grd',
     +  STATUS='OLD')
c        OPEN(UNIT=IFILE+1,FILE='polpdf-gridfiles/grad_posi.grd',
        lunplx(2)=NextUn()
        OPEN(lunplx(2),FILE='polpdf-gridfiles/grad_posi.grd',
     +  STATUS='OLD')
        NPAR=11

      ELSE IF(ISET.EQ.2) THEN
c        OPEN(UNIT=IFILE,FILE='polpdf-gridfiles/aacnlo08_node.grd',
        lunplx(1)=NextUn()
        OPEN(lunplx(1),FILE='polpdf-gridfiles/aacnlo08_node.grd',
     +  STATUS='OLD')
c        OPEN(UNIT=IFILE+1,FILE='polpdf-gridfiles/grad_node.grd',
        lunplx(2)=NextUn()
        OPEN(lunplx(2),FILE='polpdf-gridfiles/grad_node.grd',
     +  STATUS='OLD')
        NPAR=11

      ELSE
        WRITE(*,1010) ISET
 1010   FORMAT(' ','AAC08PDF ERROR: ISET =', I3)
        STOP
      END IF

      DO I=1,11
        DO J=1,11
          EM(I,J)=EMI(I,J,ISET)
        ENDDO
      ENDDO

      DO J=1,NQ
        DO K=1,NX-1
c          READ(IFILE,1025) (PDFG(K,J,I), I=1,6)
          READ(lunplx(1),1025) (PDFG(K,J,I), I=1,6)
          DO NR=1,NPAR
            NI=8+NPDF*(NR-1)
c            READ(IFILE+1,1025) (PDFG(K,J,I),I=NI,NI+NPDF-1)
            READ(lunplx(2),1025) (PDFG(K,J,I),I=NI,NI+NPDF-1)
          ENDDO
        ENDDO
      ENDDO

 1025 FORMAT(1X,6(1PE14.5))
c      CLOSE(IFILE)
c      CLOSE(IFILE+1)
      close(lunplx(1))
      close(lunplx(2))

      DO I=1,ND
        DO J=1,NQ
          PDFG(NX,J,I)=0.D0 ! x=1 XPPDF=0.D0
          CALL LSPLINE1(NX,XG,PDFG,BXG,CXG,DXG,I,J)
        ENDDO
      ENDDO

      IREAD(1)=2
      IREAD(2)=ISET
   20 CONTINUE

      DO I=1,7
        XPPDF(I-4)=0.D0
      END DO
      DO J=1,11
        DO I=1,6
          GRADFF(I,J)=0.D0
        ENDDO
      ENDDO

C CHECK X AND Q2 VALUES.  
      IF (FIRST) THEN
        ixoutmn=0
        ixoutmx=0
        iqoutmn=0
        iqoutmx=0
        FIRST=.FALSE.
      ENDIF    
      X=Xin
      Q2=Q2in
      XMIND=1.0D-9
      XMAXD=1.0D0
      Q2MIND=1.0D0
      Q2MAXD=1.D8
      if (Xin.lt.XMIND) then 
        X=XMIND
        ixoutmn=ixoutmn+1
      endif
      if (Xin.gt.XMAXD) then
        X=XMAXD
        ixoutmx=ixoutmx+1
      endif
      if (Q2in.lt.Q2MIND) then 
        Q2=Q2MIND
        iqoutmn=iqoutmn+1
      endif
      if (Q2in.gt.Q2MAXD) then
        Q2=Q2MAXD
        iqoutmx=iqoutmx+1
      endif

C INTERPOLATION.
C X: CUBIC SPLINE INTERPOLATION, LOG(Q2): LINEAR INTERPOLATION.
chs      J=ISERCH1(NQ,QG,Q2)
      J=ISRCH1Q(NQ,QG,Q2)
      IF(J.EQ.NQ) J=NQ-1
chs      K=ISERCH1(NX,XG,X)
      K=ISRCH1X(NX,XG,X)
      DO I=1,ND
        DX=X-XG(K)
        PDFJ1(I)=PDFG(K,J,I)
     >       +DX*(BXG(K,J,I)+DX*(CXG(K,J,I)+DX*DXG(K,J,I)))
        PDFJ2(I)=PDFG(K,J+1,I)
     >       +DX*(BXG(K,J+1,I)+DX*(CXG(K,J+1,I)+DX*DXG(K,J+1,I)))
      ENDDO
 
C -- Polarized PDFs --
      T=(DLOG(Q2)-DLOG(QG(J)))/(DLOG(QG(J+1))-DLOG(QG(J)))
      DO I=1,3
        XPPDF(I-1)=(1.D0-T)*PDFJ1(I)+T*PDFJ2(I)     ! g, u, d
        XPPDF(-I)=(1.D0-T)*PDFJ1(I+3)+T*PDFJ2(I+3)  ! ub, db, sb
      ENDDO
      XPPDF(3)=XPPDF(-3)

C -- Gradient terms of parameters for polarized PDFs --
      DO J=1,NPAR
        DO I=1,NPDF
          NI=7+NPDF*(J-1)
          GRADFF(I,J)=(1.D0-T)*PDFJ1(NI+I)+T*PDFJ2(NI+I)
        ENDDO
      ENDDO

      DO J=1,NPAR
        DO I=1,3
          GRAD(I-1,J)=GRADFF(I,J)   ! g, u, d
          GRAD(-I,J)=GRADFF(I+3,J)  ! ub, db, sb
        ENDDO
        GRAD(3,J)=GRAD(-3,J) ! s=sb
      ENDDO

      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE LSPLINE1(N,X,Y,B,C,D,I,J)
C ---------------------------------------------------------------------
C CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.
C INTERPOLATION SUBROUTINES ARE TAKEN FROM
C G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,
C COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977).
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NQ=33, NX=117, ND=74)
      DIMENSION Y(NX,NQ,ND),B(NX,NQ,ND),C(NX,NQ,ND),D(NX,NQ,ND)
     1         ,X(NX) 
      NM1=N-1
      IF(N.LT.2) RETURN
      IF(N.LT.3) GO TO 250
      D(1,J,I)=X(2)-X(1)
      C(2,J,I)=(Y(2,J,I)-Y(1,J,I))/D(1,J,I)
      DO 210 K=2,NM1
        D(K,J,I)=X(K+1)-X(K)
        B(K,J,I)=2.0D0*(D(K-1,J,I)+D(K,J,I))
        C(K+1,J,I)=(Y(K+1,J,I)-Y(K,J,I))/D(K,J,I)
        C(K,J,I)=C(K+1,J,I)-C(K,J,I)
  210 CONTINUE
      B(1,J,I)=-D(1,J,I)
      B(N,J,I)=-D(N-1,J,I)
      C(1,J,I)=0.0D0
      C(N,J,I)=0.0D0
      IF(N.EQ.3) GO TO 215
      C(1,J,I)=C(3,J,I)/(X(4)-X(2))-C(2,J,I)/(X(3)-X(1))
      C(N,J,I)=C(N-1,J,I)/(X(N)-X(N-2))-C(N-2,J,I)/(X(N-1)-X(N-3))
      C(1,J,I)=C(1,J,I)*D(1,J,I)**2.0D0/(X(4)-X(1))
      C(N,J,I)=-C(N,J,I)*D(N-1,J,I)**2.0D0/(X(N)-X(N-3))
  215 CONTINUE
      DO 220 K=2,N
        T=D(K-1,J,I)/B(K-1,J,I)
        B(K,J,I)=B(K,J,I)-T*D(K-1,J,I)
        C(K,J,I)=C(K,J,I)-T*C(K-1,J,I)
  220 CONTINUE
      C(N,J,I)=C(N,J,I)/B(N,J,I)
      DO 230 IB=1,NM1
        K=N-IB
        C(K,J,I)=(C(K,J,I)-D(K,J,I)*C(K+1,J,I))/B(K,J,I)
  230 CONTINUE
      B(N,J,I)=(Y(N,J,I)-Y(NM1,J,I))/D(NM1,J,I)
     1        +D(NM1,J,I)*(C(NM1,J,I)+2.0D0*C(N,J,I))
      DO 240 K=1,NM1
        B(K,J,I)=(Y(K+1,J,I)-Y(K,J,I))/D(K,J,I)
     1          -D(K,J,I)*(C(K+1,J,I)+2.0D0*C(K,J,I))
        D(K,J,I)=(C(K+1,J,I)-C(K,J,I))/D(K,J,I)
        C(K,J,I)=3.0D0*C(K,J,I)
  240 CONTINUE
      C(N,J,I)=3.0D0*C(N,J,I)
      D(N,J,I)=D(N-1,J,I)
      RETURN
  250 CONTINUE
      B(1,J,I)=(Y(2,J,I)-Y(1,J,I))/(X(2)-X(1))
      C(1,J,I)=0.0D0
      D(1,J,I)=0.0D0
      B(2,J,I)=B(1,J,I)
      C(2,J,I)=0.0D0
      D(2,J,I)=0.0D0
      RETURN
      END
C ---------------------------------------------------------------------
      INTEGER FUNCTION ISRCH1X(N,X,Y)
C ---------------------------------------------------------------------
C THIS FUNCTION SEARCHES "I" WHICH SATISFIES THE RELATION
C X(I) <= Y < X(I+1) BY USING A BINARY SEARCH.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
chs      DIMENSION X(118)
      DIMENSION X(117)

      MIN=1
      MAX=N+1

   10 CONTINUE
      MID=(MIN+MAX)/2
      IF(Y.LT.X(MID)) THEN
        MAX=MID
      ELSE
        MIN=MID
      END IF
      IF((MAX-MIN).GT.1) GO TO 10

chs      ISERCH1=MIN
      ISRCH1X=MIN

      RETURN
      END
C ---------------------------------------------------------------------
      INTEGER FUNCTION ISRCH1Q(N,X,Y)
C ---------------------------------------------------------------------
C THIS FUNCTION SEARCHES "I" WHICH SATISFIES THE RELATION
C X(I) <= Y < X(I+1) BY USING A BINARY SEARCH.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
chs      DIMENSION X(118)
      DIMENSION X(33)

      MIN=1
      MAX=N+1

   10 CONTINUE
      MID=(MIN+MAX)/2
      IF(Y.LT.X(MID)) THEN
        MAX=MID
      ELSE
        MIN=MID
      END IF
      IF((MAX-MIN).GT.1) GO TO 10

chs      ISERCH1=MIN
      ISRCH1Q=MIN

      RETURN
      END
C *********************************************************************
C THE END OF THE PROGRAM.
C *********************************************************************


C=============================================================================
C=============================================================================


*********************************************************************
*                                                                   *
*    POLARIZED RADIATIVELY GENERATED LO AND NLO PARTON DENSITIES    *
*                                                                   *
*         M. GLUCK, E. REYA, M. STRATMANN AND W. VOGELSANG,         *
*              Phys.Rev.D63:094005,2001,    hep-ph/0011215          *
*                                                                   *
*               PROBLEMS/QUESTIONS TO wvogelsang@bnl.gov            *
*            OR TO marco.stratmann@physik.uni-regensburg.de         *
*                                                                   *
*   INPUT:   ISET = number of the parton set :                      *
*            ISET = 1  'STANDARD' SCENARIO, NEXT-TO-LEADING ORDER   *
*                      (MS-bar)                                     * 
*                      (DATA FILE 'std2000_nlo.grid' UNIT=11, TO BE *
*                       DEFINED BY THE USER )                       *
*            ISET = 2  'VALENCE' SCENARIO,  NEXT-TO-LEADING ORDER   *
*                      (MS-bar)                                     *   
*                      (DATA FILE 'val2000_nlo.grid' UNIT=22, TO BE *
*                       DEFINED BY THE USER )                       *
*            ISET = 3  'STANDARD' SCENARIO, LEADING ORDER           *
*                      (DATA FILE 'std2000_lo.grid' UNIT=33, TO BE  *
*                       DEFINED BY THE USER )                       *
*            ISET = 4  'VALENCE' SCENARIO,  LEADING ORDER           *
*                      (DATA FILE 'val2000_lo.grid' UNIT=44, TO BE  *
*                       DEFINED BY THE USER )                       *
*                                                                   *
*            X  = Bjorken-x       (between  1.E-4  and  1)          *
*            Q2 = scale in GeV**2 (between  0.8  and   1.E6)        *
*                                                                   *
*   OUTPUT:  U = x * DELTA u                                        *
*            D = x * DELTA d                                        *        
*            UB = x * DELTA ubar                                    *   
*            DB = x * DELTA dbar                                    * 
*            ST = x * DELTA STRANGE                                 *     
*            GL = x * DELTA GLUON                                   *
*            G1P = g_1^proton                                       *
*            G1N = g_1^neutron                                      * 
*                                                                   *
*          (  For the parton distributions always x times           *
*                   the distribution is returned .                  *
*                 This is NOT the case for g1(p,n)  )               *
*                                                                   *
*            The sets are the result of a combined fit to           *
*            data for the spin asymmetries A_1 (p,n,d)              *
*                                                                   *
*            Note: No charm is included                             *
*                                                                   *
*   COMMON:  The main program or the calling routine has to have    *
*            a common block  COMMON / INTINI / IINI , and  IINI     *
*            has always to be zero when PARPOL is called for the    *
*            first time or when 'ISET' has been changed.            *
*                                                                   *
*********************************************************************
*
      SUBROUTINE PARPOL (ISET,Xin,Q2in,U,D,UB,DB,ST,GL,G1P,G1N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPART=8, NX=42, NQ=30, NARG=2)
      DIMENSION XUF(NX,NQ), XDF(NX,NQ), XUBF(NX,NQ), XDBF(NX,NQ), 
     1          XSF(NX,NQ), XGF(NX,NQ), XG1P(NX,NQ), XG1N(NX,NQ),
     2          PARTON (NPART,NQ,NX-1), QS(NQ), XB(NX), XT(NARG), 
     3          NA(NARG), ARRF(100)
chs     3          NA(NARG), ARRF(NX+NQ)
      COMMON / INTINI / IINI
      SAVE XUF, XDF, XUBF, XDBF, XSF, XGF, XG1P, XG1N, NA, ARRF
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 0.8D0, 1.0D0, 1.25d0, 1.5D0, 2.d0, 2.5D0, 
     1           4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2           1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3           3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4, 
     4           1.0D5, 1.8D5, 3.2D5, 5.8D5, 1.0D6  /
       DATA XB / 
     1           1.D-4, 1.5D-4, 2.2D-4, 3.2D-4, 4.8D-4, 7.D-4,
     2           1.D-3, 1.5D-3, 2.2D-3, 3.2D-3, 4.8D-3, 7.D-3,
     3           1.D-2, 1.5D-2, 2.2D-2, 3.2D-2, 5.0D-2, 7.5D-2,
     4           0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5           0.3, 0.325, 0.35, 0.375, 0.4, 0.45,  0.5, 0.55,
     6           0.6, 0.65,  0.7,  0.75,  0.8, 0.85,  0.9, 1.0 /

C...Added for interface to DJANGOH (TM&HS, 10.01.2013)
      DOUBLE PRECISION Xin, Q2in
      integer ixoutmn,iqoutmn,ixoutmx,iqoutmx
      DOUBLE PRECISION XMIND, XMAXD, Q2MIND, Q2MAXD
      COMMON /HSDPEX/ ixoutmn, ixoutmx, iqoutmn, iqoutmx,
     *                XMIND, XMAXD, Q2MIND, Q2MAXD
      logical FIRST
      data FIRST /.true./

*...CHECK OF X AND Q2 VALUES :
      IF (FIRST) THEN
        ixoutmn=0
        ixoutmx=0
        iqoutmn=0
        iqoutmx=0
        FIRST=.FALSE.
      ENDIF
      X=Xin
      Q2=Q2in
      XMIND=1.0D-4
      XMAXD=1.0D0
      Q2MIND=0.8D0
      Q2MAXD=1.D6
      if (Xin.lt.XMIND) then 
        X=XMIND
        ixoutmn=ixoutmn+1
      endif
      if (Xin.gt.XMAXD) then
        X=XMAXD
        ixoutmx=ixoutmx+1
      endif
      if (Q2in.lt.Q2MIND) then 
        Q2=Q2MIND
        iqoutmn=iqoutmn+1
      endif
      if (Q2in.gt.Q2MAXD) then
        Q2=Q2MAXD
        iqoutmx=iqoutmx+1
      endif
*...INITIALIZATION :
*    SELECTION AND READING OF THE GRID :
*    FILE - NO. = 11 FOR NLO 'STANDARD' SCENARIO ( FIRST NUMBER IN THE 
*                                                  GRID: 1.3478E-03 )
*    FILE - NO. = 22 FOR NLO 'VALENCE'  SCENARIO ( FIRST NUMBER IN THE 
*                                                  GRID: 1.5146E-05 )
*    FILE - NO. = 33 FOR  LO 'STANDARD' SCENARIO ( FIRST NUMBER IN THE 
*                                                  GRID: 3.4686E-03 )     
*    FILE - NO. = 44 FOR  LO 'VALENCE'  SCENARIO ( FIRST NUMBER IN THE 
*                                                  GRID: 2.4395E-04 )
      IF (IINI.NE.0) GOTO 16
      IF (ISET.EQ.1) THEN
c       IIREAD=11
       IIREAD=NextUn() 
c       OPEN(UNIT=11,FILE='polpdf-gridfiles/std2000_nlo_g1.grid',
       OPEN(IIREAD,FILE='polpdf-gridfiles/std2000_nlo_g1.grid',
     * STATUS='OLD')
      ELSE IF (ISET.EQ.2) THEN
c       IIREAD=22
       IIREAD=NextUn() 
c       OPEN(UNIT=22,FILE='polpdf-gridfiles/val2000_nlo_g1.grid',
       OPEN(IIREAD,FILE='polpdf-gridfiles/val2000_nlo_g1.grid',
     * STATUS='OLD')
      ELSE IF (ISET.EQ.3) THEN
c       IIREAD=33
       IIREAD=NextUn() 
c       OPEN(UNIT=33,FILE='polpdf-gridfiles/std2000_lo_g1.grid',
       OPEN(IIREAD,FILE='polpdf-gridfiles/std2000_lo_g1.grid',
     * STATUS='OLD')
      ELSE IF (ISET.EQ.4) THEN
c       IIREAD=44
       IIREAD=NextUn() 
c       OPEN(UNIT=44,FILE='polpdf-gridfiles/val2000_lo_g1.grid',
       OPEN(IIREAD,FILE='polpdf-gridfiles/val2000_lo_g1.grid',
     * STATUS='OLD')
      ELSE
        WRITE(6,93)
  93    FORMAT (2X,'PARTON INTERPOLATION: ISET OUT OF RANGE')
        GOTO 60
      END IF
C
       DO 151 M = 1, NX-1
       DO 152 N = 1, NQ
       READ(IIREAD,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1                 PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M),
     2                 PARTON(7,N,M), PARTON(8,N,M)
  90   FORMAT (8(1PE12.4))
 152  CONTINUE
 151  CONTINUE

      close(IIREAD)
C
      IINI = 1
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0 = XB(IX) 
        XB1 = 1.D0-XB(IX)
        XUF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0)
        XDF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0)
        XUBF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**8 * XB0**0.5) 
        XDBF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**8 * XB0**0.5) 
        XSF(IX,IQ)  = PARTON(5,IQ,IX) / (XB1**8 * XB0**0.5) 
        XGF(IX,IQ)  = PARTON(6,IQ,IX) / (XB1**5 * XB0**0.5)
        XG1P(IX,IQ)  = PARTON(7,IQ,IX) / XB1**3
        XG1N(IX,IQ)  = PARTON(8,IQ,IX) / XB1**3
  20  CONTINUE
        XUF(NX,IQ) = 0.D0
        XDF(NX,IQ) = 0.D0
        XUBF(NX,IQ) = 0.D0
        XDBF(NX,IQ) = 0.D0
        XSF(NX,IQ)  = 0.D0
        XGF(NX,IQ)  = 0.D0
        XG1P(NX,IQ)  = 0.D0
        XG1N(NX,IQ)  = 0.D0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
  16  CONTINUE
*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      U = DFINT(NARG,XT,NA,ARRF,XUF) * (1.D0-X)**3 * X
      D = DFINT(NARG,XT,NA,ARRF,XDF) * (1.D0-X)**4 * X
      UB = DFINT(NARG,XT,NA,ARRF,XUBF) * (1.D0-X)**8 * X**0.5
      DB = DFINT(NARG,XT,NA,ARRF,XDBF) * (1.D0-X)**8 * X**0.5
      ST = DFINT(NARG,XT,NA,ARRF,XSF)  * (1.D0-X)**8 * X**0.5
      GL = DFINT(NARG,XT,NA,ARRF,XGF)  * (1.D0-X)**5 * X**0.5
      G1P = DFINT(NARG,XT,NA,ARRF,XG1P)  * (1.D0-X)**3
      G1N = DFINT(NARG,XT,NA,ARRF,XG1N)  * (1.D0-X)**3
 60   RETURN
      END


C=============================================================================
C=============================================================================


**********************************************************************
*                                                                    *
*         POLARIZED NLO QCD PARTON DENSITIES                         *
*                                                                    *
*         E. LEADER, A.V. SIDOROV AND D.B. STAMENOV                  *
*                                                                    *
*         Phys. Rev. D82 (2010) 114018 [arXiv:1010.0574]             *
*                                                                    *
*         PROBLEMS/QUESTIONS TO sidorov@theor.jinr.ru                *
*         OR TO stamenov@inrne.bas.bg                                *
*                                                                    *
*   Two sets of polarized NLO QCD parton densities corresponding to  *
*   positive and sign changing gluon densities are presented in      * 
*   MS-bar scheme.                                                   *
*                                                                    *
*   The sets of PDFs are obtained from a combined NLO QCD analysis   *
*   of the world polarized inclusive and semi-inclusive DIS data.    *
*                                                                    *
*   Heavy quark thresholds Q(H)=M(H) in the BETA function:           *
*              M(c) = 1.43 GeV,   M(b) = 4.3 GeV.                    *
*                                                                    *
*      NLO:  LAMBDA(3) = 0.366,     LAMBDA(4) = 0.311,               *
*            LAMBDA(5) = 0.224                                       *
*   in the BETA function (for details of calculation of strong       *
*   running coupling alpha_s see the paper).                         *
*                                                                    *
*   INPUT:   ISET = number of the parton set                         *
*             (TO BE DEFINED BY THE USER ):                          *
*            ISET = 1   NEXT-TO-LEADING ORDER (xDelta G > 0)         *
*                      (DATA FILE 'NLO_MS_delGpos.grid' UNIT=11)     *
*                                                                    *
*            ISET = 2   NEXT-TO-LEADING ORDER (sign-changing         *
*                       xDelta G)                                    *
*                      (DATA FILE 'NLO_MS_chsign_delG.grid' UNIT=22) *
*                                                                    *
*            X  = Bjorken-x       (between  1.E-5  and  1)           *
*            Q2 = scale in GeV**2 (between  1.0 and 0.58E6)          *
*                                                                    *
*   OUTPUT:  UUB = x *(DELTA u + DELTA ubar)                         *
*            DDB = x *(DELTA d + DELTA dbar)                         *
*            U   = x * DELTA u                                       *
*            D   = x * DELTA d                                       *
*            UB  = x * DELTA ubar                                    *
*            DB  = x * DELTA dbar                                    *
*            ST  = x * DELTA sbar                                    *
*            GL  = x * DELTA GLUON                                   *
*                                                                    *
*          NOTE: The assumption DELTA s = DELTA sbar is used.        *
*                                                                    *
*                                                                    *
*   COMMON:  The main program or the calling routine has to have     *
*            a common block  COMMON / INTINI / IINI , and  IINI      *
*            has always to be zero when LSS2010 is called for the    *
*            first time or when 'ISET' has been changed.             *
*                                                                    *
**********************************************************************

      SUBROUTINE
     1  LSS2010(ISET,Xin,Q2in,UUB,DDB,U,D,UB,DB,ST,GL)


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPART=8, NX=48, NQ=28, NARG=2)
      DIMENSION XUUBF(NX,NQ), XDDBF(NX,NQ),
     1  XUF(NX,NQ), XDF(NX,NQ), XUBF(NX,NQ), XDBF(NX,NQ),
     1  XSF(NX,NQ), XGF(NX,NQ)
     1 ,PARTON (NPART,NQ,NX),
     2  QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(100)
chs     2  QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ)
      COMMON / INTINI / IINI
      CHARACTER  STAR*2
      SAVE XUUBF,XDDBF, XUF, XDF, XUBF, XDBF, XSF, XGF,
     1     NA, ARRF

*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 1.0D0, 1.25D0, 1.5D0, 2.D0, 2.5D0,
     1           4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2           1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3           3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4,
     4           1.0D5, 1.8D5, 3.2D5, 5.8D5 /
       DATA XB /
     1           1.D-5, 1.5D-5, 2.2D-5, 3.2D-5, 4.8D-5, 7.D-5,
     2           1.D-4, 1.5D-4, 2.2D-4, 3.2D-4, 4.8D-4, 7.D-4,
     3           1.D-3, 1.5D-3, 2.2D-3, 3.2D-3, 4.8D-3, 7.D-3,
     4           1.D-2, 1.5D-2, 2.2D-2, 3.2D-2, 5.0D-2, 7.5D-2,
     5           0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     6           0.3, 0.325, 0.35, 0.375, 0.4, 0.45,  0.5, 0.55,
     7           0.6, 0.65,  0.7,  0.75,  0.8, 0.85,  0.9, 1.0 /

C...Added for interface to DJANGOH (TM&HS, 10.01.2013)
      DOUBLE PRECISION Xin, Q2in
      integer ixoutmn,iqoutmn,ixoutmx,iqoutmx
      DOUBLE PRECISION XMIND, XMAXD, Q2MIND, Q2MAXD
      COMMON /HSDPEX/ ixoutmn, ixoutmx, iqoutmn, iqoutmx,
     *                XMIND, XMAXD, Q2MIND, Q2MAXD
      logical FIRST
      data FIRST /.true./

*...CHECK OF X AND Q2 VALUES :
      IF (FIRST) THEN
        ixoutmn=0
        ixoutmx=0
        iqoutmn=0
        iqoutmx=0
        FIRST=.FALSE.
      ENDIF
      X=Xin
      Q2=Q2in
      XMIND=1.0D-5
      XMAXD=1.0D0
      Q2MIND=1D0
      Q2MAXD=0.58D6
      if (Xin.lt.XMIND) then 
        X=XMIND
        ixoutmn=ixoutmn+1
      endif
      if (Xin.gt.XMAXD) then
        X=XMAXD
        ixoutmx=ixoutmx+1
      endif
      if (Q2in.lt.Q2MIND) then 
        Q2=Q2MIND
        iqoutmn=iqoutmn+1
      endif
      if (Q2in.gt.Q2MAXD) then
        Q2=Q2MAXD
        iqoutmx=iqoutmx+1
      endif
*...INITIALIZATION :
*    SELECTION AND READING OF THE GRID :
*                                  ( third NUMBER IN THE grid)
*   INPUT:   ISET = number of the parton set
*            ISET = 1 NEXT-TO-LEADING ORDER (xDelta G > 0)
*   FILE - NO. = 11                                     8.1044E-05
*
*            ISET = 2 NEXT-TO-LEADING ORDER (sign changing xDelta G)
*   FILE - NO. = 22                                     1.3483E-04


      IF (IINI.NE.0) GOTO 16
      IF (ISET.EQ.1) THEN

c       IIREAD=11
       IIREAD=NextUn()
c       OPEN(UNIT=11,FILE='polpdf-gridfiles/NLO_MS_delGpos.grid',
       OPEN(IIREAD,FILE='polpdf-gridfiles/NLO_MS_delGpos.grid',
     * STATUS='OLD')
      ELSE IF (ISET.EQ.2) THEN
c       IIREAD=22
       IIREAD=NextUn()
c       OPEN(UNIT=22,FILE='polpdf-gridfiles/NLO_MS_chsign_delG.grid',
       OPEN(IIREAD,FILE='polpdf-gridfiles/NLO_MS_chsign_delG.grid',
     * STATUS='OLD')

      ELSE
        WRITE(6,93)
  93    FORMAT (2X,'PARTON INTERPOLATION: ISET OUT OF RANGE')
        GOTO 60
      END IF
C

      READ(IIREAD,2004) STAR
 2004 FORMAT (A2)
       DO 151 N = 1, NQ
       DO 152 M = 1, NX

      if(Iset.eq.1.or.Iset.eq.2) then
       READ(IIREAD,190) q2gri, xgri,
     1                 PARTON( 1,N,M), PARTON( 2,N,M), PARTON( 3,N,M),
     1                 PARTON( 4,N,M), PARTON( 5,N,M), PARTON( 6,N,M),
     1                 PARTON( 7,N,M), PARTON( 8,N,M)

 190   FORMAT (2d9.3,8(1pd12.4))

      endif
 152  CONTINUE
 151  CONTINUE

      close(IIREAD)
C
      IINI = 1
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1


        XB0 = XB(IX)
        XB1 = 1.D0-XB(IX)


        XUUBF(iX,iQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0**0.5)
        XDDBF(iX,iQ) = PARTON(2,IQ,IX) / (XB1**3 * XB0**0.5)
        XUF(iX,iQ) =  PARTON(3,IQ,IX) / (XB1**3 * XB0**0.5)

        XDF(IX,IQ) =  PARTON(4,IQ,IX) / (XB1**3 * XB0**0.5)
        XUBF(IX,IQ) =  PARTON(5,IQ,IX) / (XB1**7 * XB0**0.5)
        XDBF(IX,IQ) =  PARTON(6,IQ,IX) / (XB1**7 * XB0**0.5)
        XSF(IX,IQ)  =  PARTON(7,IQ,IX) / (XB1**7 * XB0**0.5)
        XGF(IX,IQ)  =  PARTON(8,IQ,IX) / (XB1**6 * XB0**3.)


chs 2001   FORMAT (2e9.3,13(1pe12.4))
              


  20  CONTINUE
        XUUBF(nX,iQ) =0.d0
        XDDBF(nX,iQ) =0.d0
        XUF(NX,IQ) = 0.D0
        XDF(NX,IQ) = 0.D0
        XUBF(NX,IQ) = 0.D0
        XDBF(NX,IQ) = 0.D0
        XSF(NX,IQ)  = 0.D0
        XGF(NX,IQ)  = 0.D0


  10  CONTINUE
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
  16  CONTINUE
*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      UUB = DFINT(NARG,XT,NA,ARRF,XUUBF)  * (1.D0-X)**3 * X**0.5
      DDB = DFINT(NARG,XT,NA,ARRF,XDDBF)  * (1.D0-X)**3 * X**0.5
      U  = DFINT(NARG,XT,NA,ARRF,XUF)   * (1.D0-X)**3 * X**0.5
      D  = DFINT(NARG,XT,NA,ARRF,XDF)   * (1.D0-X)**3 * X**0.5
      UB  = DFINT(NARG,XT,NA,ARRF,XUBF)   * (1.D0-X)**7 * X**0.5
      DB  = DFINT(NARG,XT,NA,ARRF,XDBF)   * (1.D0-X)**7 * X**0.5
      ST  = DFINT(NARG,XT,NA,ARRF,XSF)    * (1.D0-X)**7 * X**0.5
      GL = DFINT(NARG,XT,NA,ARRF,XGF)    * (1.D0-X)**6 * X**3.

 60   RETURN
      END


C=============================================================================
C=============================================================================


      subroutine polnlo(iflag,x,q2,uval,dval,glue,ubar,dbar,str)
      implicit double precision (a-h,o-z)
      dimension aux(6)

c ---- NLO (MSbar) polarized parton distributions as described in
c ---- T. Gehrmann and W.J. Stirling: "Polarized Parton Distributions 
c ---- of the Nucleon", Phys.Rev. D53 (1996) 6100.
c ----
c ---- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c ---- before calling this subroutine for the first time, the
c ---- LO grids need to be installed by CALL NLOINI 
c ---- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c ----
c ---- INPUT parameters:
c ---- iflag = selects gluon set    0:A, 1:B, 2:C 
c ---- x     = xbjorken
c ---- q2    = Q^2
c ---- OUTPUT parameters (note that always x*distribution is returned)
c ---- uval  = u-ubar
c ---- dval  = d-dbar 
c ---- glue 
c ---- ubar  = 1/2 usea 
c ---- dbar = 1/2 dsea
c ---- str = sbar = 1/2 strsea
c ---- 
c ---- please report any problems to: gehrt@mail.desy.de

      call rdarry(x,q2,aux,iflag)
      uval=aux(1)*(x**0.6d0*(1.d0-x)**3)
      dval=aux(2)*(x**0.75d0*(1.d0-x)**4)
      glue=aux(3)*(x**0.5d0*(1.d0-x)**5)
      ubar=aux(4)*(x**0.5d0*(1.d0-x)**6)
      dbar=aux(5)*(x**0.5d0*(1.d0-x)**6)
      str=aux(6)*(x**0.5d0*(1.d0-x)**6)
      
      return
      end


      subroutine nloini
      implicit double precision (a-h,o-z)
      common/pdist/arraya(151,20,6),arrayb(151,20,6),arrayc(151,20,6)

      lunpl1=NextUn()
      open(lunpl1,file='polpdf-gridfiles/polnloA.dat')
      lunpl2=NextUn()
      open(lunpl2,file='polpdf-gridfiles/polnloB.dat')
      lunpl3=NextUn()
      open(lunpl3,file='polpdf-gridfiles/polnloC.dat')

      do i=1,20
         do j=1,151
            read(lunpl1,901) arraya(j,i,1),arraya(j,i,2),arraya(j,i,3),
     .                       arraya(j,i,4),arraya(j,i,5),arraya(j,i,6)
            read(lunpl2,901) arrayb(j,i,1),arrayb(j,i,2),arrayb(j,i,3),
     .                       arrayb(j,i,4),arrayb(j,i,5),arrayb(j,i,6)
            read(lunpl3,901) arrayc(j,i,1),arrayc(j,i,2),arrayc(j,i,3),
     .                       arrayc(j,i,4),arrayc(j,i,5),arrayc(j,i,6)
         enddo
      enddo
      close(lunpl3)
      close(lunpl2)
      close(lunpl1)

  901 format(6f14.9)

      return
      end


C=============================================================================
C=============================================================================


      subroutine polpar(iflag,x,q2,uval,dval,glue,qbar,str)
      implicit double precision (a-h,o-z)
      dimension aux(5)

c ---- LO polarized parton distributions as described in
c ---- T. Gehrmann and W.J. Stirling: "Polarized Parton Distributions 
c ---- of the Nucleon", Phys.Rev. D53 (1996) 6100.
c ----
c ---- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c ---- before calling this subroutine for the first time, the
c ---- LO grids need to be installed by CALL POLINI 
c ---- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c ----
c ---- INPUT parameters:
c ---- iflag = selects gluon set    0:A, 1:B, 2:C 
c ---- x     = xbjorken
c ---- q2    = Q^2
c ---- OUTPUT parameters (note that always x*distribution is returned)
c ---- uval  = u-ubar
c ---- dval  = d-dbar 
c ---- glue 
c ---- qbar  = ubar = dbar = 1/2 usea = 1/2 dsea
c ---- str   = sbar = 1/2 strsea
c ---- 
c ---- please report any problems to: gehrt@mail.desy.de

      call rdarry(x,q2,aux,iflag)
      uval=aux(1)*(x**0.6d0*(1.d0-x)**3)
      dval=aux(2)*(x**0.75d0*(1.d0-x)**4)
      glue=aux(3)*(x**0.5d0*(1.d0-x)**5)
      qbar=aux(4)*(x**0.5d0*(1.d0-x)**6)
      str=aux(5)*(x**0.5d0*(1.d0-x)**6)
            
      return
      end


      subroutine rdarry(x,q2,aux,iflag)
      implicit double precision (a-h,o-z)
      implicit integer(i-n)
      dimension aux(5)
      common/pdist/arraya(151,20,6),arrayb(151,20,6),arrayc(151,20,6)
chs      common/pdist/arraya(151,20,5),arrayb(151,20,5),arrayc(151,20,5)
     
C...Added for interface to DJANGOH (TM&HS, 10.01.2013)
      integer ixoutmn,iqoutmn,ixoutmx,iqoutmx
      DOUBLE PRECISION XMIND, XMAXD, Q2MIND, Q2MAXD
      COMMON /HSDPEX/ ixoutmn, ixoutmx, iqoutmn, iqoutmx,
     *                XMIND, XMAXD, Q2MIND, Q2MAXD
      logical FIRST
      data FIRST /.true./

      nx=151
      ndata=nx-1
      nq2pts=20
      nq2inv=nq2pts-1
      q2sta=1.d0
      q2fin=1.d6
      ymin=5.d0
      XMIND=10.d0**(-ymin)
      XMAXD=1.d0
      Q2MIND=q2sta
      Q2MAXD=q2fin

      IF (FIRST) THEN
        ixoutmn=0
        ixoutmx=0
        iqoutmn=0
        iqoutmx=0
        FIRST=.FALSE.
      ENDIF
      if (q2.lt.q2sta) then
        iqoutmn=iqoutmn+1
C         write (6,*) '*** Q^2-value out of range ***'
C         write (6,*) '*** Q^2 set to minimal value ***'
        q2=q2sta
C         write (6,*) q2
      endif
      if (q2.gt.q2fin) then
        iqoutmx=iqoutmx+1
C         write (6,*) '*** Q^2-value out of range ***'
C         write (6,*) '*** Q^2 set to maximal value ***'
        q2=q2fin
C         write (6,*) q2
      endif
      if (x.lt.XMIND) then
        ixoutmn=ixoutmn+1
C         write (6,*) '*** x-value out of range ***'
C         write (6,*) '*** x set to minimal value ***'
        x=XMIND
C         write (6,*) x
      endif
      if (x.gt.XMAXD) then
        ixoutmx=ixoutmx+1
C         write (6,*) '*** x-value out of range ***'
C         write (6,*) '*** x set to minimal value ***'
        x=XMAXD
C         write (6,*) x
      endif

      y=dlog10(x)
      ram=(y+ymin)*ndata/ymin+1.d0
      iram=int(ram)
      fraci=ram-dble(iram)
      ram=dlog(q2/q2sta)*nq2inv/dlog(q2fin/q2sta)+1.d0
      jram=int(ram)
      fracj=ram-dble(jram)
      
      if (iflag.eq.0) then
         do 100 i=1,5
            aux(i)=(arraya(iram,jram,i)*(1.d0-fraci)
     .             +arraya(iram+1,jram,i)*fraci)*(1.d0-fracj)+
     .             (arraya(iram,jram+1,i)*(1.d0-fraci)
     .             +arraya(iram+1,jram+1,i)*fraci)*fracj
  100    continue
      else
      if (iflag.eq.1) then
         do 200 i=1,5
            aux(i)=(arrayb(iram,jram,i)*(1.d0-fraci)
     .             +arrayb(iram+1,jram,i)*fraci)*(1.d0-fracj)+
     .             (arrayb(iram,jram+1,i)*(1.d0-fraci)
     .             +arrayb(iram+1,jram+1,i)*fraci)*fracj
  200    continue
      else
         do 300 i=1,5
           aux(i)=(arrayc(iram,jram,i)*(1.d0-fraci)
     .            +arrayc(iram+1,jram,i)*fraci)*(1.d0-fracj)+
     .            (arrayc(iram,jram+1,i)*(1.d0-fraci)
     .            +arrayc(iram+1,jram+1,i)*fraci)*fracj
  300    continue
      endif
      endif

      return
      end


      subroutine polini
      implicit double precision (a-h,o-z)
      common/pdist/arraya(151,20,6),arrayb(151,20,6),arrayc(151,20,6)
chs      common/pdist/arraya(151,20,5),arrayb(151,20,5),arrayc(151,20,5)

      lunpl1=NextUn()
      open(lunpl1,file='polpdf-gridfiles/polparA.dat')
      lunpl2=NextUn()
      open(lunpl2,file='polpdf-gridfiles/polparB.dat')
      lunpl3=NextUn()
      open(lunpl3,file='polpdf-gridfiles/polparC.dat')

      do i=1,20
         do j=1,151
            read(lunpl1,901) arraya(j,i,1),arraya(j,i,2),
     .                   arraya(j,i,3),arraya(j,i,4),arraya(j,i,5)
            read(lunpl2,901) arrayb(j,i,1),arrayb(j,i,2),arrayb(j,i,3),
     .                   arrayb(j,i,4),arrayb(j,i,5)
            read(lunpl3,901) arrayc(j,i,1),arrayc(j,i,2),arrayc(j,i,3),
     .                   arrayc(j,i,4),arrayc(j,i,5)
         enddo
      enddo
      close(lunpl3)
      close(lunpl2)
      close(lunpl1)

  901 format(5f14.9)

      return
      end


C=============================================================================
C=============================================================================


*********************************************************************
*                                                                   *
*        POLARIZED PARTON DISTRIBUTION FUNCTION WITH ERRORS         *
*                        FOR LO AND NLO                             *
*                                                                   *
*            Johannes Bluemlein and Helmut Boettcher                *
*                        hep-ph/0203155                             *
*                                                                   *
*   PROBLEMS/QUESTIONS TO blumlein@ifh.de OR hboett@ifh.de          *
*                                                                   *
*   INPUT:   iset = number of the parton set :                      *
*                                                                   *
*            iset = 1  LEADING ORDER - Scenario 1                   *
*                      data file: bb01_xpdf_lo_8pN.grid             * 
*                                                                   *
*            iset = 2  LEADING ORDER - Scenario 2                   *
*                      data file: bb01_xpdf_lo_8pS.grid             * 
*                                                                   *
*            iset = 3  NEXT-TO-LEADING ORDER (MS-bar) - Scenario 1  *
*                      data file: bb01_xpdf_nlo_8pN.grid            * 
*                                                                   *
*            iset = 4  NEXT-TO-LEADING ORDER (MS-bar) - Scenario 2  *
*                      data file: bb01_xpdf_nlo_8pS.grid            * 
*                                                                   *
*            x  = Bjorken-x       (between  1.0E-9 and 1.0E+0)      *
*            q2 = scale in GeV**2 (between  0.1E+1 and 1.0E+6)      *
*                                                                   *
*   OUTPUT:  UV = DELTA u_v    & error: DUV                         *
*            DV = DELTA d_v    & error: DDV                         *
*            GL = DELTA gluon  & error: DGL                         *
*            QB = DELTA qbar   & error: DQB                         *   
*            G1P = g_1^proton  & error: DG1P                        *
*            G1N = g_1^neutron & error: DG1N                        * 
*                                                                   *
*            NOTE:                                                  *
*            -----                                                  *
*            For the parton distributions and for g1(p,n)           *
*            always x times the distribution is returned.           *
*                                                                   *
*            The sets are the result of a combined fit to the       *
*            world data for the spin asymmetries, i.e. A_1(p,n,d)   *
*            or g_1/F_1(p,n,d).                                     *
*                                                                   *
*            The subroutine PPDF returns the BB polarized parton    *
*            distribution function values and g1(p,n) at the given  *
*            point in Q**2 and XB by interpolating the data on the  *
*            specified grid.                                        *
*                                                                   *
*            Note: No charm is included                             *
*                                                                   *
*   COMMON:  The main program or the calling routine has to have    *
*            a common block  COMMON / INTINI / IINI , and  IINI     *
*            has always to be ZERO when POLPDF is called for the    *
*            first time or when 'ISET' has been changed.            *
*                                                                   *
*********************************************************************
*
      SUBROUTINE PPDF(ISET, Xin, Q2in, UV, DUV, DV, DDV, GL, DGL,
     1                QB, DQB, G1P, DG1P, G1N, DG1N)
*     ---------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*     
      PARAMETER (NPDF=6, NQ=51, NX=104 )
*
      CHARACTER  STAR*2
*
      COMMON / INTINI / IINI
*
      DIMENSION QS(NQ),XB(NX),XPDF(NX,NQ,NPDF),XDPDF(NX,NQ,NPDF)
      DIMENSION AXB(NX,NQ,NPDF),BXB(NX,NQ,NPDF),CXB(NX,NQ,NPDF)
      DIMENSION AXBE(NX,NQ,NPDF),BXBE(NX,NQ,NPDF),CXBE(NX,NQ,NPDF)
*
      DIMENSION PDF1(NPDF),PDF2(NPDF),DPDF1(NPDF),DPDF2(NPDF)
      DIMENSION PDF(NPDF),DPDF(NPDF)
*
      SAVE QS, XB, XPDF, XDPDF
      SAVE AXB, BXB, CXB, AXBE, BXBE, CXBE

C...Added for interface to DJANGOH (TM&HS, 10.01.2013)
      DOUBLE PRECISION Xin, Q2in
      integer ixoutmn,iqoutmn,ixoutmx,iqoutmx
      DOUBLE PRECISION XMIND, XMAXD, Q2MIND, Q2MAXD
      COMMON /HSDPEX/ ixoutmn, ixoutmx, iqoutmn, iqoutmx,
     *                XMIND, XMAXD, Q2MIND, Q2MAXD
      logical FIRST
      data FIRST /.true./

*
*     BJORKEN-X AND Q**2 VALUES OF THE GRID WHICH WILL BE READ IN 
*     TOGETHER WITH THE OTHER DATA :
*
c      DATA QS /
c     1 0.10000D+01, 0.13183D+01, 0.17378D+01, 0.22909D+01,0.30200D+01,
c     2 0.39811D+01, 0.52481D+01, 0.69183D+01, 0.91201D+01,0.12023D+02,
c     3 0.15849D+02, 0.20893D+02, 0.27542D+02, 0.36308D+02,0.47863D+02,
c     4 0.63096D+02, 0.83176D+02, 0.10965D+03, 0.14454D+03,0.19055D+03,
c     5 0.25119D+03, 0.33113D+03, 0.43652D+03, 0.57544D+03,0.75858D+03,
c     6 0.10000D+04, 0.13183D+04, 0.17378D+04, 0.22909D+04,0.30200D+04,
c     7 0.39811D+04, 0.52481D+04, 0.69183D+04, 0.91201D+04,0.12023D+05,
c     8 0.15849D+05, 0.20893D+05, 0.27542D+05, 0.36308D+05,0.47863D+05,
c     9 0.63096D+05, 0.83176D+05, 0.10965D+06, 0.14454D+06,0.19055D+06,
c     A 0.25119D+06, 0.33113D+06, 0.43652D+06, 0.57544D+06,0.75858D+06,
c     B 0.10000D+07 /
c       DATA XB /
c     1 0.10000D-08, 0.13183D-08, 0.17378D-08, 0.22909D-08, 0.30200D-08,
c     2 0.39811D-08, 0.52481D-08, 0.69183D-08, 0.91201D-08, 0.12023D-07,
c     3 0.15849D-07, 0.20893D-07, 0.27542D-07, 0.36308D-07, 0.47863D-07,
c     4 0.63096D-07, 0.83176D-07, 0.10965D-06, 0.14454D-06, 0.19055D-06,
c     5 0.25119D-06, 0.33113D-06, 0.43652D-06, 0.57544D-06, 0.75858D-06,
c     1 0.10000D-05, 0.13183D-05, 0.17378D-05, 0.22909D-05, 0.30200D-05,
c     2 0.39811D-05, 0.52481D-05, 0.69183D-05, 0.91201D-05, 0.12023D-04,
c     3 0.15849D-04, 0.20893D-04, 0.27542D-04, 0.36308D-04, 0.47863D-04,
c     4 0.63096D-04, 0.83176D-04, 0.10965D-03, 0.14454D-03, 0.19055D-03,
c     5 0.25119D-03, 0.33113D-03, 0.43652D-03, 0.57544D-03, 0.75858D-03,
c     1 0.10000D-02, 0.13183D-02, 0.17378D-02, 0.22909D-02, 0.30200D-02,
c     2 0.39811D-02, 0.52481D-02, 0.69183D-02, 0.91201D-02, 0.12023D-01,
c     3 0.15849D-01, 0.20893D-01, 0.27542D-01, 0.36308D-01, 0.47863D-01,
c     4 0.63096D-01, 0.83176D-01,
c     6 0.10000D+00, 0.12500D+00, 0.15000D+00, 0.17500D+00, 0.20000D+00,
c     7 0.22500D+00, 0.25000D+00, 0.27500D+00, 0.30000D+00, 0.32500D+00,
c     8 0.35000D+00, 0.37500D+00, 0.40000D+00, 0.42500D+00, 0.45000D+00,
c     9 0.47500D+00, 0.50000D+00, 0.52500D+00, 0.55000D+00, 0.57500D+00,
c     A 0.60000D+00, 0.62500D+00, 0.65000D+00, 0.67500D+00, 0.70000D+00,
c     B 0.72500D+00, 0.75000D+00, 0.77500D+00, 0.80000D+00, 0.82500D+00,
c     C 0.85000D+00, 0.87500D+00, 0.90000D+00, 0.92500D+00, 0.95000D+00,
c     D 0.97500D+00, 0.10000D+01 /
*
*     CHECK OF X AND Q2 VALUES : 
      IF (FIRST) THEN
        ixoutmn=0
        ixoutmx=0
        iqoutmn=0
        iqoutmx=0
        FIRST=.FALSE.
      ENDIF
*
      X=Xin
      Q2=Q2in
      XMIND=1.0D-9
      XMAXD=1.0D0
      Q2MIND=1.0D0
      Q2MAXD=1.D6
      if (Xin.lt.XMIND) then 
        X=XMIND
        ixoutmn=ixoutmn+1
      endif
      if (Xin.gt.XMAXD) then
        X=XMAXD
        ixoutmx=ixoutmx+1
      endif
      if (Q2in.lt.Q2MIND) then 
        Q2=Q2MIND
        iqoutmn=iqoutmn+1
      endif
      if (Q2in.gt.Q2MAXD) then
        Q2=Q2MAXD
        iqoutmx=iqoutmx+1
      endif
*
      IF (IINI .NE. 0) GOTO 10
*
*     READING OF THE DATA SET :
*     
c     LIN = 10
      LIN=NextUn()
*
      IF (ISET .EQ. 1) THEN
*
*     open input unit for LO data set N (Scenario 1)
*
         OPEN(UNIT=LIN,FILE='polpdf-gridfiles//bb01_xpdf_lo_8pN.grid',
     1              STATUS='UNKNOWN')
*
      ELSEIF (ISET .EQ. 2) THEN
*
*     open input unit for LO data set S (Scenario 2)
*
         OPEN(UNIT=LIN,FILE='polpdf-gridfiles//bb01_xpdf_lo_8pS.grid',
     1              STATUS='UNKNOWN')
*
      ELSEIF (ISET .EQ. 3) THEN
*
*     open input unit for NLO data set N (Scenario 1)
*
         OPEN(UNIT=LIN,FILE='polpdf-gridfiles//bb01_xpdf_nlo_8pN.grid',
     1              STATUS='UNKNOWN')
*
      ELSEIF (ISET .EQ. 4) THEN
*
*     open input unit for NLO data set S (Scenario 2)
*
         OPEN(UNIT=LIN,FILE='polpdf-gridfiles//bb01_xpdf_nlo_8pS.grid',
     1              STATUS='UNKNOWN')
*
      ELSE
*
         WRITE(6,2003)
 2003    FORMAT (2X,'PARTON INTERPOLATION: ISET OUT OF RANGE -- '
     1        //'STOP')
         STOP
*
      ENDIF
*
      READ(LIN,2004) STAR
 2004 FORMAT (A2)
      DO IQ = 1,NQ
         DO IX = 1,NX
            READ(LIN,2005) QS(IQ),XB(IX),
     1           (XPDF(IX,IQ,IPDF),XDPDF(IX,IQ,IPDF),IPDF=1,NPDF)
 2005       FORMAT (14(1P,E13.5))
         ENDDO
      ENDDO
*     
      CLOSE(LIN)
*
*     CALCULATE SPLINE COEFFIFIENTS FOR THE INTERPOLATION IN X :
*
      DO IPDF = 1,NPDF
         DO IQ = 1,NQ
            CALL SPLINE(NX,XB,XPDF,AXB,BXB,CXB,IPDF,IQ)
            CALL SPLINE(NX,XB,XDPDF,AXBE,BXBE,CXBE,IPDF,IQ)
         ENDDO
      ENDDO
*
      IINI = 1
*
 10   CONTINUE
*
*     INTERPOLATION :
*     X: CUBIC SPLINE INTERPOLATION, LOG(Q2): LINEAR INTERPOLATION
*
chs      IQ = ISEARCH(NQ,QS,Q2)
      IQ = ISRCH2Q(NQ,QS,Q2)
      IF (IQ .EQ. NQ) IQ = NQ-1
chs      IX = ISEARCH(NX,XB,X)
      IX = ISRCH2X(NX,XB,X)
*
      DX = X - XB(IX)
*
      DO IPDF = 1,NPDF
         PDF1(IPDF) = XPDF(IX,IQ,IPDF)
     1        + DX*(AXB(IX,IQ,IPDF) + DX*(BXB(IX,IQ,IPDF) 
     2        + DX*CXB(IX,IQ,IPDF)))
         PDF2(IPDF) = XPDF(IX,IQ+1,IPDF)
     1        + DX*(AXB(IX,IQ+1,IPDF) + DX*(BXB(IX,IQ+1,IPDF) 
     2        + DX*CXB(IX,IQ+1,IPDF)))
         DPDF1(IPDF) = XDPDF(IX,IQ,IPDF)
     1        + DX*(AXBE(IX,IQ,IPDF) + DX*(BXBE(IX,IQ,IPDF) 
     2        + DX*CXBE(IX,IQ,IPDF)))
         DPDF2(IPDF) = XDPDF(IX,IQ+1,IPDF)
     1        + DX*(AXBE(IX,IQ+1,IPDF) + DX*(BXBE(IX,IQ+1,IPDF) 
     2        + DX*CXBE(IX,IQ+1,IPDF)))
      ENDDO
*
      TQ = (DLOG(Q2)-DLOG(QS(IQ))) / (DLOG(QS(IQ+1))-DLOG(QS(IQ)))
*
      DO IPDF = 1,NPDF
         PDF(IPDF)  = (1.0D0-TQ)*PDF1(IPDF) + TQ*PDF2(IPDF)
         DPDF(IPDF) = (1.0D0-TQ)*DPDF1(IPDF) + TQ*DPDF2(IPDF)
      ENDDO
*
      UV    = PDF(1)
      DUV   = DPDF(1)
      DV    = PDF(2)
      DDV   = DPDF(2)
      GL    = PDF(3)
      DGL   = DPDF(3)
      QB    = PDF(4) / 6.0d0
      DQB   = DPDF(4) / 6.0d0
      G1P   = PDF(5)
      DG1P  = DPDF(5)
      G1N   = PDF(6)
      DG1N  = DPDF(6)
*
      RETURN
      END
*
* ---------------------------------------------------------------------
      SUBROUTINE SPLINE(N,X,Y,B,C,D,I,J)
* ---------------------------------------------------------------------
* CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.
* INTERPOLATION SUBROUTINES ARE TAKEN FROM
* G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,
* COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977).
*
* SUBROUTINE TAKEN FROM AAC GROUP (KUMANO et al.)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER (NPDF=6, NQ=51, NX=104 )
*
      DIMENSION X(NX), Y(NX,NQ,NPDF),
     1          B(NX,NQ,NPDF), C(NX,NQ,NPDF), D(NX,NQ,NPDF)
*
      NM1=N-1
      IF(N.LT.2) RETURN
      IF(N.LT.3) GO TO 250
      D(1,J,I)=X(2)-X(1)
      C(2,J,I)=(Y(2,J,I)-Y(1,J,I))/D(1,J,I)
      DO 210 K=2,NM1
         D(K,J,I)=X(K+1)-X(K)
         B(K,J,I)=2.0D0*(D(K-1,J,I)+D(K,J,I))
         C(K+1,J,I)=(Y(K+1,J,I)-Y(K,J,I))/D(K,J,I)
         C(K,J,I)=C(K+1,J,I)-C(K,J,I)
  210 CONTINUE
      B(1,J,I)=-D(1,J,I)
      B(N,J,I)=-D(N-1,J,I)
      C(1,J,I)=0.0D0
      C(N,J,I)=0.0D0
      IF(N.EQ.3) GO TO 215
      C(1,J,I)=C(3,J,I)/(X(4)-X(2))-C(2,J,I)/(X(3)-X(1))
      C(N,J,I)=C(N-1,J,I)/(X(N)-X(N-2))-C(N-2,J,I)/(X(N-1)-X(N-3))
      C(1,J,I)=C(1,J,I)*D(1,J,I)**2.0D0/(X(4)-X(1))
      C(N,J,I)=-C(N,J,I)*D(N-1,J,I)**2.0D0/(X(N)-X(N-3))
 215  CONTINUE
      DO 220 K=2,N
         T=D(K-1,J,I)/B(K-1,J,I)
         B(K,J,I)=B(K,J,I)-T*D(K-1,J,I)
         C(K,J,I)=C(K,J,I)-T*C(K-1,J,I)
 220  CONTINUE
      C(N,J,I)=C(N,J,I)/B(N,J,I)
      DO 230 IB=1,NM1
         K=N-IB
         C(K,J,I)=(C(K,J,I)-D(K,J,I)*C(K+1,J,I))/B(K,J,I)
 230  CONTINUE
      B(N,J,I)=(Y(N,J,I)-Y(NM1,J,I))/D(NM1,J,I)
     1     +D(NM1,J,I)*(C(NM1,J,I)+2.0D0*C(N,J,I))
      DO 240 K=1,NM1
         B(K,J,I)=(Y(K+1,J,I)-Y(K,J,I))/D(K,J,I)
     1        -D(K,J,I)*(C(K+1,J,I)+2.0D0*C(K,J,I))
         D(K,J,I)=(C(K+1,J,I)-C(K,J,I))/D(K,J,I)
         C(K,J,I)=3.0D0*C(K,J,I)
 240  CONTINUE
      C(N,J,I)=3.0D0*C(N,J,I)
      D(N,J,I)=D(N-1,J,I)
      RETURN
 250  CONTINUE
      B(1,J,I)=(Y(2,J,I)-Y(1,J,I))/(X(2)-X(1))
      C(1,J,I)=0.0D0
      D(1,J,I)=0.0D0
      B(2,J,I)=B(1,J,I)
      C(2,J,I)=0.0D0
      D(2,J,I)=0.0D0
      RETURN
      END
*
* ---------------------------------------------------------------------
      INTEGER FUNCTION ISRCH2Q(N,X,Y)
* ---------------------------------------------------------------------
* THIS FUNCTION SEARCHES "I" WHICH SATISFIES THE RELATION
* X(I) <= Y < X(I+1) BY USING A BINARY SEARCH.
*
* FUNCTION TAKEN FROM AAC GROUP (KUMANO et al.)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
chs      PARAMETER (NPDF=6, NQ=51, NX=104 )
*
chs      DIMENSION X(NX)
      DIMENSION X(51)
*
      MIN=1
      MAX=N+1
*
   10 CONTINUE
      MID=(MIN+MAX)/2
      IF(Y.LT.X(MID)) THEN
        MAX=MID
      ELSE
        MIN=MID
      END IF
      IF((MAX-MIN).GT.1) GO TO 10
*
chs      ISEARCH=MIN
      ISRCH2Q=MIN
*
      RETURN
      END
*
* *********************************************************************
*
* ---------------------------------------------------------------------
      INTEGER FUNCTION ISRCH2X(N,X,Y)
* ---------------------------------------------------------------------
* THIS FUNCTION SEARCHES "I" WHICH SATISFIES THE RELATION
* X(I) <= Y < X(I+1) BY USING A BINARY SEARCH.
*
* FUNCTION TAKEN FROM AAC GROUP (KUMANO et al.)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
chs      PARAMETER (NPDF=6, NQ=51, NX=104 )
*
chs      DIMENSION X(NX)
      DIMENSION X(104)
*
      MIN=1
      MAX=N+1
*
   10 CONTINUE
      MID=(MIN+MAX)/2
      IF(Y.LT.X(MID)) THEN
        MAX=MID
      ELSE
        MIN=MID
      END IF
      IF((MAX-MIN).GT.1) GO TO 10
*
chs      ISEARCH=MIN
      ISRCH2X=MIN
*
      RETURN
      END
*
* *********************************************************************



C=============================================================================
C=============================================================================


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                            C
C       POLARIZED PARTON DISTRIBUTIONS (NLO)                                 C
C       FROM  PRD71, 094018 (2005)                                           C
C       hep-ph/0504155                                                       C
C       D. de Florian G.A. Navarro and R. Sassot.                            C
C                                                                            C
C       MODE: 1 SET 1 NLO (MSbar) (based on KRETZER ffs)                     C
C             2 SET 2 NLO (MSbar) (based on KKP     ffs)                     C
C             3 SET 1  LO         (based on KRETZER ffs)                     C
C             4 SET 2  LO         (based on KKP     ffs)                     C
C                                                                            C
C                    Q2=Q^2                                                  C
C                    DUV :    X * U VALENCE DISTRIBUTION                     C
C                    DDV :    X * D VALENCE DISTRIBUTION                     C
C                    DUBAR :  X * UBAR DISTRIBUTION                          C
C                    DDBAR :  X * DBAR DISTRIBUTION                          C
C                    DSTR :   X * STRANGE DISTRIBUTION                       C
C                    DGLU :   X * GLUON DISTRIBUTION                         C
C                    G1P :    X * POLARIZED STRUCTURE FUNCTION (PROTON)      C
C                    G1N :    X * POLARIZED STRUCTURE FUNCTION(NEUTRON)      C
C                                                                            C
C       REMEMBER: ALWAYS X*DISTRIBUTION !!!                                  C
C       BEFORE CALLING THE SUBRUTINE `POLFIT` FOR THE FIRST TIME, THE        C
C       SUBROUTINE `INIDNS` MUST BE CALLED (ONLY ONCE) TO READ THE GRIDS.    C
C              (CALL INIDNS)                                                 C
C       RANGE OF VALIDITY OF THE INTERPOLATION:                              C
C       10**(-4)< X < 0.9                                                    C
C       1 < Q**2 < 5*10**4                                                   C
C                                                                            C
C       IN CASE OF PROBLEMS, DOUBTS, ETC, PLEASE REPORT TO                   C
C        dflorian@df.uba.ar                                                  C
C        sassot@df.uba.ar                                                    C
C                                                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       SUBROUTINE POLFIT(MODE,Xin,Q2in,DUV,DDV,DUBAR,DDBAR,DSTR,DGLU,
     #  G1P,G1N)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       	DIMENSION GU(6,26,76),GD(6,26,76),GUB(6,26,76),GDB(6,26,76),
     #            GS(6,26,76), GG(6,26,76), GP(6,26,76), GN(6,26,76)
        DIMENSION DUVG(76,26), DDVG(76,26), DUBARG(76,26),DDBARG(76,26),        
     #            DSTRG(76,26),DGLUG(76,26),GPG(76,26),GNG(76,26),XQ(2)
        COMMON/ GRIDDNS / GU,GD,GUB,GDB,GS,GG,GP,GN         

C...Added for interface to DJANGOH (TM&HS, 10.01.2013)
      DOUBLE PRECISION Xin, Q2in
      integer ixoutmn,iqoutmn,ixoutmx,iqoutmx
      DOUBLE PRECISION XMIND, XMAXD, Q2MIND, Q2MAXD
      COMMON /HSDPEX/ ixoutmn, ixoutmx, iqoutmn, iqoutmx,
     *                XMIND, XMAXD, Q2MIND, Q2MAXD
      logical FIRST
      data FIRST /.true./

        DO  K = 1, 26
        DO  J = 1, 76
        DUVG(J,K) = GU(MODE,K,J)  
        DDVG(J,K) = GD(MODE,K,J)  
        DUBARG(J,K) = GUB(MODE,K,J)   
        DDBARG(J,K) = GDB(MODE,K,J)   
        DSTRG(J,K)  = GS(MODE,K,J) 
        DGLUG(J,K)  = GG(MODE,K,J)  
        GPG(J,K)  = GP(MODE,K,J)  
        GNG(J,K)  = GN(MODE,K,J) 
       	END DO
        END DO 

C CHECK X AND Q2 VALUES.
      IF (FIRST) THEN
        ixoutmn=0
        ixoutmx=0
        iqoutmn=0
        iqoutmx=0
        FIRST=.FALSE.
      ENDIF
      X=Xin
      Q2=Q2in
      XMIND=1.0D-4
      XMAXD=0.9D0
      Q2MIND=1.0D1
      Q2MAXD=5.0D4
      if (Xin.lt.XMIND) then 
        X=XMIND
        ixoutmn=ixoutmn+1
      endif
      if (Xin.gt.XMAXD) then
        X=XMAXD
        ixoutmx=ixoutmx+1
      endif
      if (Q2in.lt.Q2MIND) then 
        Q2=Q2MIND
        iqoutmn=iqoutmn+1
      endif
      if (Q2in.gt.Q2MAXD) then
        Q2=Q2MAXD
        iqoutmx=iqoutmx+1
      endif

  
        XQ(1) = DLOG(X)
        XQ(2) = DLOG(Q2)
        X3=(1.D0-X)**3.D0
        X4=(1.D0-X)**4.D0 
        X5=X**0.5D0
        X6=X**0.6D0
        X7=X**0.7D0
        DUV = PERINOLA(XQ,DUVG) * X3* X6*X
        DDV = PERINOLA(XQ,DDVG) * X4 * X7*X 
        DUBAR = PERINOLA(XQ,DUBARG) * X3 * X5*X
        DDBAR = PERINOLA(XQ,DDBARG) * X3 * X5*X
        DSTR = PERINOLA(XQ,DSTRG)  * X3 * X5*X
        DGLU = PERINOLA(XQ,DGLUG)  * X3 * X5*X
        G1P = PERINOLA(XQ,GPG)  * X3 * X5/X
        G1N = PERINOLA(XQ,GNG)  * X4 * X5/X
c        G1D=(G1P+G1N)*0.5D0*(1.D0-1.5D0*0.058D0)
        RETURN
        END


      SUBROUTINE INIDNS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XARRAY(102), GU(6,26,76), GD(6,26,76), GUB(6,26,76),
     # GDB(6,26,76),GS(6,26,76),GG(6,26,76),GP(6,26,76),GN(6,26,76) 
      COMMON/ XARRAY / XARRAY
      COMMON/ GRIDDNS / GU,GD,GUB,GDB,GS,GG,GP,GN
      integer lunplx
      dimension lunplx(4)

c       OPEN(UNIT=10,FILE='polpdf-gridfiles/DNS1NLO.GRID',STATUS='OLD')
c       OPEN(UNIT=11,FILE='polpdf-gridfiles/DNS2NLO.GRID',STATUS='OLD')
c       OPEN(UNIT=12,FILE='polpdf-gridfiles/DNS1LO.GRID',STATUS='OLD')
c       OPEN(UNIT=13,FILE='polpdf-gridfiles/DNS2LO.GRID',STATUS='OLD')
      lunplx(1)=NextUn()
      OPEN(lunplx(1),FILE='polpdf-gridfiles/DNS1NLO.GRID',STATUS='OLD')
      lunplx(2)=NextUn()
      OPEN(lunplx(2),FILE='polpdf-gridfiles/DNS2NLO.GRID',STATUS='OLD')
      lunplx(3)=NextUn()
      OPEN(lunplx(3),FILE='polpdf-gridfiles/DNS1LO.GRID',STATUS='OLD')
      lunplx(4)=NextUn()
      OPEN(lunplx(4),FILE='polpdf-gridfiles/DNS2LO.GRID',STATUS='OLD')

c      	DO IFILE=10,13
        DO if=1,4
        IFILE=lunplx(if)
        DO  K = 1, 76 
        DO  J = 1, 26
c        I=IFILE-9
        I=if
        READ(IFILE,40) GU(I,J,K), GD(I,J,K), GUB(I,J,K),GDB(I,J,K),  
     #                 GS(I,J,K), GG(I,J,K), GP(I,J,K), GN(I,J,K)
        END DO
        END DO
        CLOSE(IFILE)
        END DO
  40   FORMAT (8(1PE15.7))

        DO LX1=0,25,1
        XX=10.D0**(-4.D0+LX1/25.D0*3.D0)
        XARRAY(LX1+1) = DLOG(XX)
        END DO
        DO LX1=1,50,1
        XX=0.1D0+LX1/50.D0*0.9D0
        XARRAY(LX1+26) = DLOG(XX)
        END DO
        DO LQ=0,25,1
        QQ=0.6D0*10**(5.D0*LQ/25.D0)
        XARRAY(LQ+77) = DLOG(QQ)
        END DO
 
        RETURN
        END 


