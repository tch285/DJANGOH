Ch 08/08/05
Ch mod 05/06/21
C
      SUBROUTINE HSUSER(ICALL,X,Y,Q2)
C...User analysis routine:
C   ICALL=1 initialization before event generation
C        =2 analysis for each generated event 
C        =3 final call after completed event generation
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
C...Declarations for heracles
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSOPTN/ INT2(5),INT3(15),ISAM2(5),ISAM3(15),
     *                IOPLOT,IPRINT,ICUT
      COMMON /HSNUME/ SIGTOT,SIGTRR,SIGG(20),SIGGRR(20),NEVENT,NEVE(20)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSIKP/  S,T,U,SS,TS,US,DKP,DKPS,DKQ,DKQS
      COMMON /HSLABP/ EH,PH,EQH,PQH,ESH,PSH,COSEH,SINEH
      PARAMETER (NMXHEP=2000)
      COMMON /HEPEVT/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &                JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),
     &                PHEP(5,NMXHEP),VHKK(4,NMXHEP)
      COMMON /HSCHNN/ ICHNN
      COMMON /HSONLY/ IHSONL
      COMMON /HSPASS/ NREJCW
      CHARACTER OUTFILENAM*80
      COMMON /HSOUTF/ OUTFILENAM,ICH
      REAL RTIME
C...Declarations for lepto65
      COMMON /DJPASS/ NTOT,NPASS,NQELAS,NFAILL,NFAILQ
      COMMON /DJFAIL/ NFAILI(10)
      COMMON/LEPTOU/CUT(14),LST(40),PARL(30),XSCH,YSCH,W2SCH,Q2SCH,USCH
      REAL          CUT            ,PARL    ,XSCH,YSCH,W2SCH,Q2SCH,USCH
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      REAL                      P        ,V
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      REAL                    PARU               ,PARJ
      INTEGER MSTU,MSTJ
      REAL (8) NUHAD,NU
C...Number of events from sophia
      COMMON /SPPASS/ NSOPH,NSPOUT,NFAILP,NSPACC
C...Structure  functions and asymmetries
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
      COMMON /HSFGNC/ F1NC,F3NC,G1NC,G3NC
      COMMON /HSFGCC/ F1CC,F3CC,G1CC,G5CC
      COMMON /HSASYM/ A1NC
C
      DIMENSION IFLCNT(-6:6)
      DIMENSION NFAILC(10)
      DIMENSION NMIS(0:12)
      DIMENSION PSUM(4)
      LOGICAL LFIRST
      DATA LFIRST /.TRUE./
ctest      DATA NEVMOD/1000/
      DATA NEVMOD/100/
      SAVE IFLCNT
C
      IF(LFIRST) THEN
        LFIRST=.FALSE.
        IEIN=1
        DO 10 I=-6,6
          IFLCNT(I)=0
 10     CONTINUE
        DO 11 I=1,10
          NFAILC(I)=0
 11     CONTINUE
        DO 12 I=0,12
          NMIS(I)=0
 12     CONTINUE
        CALL TIMEX(RTIME)
        TIMINI=RTIME
        N10CNT=0
        GSP=SP-MPRO2-MEI2
C        LUNEVT=NextUn()
C        OPEN(LUNEVT,FILE=OUTFILENAM(1:ICH)//'_evt.dat',STATUS='NEW')
      ENDIF

C-----------------------------------------------------------------------
C...Initialization of user action
      IF (ICALL.EQ.1) THEN
         GOTO 100
      ELSEIF (ICALL.EQ.2) THEN
         GOTO 200
      ELSEIF (ICALL.EQ.3) THEN
         GOTO 300
      ENDIF
 
 100  CONTINUE
      IF (ICALL.EQ.1) THEN
c ---------------------------------------------------------------------
c     Open ascii output file
c ---------------------------------------------------------------------
C       EXTENSION='_evt.dat'
        LUNEVT=NextUn()
        OPEN(LUNEVT,FILE=OUTFILENAM(1:ICH)//'_evt.dat',STATUS='NEW')
        write(LUNOUT,*) 'the outputfile will be named: '
     &            ,OUTFILENAM(1:ICH)//'_evt.dat'

C   This is what we write in the ascii-file

C   This is what we write in the ascii-file

        write(LUNEVT,*)' DJANGOH EVENT FILE '
        write(LUNEVT,*)'============================================'
        write(LUNEVT,30)
30      format('I, ievent, IChannel, process, subprocess, nucleon, ',
     &  'struckparton, partontrck, Y, Q2, X, W2, NU, trueY, trueQ2, ',
     &  'trueX, trueW2,trueNu, SIGtot, errSIGtot, D, F1NC, F3NC, ',
     &  'G1NC,G3NC, A1NC, F1CC, F3CC, G1CC, G5CC, nrTracks, status')
        write(LUNEVT,*)'============================================'

        write(LUNEVT,*)' I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5) ',
     &  'P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)',
     &  '  V(I,3)'
        write(LUNEVT,*)'============================================'

      ENDIF
C-----------------------------------------------------------------------
C...Initialize flavor counters
      DO 13 I=-6,6
        IFLCNT(I)=0
 13   CONTINUE
      RETURN

 200  CONTINUE
C-----------------------------------------------------------------------
C...User's action with generated event
      IF (ICHNN.EQ.6.OR.ICHNN.EQ.7.OR.ICHNN.EQ.8.OR.ICHNN.EQ.12) THEN
        DKPRO=EPRO*PHEP(4,3)+PPRO*PHEP(3,3)
        Q2HAD=-TS
        XHAD=Q2HAD/(Y*GSP-2D0*DKPRO)
        YHAD=Q2HAD/XHAD/GSP
      ELSE
        Q2HAD=Q2SCH
        XHAD=XSCH
        YHAD=YSCH
      ENDIF
      W2HAD=YHAD*(1D0-XHAD)*GSP+MPRO2
      WHAD=SQRT(W2HAD)
 
C...Count flavors
      IF (NSPACC.EQ.0) THEN
        DO 20 I=-6,6
   20   IF (LST(25).EQ.I) IFLCNT(I)=IFLCNT(I)+1
      ENDIF

C...count errors in DJGEVT
      DO 21 I=1,10
   21 IF (NFAILI(I).NE.0) NFAILC(I)=NFAILC(I)+1
 
      IF (IHSONL.EQ.0) THEN
C...commented out by eca
C...Check energy-momentum conservation
      IF (LST(21).EQ.0) THEN
C      CALL LUEDIT(1)
C...Print every 'nevmod'th event
C      IF (MOD(NEVHEP,NEVMOD).EQ.0) THEN
C        CALL LULIST(1)
C        write(LUNOUT,*) ' LST(21) = ',LST(21),'    MSTU(24) = ',MSTU(24),
C     .                                   '      NSPACC = ',NSPACC
C      ENDIF
C
      DO 22 J=1,4
   22 PSUM(J)=0D0
      DO 23 I=1,N
      DO 23 J=1,4
   23 IF (K(I,1).GT.0.AND.K(I,1).LE.10) PSUM(J)=PSUM(J)+P(I,J)
      PSUMT=SQRT(PSUM(1)**2+PSUM(2)**2)
      IF (PSUMT.GT.0.1D0) NMIS(1)=NMIS(1)+1
      IF (PSUMT.GT.1.0D0) NMIS(2)=NMIS(2)+1
      IF (PSUMT.GT.2.0D0) NMIS(3)=NMIS(3)+1
      IF (PSUMT.GT.10.0D0) NMIS(4)=NMIS(4)+1
      IF (ABS(PSUM(3)-PELE+PPRO).GT.0.1D0) NMIS(5)=NMIS(5)+1
      IF (ABS(PSUM(3)-PELE+PPRO).GT.1.0D0) NMIS(6)=NMIS(6)+1
      IF (ABS(PSUM(3)-PELE+PPRO).GT.2.0D0) NMIS(7)=NMIS(7)+1
      IF (ABS(PSUM(3)-PELE+PPRO).GT.10.0D0) NMIS(8)=NMIS(8)+1
      IF (ABS(PSUM(4)-EELE-EPRO).GT.0.1D0) NMIS(9)=NMIS(9)+1
      IF (ABS(PSUM(4)-EELE-EPRO).GT.1.0D0) NMIS(10)=NMIS(10)+1
      IF (ABS(PSUM(4)-EELE-EPRO).GT.2.0D0) NMIS(11)=NMIS(11)+1
      IF (ABS(PSUM(4)-EELE-EPRO).GT.10.0D0) NMIS(12)=NMIS(12)+1
      PSUML=PSUM(3)-PELE+PPRO
      ESUM=PSUM(4)-EELE-EPRO
C      DSIGH1=DSIGHB/BINW(21)
C      CALL HFF1(21,NID(21),SNGL(PSUMT),SNGL(DSIGH1))
C      DSIGH1=DSIGHB/BINW(22)
C      CALL HFF1(22,NID(22),SNGL(PSUML),SNGL(DSIGH1))
C      DSIGH1=DSIGHB/BINW(23)
C      CALL HFF1(23,NID(23),SNGL(ESUM),SNGL(DSIGH1))
C      DSIGH1=DSIGHB/BINW(24)
C      CALL HFF2(24,NID(24),SNGL(ESUM),SNGL(WHAD),SNGL(DSIGH1))
      IF ((PSUMT.GT.10D0.OR.ABS(PSUM(3)-PELE+PPRO).GT.10D0
     &   .OR.ABS(PSUM(4)-EELE-EPRO).GT.10D0).AND.N10CNT.LT.3) THEN
        N10CNT=N10CNT+1
        write(LUNOUT,*) ' '
        write(LUNOUT,*) ' '
        write(LUNOUT,*) ' '
        write(LUNOUT,*) ' DJUSER: E-P-mismatch > 10GeV: '
        write(LUNOUT,*) ' NTOT = ',NTOT,'    NEVHEP = ',NEVHEP
        write(LUNOUT,*) ' ICHNN = ',ICHNN
        write(LUNOUT,*) ' PSUM = ',PSUM(1),PSUM(2)
        write(LUNOUT,*) '        ',PSUM(3),PSUM(4)
        write(LUNOUT,*) ' LST(21) = ',LST(21),'    MSTU(24) = ',MSTU(24)
        CALL LULIST(1)
        write(LUNOUT,*) ' '
        write(LUNOUT,*) ' X = ',X
        write(LUNOUT,*) ' Y = ',Y
        write(LUNOUT,*) ' Q2 = ',Q2
        write(LUNOUT,*) ' XHAD = ',XHAD
        write(LUNOUT,*) ' YHAD = ',YHAD
        write(LUNOUT,*) ' Q2HAD = ',Q2HAD
        write(LUNOUT,*) ' WHAD = ',WHAD
        WMOM=(P(2,4)+P(3,4))**2-(P(2,1)+P(3,1))**2
     &      -(P(2,2)+P(3,2))**2-(P(2,3)+P(3,3))**2
        WMOM=SQRT(WMOM)
        write(LUNOUT,*) ' WHAD from momenta on LUJETS = ',WMOM
        WMOM=P(2,4)*P(3,4)-P(2,1)*P(3,1)-P(2,2)*P(3,2)-P(2,3)*P(3,3)
        WMOM=(PHEP(4,6)+PHEP(4,5)-PHEP(4,1)-PHEP(4,3))**2
     &      -(PHEP(1,6)+PHEP(1,5)-PHEP(1,1)-PHEP(1,3))**2
     &      -(PHEP(2,6)+PHEP(2,5)-PHEP(2,1)-PHEP(2,3))**2
     &      -(PHEP(3,6)+PHEP(3,5)-PHEP(3,1)-PHEP(3,3))**2
        WMOM=SQRT(WMOM)
        write(LUNOUT,*) ' WHAD from momenta on HEPEVT = ',WMOM
        IF (ICHNN.GT.2) THEN
          OMEGA=(PH*DKQ+PQH*DKP)/(PH*EQH+PQH*EH)
          DCTHGA=(DKP/OMEGA-ME2/2D0/EH)/PH
          CTHGA=1D0-DCTHGA
          DKPRO=OMEGA*(EPRO+CTHGA*PPRO)
          WHSA=Y*(1D0-X)*(SP-MPRO2)+MPRO2-2D0*(DKP-DKPS+DKPRO)
          WHSA=SQRT(WHSA)
          write(LUNOUT,*) ' WHAD as in HSACPT = ',WHSA
        ENDIF 
      ENDIF
      ENDIF
      ENDIF
 
C...average time per event
  299 CONTINUE
      IF(MOD(NEVHEP,NEVMOD).EQ.0) THEN
        WRITE(LUNOUT,2200) NEVHEP
        WRITE(LUNOUT,2201) ICHNN
        WRITE(LUNOUT,2202) 'leptonic: ',X,Y,Q2
        WRITE(LUNOUT,2202) 'hadronic: ',XHAD,YHAD,Q2HAD
        WRITE(LUNOUT,2202) 'scaled:   ',XSCH,YSCH,Q2SCH
        WRITE(LUNOUT,2203) PHEP(1,1),PHEP(2,1),PHEP(3,1),PHEP(4,1)
     &                    ,PHEP(5,1)
        WRITE(LUNOUT,2204) PHEP(1,2),PHEP(2,2),PHEP(3,2),PHEP(4,2)
     &                    ,PHEP(5,2)
        WRITE(LUNOUT,2205) PHEP(1,3),PHEP(2,3),PHEP(3,3),PHEP(4,3)
     &                    ,PHEP(5,3)
 2200 FORMAT(//,' ***** Event No ',I12)
 2201 FORMAT(/,' CHANNEL = ',I3,9x,'X',12x,'Y',12x,'Q2')
 2202 FORMAT(1x,A10,9x,2(E9.3,4x),E9.3)
 2203 FORMAT('-----------',
     &      9x,'Px',10x,'Py',10x,'Pz',10x,'E ',10x,'M'/,
     &       ' HS e-scat: ',5F12.4)
 2204 FORMAT(' HS q-scat: ',5F12.4)
 2205 FORMAT(' HS ga-rad: ',5F12.4)
 2206 FORMAT(/,' LST(21) = ',I3)
        WRITE(LUNOUT,2206) LST(21)
        CALL LULIST(1)
        CALL TIMEX(RTIME)
        TIMFIN=RTIME
        RTIME=TIMFIN-TIMINI
        TIMINI=TIMFIN
        TIMEVT=RTIME/NEVMOD
        WRITE(LUNOUT,2001) NEVHEP,TIMEVT
      ENDIF
C...write events to file
         nrtrack=N
         NUHAD=(W2HAD+YHAD-MPRO2)/(2.*MPRO2)
         W2=Y*(1D0-X)*GSP+MPRO2
         NU=(W2+Y-MPRO2)/(2.*MPRO2)
         D=(S*S-U*U)/(S*S+U*U)
         CALL HSFG(X,Q2,LLEPT)
         write(LUNEVT,32) 0, NEVHEP, ICHNN, LST(23), LST(24), LST(22),
     &   LST(25), LST(26), Y, Q2, X, W2, NU, YHAD, Q2HAD, XHAD, W2HAD,
     &   NUHAD, SIGTOT,SIGTRR, D,F1NC,F3NC,G1NC,G3NC,A1NC,
     &   F1CC,F3CC,G1CC,G5CC, nrtrack, LST(21)
 32      format((I4,1x),(I10,1x),6(I6,1x),22(E18.10,1x),I12,I6)
         write(LUNEVT,*)'============================================'

         DO I=1,N
           if (K(I,3).le.N) then
             write(LUNEVT,34) (I),K(I,1),K(I,2),K(I,3),K(I,4),K(I,5)
     &                    ,P(I,1),P(I,2),-1.*P(I,3),P(I,4),P(I,5)
     &                    ,V(I,1),V(I,2),V(I,3)
           endif
         ENDDO
 34      format(6(I10,1x),8(f15.6,1x))
         write(LUNEVT,*)'=============== Event finished ==============='

C       write(LUNEVT,2301) ICHNN
C       write(LUNEVT,2302) IDHEP(1)
C       write(LUNEVT,2303) PHEP(1,1),PHEP(2,1),PHEP(3,1),PHEP(4,1),PHEP(5,1)
C       write(LUNEVT,2302) IDHEP(2)
C       write(LUNEVT,2303) PHEP(1,2),PHEP(2,2),PHEP(3,2),PHEP(4,2),PHEP(5,2)
C       write(LUNEVT,2302) IDHEP(3)
C       write(LUNEVT,2303) PHEP(1,3),PHEP(2,3),PHEP(3,3),PHEP(4,3),PHEP(5,3)
C 2301 FORMAT(I4)
C 2302 FORMAT(I4)
C 2303 FORMAT(5E16.8)

      RETURN 
 300  CONTINUE
C-----------------------------------------------------------------------
C...Final call, overall output, program performance
 
      IF (IHSONL.EQ.0) THEN
        WRITE(LUNOUT,3001) NTOT,NREJCW,NPASS,NFAILQ,NFAILL,NSOPH,NFAILP
        IF (NTOT.NE.NPASS) WRITE(LUNOUT,3003) PARL(24)
C...Flavor distribution
        WRITE(LUNOUT,3004) IFLCNT(1),IFLCNT(2),IFLCNT(3)
     &                    ,IFLCNT(4),IFLCNT(5),IFLCNT(6)
     &                    ,IFLCNT(-1),IFLCNT(-2),IFLCNT(-3)
     &                    ,IFLCNT(-4),IFLCNT(-5),IFLCNT(-6)
C...Failures in hadronization
        WRITE(LUNOUT,3005) NFAILC(1),NFAILC(2),NFAILC(3)
     &                    ,NFAILC(4),NFAILC(5),NFAILC(6)
     &                    ,NFAILC(7),NFAILC(8),NFAILC(9)
     &                    ,NFAILC(10)
C...Energy-momentum mismatch
        WRITE(LUNOUT,3006) NMIS(1),NMIS(2),NMIS(3),NMIS(4)
     &                    ,NMIS(5),NMIS(6),NMIS(7),NMIS(8),NMIS(9)
     &                    ,NMIS(10),NMIS(11),NMIS(12)
      ENDIF

      RETURN
 
2001  FORMAT(/,' AVERAGE TIME PER EVENT AFTER ',I10,' EVENTS: ',
     &            1PE11.2,' SEC')
C
3001  FORMAT(/,' ******************************************************'
     F,'************************',/
     F      ,' Program performance ',/
     F      ,1X,I12,' Events were accepted by HERACLES',/
     F      ,1X,I12,' Events do not have min W-remnant in HERACLES',/
     F      ,1X,I12,' Events passed fragmentation in LEPTO',/
     F      ,1X,I12,' Events not accepted for fragmentation in LEPTO',/
     F      ,1X,I12,' Events failed fragmentation in LEPTO',/
     F      ,1X,I12,' Events passed fragmentation in SOPHIA',/
     F      ,1X,I12,' Events failed fragmentation in SOPHIA') 
3003  FORMAT(' Cross section was corrected: ',/
     F      ,' Total cross section is now    SIGTOT = ',E12.5,' nb',/)
3004  FORMAT(/,' Distribution of flavors (DJ):',/
     F        ,10X,'d',9X,'u',9X,'s',9X,'c',9X,'b',9X,'t',/,1X,6I10,/
     F        ,7X,'dbar',6X,'ubar',6X,'sbar',6X,'cbar',6X,'bbar',6X
     F        ,'tbar',/,1X,6I10,/)
3005  FORMAT(/,' Errors in DJGEVT: Hadronization failed:',/
     F        ,' Nfail(1) = ',I8,'     Nfail(2)  = ',I8,/
     F        ,' Nfail(3) = ',I8,'     Nfail(4)  = ',I8,/
     F        ,' Nfail(5) = ',I8,'     Nfail(6)  = ',I8,/
     F        ,' Nfail(7) = ',I8,'     Nfail(8)  = ',I8,/
     F        ,' Nfail(9) = ',I8,'     Nfail(10) = ',I8,/)
3006  FORMAT(/,' Energy-momentum mismatch: ',/
     F        ,' PmisT >  0.1: NMIS(1)  = ',I8,/
     F        ,' PmisT >  1.0: NMIS(2)  = ',I8,/
     F        ,' PmisT >  2.0: NMIS(3)  = ',I8,/
     F        ,' PmisT > 10.0: NMIS(4)  = ',I8,/
     F        ,' PmisL >  0.1: NMIS(5)  = ',I8,/
     F        ,' PmisL >  1.0: NMIS(6)  = ',I8,/
     F        ,' PmisL >  2.0: NMIS(7)  = ',I8,/
     F        ,' PmisL > 10.0: NMIS(8)  = ',I8,/
     F        ,'     E >  0.1: NMIS(9)  = ',I8,/
     F        ,'     E >  1.0: NMIS(10) = ',I8,/
     F        ,'     E >  2.0: NMIS(11) = ',I8,/
     F        ,'     E > 10.0: NMIS(12) = ',I8,/)
 
      END
C 
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C...dummy routines
C...when ARIADNE is not linked: uncomment the following two routines
C 
      SUBROUTINE ARINIT(MODE)
      CHARACTER MODE*(*)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON/LEPTOU/CUT(14),LST(40),PARL(30),X,Y,W2,Q2,U
       WRITE(LUNOUT,101)
  101  FORMAT(/'*************************************************',
     &   /'WARNING: Ariadne is inactive because:',
     &   /'       1/ ARIADNE program is not linked or',
     &   /'       2/ dummy routines in DJUSER are active',
     &   /'         (should be commented out when ARIADNE is linked',
     &   /'  CASCADE flag set to 12, i.e. ME+PS option',
     &   /'*************************************************'/)
        LST(8)=12
      RETURN
      END
cck---------------------------------------------------------------------
      SUBROUTINE AREXEC
cck..dummy
      RETURN
      END
cck---------------------------------------------------------------------
      SUBROUTINE ARPRDA
cck..dummy
      RETURN
      END
C...requisite HSFG subroutine
      SUBROUTINE HSFG(X,Q2,LL)
C-----Computation of electroweak structure functions 
C-----in terms of PDFs for the point X,Y,Q2
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      DOUBLE PRECISION VAFI(2,3,2),FLIND(2,3,2,2),VAEL(2,2)
      COMMON /HSUNTS/ LUNTES,LUNDAT,LUNIN,LUNOUT,LUNRND
      COMMON /HSPDFQ/ QU,QBU,QD,QBD,QS,QBS,QC,QBC,QB,QBB,QT,QBT
      COMMON /HSDPDF/ DQU,DQBU,DQD,DQBD,DQS,DQBS,DQC,DQBC,DQB,DQBB,
     *                DQT,DQBT
      COMMON /HSFGNC/ F1NC,F3NC,G1NC,G3NC
      COMMON /HSFGCC/ F1CC,F3CC,G1CC,G5CC
      COMMON /HSASYM/ A1NC
      COMMON /HSNUCL/ HNA,HNZ,INUMOD
      COMMON /HSPARM/ POLARI,HPOLAR,LLEPT,LQUA
C-----Get PDFs

       
      CALL HSPVER(X,Q2)
      IF (HPOLAR.NE.0D0) THEN
         CALL HSDPVR(X,Q2)
      ENDIF
      CALL HSIROT
      

C-----Constants and couplings     
      PI=4D0*DATAN(1D0) 
      ALPHA=1D0/137.0359895D0
      GF=1.166389D-5
      MZ=91.1867D0
      MW=80.4639D0
      CW=MW/MZ
      CW2=CW*CW
      SW2=1D0-CW2
      VAFI(2,1,1)=0D0
      VAFI(2,2,1)=0D0
      VAFI(2,3,1)=0D0
      VAFI(2,1,2)=-1D0/2D0
      VAFI(2,2,2)=1D0/2D0
      VAFI(2,3,2)=-1D0/2D0
      VAFI(1,1,1)=-1D0
      VAFI(1,2,1)=2D0/3D0
      VAFI(1,3,1)=-1D0/3D0
      VAFI(1,1,2)=(4D0*SW2-1D0)/2D0
      VAFI(1,2,2)=(1D0-8D0*SW2/3D0)/2D0
      VAFI(1,3,2)=(4D0*SW2/3D0-1D0)/2D0

      IGAMMA = 1
      IZ = 2
      IEL = 1
      IFU = 2
      IFD = 3
      INDV = 1
      INDA = 2    
 
      B=MZ*MZ*GF/PI/ALPHA/2D0/DSQRT(2D0)*(Q2/(Q2+MZ*MZ))
      VAEL(1,1)=1D0
      VAEL(1,2)=B*(1D0-2D0*SW2)
      VAEL(2,1)=VAEL(IGAMMA,IZ)
      VAEL(2,2)=VAEL(IGAMMA,IZ)*VAEL(IGAMMA,IZ)
C
      DO 1 IF=IEL,IFD
        DO 1 IB1=IGAMMA,IZ
          DO 1 IB2=IGAMMA,IZ
          FLIND(INDV,IF,IB1,IB2)=
     *     2D0*(VAFI(INDV,IF,IB1)*VAFI(INDV,IF,IB2)
     *         +VAFI(INDA,IF,IB1)*VAFI(INDA,IF,IB2))
    1 CONTINUE
      DO 2 IF=IEL,IFD
        DO 2 IB1=IGAMMA,IZ
          DO 2 IB2=IGAMMA,IZ
          FLIND(INDA,IF,IB1,IB2)=
     *     2D0*(VAFI(INDV,IF,IB1)*VAFI(INDA,IF,IB2)
     *         +VAFI(INDA,IF,IB1)*VAFI(INDV,IF,IB2))
    2 CONTINUE

      F1NC=0D0
      G1NC=0D0
      F3NC=0D0
      G3NC=0D0

C-----compute NC structure functions
      DO 12 IB1=IGAMMA,IZ
          DO 12 IB2=IGAMMA,IZ
           F1NC=F1NC
     *    +VAEL(IB1,IB2)*FLIND(INDV,IFU,IB1,IB2)*(QU+QBU+QC+QBC+QT+QBT)
     *    +VAEL(IB1,IB2)*FLIND(INDV,IFD,IB1,IB2)*(QD+QBD+QS+QBS+QB+QBB)
           F3NC=F3NC
     *    +VAEL(IB1,IB2)*FLIND(INDA,IFU,IB1,IB2)*(QU-QBU+QC-QBC+QT-QBT)
     *    +VAEL(IB1,IB2)*FLIND(INDA,IFD,IB1,IB2)*(QD-QBD+QS-QBS+QB-QBB)
   12 CONTINUE
        DO 13 IB1=IGAMMA,IZ
          DO 13 IB2=IGAMMA,IZ
           G1NC=G1NC
     *    +VAEL(IB1,IB2)*FLIND(INDV,IFU,IB1,IB2)*
     *     (DQU+DQBU+DQC+DQBC+DQT+DQBT)
     *    +VAEL(IB1,IB2)*FLIND(INDV,IFD,IB1,IB2)*
     *     (DQD+DQBD+DQS+DQBS+DQB+DQBB)
           G3NC=G3NC
     *    +VAEL(IB1,IB2)*FLIND(INDA,IFU,IB1,IB2)*
     *     (DQU-DQBU+DQC-DQBC+DQT-DQBT)
     *    +VAEL(IB1,IB2)*FLIND(INDA,IFD,IB1,IB2)*
     *     (DQD-DQBD+DQS-DQBS+DQB-DQBB)
   13 CONTINUE

C      F1NC=(A1U+POL*B1U)*(QU+QBU+QC+QBC+QT+QBT)
C     *    +(A1D+POL*B1D)*(QD+QBD+QS+QBS+QB+QBB)

C      F3NC=-LL*(B3U+POL*A3U)*(QU-QBU+QC-QBC+QT-QBT)
C     *     -LL*(B3D+POL*A3D)*(QD-QBD+QS-QBS+QB-QBB)

C      G3NC=(A3U+POL*B3U)*(DQU-DQBU+DQC-DQBC+DQT-DQBT)
C     *    +(A3D+POL*B3D)*(DQD-DQBD+DQS-DQBS+DQB-DQBB)

C      G1NC=-LL*(B1U+POL*A1U)*(DQU+DQBU+DQC+DQBC+DQT+DQBT)
C     *     -LL*(B1D+POL*A1D)*(DQD+DQBD+DQS+DQBS+DQB+DQBB)
C-----construct asymmetry
      A1NC=G1NC/F1NC

C-----compute CC stucture function
      IF (LL.EQ.-1) THEN
       F1CC=QU+QC+QBD+QBS
       F3CC=2D0*(QU+QC-QBD-QBS)
       G1CC=DQU+DQC+DQBD+DQBS
       G5CC=-DQU-DQC+DQBD+DQBS
      ELSEIF (LL.EQ.1) THEN
       F1CC=QBU+QBC+QD+QS
       F3CC=2D0*(-QBU-QBC+QD+QS)
       G1CC=DQBU+DQBC+DQD+DQS
       G5CC=DQBU+DQBC-DQD-DQS
      ELSE
       WRITE(LUNOUT,*) ' ERROR IN HSFG: WRONG LEPTON CHARGE = ',LL
       STOP
      ENDIF
      RETURN
      END
