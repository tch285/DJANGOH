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
C...Number of events from sophia
      COMMON /SPPASS/ NSOPH,NSPOUT,NFAILP,NSPACC
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
        LUNEVT=NextUn()
        OPEN(LUNEVT,FILE=OUTFILENAM(1:ICH)//'_evt.dat',STATUS='NEW')
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
C-----------------------------------------------------------------------
C...Initialize flavor counters
      DO 13 I=-6,6
        IFLCNT(I)=0
 13   CONTINUE
      RETURN

 200  CONTINUE
C-----------------------------------------------------------------------
C...User's action with generated event

C...Get hadronic kinematic variables (DIS)
      Q2HAD=-TS
      IF (ICHNN.EQ.6.OR.ICHNN.EQ.7.OR.ICHNN.EQ.8.OR.ICHNN.EQ.12) THEN
        DKPRO=EPRO*PHEP(4,3)+PPRO*PHEP(3,3)
        XHAD=Q2HAD/(Y*GSP-2D0*DKPRO)
        YHAD=Q2HAD/XHAD/GSP
      ELSE
        XHAD=X
        YHAD=Y
      ENDIF
      W2HAD=YHAD*(1D0-XHAD)*GSP+MPRO2
      WHAD=SQRT(W2HAD)
 
C...Count flavors
      IF (NSPACC.EQ.0) THEN
        DO 20 I=-6,6
          IF (LST(25).EQ.I) IFLCNT(I)=IFLCNT(I)+1
 20     CONTINUE
      ENDIF

C...Count errors in DJGEVT
      DO 21 I=1,10
        IF (NFAILI(I).NE.0) NFAILC(I)=NFAILC(I)+1
 21   CONTINUE

      IF (IHSONL.EQ.0) THEN
C...Event characteristics, only after successful hadronization
         IF (LST(21).EQ.0) THEN

C...Print every 'nevmod'th event
            IF (MOD(NEVHEP,NEVMOD).EQ.0) THEN
              WRITE(LUNOUT,2200) NEVHEP
              WRITE(LUNOUT,2201) ICHNN
              WRITE(LUNOUT,2206) LST(21),MSTU(24),NSPACC
              WRITE(LUNOUT,2202) 'lepton x,y,Q2:   ',X,Y,Q2
              WRITE(LUNOUT,2202) 'hadron x,y,Q2:   ',XHAD,YHAD,Q2HAD
              WRITE(LUNOUT,2202) 'hadron mass W:   ',WHAD
              WRITE(LUNOUT,2203) PHEP(1,1),PHEP(2,1),PHEP(3,1),PHEP(4,1)
     &                          ,PHEP(5,1)
              WRITE(LUNOUT,2204) PHEP(1,2),PHEP(2,2),PHEP(3,2),PHEP(4,2)
     &                          ,PHEP(5,2)
              WRITE(LUNOUT,2205) PHEP(1,3),PHEP(2,3),PHEP(3,3),PHEP(4,3)
     &                          ,PHEP(5,3)
              CALL LULIST(1)
            ENDIF
 
C...Remove inactive elements from the event record
            CALL LUEDIT(1)

C...Check energy-momentum conservation
            DO 22 J=1,4
              PSUM(J)=0D0
 22         CONTINUE
            DO 231 I=1,N
            DO 232 J=1,4
            IF (K(I,1).GT.0.AND.K(I,1).LE.10) PSUM(J)=PSUM(J)+P(I,J)
 232        CONTINUE
 231        CONTINUE
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
            IF ((PSUMT.GT.10D0.OR.ABS(PSUM(3)-PELE+PPRO).GT.10D0
     &      .OR.ABS(PSUM(4)-EELE-EPRO).GT.10D0).AND.N10CNT.LT.3) THEN
            N10CNT=N10CNT+1
            WRITE(LUNOUT,*) ' '
            WRITE(LUNOUT,*) ' '
            WRITE(LUNOUT,*) ' DJUSER: E-P-mismatch > 10GeV: '
            WRITE(LUNOUT,*) ' NTOT = ',NTOT,'    NEVHEP = ',NEVHEP
            WRITE(LUNOUT,*) ' ICHNN = ',ICHNN
            WRITE(LUNOUT,*) ' PSUM = ',PSUM(1),PSUM(2)
            WRITE(LUNOUT,*) '        ',PSUM(3),PSUM(4)
            WRITE(LUNOUT,*) 
     &      ' LST(21) = ',LST(21),'    MSTU(24) = ',MSTU(24)
            CALL LULIST(1)
            WRITE(LUNOUT,*) ' '
            WRITE(LUNOUT,*) ' X = ',X
            WRITE(LUNOUT,*) ' Y = ',Y
            WRITE(LUNOUT,*) ' Q2 = ',Q2
            WRITE(LUNOUT,*) ' XHAD = ',XHAD
            WRITE(LUNOUT,*) ' YHAD = ',YHAD
            WRITE(LUNOUT,*) ' Q2HAD = ',Q2HAD
            WRITE(LUNOUT,*) ' WHAD = ',WHAD
            WMOM=(P(2,4)+P(3,4))**2-(P(2,1)+P(3,1))**2
     &          -(P(2,2)+P(3,2))**2-(P(2,3)+P(3,3))**2
            WMOM=SQRT(WMOM)
            WRITE(LUNOUT,*) ' WHAD from momenta on LUJETS = ',WMOM
            WMOM=P(2,4)*P(3,4)-P(2,1)*P(3,1)-P(2,2)*P(3,2)-P(2,3)*P(3,3)
            WMOM=(PHEP(4,6)+PHEP(4,5)-PHEP(4,1)-PHEP(4,3))**2
     &          -(PHEP(1,6)+PHEP(1,5)-PHEP(1,1)-PHEP(1,3))**2
     &          -(PHEP(2,6)+PHEP(2,5)-PHEP(2,1)-PHEP(2,3))**2
     &          -(PHEP(3,6)+PHEP(3,5)-PHEP(3,1)-PHEP(3,3))**2
            WMOM=SQRT(WMOM)
            WRITE(LUNOUT,*) ' WHAD from momenta on HEPEVT = ',WMOM
            IF (ICHNN.GT.2) THEN
              OMEGA=(PH*DKQ+PQH*DKP)/(PH*EQH+PQH*EH)
              DCTHGA=(DKP/OMEGA-MEF2/2D0/EH)/PH
              CTHGA=1D0-DCTHGA
              DKPRO=OMEGA*(EPRO+CTHGA*PPRO)
              WHSA=Y*(1D0-X)*(SP-MPRO2)+MPRO2-2D0*(DKP-DKPS+DKPRO)
              WHSA=SQRT(WHSA)
              WRITE(LUNOUT,*) ' WHAD as in HSACPT = ',WHSA
            ENDIF 
          ENDIF
        ENDIF
      ENDIF
 
C...Print clean event record
      IF(MOD(NEVHEP,NEVMOD).EQ.0) THEN
        CALL LULIST(1)
C...Average time per event
        CALL TIMEX(RTIME)
        TIMFIN=RTIME
        RTIME=REAL(TIMFIN-TIMINI)
        TIMINI=TIMFIN
        TIMEVT=RTIME/NEVMOD
        WRITE(LUNOUT,2001) NEVHEP,TIMEVT
      ENDIF

C...Write events to a file
c      WRITE(LUNEVT,2301) ICHNN
c      WRITE(LUNEVT,2302) IDHEP(1)
c      WRITE(LUNEVT,2303)
c     & PHEP(1,1),PHEP(2,1),PHEP(3,1),PHEP(4,1),PHEP(5,1)
c      WRITE(LUNEVT,2302) IDHEP(2)
c      WRITE(LUNEVT,2303)
c     & PHEP(1,2),PHEP(2,2),PHEP(3,2),PHEP(4,2),PHEP(5,2)
c      WRITE(LUNEVT,2302) IDHEP(3)
c      WRITE(LUNEVT,2303)
c     & PHEP(1,3),PHEP(2,3),PHEP(3,3),PHEP(4,3),PHEP(5,3)

 2200 FORMAT(///,42(' *'),/,' ***** Event No ',I12)
 2201 FORMAT(/,' CHANNEL = ',I3)
 2202 FORMAT(1x,A17,9x,2(E9.3,4x),E9.3)
 2203 FORMAT('-----------',
     &      9x,'Px',10x,'Py',10x,'Pz',10x,'E ',10x,'M'/,
     &       ' HS e-scat: ',5F12.4)
 2204 FORMAT(' HS q-scat: ',5F12.4)
 2205 FORMAT(' HS ga-rad: ',5F12.4)
 2206 FORMAT(' LST(21) = ',I3,'  MSTU(24) = ',I3,'  NSPACC = ',I3,/)

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
