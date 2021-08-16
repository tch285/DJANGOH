C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C...  For test run (23.06.2021) 
C     This routine is called after reading input and after
C     checking input for consistency, but before integration
C     and event sampling. It can be used to perform calculations
C     where integrated cross sections are not needed.
C     Use with care!

      subroutine tstrrr

C...Print a table of structure function values

      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      COMMON /HSELAB/ SP,EELE,PELE,EPRO,PPRO
      COMMON /HSGSW1/ MEI,MEF,MQI,MQF,MEI2,MEF2,MQI2,MQF2,MPRO,MPRO2
      COMMON /HSCUTS/ XMIN,XMAX,Q2MIN,Q2MAX,YMIN,YMAX,WMIN,GMIN
ckc..from L61
      COMMON /LEPTOU/CUT(14),LST(40),PARL(30),XHAD,YHAD,W2HAD,Q2HAD,UHAD
      REAL CUT,PARL,ULMASS,XHAD,YHAD,W2HAD,Q2HAD,UHAD,WR
C...L61 (L52 also) uses old PYTHIA common/PYPARA/
      COMMON /LYPARA/ IPY(80),PYPAR(80),PYVAR(80)
      REAL                    PYPAR    ,PYVAR

      integer i,ifail 
      double precision w

      write(*,*) ' '
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) ' '
      write(*,*) ' TSTRRR: '
      write(*,*) ' Test HSWCUT '
      write(*,*) ' '

      write(*,*) ' reset PYPAR(12) to 1.0 '
      pypar(12)=1.0d0
      write(*,*) ' reset PARL(3) to 0.2 '
      parl(3)=0.2d0

      write(*,*) ' pypar(12)',pypar(12)
      write(*,*) ' ulmass(1) = ',ulmass(1)
      write(*,*) ' ulmass(2) = ',ulmass(2)
      write(*,*) ' ulmass(12) = ',ulmass(12)
      write(*,*) ' ulmass(2212) = ',ulmass(2212)

      do 10 i=1,30
      w=1.0d0+dfloat(i)/1d1
      call HSWCUT(w,1,12,ifail) 
      write(*,*) ' w = ',w,' ifail = ',ifail
 10   continue

      return
      end
