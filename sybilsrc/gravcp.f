      SUBROUTINE GRAVCP(NCOMP,IR,DELRHO,ZMIN,ZMAX,BOUND,NL,
     1  GSUM,RZERO,IVALS,NG,RMAX,PROFLN,ZMEAS,DSCALE,RHSCALE,IERR)
C
C  routine to calculate gravity profile, placed in GSUM, from
C  surface integral defined by boundary segments in BOUND(4,NL)
C  NL is number of boundary segments, listed in order X1,Y1,X2,Y2
C  with points 1 and 2 in anticlockwise order around body of 
C  anomalous density DELRHO.  Radial coordinate of measurement
C  position is stored in RZERO, calculated from profile length
C  PROFLN, and number of measurment points NG.
C  For the purpose of computing the gravity anomaly, the structure
C  at r = RMAX is assumed to continue with the same stratification
C  out to infinite r.
C
      REAL RZERO,DRG,BIGG,CONST1,DIFF,DELRHO,RADB,RADE
      REAL*8 GSUMS,GSUM2
      DIMENSION RZERO(NG),GSUM(NG),IVALS(NG),BOUND(4,NL)
C
C   Anomaly is zeroed, only for first call with IR = 1
C   Subsequent contributions are added in
C
      DRG=PROFLN/FLOAT(NG-1)
      IF(IR.EQ.1)THEN
        DO K=1,NG
          RZERO(K)=FLOAT(K-1)*DRG
          IVALS(K)=1
          GSUM(K)=0.0
        ENDDO
      END IF
      IF(NL.EQ.0)RETURN
C
C   The constants required here are transmitted from plot_gravity,
C    except for BIGG and the scale factor of 1.0e+5 required for
C    conversion from MKS to mgal.
C     DSCALE is the assumed length scale (default 10^5 m)
C     RHSCALE is the assumed density scale (default 30 kgm^-3)
C     ZMEAS is the measurement height (default of 1.05 based
C              on .05*DSCALE above upper surface presumed to be at Z =1)
C
C     ZMEAS=1.05
C     DSCALE=1.0E+5
C     RHSCALE=30.0
      PI=3.141592654
      PHI=0.5*PI
      BIGG=1.E5*6.6726E-11
      CONST1=-4.0*BIGG*DELRHO*DSCALE*RHSCALE
      EPS=1.0E-5
      EPS1=1.0E-10
      NLR=NL/2
      IVERB=1
C
C   First determine Rmax
C
      RMAXL=BOUND(1,1)
      RMIN=BOUND(1,1)
      DO 70 J=1,NL
      IF(BOUND(1,J).GT.RMAXL)RMAXL=BOUND(1,J)
      IF(BOUND(3,J).GT.RMAXL)RMAXL=BOUND(3,J)
      IF(BOUND(1,J).LT.RMIN)RMIN=BOUND(1,J)
      IF(BOUND(3,J).LT.RMIN)RMIN=BOUND(3,J)
   70 CONTINUE
C
C   Next determine Zmin and Zmax, extremes on r = Rmax
C
      ZMIN=999.999
      ZMAX=-999.999
      INFLAY=0
      DO 72 J=1,NL
      IF((ABS(BOUND(1,J)-RMAX).LT.EPS).AND.
     1   (ABS(BOUND(3,J)-RMAX).LT.EPS))THEN
        INFLAY=1
        IF(BOUND(2,J).LT.ZMIN)ZMIN=BOUND(2,J)
        IF(BOUND(4,J).GT.ZMAX)ZMAX=BOUND(4,J)
      END IF
   72 CONTINUE
C
C    output some information about calculation to background window
C
      IF(NCOMP.NE.2)THEN
        WRITE(6,10001)
10001   FORMAT('Warning: the gravity calculation is only valid ',
     1         'for NCOMP = 2 (cylindical symmetry about r = 0)')
        IERR=1
        RETURN
      END IF
      WRITE(6,10002)IR,NL,DELRHO,ZMIN,ZMAX
10002 FORMAT('Gravity profile calculation called for region',I3,
     1' with',I5,' boundary segments and',/,'density = ',F8.3,
     2' with limits in z =',2F9.4)
      IF(ABS(RMIN).GT.EPS)THEN
        WRITE(*,*)'Warning: minimum radius for this region is RMIN = ',
     1  RMIN
      END IF
      IF(INFLAY.NE.1)THEN
        WRITE(6,*)'This region is closed at outer radius r = ',RMAXL
      END IF
C
C   Commence integration for each measurement point I, with NSUB
C   subdivisions of the boundary segments.  You may experiment with increasing
C   NSUB in order to reduce discretisation error, but at some level round-off
C   error will become a problem if NSUB is too large.
C
      NSUB=8
      DO 80 I=1,NG
        RZER=FLOAT(I-1)*DRG
        GSUMS=0.0
        DO 60 J=1,NL
          DR=BOUND(3,J)-BOUND(1,J)
          DZ=BOUND(4,J)-BOUND(2,J)
          DL=SQRT(DR*DR+DZ*DZ)
          CALP=DR/DL
          RADBF=BOUND(1,J)
          RADEF=BOUND(3,J)
C
C   Ignore the path segment if both endpoints lie on r=0 or r=Rmax
C
          IF(((RADBF.GE.EPS).OR.(RADEF.GE.EPS)).AND.
     1       ((ABS(RADBF-RMAX).GE.EPS).OR.(ABS(RADEF-RMAX).GE.EPS)))THEN
C
C   First add the path contributions from the external boundary
C
            ZBF=BOUND(2,J)-ZMEAS
            ZEF=BOUND(4,J)-ZMEAS
C
C  each segment may be subdivided NSUB times to improve accuracy
C   (significant when RAD is close to RZER)
C
            DL=DL/FLOAT(NSUB)
            DO 59 JSUB=1,NSUB
              X1=FLOAT(JSUB-1)/FLOAT(NSUB)
              X2=FLOAT(JSUB)/FLOAT(NSUB)
              RADB=RADBF*(1.0-X1)+RADEF*X1
              ZB=ZBF*(1.0-X1)+ZEF*X1
              RADE=RADBF*(1.0-X2)+RADEF*X2
              ZE=ZBF*(1.0-X2)+ZEF*X2
              AVALB=SQRT((RADB+RZER)*(RADB+RZER)+ZB*ZB)
              AVALE=SQRT((RADE+RZER)*(RADE+RZER)+ZE*ZE)
              IF(AVALB.GT.EPS1)THEN
                AKVALB=(2.0*(SQRT(ABS(RADB*RZER))))/AVALB
                DIFFB=RADB*ELLF(PHI,AKVALB)/AVALB
              ELSE
                DIFFB=0.0
              END IF
              IF(AVALE.GT.EPS1)THEN
                AKVALE=(2.0*(SQRT(ABS(RADE*RZER))))/AVALE
                DIFFE=RADE*ELLF(PHI,AKVALE)/AVALE
              ELSE
                DIFFE=0.0
              END IF
              GSUMS=GSUMS+CALP*DL*0.5*(DIFFB+DIFFE)
   59       CONTINUE
          END IF
   60   CONTINUE
C
C   For an infinite layer add the return path contributions along the horizontal 
C   surfaces at Zmin and Zmax: inward path along Zmin (K=2), outward along Zmax (K=1).
C
        GSUM2=0.0
        IF(INFLAY.EQ.1)THEN
        DR=RMAX/FLOAT(NLR)
        DO 65 K=1,2
          IF(K.EQ.1)THEN
            ZM=ZMAX-ZMEAS
          ELSE
            ZM=ZMIN-ZMEAS
          END IF
        DO J=1,NLR
          IF(K.EQ.1)THEN
            RADBF=DR*FLOAT(J-1)
            RADEF=DR*FLOAT(J)
          ELSE
            RADBF=DR*FLOAT(J)
            RADEF=DR*FLOAT(J-1)
          END IF
          SIGNDR=FLOAT(3-2*K)*DR
          SIGNDR=SIGNDR/FLOAT(NSUB)
          DO 64 JSUB=1,NSUB
            X1=FLOAT(JSUB-1)/FLOAT(NSUB)
            X2=FLOAT(JSUB)/FLOAT(NSUB)
            RADB=RADBF*(1.0-X1)+RADEF*X1
            RADE=RADBF*(1.0-X2)+RADEF*X2
            AVALB=SQRT((RADB+RZER)*(RADB+RZER)+ZM*ZM)
            AVALE=SQRT((RADE+RZER)*(RADE+RZER)+ZM*ZM)
            IF(AVALB.GT.EPS1)THEN
              AKVALB=(2.0*(SQRT(ABS(RADB*RZER))))/AVALB
              DIFFB=RADB*ELLF(PHI,AKVALB)/AVALB
            ELSE
              DIFFB=0.0
            END IF
            IF(AVALE.GT.EPS1)THEN
              AKVALE=(2.0*(SQRT(ABS(RADE*RZER))))/AVALE
              DIFFE=RADE*ELLF(PHI,AKVALE)/AVALE
            ELSE
              DIFFE=0.0
            END IF
            GSUM2=GSUM2+SIGNDR*0.5*(DIFFB+DIFFE)
   64     CONTINUE
        ENDDO         ! on J
   65   CONTINUE      ! on K
        END IF
        GSUM(I)=GSUM(I)+CONST1*(GSUMS+GSUM2)
C       WRITE(6,*)'in gravcp: I,GSUMS,GSUM2,GSUM(I) =',I,GSUMS,
C    :             GSUM2,GSUM(I)
   80 CONTINUE
      IERR=0
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   The next two subroutines are taken from Numerical Recipes
C
      FUNCTION ellf(phi,ak)
      REAL ellf,ak,phi
C  USES rf
C
C  Legendre elliptical integral of the 1st kind F(phi,ak) (refer
C  EQUATION 11 OR EQUATION (6.11.17) in numerical recipies.
C  Integral is evaluated using Carlsons function RF.
C  The argument ranges are 0<=phi<=90,0<=k sin phi<=1
C
      REAL s,rf
      s=sin(phi)
      ellf=s*rf(cos(phi)**2,(1.-s*ak)*(1.+s*ak),1.)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION rf(x,y,z)
C
C  Computes carlsons elliptical integral of the first kind, RF(x,y,z).
C  x,y and Z must be non negative and at most one can be zero.  TINY
C  must be at least 5 times the machines underflow limit, BIG at most
C  one fifth the machine overflow limit.
C
      REAL rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.08,TINY=1.5e-38,BIG=3.E37,THIRD=1./3.,
     *C1=1./24.,C2=.1,C3=3./44.,C4=1./14.)
      REAL alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,
     *z).gt.BIG)WRITE(*,*)'*** WARNING *** invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25*(xt+alamb)
        yt=.25*(yt+alamb)
        zt=.25*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return
      END
C (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.

