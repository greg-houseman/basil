C*--------------------------------------------------------------------
C*    Basil / Sybil:   topog.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------

      SUBROUTINE TOPO(S,NX3,NY3,SS0,TL0,ELEV0,BGAM0,BGAM1,
     1                    EX,EY,NOR,LEM,NUP,IELFIX,IHELP,
     2                    NE,NP3)
C
C     Scales the interpolated crustal thickness matrix to
C     topography by the relation  T = ELEV0 + BGAM*(S-SS0)
C     where BGAM is either BGAM0 or BGAM1 depending on the
C     value in the interpolation point matrix IELFIX.  The
C     criterion for choosing BGAM0 or BGAM1 is to look at the
C     integration point which is within the element, and spatially
C     closest to the interpolation point.  IHELP contains the
C     element numbers for each interpolation point and is set up
C     in NTRPLT.
C
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION NOR(NUP),IELFIX(NUP)
      DIMENSION IHELP(NP3),LEM(6,NE)
      DIMENSION S(NX3,NY3)
      DIMENSION ALFA(3,7),D(7),XV(3),YV(3)
      DATA ALFA/0.79742699,0.10128651,0.10128651,
     1          0.10128651,0.79742699,0.10128651,
     2          0.10128651,0.10128651,0.79742699,
     3          0.47014206,0.05971587,0.47014206,
     4          0.47014206,0.47014206,0.05971587,
     5          0.05971587,0.47014206,0.47014206,
     6          0.33333333,0.33333333,0.33333333/
C
C    For every interpolation point
C
      NX2=NX3-1
      NY2=NY3-1
      NX=NX3-3
      HS=1.0/FLOAT(NX)
      DO 100 J=2,NY2
        Y=FLOAT(J-2)*HS
        DO I=2,NX2
          X=FLOAT(I-2)*HS
          N=IHELP((J-1)*NX3+I)
C
C    if interpolation point is not within mesh
C
          IF(N.EQ.0)THEN
            S(I,J)=-5.0+1.0E-4*MOD(I+J,9)
          ELSE
C
C      get location of nodal vertices
C
          DO M=1,3
            LN=NOR(LEM(M,N))
            XV(M)=EX(LN)
            YV(M)=EY(LN)
          ENDDO
C
C    get location of 7 integration points for element N
C
          DO K=1,7
            XIP=XV(1)*ALFA(1,K) + XV(2)*ALFA(2,K) + XV(3)*ALFA(3,K)
            YIP=YV(1)*ALFA(1,K) + YV(2)*ALFA(2,K) + YV(3)*ALFA(3,K)
            XD=X-XIP
            YD=Y-YIP
            D(K)=SQRT(XD*XD + YD*YD)
          ENDDO
C
C     find integration point closest to interpolation point
C
          DMIN=D(1)
          KMIN=1
          DO K=2,7
            IF(D(K).LT.DMIN)THEN
              DMIN=D(K)
              KMIN=K
            ENDIF
          ENDDO
C
C     Now scale the crustal thickness
C
          IF(IELFIX(N).EQ.0)THEN
            S(I,J)=ELEV0 + BGAM0*TL0*(S(I,J)-SS0)
          ELSE
            S(I,J)=ELEV0 + BGAM1*TL0*(S(I,J)-SS0)
          ENDIF
        ENDIF           !  value of N zero or non-zero
        ENDDO           !  on I
  100 CONTINUE          !  on J
      RETURN
      END

      SUBROUTINE SNORM(S,NX3,NY3,SFAC,SZERO)
C
C    The logarithm of crustal thickness is converted to a linear
C      number, scaled, and zero-shifted
C
      DIMENSION S(NX3,NY3)
      DO J=1,NY3
        DO I=1,NX3
          S(I,J)=SFAC*EXP(S(I,J))+SZERO
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE VELOVECTORS(UX,UY,XCMIN,YCMIN,DX,DY,
     :                           M1,M2,N1,N2,NXY,NX3)
C
C    Routine to calculate velocity vectors in the x,y plane for
C    spherical cases (NCOMP=-1)
C    Radius is assumed to be 1
C
      DIMENSION UX(NXY),UY(NXY)

      DO I=M1,M2
        DO J=N1,N2
          K=(J-1)*NX3+I
          XPT=XCMIN+FLOAT(I-M1)*DX
          YPT=YCMIN+FLOAT(J-N1)*DY
          YCOT = 1.0/TAN(YPT)
          TMP = SQRT(1.0+XPT*XPT*YCOT*YCOT)
          UX(K)= (UX(K)*TMP + UY(K)*XPT*YCOT)/TMP
          UY(K)= UY(K)/TMP
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE PROJECTXY(EX,EY,XMID,YMID,NUP,NCOMP,IERROR)
C
C    Routine to project data in longitude, latitude
C    into x, y space
C    Radius is assumed to be 1
C    Latitude values adjusted to the range 0->180 (0 at S.Pole)
C    Longitude values centred on the longitude at centre of projection
C    Convert to radians (DTOR=PI/180)
C
      INCLUDE 'limits.parameters'
      DIMENSION EX(NUP),EY(NUP)
C     WRITE(6,1000)EX(1),EY(1)
C
C    this block for sinusoidal equal area projection
C
      IF (NCOMP.LT.0) THEN
        DO 10 I=1,NUP
          EY(I)=EY(I)*DTOR+R_PION2
          EX(I)=(EX(I)-XMID)*DTOR*SIN(EY(I))
  10    CONTINUE
C
C    this block uses invariant X-scale; actually used ?
C
      ELSE
        YVAL=SIN(YMID*DTOR+R_PION2)
        DO 20 I=1,NUP
          EY(I)=EY(I)*DTOR+R_PION2
          EX(I)=(EX(I)-XMID)*DTOR*YVAL
  20    CONTINUE
      END IF
C     WRITE(6,1001)EX(1),EY(1)
      RETURN
 1000 FORMAT('Projecting lat, long (',F6.2,', ',F6.2,')')
 1001 FORMAT('to mesh coordinates (x,y) (',F6.3,', ',F6.3,')')
      END

      SUBROUTINE PROJECTDEG(EX,EY,XMID,YMID,NUP,NCOMP,IERROR)
C
C    Routine to project data in x, y
C    into longitude, latitude space
C    Latitude values adjusted to the range 0->180 (0 at S.Pole)
C    Longitude values centred on the longitude at centre of projection
C    Convert to degrees (RTOD=180/PI)
C
      INCLUDE 'limits.parameters'
      DIMENSION EX(NUP),EY(NUP)

      IERROR=0
C     WRITE(6,1000)EX(1),EY(1)
C
C    this block for sinusoidal equal area projection
C
      IF (NCOMP.LT.0) THEN
        DO 10 I=1,NUP
          EX(I)=EX(I)/SIN(EY(I))/DTOR+XMID
          EY(I)=(EY(I)-R_PION2)/DTOR
  10    CONTINUE
C
C    this block uses invariant X-scale; actually used ?
C
      ELSE
        YVAL=ASIN(YMID-R_PION2)/DTOR
        DO 20 I=1,NUP
          EY(I)=(EY(I)-R_PION2)/DTOR
          EX(I)=EX(I)/SIN(YVAL)/DTOR+XMID
  20    CONTINUE
      END IF
C     WRITE(6,1001)EX(1),EY(1)
      RETURN
 1000 FORMAT('Projecting x,y (',F6.2,', ',F6.2,')')
 1001 FORMAT('to lat, long (',F6.3,', ',F6.3,')')
      END


C
C   Stuff without a home
C
C   2. Crustal thickness is integrated along the x-direction to
C      estimate shortening as a function of y
C
C    For a given y-coordinate, estimate the x-coordinates of the
C     limiting x-boundaries: Strategy - go around boundary using
C     array NEDGES.  When YCOORD is crossed, interpolate to find
C     XMIN and XMAX.
C
C     WRITE(LUW,*)'   J    YCOORD    XMIN   XMAX'
C     WRITE(LUW,*)'   J     YCOORD      EIND      XMIN'
C     SS0C=EXP(SS0)
C     NY2=NY3-1
C     NX2=NX3-1
C     SUMA=0.0
C     SUMI=0.0
C     DO 51 J=2,NY2
C     YCOORD=FLOAT(J-2)*DY+YCMIN
C     NNE=NOR(NEDGES(1))
C     XX1=EX(NNE)
C     YY1=EY(NNE)
C     XMIN=99999.9
C     XMAX=-99999.9
C     DO 50 K=2,NED
C     NNE=NOR(NEDGES(K))
C     XX2=EX(NNE)
C     YY2=EY(NNE)
C     IF((ABS((YY2-YY1)).LT.1.E-4).AND.
C    1                     (ABS((YY2-YCOORD)).LT.1.E-4))THEN
C       XMIN=MIN(XMIN,XX2)
C       XMAX=MAX(XMAX,XX2)
C     ELSE IF((YY2-YCOORD)*(YY1-YCOORD).LE.1.E-6)THEN
C       XLIM=XX1+(XX2-XX1)*(YCOORD-YY1)/(YY2-YY1)
C       IF(YY2.GT.YY1)THEN
C         XMAX=XLIM
C       ELSE
C         XMIN=XLIM
C       END IF
C     END IF
C     XX1=XX2
C     YY1=YY2
C  50 CONTINUE
CC    WRITE(LUW,*)J,YCOORD,XMIN,XMAX
C     IF((XMIN.GT.99999.0).OR.(XMAX.LT.-99999.0))THEN
C     WRITE(LUW,*)'Problem in SNORM, Loop 50'
C     STOP
C     END IF
C
CC   Integrate between XMIN and XMAX using trapezium rule
C
C     SUM=0.0
C     DO 20 I=2,NX2
C     XCOORD=FLOAT(I-2)*DX+XCMIN
C     IF(XCOORD.LE.(XMIN-DX*0.5))THEN
C       WEIGHT=0.0
C     ELSE IF(XCOORD.GE.(XMIN+DX*0.5))THEN
C       WEIGHT=DX
C       IF(XCOORD.GE.(XMAX-DX*0.5))WEIGHT=XMAX-XCOORD+DX*0.5
C     ELSE
C       WEIGHT=XCOORD-XMIN+DX*0.5
C     END IF
C     SUM=SUM+WEIGHT*S(I,J)
C  20 CONTINUE
CC    WRITE(LUW,*)XMIN,XMAX,SFAC,SS0C,SUM,DX,S(20,20)
C     EIND=SUM/(SFAC*SS0C) - (XMAX-XMIN)
C     DYF=DY
C     IF((J.EQ.2).OR.(J.EQ.NY2))DYF=DY*0.5
C     SUMA=SUMA+EIND*DYF
C
CC   XZERO depends on innitial boundary geometry:
CC   XZERO=0.0 for D series, =-0.83910 for E series
CC   XZERO=-1.67820+YCOORD*0.83910 for F series
C
C     XZERO=-1.67820+YCOORD*0.83910
C     XMIN=XMIN-XZERO
C     SUMI=SUMI+XMIN*DYF
C     WRITE(LUW,10001)J,YCOORD,EIND,XMIN
C10001 FORMAT(I6,3F12.6)
C  51 CONTINUE
C     WRITE(LUW,*)' '
C     PCERR=(SUMI-SUMA)*100.0/SUMI
C     WRITE(LUW,10002)SUMI,SUMA,PCERR
C10002 FORMAT('Indentation area: ',G12.5,/,'Estimate from crustal',
C    1' thickness distribution: ',G12.5,'   % error: ',F8.3,/)
C     RETURN
