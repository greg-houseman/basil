C*--------------------------------------------------------------------
C*    Basil / Sybil:   profil.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE NTRPLN(S,RLINE,NL,SS0,XP1,YP1,XP2,YP2,
     :                  TBXOFF,TBYOFF,VOFFSET,EX,EY,
     :                  XREFM,YREFM,LEM,NOR,ILINE,NE,NUP,NP3,
     :                  NCOMP,IVERBOSE,LABEL,LEN)
C
C    Routine to interpolate the nodal values from the finite
C     element mesh onto a line traversing the solution region
C     NTRPLN uses quadratic interpolation (e.g. velocity)
C     NTRPLN2 uses linear interpolation (e.g. pressure)
C
      CHARACTER FNAME*10,LABEL*10
      LOGICAL XIST
      DIMENSION EX(NUP),EY(NUP),LEM(6,NE),NOR(NUP)
      DIMENSION S(NUP),RLINE(NL),ILINE(NP3)
      DOUBLE PRECISION BY(3),CX(3),CNL,COOL(3),TRI,X1,X2,X3,Y1,Y2,Y3
      DATA EPSL/1.0E-3/
      SAVE EPSL
C
C   Zero the arrays
C   RLINE(NL) contains the intrpolated values on the line
C   ILINE(NL) contains a switch indicating a value in RLINE
C   Line runs between (XP1,YP1) and (XP2,YP2)
C
      DO I=1,NL
        ILINE(I)=0
        RLINE(I)=SS0
      ENDDO
      DXP=(XP2-XP1)/FLOAT(NL-1)
      DYP=(YP2-YP1)/FLOAT(NL-1)
      EPS=0.01*SQRT(DXP*DXP+DYP*DYP)
C
C    Look at each element in turn
C
      DO 60 N=1,NE
C
C     Calculate the geometrical coefficients
C
        XMIN=999.
        YMIN=999.
        XMAX=-999.
        YMAX=-999.
        DO 20 K=1,3
          K2=MOD(K,3)+1
          K3=MOD(K+1,3)+1
          LK=NOR(LEM(K,N))
          LK2=NOR(LEM(K2,N))
          LK3=NOR(LEM(K3,N))
          X1=EX(LK)
          X2=EX(LK2)
          X3=EX(LK3)
          Y1=EY(LK)
          Y2=EY(LK2)
          Y3=EY(LK3)
          XMIN=MIN(XMIN,X1)
          XMAX=MAX(XMAX,X1)
          YMIN=MIN(YMIN,Y1)
          YMAX=MAX(YMAX,Y1)
          BY(K)=Y2    - Y3
          CX(K)=X3    - X2
   20   CONTINUE
C
C    TRI is twice the area of the triangle element
C
      TRI=CX(3)*BY(2)-CX(2)*BY(3)
C
C     Now look for the points on the line that are inside the tri.
C
        XP = XP1-DXP
        YP = YP1-DYP
        DO 80 J=1,NL
          XP=XP+DXP
          YP=YP+DYP
          IF((YP.GT.YMAX+EPS).OR.(YP.LT.YMIN-EPS).OR.
     1       (XP.GT.XMAX+EPS).OR.(XP.LT.XMIN-EPS).OR.
     2       (ILINE(J).NE.0)) GO TO 80
C
C     Calculate the natural coordinates
C
          DO 45 K1=1,3
            K2=MOD(K1,3)+1
            LK2=NOR(LEM(K2,N))
            CNL=(CX(K1)*(YP-EY(LK2)) + BY(K1)*(XP-EX(LK2)))/TRI
            IF((CNL.GT.1.0+EPSL).OR.(CNL.LT.-EPSL))GO TO 80
            COOL(K1)=CNL
   45     CONTINUE
C
C     If we reach here the point is within the triangle, so
C         Interpolate !
C
          SI=0.0
          DO K=1,6
            LK=LEM(K,N)
            IF(LK.LT.0)LK=-LK
            SK=S(LK)
            IF(K.LE.3)THEN
              K2=MOD(K,3)+1
              K3=MOD(K+1,3)+1
              CNK=COOL(K)*(COOL(K) - COOL(K2) - COOL(K3))
            ELSE
              K2=K-3
              K3=MOD(K+1,3)+1
              CNK=4.0*COOL(K2)*COOL(K3)
            ENDIF
            SI=SI + SK*CNK
          ENDDO
          RLINE(J)=SI
          ILINE(J)=N
   80   CONTINUE
   60 CONTINUE
C
C   Output the interpolated (x,y,value) if required
C
C      IVERBOSE=0
C      IF(IVERBOSE.NE.0)THEN
C          FNAME(1:10) = ' '
C          FNAME='sybil.'//LABEL(1:LEN)
C          INQUIRE(FILE=FNAME,EXIST=XIST)
C          IF (XIST) THEN
C              OPEN(10,FILE=FNAME,STATUS='old',
C     :                POSITION='append',IOSTAT=IOS)
C          ELSE
C              OPEN(10,FILE=FNAME)
C          ENDIF
C          WRITE(6,10002)LABEL(1:LEN),FNAME,LABEL(1:LEN)
C          WRITE(10,10004)LABEL(1:LEN) ! header
C          DXP=(XP2-XP1)/FLOAT(NL-1)
C          DYP=(YP2-YP1)/FLOAT(NL-1)
C          DO 212 I=1,NL
C              X=XP1+DXP*FLOAT(I-1)
C              Y=YP1+DYP*FLOAT(I-1)
C              IF(NCOMP.EQ.-1)THEN
C                CALL PROJECTDEG(X,Y,XREFM,YREFM,1,NCOMP,IERR)
C              ENDIF
C              WRITE(10,10003)X,Y,RLINE(I)
C  212     CONTINUE
C          CLOSE(10)
C      ENDIF
C10002 FORMAT('interpolated values of ',A4,' written to ',A10,
C     1 '; 3 columns: X, Y, ',A4)
C10004 FORMAT('X, Y, ',A4)
C10003 FORMAT(2F12.5,G15.6)
C
C    check that there are some data points on the nominated profile
C
      ICNT = 0
      DO 100 I=1,NL
        IF (ILINE(I).NE.0) ICNT = ICNT+1
 100  CONTINUE
      IF (ICNT.EQ.0) WRITE(6,10001)ICNT
10001 FORMAT('Warning from NTRPLN: no interpolation points',
     :' lie within the solution domain')
      RETURN
      END

      SUBROUTINE NTRPLN2(S,RLINE,NL,SS0,XP1,YP1,XP2,YP2,
     :                   TBXOFF,TBYOFF,VOFFSET,
     1                   EX,EY,LEM,NOR,NE,NUP,NN,NP3,IVERBOSE)
C
C
C    Routine to interpolate the nodal values from the finite
C     element mesh onto a line traversing the solution region
C     NTRPLN uses quadratic interpolation (e.g. velocity)
C     NTRPLN2 uses linear interpolation (e.g. pressure)
C
      DIMENSION EX(NUP),EY(NUP),LEM(6,NE),NOR(NUP)
      DIMENSION S(NN),RLINE(NL)
      DOUBLE PRECISION BY(3),CX(3),CNL,COOL(3),TRI,X1,X2,X3,Y1,Y2,Y3
C     DOUBLE PRECISION A0(3),DPX1,DPX2,DPX3,DPY1,DPY2,DPY3
      DATA EPSL/1.0E-3/
      SAVE EPSL
C
C   Zero the arrays
C   RLINE(NL) contains the intrpolated values on the line
C   ILINE(NL) not used in this routine
C   Line runs between (XP1,YP1) and (XP2,YP2)
C
      DO I=1,NL
        RLINE(I)=SS0
      ENDDO
      DXP=(XP2-XP1)/FLOAT(NL-1)
      DYP=(YP2-YP1)/FLOAT(NL-1)
      EPS=0.01*SQRT(DXP*DXP+DYP*DYP)
C
C    Look at each element in turn
C
      DO 60 N=1,NE
C
C     Calculate the geometrical coefficients
C
        XMIN=999.
        YMIN=999.
        XMAX=-999.
        YMAX=-999.
        DO 20 K=1,3
          K2=MOD(K,3)+1
          K3=MOD(K+1,3)+1
          LK=NOR(LEM(K,N))
          LK2=NOR(LEM(K2,N))
          LK3=NOR(LEM(K3,N))
          X1=EX(LK)
          X2=EX(LK2)
          X3=EX(LK3)
          Y1=EY(LK)
          Y2=EY(LK2)
          Y3=EY(LK3)
          XMIN=MIN(XMIN,X1)
          XMAX=MAX(XMAX,X1)
          YMIN=MIN(YMIN,Y1)
          YMAX=MAX(YMAX,Y1)
C     DPX1=X1
C     DPX2=X2
C     DPX3=X3
C     DPY1=Y1
C     DPY2=Y2
C     DPY3=Y3
C     A0(K)=DPX2*DPY3 - DPX3*DPY2
          BY(K)=Y2    - Y3
          CX(K)=X3    - X2
   20   CONTINUE
C
C    TRI is twice the area of the triangle element
C
C     TRI=DPX2*DPY3-DPX3*DPY2+DPX3*DPY1-DPX1*DPY3+DPX1*DPY2-DPX2*DPY1
      TRI=CX(3)*BY(2)-CX(2)*BY(3)
C
C     Now look for the points on the line that are inside the tri.
C
        DO 80 J=1,NL
          XP=XP1+DXP*FLOAT(J-1)
          YP=YP1+DYP*FLOAT(J-1)
          IF((YP.GT.YMAX+EPS).OR.(YP.LT.YMIN-EPS).OR.
     1       (XP.GT.XMAX+EPS).OR.(XP.LT.XMIN-EPS))GO TO 80
C
C     Calculate the natural coordinates
C
          DO 45 K1=1,3
C           CNL=(A0(K1) + XP*BY(K1) + YP*CX(K1))/TRI
            K2=MOD(K,3)+1
            LK2=NOR(LEM(K2,N))
            CNL=(CX(K1)*(YP-EY(LK2)) + BY(K1)*(XP-EX(LK2)))/TRI
            IF((CNL.GT.1.0+EPSL).OR.(CNL.LT.-EPSL))GO TO 80
            COOL(K1)=CNL
   45     CONTINUE
C
C     If we reach here the point is within the triangle; Interpolate !
C
          SI=0.0
          DO K=1,3
            LK=NOR(LEM(K,N))
            SK=S(LK)
            SI=SI + SK*COOL(K)
          ENDDO
          RLINE(J)=SI
   80   CONTINUE
   60 CONTINUE
      IF (IVERBOSE.NE.0) THEN
        WRITE(6,*)' data values from NTRPLN2, between '
        WRITE(6,*)'    (',XP1,',',YP1,') and (',XP2,',',YP2,')'
      END IF
      RETURN
      END

      SUBROUTINE NTRPLND(DENS,RLINE,NL,BGD,XP1,YP1,XP2,YP2,
     :                   TBXOFF,TBYOFF,VOFFSET,EX,EY,LEM,
     :                   NOR,ILINE,NE,NUP,NP3,IVERBOSE)
C
C    Routine to interpolate the nodal values from the finite
C     element mesh onto a line traversing the solution region
C     NTRPLN uses quadratic interpolation (e.g. velocity)
C     NTRPLN2 uses linear interpolation (e.g. pressure)
C     NTRPLND interpolates density array from integration pts
C
      DIMENSION EX(NUP),EY(NUP),LEM(6,NE),NOR(NUP)
      DIMENSION DENS(7,NE),RLINE(NL),ILINE(NP3),BNODEV(6)
      DOUBLE PRECISION BY(3),CX(3),CNL,COOL(3),TRI,X1,X2,X3,Y1,Y2,Y3
C     DOUBLE PRECISION A0(3),DPX1,DPX2,DPX3,DPY1,DPY2,DPY3
      DATA EPSL/1.0E-3/
      SAVE EPSL
C
C   Zero the arrays
C   RLINE(NL) contains the intrpolated values on the line
C   ILINE(NL) contains a switch indicating a value in RLINE
C   Line runs between (XP1,YP1) and (XP2,YP2)
C
      DO I=1,NL
        ILINE(I)=0
        RLINE(I)=BGD
      ENDDO
      DXP=(XP2-XP1)/FLOAT(NL-1)
      DYP=(YP2-YP1)/FLOAT(NL-1)
      EPS=0.01*SQRT(DXP*DXP+DYP*DYP)
C
C    Look at each element in turn
C
      DO 60 N=1,NE
C
C     Calculate the geometrical coefficients
C
        XMIN=999.
        YMIN=999.
        XMAX=-999.
        YMAX=-999.
        DO 20 K=1,3
          K2=MOD(K,3)+1
          K3=MOD(K+1,3)+1
          LK=NOR(LEM(K,N))
          LK2=NOR(LEM(K2,N))
          LK3=NOR(LEM(K3,N))
          X1=EX(LK)
          X2=EX(LK2)
          X3=EX(LK3)
          Y1=EY(LK)
          Y2=EY(LK2)
          Y3=EY(LK3)
          XMIN=MIN(XMIN,X1)
          XMAX=MAX(XMAX,X1)
          YMIN=MIN(YMIN,Y1)
          YMAX=MAX(YMAX,Y1)
C         DPX1=X1
C         DPX2=X2
C         DPX3=X3
C         DPY1=Y1
C         DPY2=Y2
C         DPY3=Y3
C         A0(K)=DPX2*DPY3 - DPX3*DPY2
          BY(K)=Y2    - Y3
          CX(K)=X3    - X2
   20   CONTINUE
C
C    TRI is twice the area of the triangle element
C
C     TRI=DPX2*DPY3-DPX3*DPY2+DPX3*DPY1-DPX1*DPY3+DPX1*DPY2-DPX2*DPY1
        TRI=CX(3)*BY(2)-CX(2)*BY(3)
C
C     Now look for the points on the line that are inside the tri.
C
        XP = XP1-DXP
        YP = YP1-DYP
        XOFF=0.0
        YOFF=0.0
        VOFF=0.0
        DO 80 J=1,NL
          XP=XP+DXP-XOFF
          YP=YP+DYP-YOFF
          IF((YP.GT.YMAX+EPS).OR.(YP.LT.YMIN-EPS).OR.
     1       (XP.GT.XMAX+EPS).OR.(XP.LT.XMIN-EPS))GO TO 80
C
C     Calculate the natural coordinates
C
          DO 45 K1=1,3
C           CNL=(A0(K) + XP*BY(K) + YP*CX(K))/TRI
            K2=MOD(K1,3)+1
            LK2=NOR(LEM(K2,N))
            CNL=(CX(K1)*(YP-EY(LK2)) + BY(K1)*(XP-EX(LK2)))/TRI
            IF((CNL.GT.1.0+EPSL).OR.(CNL.LT.-EPSL))GO TO 80
            COOL(K1)=CNL
   45     CONTINUE
C
C    If we reach here the point is within the triangle; Interpolate !
C
          DO 25 K=1,3
            KM=MOD(K+1,3)+4
            BNODEV(K)=1.83095*DENS(K,N)-1.5*DENS(7,N)+
     1                                    0.669052*DENS(KM,N)
            BNODEV(KM)=0.0581402*DENS(K,N)-0.375*DENS(7,N)+
     1                                    1.31686*DENS(KM,N)
   25     CONTINUE
          SI=0.0
          DO K=1,6
            SK=BNODEV(K)+VOFF
            IF(K.LE.3)THEN
              K2=MOD(K,3)+1
              K3=MOD(K+1,3)+1
              CNK=COOL(K)*(COOL(K) - COOL(K2) - COOL(K3))
            ELSE
              K2=K-3
              K3=MOD(K+1,3)+1
              CNK=4.0*COOL(K2)*COOL(K3)
            END IF
            SI=SI + SK*CNK
          ENDDO
          RLINE(J)=SI
          ILINE(J)=N
   80   CONTINUE
   60 CONTINUE
C
      ICNT = 0
      DO 100 I=1,NL
        IF (ILINE(I).NE.0) ICNT = ICNT+1
 100  CONTINUE
      IF(ICNT.EQ.0)
     :  WRITE(6,*)'Warning: no data found on interpolation line'
C
C  Output moved to profile operation instead of interpolate operation
C
C     IF (IVERBOSE.NE.0) THEN
C       WRITE(6,1001)ICNT,XP1,YP1,XP2,YP2
C1001   FORMAT(I5,' data between (',F10.4,','F10.4,') and (',
C    :        F10.4,',',F10.4,')')
C     END IF
C
      RETURN
      END

      SUBROUTINE CALCPROFIL(XDATA,YDATA,IDATA,NL,XDATMIN,XDATMAX,
     1          YDATMIN,YDATMAX,HMIN,VMIN,HMAX,VMAX,IVERBOSE) 
      DIMENSION XDATA(NL),YDATA(NL),IDATA(NL)

C
C    determine length and scale factors
C
      AH=(HMAX-HMIN)/(XDATMAX-XDATMIN)
      BH=HMAX-AH*XDATMAX
      AV=(VMAX-VMIN)/(YDATMAX-YDATMIN)
      BV=VMAX-AV*YDATMAX
C
C     Find the start of the line
C
C     IF(IDASH.NE.0)CALL DASHLN(IDASH,1)
      DO 50 K=1,NL
      IF(IDATA(K).NE.0)THEN
      HCO=AH*XDATA(K)+BH
      VCO=AV*YDATA(K)+BV
      CALL PLOTU(HCO,VCO,3)
      KSTART=K+1
      GO TO 60
      END IF
   50 CONTINUE
   60 CONTINUE
      DO 100 J=KSTART,NL
      HCO=AH*XDATA(J)+BH
      VCO=AV*YDATA(J)+BV
      IF(IDATA(J).NE.0)THEN
      CALL PLOTU(HCO,VCO,2)
      ELSE
      CALL PLOTU(HCO,VCO,3)
      END IF
  100 CONTINUE
      RETURN
      END

      SUBROUTINE EDGES(EDGEPTS,YLWR,YUPR,EX,EY,PTS,IBNGH,IBC,
     1                                      NOR,NUP,NE,NBP,NL)
C
C   this routine finds the endpoints (in X) of a set of NL
C   lines which are needed to construct a 2D profile in which
C   integrals wrt to X will be plotted against Y.  
C
C   X and Y may be swapped in the calling routine.
C
C   The routine finds the min and max coordinates of the places
C   where each line intersects the perimeter of the region and
C   determines an average value along each line to construct a
C   profile in the orthogonal direction.
C
      DIMENSION EDGEPTS(4,NL)
      DIMENSION EX(NUP),EY(NUP)
      DIMENSION NOR(NUP),IBC(NBP),IBNGH(NBP*2)
      DIMENSION PTS(NBP)
      PARAMETER (IX1INDX=1,IY1INDX=2,IX2INDX=3,IY2INDX=4)

      EPS=1.0E-5
      DY=(YUPR-YLWR)/FLOAT(NL-1)
C
C     for each line, set the Y-value and find intersections with boundary
C
      DO 10 J = 1,NL
        YVAL=YLWR+FLOAT(J-1)*DY
        XMIN=1.E10
        XMAX=-1.E10
        NUMPTS = 0
C
C     check where the perimeter of the solution region crosses the line
C
        DO 20 I = 1,NBP
C
C    FIX - check for IBCTYPE of 11 or 13 (internal fault nodes)
C    this comment suggests unresolved problem if fault present - not clear
C
C         IF (IBCTYP(J).EQ.11.OR.IBCTYP(J).EQ.13) GOTO 20
          NODE = NOR(IBC(I))
          NODENEXT = NOR(IBNGH(I+NBP))
          X1=EX(NODE)
          Y1=EY(NODE)
          X2=EX(NODENEXT)
          Y2=EY(NODENEXT)
          IF((Y1-YVAL)*(Y2-YVAL).LE.0.0)THEN
C
C    we have crossed, or at least hit the perimeter, store the intersection
C    coordinates
C
            IF(ABS(Y2-Y1).GT.1.E-6)THEN
              XVAL=X1+(YVAL-Y1)*(X2-X1)/(Y2-Y1)
            ELSE
              XVAL=X1
            ENDIF
            IF(XVAL.LT.XMIN)XMIN=XVAL
            IF(XVAL.GT.XMAX)XMAX=XVAL
            NUMPTS = NUMPTS+1 ! multiple crossings accepted
          END IF
  20    CONTINUE
C       WRITE(*,*)'in EDGES, J, YVAL, XMIN, XMAX =',
C    :J, YVAL, XMIN, XMAX
C
C    if the line didn't intersect any part of the solution region
C
        IF(XMAX.LT.XMIN)THEN
          EDGEPTS(IX1INDX,J) = 0.0
          EDGEPTS(IX2INDX,J) = 0.0
C
C    store coordinates of line beginning and end
C
        ELSE
          EDGEPTS(IX1INDX,J) = XMIN
          EDGEPTS(IX2INDX,J) = XMAX
        ENDIF
        EDGEPTS(IY1INDX,J) = YVAL
        EDGEPTS(IY2INDX,J) = YVAL
  10  CONTINUE
C
      RETURN
      END

      SUBROUTINE INTEGRAT(RLINE,ILINE,NL,HLEN,AREA,LABEL,LEN)
      CHARACTER LABEL*10,LLBL*10
      DIMENSION RLINE(NL),ILINE(NL)
C
      HINC=HLEN/(NL-1)
      DO 50 K=1,NL
        IF(ILINE(K).NE.0)THEN
          KSTART=K
          GO TO 60
        END IF
   50 CONTINUE
   60 CONTINUE
      KEND=KSTART
      DO 70 K=KSTART,NL
        IF(ILINE(K).NE.0)THEN
          KEND=K
        ELSE
          GO TO 80
        END IF
   70 CONTINUE
   80 CONTINUE

      NUMSEG = KEND-KSTART
      AREA=0.0
      AREA2=0.0
      IF (NUMSEG.GE.3) THEN
        IF (MOD(NUMSEG,2).NE.0) THEN
          AREA1 = 3.0/8.0*HINC*(RLINE(KEND-3) + 3.0*RLINE(KEND-2)
     1              + 3.0*RLINE(KEND-1) + RLINE(KEND))
          KEND = KEND-3
        ELSE  
          AREA1 = 0.0
        ENDIF
        IF (NUMSEG.NE.3) THEN
          DO 90 K=KSTART,KEND-2,2
            AREA2 = AREA2 + 1.0/3.0*HINC*(RLINE(K) +
     1                   4.0*RLINE(K+1) + RLINE(K+2))
  90      CONTINUE
        ENDIF
        AREA = AREA1 + AREA2
      ELSE
C       WRITE(LUW,*)' only ',NUMSEG,' data values in profile'
        WRITE(6,*)' only ',NUMSEG,' data values in profile'
        RETURN
      ENDIF
      LLBL(1:10) = '          '
      IF (LEN.GT.0) THEN
        IF (LEN.GT.10) LEN = 10
        LLBL(1:LEN) = LABEL(1:LEN)
C       WRITE(LUW,10001)LLBL,AREA
        WRITE(6,10001)LLBL,AREA
      ENDIF
10001 FORMAT(A10,' Integral under profile: ',G12.6)
      RETURN
      END

      SUBROUTINE INTEGRATTRPZ(RLINE,ILINE,NL,HLEN,AREA,LABEL,LEN)
      CHARACTER LABEL*30,LLBL*30
      DIMENSION RLINE(NL),ILINE(NL)
C
      HINC=HLEN/FLOAT(NL-1)
      KSTART=1
      DO 50 K=1,NL
        IF(ILINE(K).NE.0)THEN
          KSTART=K
          GO TO 60
        END IF
   50 CONTINUE
   60 CONTINUE
      KEND=KSTART
      DO 70 K=KSTART,NL
        IF(ILINE(K).NE.0)THEN
          KEND=K
        ELSE
          GO TO 80
        END IF
   70 CONTINUE
   80 CONTINUE
      AREA=0.0
      IF ((KEND-KSTART).GT.0) THEN
        AREA = AREA + 0.5*(RLINE(KSTART)+RLINE(KEND))
        DO 90 K=KSTART+1,KEND-1
          AREA = AREA + RLINE(K)
  90    CONTINUE
        AREA=AREA*HINC
C       AVG=AREA/HLEN
      ELSE
        WRITE(6,*)' less than 2 values in profile'
        RETURN
      ENDIF
      LLBL(1:10) = '          '
      IF (LEN.GT.0) THEN
        IF (LEN.GT.10) LEN = 10
        LLBL(1:LEN) = LABEL(1:LEN)
C       WRITE(LUW,10001)LLBL,AREA
C       WRITE(6,10001)LLBL,AREA
      ENDIF
10001 FORMAT(A10,' Integral under profile: ',G12.6)
      RETURN
      
      END
