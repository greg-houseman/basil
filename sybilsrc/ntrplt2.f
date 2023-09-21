
      SUBROUTINE NTRPLT2(S,RMESH,SS0,XCMIN,YCMIN,XCMAX,YCMAX,
     :                          TBXOFF,TBYOFF,VOFFSET,
     :                          EX,EY,LEM,NOR,NE,NUP,NP3,
     :                          NX3,NY3,MESHNUM,IVERBOSE)
C
C    Routine to interpolate the nodal values from the finite
C     element mesh onto a regular Cartesian mesh for contouring
C
      SAVE EPS
      DIMENSION S(1),RMESH(1)
      DIMENSION EX(NUP),EY(NUP),LEM(6,NE),NOR(NUP)
C     DIMENSION BY(3),CX(3),COOL(3)
      DOUBLE PRECISION BY(3),CX(3),COOL(3)
      DOUBLE PRECISION A0(3),TRI,DPX1,DPX2,DPX3,DPY1,DPY2,DPY3
      DATA EPS/1.0E-3/
C
      HSX=(XCMAX-XCMIN)/FLOAT(NX3-3)
      HSY=(YCMAX-YCMIN)/FLOAT(NY3-3)
C
C   Zero the arrays
C   grid point is actually within the finite element mesh.  It contains
C   the relevant element number - or else zero.
C
      NP=NX3*NY3
      DO 10 I=1,NP
   10 RMESH(I)=SS0
C
C    Look at each element in turn
C
      DO 80 N=1,NE
C
C     Calculate the geometrical coefficients
C
   12 XMIN=999.
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
      DPX1=X1
      DPX2=X2
      DPX3=X3
      DPY1=Y1
      DPY2=Y2
      DPY3=Y3
      A0(K)=DPX2*DPY3 - DPX3*DPY2
      BY(K)=Y2    - Y3
      CX(K)=X3    - X2
   20 CONTINUE
C
C    TRI is twice the area of the triangle element
C
C     TRI=X2*Y3 - X3*Y2 + X3*Y1 - X1*Y3 + X1*Y2 - X2*Y1
      TRI=DPX2*DPY3-DPX3*DPY2+DPX3*DPY1-DPX1*DPY3+DPX1*DPY2-DPX2*DPY1
C
C     Now look for the points in the mesh that are inside the tri.
C
      XOFF=0.0
      YOFF=0.0
      VOFF=0.0
      DO 75 IDUP=1,MESHNUM
        IF (IDUP.GT.1) THEN
          XOFF=XOFF+TBXOFF
          YOFF=YOFF+TBYOFF
          XMIN = XMIN + XOFF
          XMAX = XMAX + XOFF
          YMIN = YMIN + YOFF
          YMAX = YMAX + YOFF
          VOFF = VOFF + VOFFSET
        END IF
C set jmin,jmax before looping j=jmin,jmax
        JMAX=(YMAX+EPS-YCMIN)/HSY +2
        IF (JMAX.LT.2) GO TO 80
        IF (JMAX.GT.NY3) JMAX=NY3
        JMIN=(YMIN-EPS-YCMIN)/HSY +2
        IF (JMIN.LT.2) JMIN=2
        IF (JMIN.GT.NY3) GO TO 80
        IMAX=(XMAX+EPS-XCMIN)/HSX +2
        IF (IMAX.LT.2) GO TO 80
        IF (IMAX.GT.NX3) IMAX=NX3
        IMIN=(XMIN-EPS-XCMIN)/HSX +2
        IF (IMIN.LT.2) IMIN=2
        IF (IMIN.GT.NX3) GO TO 80
        YP=YCMIN + (JMIN-3)*HSY - YOFF
        DO 60 J=JMIN,JMAX
        YP=YP+HSY
        XP=XCMIN + (IMIN-3)*HSX - XOFF
        DO 50 I=IMIN,IMAX
        XP=XP+HSX
C
C     Calculate the natural coordinates
C
        DO 45 K=1,3
          CNL=(A0(K) + XP*BY(K) + YP*CX(K))/TRI
          IF((CNL.GT.1.0+EPS).OR.(CNL.LT.-EPS))GO TO 50
   45   COOL(K)=CNL
C
C     If we reach here the point is within the triangle, so
C         Interpolate !
C
              SI=0.0
              DO 35 K=1,3
                LK=NOR(LEM(K,N))
                SK=S(LK) + VOFF
                SI=SI + SK*COOL(K)
   35         CONTINUE
              IJ=(J-1)*NX3 + I
              RMESH(IJ)=SI
   50       CONTINUE
   60     CONTINUE
   75   CONTINUE
   80 CONTINUE
      RETURN
      END
