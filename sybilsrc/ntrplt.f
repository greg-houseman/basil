C*--------------------------------------------------------------------
C*    Basil / Sybil:   ntrplt.f  1.1  1 October 1998
C*
C*    Copyright (c) 1997 by G.A. Houseman, T.D. Barr, & L.A. Evans
C*    See README file for copying and redistribution conditions.
C*--------------------------------------------------------------------
C
      SUBROUTINE NTRPLT(S,RMESH,SS0,XCMIN,YCMIN,XCMAX,YCMAX,
     :                          TBXOFF,TBYOFF,VOFFSET,
     :                          EX,EY,LEM,NOR,IHELP,NE,NUP,NP3,
     :                          NX3,NY3,MESHNUM,IVERBOSE)
C
C    Routine to interpolate the nodal values from the finite
C     element mesh onto a regular Cartesian mesh for contouring
C
      DIMENSION S(NUP),RMESH(NP3)
      DIMENSION EX(NUP),EY(NUP),LEM(6,NE),NOR(NUP),IHELP(NP3)
      DOUBLE PRECISION BY(3),CX(3),CNL,COOL(3),TRI,X1,X2,X3,Y1,Y2,Y3
C     DOUBLE PRECISION A0(3),DPX1,DPX2,DPX3,DPY1,DPY2,DPY3
C
      HSX=(XCMAX-XCMIN)/FLOAT(NX3-3)
      HSY=(YCMAX-YCMIN)/FLOAT(NY3-3)
      EPSL=1.e-2/FLOAT(NX3-3)
      EPSX=HSX*0.5
      EPSY=HSY*0.5
C
C   Zero the arrays
C   IHELP is an index array which records whether an interpolation
C   grid point is actually within the finite element mesh.  It contains
C   the relevant element number - or else zero.
C
      DO I=1,NP3
        IHELP(I)=0
        RMESH(I)=SS0
      ENDDO
C
C    Look at each element in turn
C
      DO 80 N=1,NE
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
          IF(X1.LT.XMIN)XMIN=X1
          IF(X1.GT.XMAX)XMAX=X1
          IF(Y1.LT.YMIN)YMIN=Y1
          IF(Y1.GT.YMAX)YMAX=Y1
C     DPX1=X1
C     DPX2=X2
C     DPX3=X3
C     DPY1=Y1
C     DPY2=Y2
C     DPY3=Y3
C     A0(K)=DPX2*DPY3 - DPX3*DPY2
          BY(K)=Y2 - Y3
          CX(K)=X3 - X2
   20   CONTINUE
C
C    TRI is twice the area of the triangle element
C
        TRI=CX(3)*BY(2)-CX(2)*BY(3)
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
C   XCMIN,XCMAX map to I = 2,NX3-1
C   YCMIN,YCMAX map to J = 2,NY3-1
          JMAX=(YMAX+EPSY-YCMIN)/HSY +2
          IF (JMAX.LT.2) GO TO 75
          IF (JMAX.GT.NY3) JMAX=NY3
          JMIN=(YMIN-EPSY-YCMIN)/HSY +2
          IF (JMIN.LT.2) JMIN=2
          IF (JMIN.GT.NY3-1) GO TO 75       !changed
C         IF ((JMAX-JMIN).LT.2) JMAX=JMAX+1 !this caused a crash
          IMAX=(XMAX+EPSX-XCMIN)/HSX +2
          IF (IMAX.LT.2) GO TO 75
          IF (IMAX.GT.NX3) IMAX=NX3
          IMIN=(XMIN-EPSX-XCMIN)/HSX +2
          IF (IMIN.LT.2) IMIN=2
          IF (IMIN.GT.NX3-1) GO TO 75       !changed
          YP=YCMIN + (JMIN-3)*HSY - YOFF
          DO 70 J=JMIN,JMAX
            YP=YP+HSY
            XP=XCMIN + (IMIN-3)*HSX - XOFF
            DO 50 I=IMIN,IMAX
              XP=XP+HSX
C
C     Calculate the natural coordinates
C
              DO 45 K1=1,3
                K2=MOD(K1,3)+1
                LK2=NOR(LEM(K2,N))
                CNL=(CX(K1)*(YP-EY(LK2)) + BY(K1)*(XP-EX(LK2)))/TRI
                IF((CNL.GT.1.0+EPSL).OR.(CNL.LT.-EPSL))GO TO 50
                COOL(K1)=CNL
   45         CONTINUE
C
C     If we reach here the point is within the triangle, so
C         Interpolate !
C
              SI=0.0
              DO K=1,6
                LK=LEM(K,N)
                IF(LK.LT.0)LK=-LK
                SK=S(LK)+VOFF
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
              ENDDO                 ! on K
              IJ=(J-1)*NX3 + I
              RMESH(IJ)=SI
              IHELP(IJ)=N
   50       CONTINUE
   70     CONTINUE
   75   CONTINUE
   80 CONTINUE
      RETURN
      END

      SUBROUTINE NTRPLD(DENS,RMESH,BGD,XCMIN,YCMIN,XCMAX,YCMAX,
     :                          TBXOFF,TBYOFF,VOFFSET,
     :                          EX,EY,LEM,NOR,IHELP,NE,NUP,NP3,
     :                          NX3,NY3,MESHNUM,IVERBOSE)
C
C    Routine to interpolate the density values from the integration
C     points of the mesh onto a regular Cartesian mesh for contouring
C
      DIMENSION DENS(7,NE),RMESH(1),BNODEV(6)
      DIMENSION EX(NUP),EY(NUP),LEM(6,NE),NOR(NUP),IHELP(NP3)
      DOUBLE PRECISION BY(3),CX(3),CNL,COOL(3),TRI,X1,X2,X3,Y1,Y2,Y3
C     DOUBLE PRECISION A0(3),DPX1,DPX2,DPX3,DPY1,DPY2,DPY3
C     DATA EPSL/1.0E-3/
C     SAVE EPSL
C
      HSX=(XCMAX-XCMIN)/FLOAT(NX3-3)
      HSY=(YCMAX-YCMIN)/FLOAT(NY3-3)
      EPSL=1.e-2/FLOAT(NX3-3)
      EPSX=HSX*0.5
      EPSY=HSY*0.5
C
C   Zero the arrays
C   IHELP is an index array which records whether an interpolation
C   grid point is actually within the finite element mesh.  It contains
C   the relevant element number - or else zero.
C
      NP=NX3*NY3
      DO I=1,NP
        IHELP(I)=0
        RMESH(I)=BGD
      ENDDO
C
C    Look at each element in turn
C
      DO 80 N=1,NE
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
          IF(X1.GT.XMAX)XMAX=X1
          IF(Y1.LT.YMIN)YMIN=Y1
          IF(Y1.GT.YMAX)YMAX=Y1
          BY(K)=Y2    - Y3
          CX(K)=X3    - X2
   20   CONTINUE
C
C    TRI is twice the area of the triangle element
C
        TRI=CX(3)*BY(2)-CX(2)*BY(3)
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
          JMAX=(YMAX+EPSY-YCMIN)/HSY +2
          IF (JMAX.LT.2) GO TO 75
          IF (JMAX.GT.NY3) JMAX=NY3
          JMIN=(YMIN-EPSY-YCMIN)/HSY +2
          IF (JMIN.LT.2) JMIN=2
          IF (JMIN.GT.NY3) GO TO 75
          IMAX=(XMAX+EPSX-XCMIN)/HSX +2
          IF (IMAX.LT.2) GO TO 75
          IF (IMAX.GT.NX3) IMAX=NX3
          IMIN=(XMIN-EPSX-XCMIN)/HSX +2
          IF (IMIN.LT.2) IMIN=2
          IF (IMIN.GT.NX3) GO TO 75
          YP=YCMIN + (JMIN-3)*HSY - YOFF
          DO 70 J=JMIN,JMAX
            YP=YP+HSY
            XP=XCMIN + (IMIN-3)*HSX - XOFF
            DO 50 I=IMIN,IMAX
              XP=XP+HSX
C
C     Calculate the natural coordinates
C
              DO 45 K1=1,3
C               CNL=(A0(K) + XP*BY(K) + YP*CX(K))/TRI
                K2=MOD(K1,3)+1
                LK2=NOR(LEM(K2,N))
                CNL=(CX(K1)*(YP-EY(LK2)) + BY(K1)*(XP-EX(LK2)))/TRI
                IF((CNL.GT.1.0+EPSL).OR.(CNL.LT.-EPSL))GO TO 50
                COOL(K1)=CNL
   45         CONTINUE
C
C     If we reach here the point is within the triangle, so: Interpolate !
C     First calculate the boundary node values from extrapolation of
C     integration point values
C
              DO 30 K=1,3
                KM=MOD(K+1,3)+4
                BNODEV(K)=1.83095*DENS(K,N)-1.5*DENS(7,N)+
     1                                       0.669052*DENS(KM,N)
                BNODEV(KM)=0.0581402*DENS(K,N)-0.375*DENS(7,N)+
     1                                       1.31686*DENS(KM,N)
   30         CONTINUE
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
              ENDDO                 ! on K
              IJ=(J-1)*NX3 + I
              RMESH(IJ)=SI
              IHELP(IJ)=N
   50       CONTINUE
   70     CONTINUE
   75   CONTINUE
   80 CONTINUE
      RETURN
      END
C
C   subroutine NTRPLT2 deleted as apparently unused, saved in 
C   file ntrplt2.f in case later needed, but remove after 12/2014
